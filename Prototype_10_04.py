#!/usr/bin/env python
import edlib as ed
import pandas as pd
import seaborn as sns
import numpy as np
from Bio import SeqIO
import re
import os
import pysam as ps
import matplotlib.pyplot as plt
import pylab as plt2
import time
from collections import defaultdict
from matplotlib.patches import Rectangle
from collections import Counter
import pybedtools
#from pybedtools import BedTool

def CIGAR_len(cigar_str):
    """Returns the actual length of the sequence based on the CIGAR int values"""
    count = 0
    for i in re.findall(r'(\d+)([A-Z]{1})', cigar_str):
        if i[-1] == "M" or i[-1] == "I" or i[-1] == "H":
            count += int(i[0])
        elif i[-1] == "D":
            count -= int(i[0])
    return count

def CIGAR_pos(cig_str,char_str):
    """Returns the position in the Cigar string for the given character string"""
    return [pos for pos, char in enumerate(cig_str) if char == char_str]

def IS_SA(read,reverse_str,mapQ):
    """returns true if the read in the bamfile is soft-clipped"""
    if read.cigar[0][0] == 4:
        if read.has_tag("SA"):
            if read.is_reverse == reverse_str:
                if read.mapping_quality == mapQ:
                    return True

def IS_supp(read,reverse_str,mapQ):
    """returns true if the read in the bamfile is supplementary"""
    if read.is_supplementary == True:
        if read.is_reverse == reverse_str:
            if read.mapping_quality == mapQ:
                return True

def Right_dir(prim_pos,supp_pos):
    """Checks if the sequencing direction is from left to right"""
    if prim_pos[0] < supp_pos[0] and prim_pos[-1] < supp_pos[-1]:
        return True

def Left_dir(prim_pos,supp_pos):
    """"Checks if the sequencing direction is from right to left"""
    if prim_pos[0] > supp_pos[0] and prim_pos[-1] > supp_pos[-1]:
        return True

def Right_circ_dir(prim_pos,supp_pos):
    """Checks if the sequencing direction is from left to right
    spanning entire circle at least once"""
    if prim_pos[0] < supp_pos[0] and prim_pos[-1] >= supp_pos[-1]:
        return True

def Left_circ_dir(prim_pos,supp_pos):
    """Checks if the sequencing direction is from right to left
    spanning entire circle at least once"""
    if prim_pos[0] <= supp_pos[0] and prim_pos[-1] > supp_pos[-1]:
        return True

def Prim_dict(bamfile,reverse_str,mapQ,reg,start,end):
    """Creates a dictionary of the primary alignments for the soft-clipped reads"""
    Prim_dict = {}
    for read in bamfile.fetch(reg,start,end,multiple_iterators=True):
        #print(read.cigarstring)
        #checks if soft-clipped
        if IS_SA(read,reverse_str,mapQ) == True:
            pos = ps.AlignedSegment.get_reference_positions(read)

            #creates dict with key being read_id and value being read ref positions
            Prim_dict[read.query_name] = [pos[0],pos[-1]]
    return Prim_dict
bamfile = ps.AlignmentFile("BC05.aln_hg19.bam","rb")
print(len(Prim_dict(bamfile,True,60,"chr1", 243928620, 243938331)))

def Supplement_dict(bamfile,reverse_str,mapQ,reg,start,end):
    """Creates a dictionary of the primary soft-clipped reads and
    their corresponding supplementary alignments"""

    Primary_dict = Prim_dict(bamfile,reverse_str,mapQ,reg,start,end)
    SA_dict = {}

    for read in bamfile.fetch(reg,start,end,multiple_iterators=True):

        #only extract those supp in primary
        if read.query_name in Primary_dict.keys():
            if IS_supp(read,reverse_str,mapQ) == True:

                #extracts positions
                supp_pos = ps.AlignedSegment.get_reference_positions(read)
                prim_pos = Primary_dict[read.query_name]

                # Adding to positions based on direction to the empty dict
                if read.query_name not in SA_dict.keys():

                    #From left to right
                    if Right_dir(prim_pos,supp_pos) == True:
                        SA_dict[read.query_name] = [prim_pos[0],supp_pos[-1]]

                    #From right to left
                    if Left_dir(prim_pos,supp_pos) == True:
                        SA_dict[read.query_name] = [supp_pos[0],prim_pos[-1]]

                    #From left to right once
                    if Right_circ_dir(prim_pos,supp_pos) == True:
                        SA_dict[read.query_name] = [prim_pos[0],prim_pos[-1]]

                    #From right to left once
                    if Left_circ_dir(prim_pos,supp_pos) == True:
                        SA_dict[read.query_name] = [prim_pos[0],prim_pos[-1]]

                # Appends for reads with several supplementary alignments their position
                elif read.query_name in SA_dict.keys():
                    if Right_dir(prim_pos,supp_pos) == True:
                        if supp_pos[-1] not in SA_dict[read.query_name]:
                            SA_dict[read.query_name].append(supp_pos[-1])

                    # From right to left
                    if Left_dir(prim_pos, supp_pos) == True:
                        if prim_pos[-1] not in SA_dict[read.query_name]:
                            SA_dict[read.query_name].append(prim_pos[-1])

                    # From left to right once
                    if Right_circ_dir(prim_pos, supp_pos) == True:
                        if prim_pos[-1] not in SA_dict[read.query_name]:
                            SA_dict[read.query_name].append(prim_pos[-1])

                    # From right to left once
                    if Left_circ_dir(prim_pos, supp_pos) == True:
                        if prim_pos[-1] not in SA_dict[read.query_name]:
                            SA_dict[read.query_name].append(prim_pos[-1])
    return SA_dict

#BREAK IT DOWN AND TRY SOME DIFFERENT COMMANDS TO MAKE SURE YOU UNDERSTAND IT
def reduce_coords(bamfile,reverse_str,mapQ,reg,start,end):
    """ Reduce the number of end coordinates for those reads with several
    supplementary alignments by comparing the position to the other end positions"""

    Coord = Supplement_dict(bamfile, reverse_str, mapQ, reg, start, end)
    # Counter of second list element for 2-element lists
    count = Counter(v[1] for v in Coord.values() if len(v) == 2)
    # Result dict
    reduce_dic = {}
    # Iterate data entries
    for k, v in Coord.items():
        # Modify lists longer than two with at least one element in the counter
        if len(v) > 2 and any(elem in count for elem in v[1:]):
            # Replace list with first element and following element with max count
            v = [v[0], max(v[1:], key=lambda elem: count.get(elem, 0))]
        # Add to result
        reduce_dic[k] = v
    return reduce_dic

def Single_coord(bamfile,reverse_str,mapQ,reg,start,end):
    """ returns merged coordinates so only for some cases"""
    Circle_dict = reduce_coords(bamfile, reverse_str, mapQ, reg, start, end)
    # intermediate dict, but keys should be values of orig dict
    Inter_dict = defaultdict(list)
    if len(Circle_dict.values()) > len(set(tuple(row) for row in Circle_dict.values())):
        for k,v in sorted(Circle_dict.items()):
            Inter_dict[tuple(v)].append(k)
        final_dict = {tuple(v): list(k) for k, v in Inter_dict.items() if len(v) > 1}
        #print("coordinates", final_dict[list(final_dict)[0]],
        #      "supported by %s soft clipped reads" % len(list(final_dict)[0]) + " with " + str(reverse_str))
        return final_dict[list(final_dict)[0]]
    else:
        #print(len(Circle_dict))
        return Circle_dict

def File_identification(file_name,reverse_str,mapQ):
    with open(file_name) as f:
        for line in f:
            line_value = line.strip().split()
            coord = line_value[0]
            start = line_value[1]
            end = line_value[2]
            cov = line_value[3]
            coord_dict = Single_coord(bamfile, reverse_str,mapQ, str(coord), int(start), int(end))
            #print(line)
            #print(bool(coord_dict))
            #print("------------")
            if bool(coord_dict) == True:
                print("dictionary",coord_dict)
            else:
                #those dict empy, possible due to complex circle formation
                continue

def Genome_Cov_identification(bamfile,bam_name,overlap,reverse_str,mapQ):
    Bedfile = pybedtools.example_bedtool(str(os.getcwd()) +"/"+bam_name)
    Cov = Bedfile.genome_coverage(bg=True)
    Merged = Cov.merge(d=overlap)

    for region in Merged:
        Chr = region[0]
        start = region[1]
        end = region[2]
        if len(Single_coord(bamfile, reverse_str, mapQ,str(Chr), int(start), int(end))) == 0:
            continue
        else:
            print(Chr + ":" + start + "-" + end)
            print(Single_coord(bamfile, reverse_str, mapQ, str(Chr), int(start), int(end)))

bamfile = ps.AlignmentFile("BC05.aln_hg19.bam","rb")

#Genome_Cov_identification(bamfile,"BC05.aln_hg19.bam",500,True,60)
#File_identification("BC05.ge_mean5.bdg",True,60)

Complex_dict = {}
#"chr1", 106597836, 106598346
#"chr1", 16626700, 16627600
# "chr1", 243928620, 243938331
"""
for read in bamfile.fetch("chr1", 243928620, 243938331, multiple_iterators=True):
    if read.cigar[0][0] == 4:
        if read.has_tag("SA"):
            print(read.query_name)
            print(read.get_tag("SA"))
"""

print("---------------------")
print(Single_coord(bamfile,False,60,"chr1",16626700,16627600))
print("---------------------")
print(Single_coord(bamfile,False,60,"chr8",51869044,51869381))
print("---------------------")

for read in bamfile.fetch("chr1", 16626700, 16627600, multiple_iterators=True):
    if read.cigar[0][0] == 4:
        if read.has_tag("SA"):
            if read.is_reverse == False:
                print("readname", read.query_name)
                # print(read.get_tag("SA"))
                for Tag in read.get_tag("SA").split(';')[:-1]:
                    # print(Tag)
                    Column_list = Tag.split(',')
                    if Column_list[2] == '-' and Column_list[4] == str(60):
                        print(Column_list)
                        length = CIGAR_len(Column_list[3])
                        #print("length", length)
                        print("start", int(Column_list[1])-1, "end", int(Column_list[1]) + length-1)
                #print("new_read")

"""
        for SA in read.get_tag("SA").split(';'):
            print("et read")
            print(SA)
            print(SA.split(','))
            """
#Strand is either ‘+’ or ‘-’, indicating forward/reverse strand
#pos is 1-based

def TAG_extrat(bamfile):
    for read in bamfile.fetch("chr1", 106597836, 106598346, multiple_iterators=True):
        if read.has_tag("SA"):
            if read.is_reverse == False:
                print("readname",read.query_name)
                #print(read.get_tag("SA"))
                for Tag in read.get_tag("SA").split(';')[:-1]:
                    #print(Tag)
                    Column_list = Tag.split(',')
                    if Column_list[2] == '-' and Column_list[4] == str(60):
                        print("TAG_INFO",Column_list)
                        length = CIGAR_len(Column_list[3])
                        print("length",length)
                        print("start",int(Column_list[1]),"end",int(Column_list[1])+length)
                print("new_read")

            #print(read.get_tag("SA").split(';')[:-1])
