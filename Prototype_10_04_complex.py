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

def IS_SA(read,mapQ):
    """returns true if the read in the bamfile is soft-clipped"""
    if read.cigar[0][0] == 4:
        if read.has_tag("SA"):
            if read.mapping_quality == mapQ:
                return True

def IS_supp(read,mapQ):
    """returns true if the read in the bamfile is supplementary"""
    if read.is_supplementary == True:
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

def Prim_dict(bamfile,mapQ,reg,start,end):
    """Creates a dictionary of the primary alignments for the soft-clipped reads"""
    Prim_dict = {}
    for read in bamfile.fetch(reg,start,end,multiple_iterators=True):
        #print(read.cigarstring)
        #checks if soft-clipped
        if IS_SA(read,mapQ) == True:
            pos = ps.AlignedSegment.get_reference_positions(read)

            #creates dict with key being read_id and value being read ref positions
            Prim_dict[read.query_name] = [pos[0],pos[-1]]
    return Prim_dict

def Supplement_dict(bamfile,mapQ,reg,start,end):
    """Creates a dictionary of the primary soft-clipped reads and
    their corresponding supplementary alignments"""

    Primary_dict = Prim_dict(bamfile,mapQ,reg,start,end)
    SA_dict = {}

    for read in bamfile.fetch(reg,start,end,multiple_iterators=True):

        #only extract those supp in primary
        if read.query_name in Primary_dict.keys():
            if IS_supp(read,mapQ) == True:

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
                    #print("test",SA_dict[read.query_name])
                    #print("test2",SA_dict[read.query_name][1:])
                    #supp_pos[-1] > all(SA_dict[read.query_name][1:])
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

bamfile = ps.AlignmentFile("BC05.aln_hg19.bam","rb")
#print(Supplement_dict(bamfile,60,"chr1", 16626700, 16627600))

#BREAK IT DOWN AND TRY SOME DIFFERENT COMMANDS TO MAKE SURE YOU UNDERSTAND IT
def reduce_coords(bamfile,mapQ,reg,start,end):
    """ Reduce the number of end coordinates for those reads with several
    supplementary alignments by comparing the position to the other end positions"""

    Coord = Supplement_dict(bamfile, mapQ, reg, start, end)
    # Counter of second list element for 2-element lists.
    # I dont have several primary so there is no need for counting them
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

def Single_coord(bamfile,mapQ,reg,start,end):
    """ returns merged coordinates so only for some cases"""
    Circle_dict = reduce_coords(bamfile, mapQ, reg, start, end)
    # intermediate dict, but keys should be values of orig dict
    Inter_dict = defaultdict(list)

    # WRITE SOME PROPER NOTE FOR THIS LINE
    if len(Circle_dict.values()) > len(set(tuple(row) for row in Circle_dict.values())):
        for k,v in sorted(Circle_dict.items()):
            Inter_dict[tuple(v)].append(k)
        final_dict = {tuple(v): list(k) for k, v in Inter_dict.items() if len(v) > 1}
        #print("coordinates", final_dict[list(final_dict)[0]],
        #      "supported by %s soft clipped reads" % len(list(final_dict)[0])))
        #return final_dict[list(final_dict)[0]]
        return final_dict
    else:
        #print(len(Circle_dict))
        return Circle_dict

bamfile = ps.AlignmentFile("BC05.aln_hg19.bam","rb")


#"chr1", 106597836, 106598346
#"chr1", 16626700, 16627600
#"chr1", 243928620, 243938331


def Complex(bamfile,mapQ,reg, start, end):
    # den overskriver igen hvis der er flere supp
    Complex_dict = {}
    for read in bamfile.fetch(reg, start, end, multiple_iterators=True):
        if read.cigar[0][0] == 4:
            if read.has_tag("SA"):
                for Tag in read.get_tag("SA").split(';')[:-1]:
                    Column_list = Tag.split(',')
                    if Column_list[4] == str(mapQ):
                        #print("readname", read.query_name)
                        #print(read.cigarstring)
                        # print(Column_list[0],"start", pos_start, "end", pos_end, "len",length)
                        length = CIGAR_len(Column_list[3])
                        pos_start = int(Column_list[1])-1 #minus 1 as its 1-based
                        pos_end = int(Column_list[1]) + length-1
                        Complex_dict[read.query_name] = [Column_list[0],pos_start,pos_end]
    return Complex_dict

def Complex2(bamfile,mapQ,reg, start, end):
    Sing_dict = Single_coord(bamfile, mapQ, reg, start, end)
    reads_IDs = list(Sing_dict.keys())[0]

    Complex_dict = {}
    test = []
    for read in bamfile.fetch(reg, start, end, multiple_iterators=True):
        if read.query_name in reads_IDs:
            #print("yes")
            #print(read.query_name)
            if IS_SA(read, mapQ) == True:
                #print(read.get_tag("SA").split(';')[:-1])
                Tag = read.get_tag("SA").split(';')[:-1]
                #print("tag",Tag)
                Start_coord = []
                end_coord = []
                total_coord = []
                for Tagelem in Tag:
                    Column_list = Tagelem.split(',')
                    chrom = Column_list[0]
                    length = CIGAR_len(Column_list[3])
                    pos_start = int(Column_list[1]) - 1  # minus 1 as its 1-based
                    pos_end = int(Column_list[1]) + length - 1
                    if chrom == reg and start < pos_start < end and start < pos_end < end:
                        continue
                    elif int(Column_list[4]) == mapQ:
                        Start_coord.append(pos_start)
                        end_coord.append(pos_end)
                        total_coord.append(pos_start)
                        total_coord.append(pos_end)
                        #print("Start", Start_coord)
                        #print(Column_list[0], "start", pos_start, "end", pos_end, "len", length)
                        #print("readname", read.query_name)
                        #print(Column_list[0],"start", pos_start, "end", pos_end, "len",length)
                        Complex_dict[read.query_name] = [Column_list[0],pos_start,pos_end]
                print("Starts",Start_coord)
                print("end", end_coord)
                print("total",total_coord)
                print("len",len(Start_coord))
                group_val = int(len(Start_coord)/len(reads_IDs))
                Total_val = group_val*2
                print("group",[Start_coord[k:k+group_val] for k in range(0,len(Start_coord),group_val)])
                print("group", [total_coord[k:k + Total_val] for k in range(0, len(total_coord), Total_val)])

                #print([Start_coord[k::k+2] for k in range(0,len(Start_coord),2)])
            print("--------------")
    return Complex_dict


Sing_dict = Single_coord(bamfile,60,"chr1", 16626700, 16627600)
print("single dict",Sing_dict)
Sing_dict = Single_coord(bamfile,60,"chr8", 51869043, 51869438)
print("single dict",Sing_dict)

#print(Complex(bamfile,60,"chr8", 51869043, 51869438))

Sing_dict = Single_coord(bamfile,60,"chr8", 51892600, 51893000)
print("single dict",Sing_dict)

#complex = Complex(bamfile,60,"chr1", 16626700, 16627600)
#print("complex",complex)
compl = Complex2(bamfile,60,"chr1", 16626700, 16627600)

test = [[51869043, 51869042, 51869043], [51892636, 51892640, 51892709]]
test2 = [[51869043, 51869561, 51869042, 51869573, 51869043, 51869438], [51892636, 51892979, 51892640, 51892972, 51892709, 51892995]]
for i in test:
    print([min(i),max(i)])

print([[min(i),max(i)] for i in test2])


#print("Complex 2",compl)
#print(Complex(bamfile,60,"chr8", 51869043, 51869573))

#print("testval",complex['3e133f13-9c3d-4245-8732-41436df88436'])
#print(complex['b1c9023c-751e-47e1-aecd-cf33204f3c75'])
"""
print(list(Sing_dict.keys())[0])
for i in range(len(list(Sing_dict.keys())[0])):
    #extracting the different read_ID for final coordinatset for a region
    read_id = list(Sing_dict.keys())[0][i]
    #extract the complex value for corresponding reads
    matching_val = complex[str(read_id)]
    print(matching_val)
    #if Sing_dict != ma

print(Single_coord(bamfile,60,"chr8", 51892648, 51892966))
"""

#for k,v in complex.items():
#    print(v)
#print("Complex_coord_chr1",Complex(bamfile,60,"chr1", 16626700, 16627600))
#print("single_coord_chr8",Single_coord(bamfile,60,"chr8", 51869000, 51869800))
#print("Complex_coord_chr8",Complex(bamfile,60,"chr8", 51869000, 51869800))
"""
print(Complex_dict)
test = Single_coord(bamfile,60,"chr1", 16626700, 16627600)
for k,v in test.items():
    print(k)
    print(len(k))
#"chr8", 51869043, 51869590
"""

#print("prim",Prim_dict(bamfile,60,"chr1", 16626700, 16627600))
#print("complex",Complex_dict)
#print("single_coord",Single_coord(bamfile,60,"chr1", 16626700, 16627600))
#Prim_dict[read.query_name] = [pos[0],pos[-1]]


#Strand is either ‘+’ or ‘-’, indicating forward/reverse strand
#pos is 1-based

"""
        for i in reads_IDs:
            print("iiii",i)
            #print("test",read.query_name)
            if read.query_name == i:
                #print("lool")
                #print(read.query_name)
                if read.cigar[0][0] == 4:
                    #print("lool2")
                    #print(read.query_name)
                    if read.has_tag("SA"):
                        for Tag in read.get_tag("SA").split(';')[:-1]:
                            Column_list = Tag.split(',')
                            chrom = Column_list[0]
                            length = CIGAR_len(Column_list[3])
                            pos_start = int(Column_list[1]) - 1  # minus 1 as its 1-based
                            pos_end = int(Column_list[1]) + length - 1
                            if chrom == reg and start < pos_start < end and start < pos_end < end:
                                continue
                            elif int(Column_list[4]) == mapQ:
                                #print("readname", read.query_name)
                                print(Column_list[0],"start", pos_start, "end", pos_end, "len",length)
                                Complex_dict[read.query_name] = [Column_list[0],pos_start,pos_end]
    return Complex_dict
"""