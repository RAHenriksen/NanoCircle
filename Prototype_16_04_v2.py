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
import math
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

def Single_coord(bamfile,mapQ,reg,start,end,coord_val):
    """ returns merged coordinates so only for some cases"""
    Circle_dict = reduce_coords(bamfile, mapQ, reg, start, end)
    # intermediate dict, but keys should be values of orig dict
    Inter_dict = defaultdict(list)

    # WRITE SOME PROPER NOTE FOR THIS LINE
    if len(Circle_dict.values()) > len(set(tuple(row) for row in Circle_dict.values())):
        #so here the values are the keys
        for k,v in sorted(Circle_dict.items()):
            Inter_dict[tuple(v)].append(k)
        final_dict = {tuple(v): list(k) for k, v in Inter_dict.items() if len(v) > 1}
        if coord_val == 0:
            return final_dict
        elif coord_val == 1:
            #print("coordinates", final_dict[list(final_dict)[0]],
            #      "supported by %s soft clipped reads" % len(list(final_dict)[0]))
            return final_dict[list(final_dict)[0]]

    else:
        return Circle_dict

def complex_reads(bamfile,mapQ,reg, start, end,reads_IDs):
    """ extract those reads which is found using the single coordinate dictionary"""
    Read_list = []
    for read in bamfile.fetch(reg, start, end, multiple_iterators=True):
        if read.query_name in reads_IDs:
            if IS_SA(read, mapQ) == True:
                Read_list.append(read)
    return Read_list

def Grouping(Coord_list,overlap_bp):
    First = [[Coord_list[0]]]
    for i in Coord_list[1:]:
        # sets the previous and current elem as the first elem
        prev = First[-1]
        current = prev[-1]
        dist = (i-current)

        # if distance between elem in coord_list and current is below 700, then its
        # appended to the previous list
        if dist < overlap_bp:
            prev.append(i)
        else:
            First.append([i])
    return First

def chr_coord_sort(chrlist,coordlist):
    coord_idx = np.argsort(coordlist)
    Coord_sort = [coordlist[i] for i in coord_idx]
    chr_sort = [chrlist[i] for i in coord_idx]
    return Coord_sort,chr_sort

def Complex(bamfile,mapQ,reg, start, end):
    Sing_dict = Single_coord(bamfile, mapQ, reg, start, end,0)
    print("single",Sing_dict)

    Complex_dict = {}

    #one set of coordinates supported by several reads
    if len(Sing_dict) == 1:
        coord1 = list(Sing_dict.values())[0][0]
        coord2 = list(Sing_dict.values())[0][-1]
        print("coord",coord1,coord2)
        reads_IDs = list(Sing_dict.keys())[0]

    #several different set of coordinates
    else:
        reads_IDs = list(Sing_dict.keys())
        print("test",list(Sing_dict.values()))

    #print("Read IDS",reads_IDs)

    #extract the reads from the single coordinate dictionary
    Read_list = complex_reads(bamfile,mapQ,reg, start, end,reads_IDs)
    for read in Read_list[:1]:
        #coord1 = Sing_dict[read.query_name]
        #print("coo",coord1)
        #print("name",read.query_name)
        Tag = read.get_tag("SA").split(';')[:-1]
        Coord_list = []
        chroms = []
        cigar_len = []
        for Tagelem in Tag:
            #splitting the individual aligment info up
            Column_list = Tagelem.split(',')
            chrom = Column_list[0]
            length = CIGAR_len(Column_list[3])
            cigar_len.append(length)
            overlap = max(cigar_len)
            pos_start = int(Column_list[1]) - 1  #1-based
            pos_end = int(Column_list[1]) + length - 1

            # if an supp align is in the circular input region then its not necessary for the complex circle
            if chrom == reg and start - overlap <= pos_start <= end + overlap and start - overlap <= pos_end <= end + overlap:
                print("CONTINUE")
                continue

            elif int(Column_list[4]) == mapQ:
                chroms.append(chrom)
                chroms.append(chrom)
                Coord_list.append(pos_start)
                Coord_list.append(pos_end)
            if Coord_list != []:
                Coord_sort, chr_sort = chr_coord_sort(chroms, Coord_list)
                Complex_dict[read.query_name] = [reg, coord1, coord2]

                if len(Coord_sort) == 2:
                    Complex_dict[read.query_name].extend(
                        [chr_sort[0], min(Coord_sort), max(Coord_sort)])
                else:
                    Grouped_list = Grouping(Coord_sort, overlap)
                    Group_size = [len(i) for i in Grouped_list]
                    for i in range(len(Grouped_list)):
                        chr_group = [chr_sort[x:x + Group_size[i]] for x in range(0, len(chr_sort),Group_size[i])]
                        Complex_dict[read.query_name].extend(
                            [chr_group[i][0], min(Grouped_list[i]), max(Grouped_list[i])])
    if Complex_dict == {}:
        print("not complex")
        return Single_coord(bamfile, mapQ, reg, start, end,1)
    else:
        print("complex")
        return Complex_dict


bamfile = ps.AlignmentFile("BC05.aln_hg19.bam","rb")

"""
print("final coordinates")
#several different circle coordinates
print(Single_coord(bamfile,60,"chr1", 14682287, 14687166,1))
print(Single_coord(bamfile,60,"chr18", 36759381, 36764356,1))
# only one circle coordinate, then just returning the coordinate and not all of the reads
print(Single_coord(bamfile,60,"chr1", 243928620, 243938331,1))
# a region which is in the .bdg file but not supported by any soft-clipped reads
print(Single_coord(bamfile,60,"chr1", 106597836, 106597836,1))
"""

print("complex circles")
#actual complex regions


print(Complex(bamfile,60,"chr1", 14682287, 14687166))
#print(Complex(bamfile,60,"chr1", 106597836, 106598346))

#print(Complex(bamfile,60,"chr1", 16626700, 16627600))
#print(Complex(bamfile,60,"chr8", 51869000, 51869500))
#print(Complex(bamfile,60,"chr7", 75713114, 75717415))

#print(Complex(bamfile,60,"chr1", 243928620, 243938331))
#print(Complex(bamfile,60,"chr8", 122226956, 122231216))

print("---------")
#bamfile = ps.AlignmentFile("BC01.aln_hg19.bam","rb")
#print(Complex(bamfile,60,"chr2", 82083496, 82087081))


"""
for k,v in lol.items():
    print("reads",k)
    print("val",v)
"""

#not a complex region, then just returning the single_coord dictionary
#print(Complex2(bamfile,60,"chr1", 243928620, 243938331))
