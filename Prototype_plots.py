#!/usr/bin/env python
import edlib as ed
import pandas as pd
import seaborn as sns
import numpy as np
import pybedtools as bt
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
            if read.mapping_quality >= mapQ:
                return True

def IS_supp(read,mapQ):
    """returns true if the read in the bamfile is supplementary"""
    if read.is_supplementary == True:
        if read.mapping_quality >= mapQ:
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

def Single_coord(bamfile,mapQ,reg,start,end):
    """ returns merged coordinates for circles w. several reads having exact same coordinates. If not
    a dictionary with all reads and potential coordinates are returned"""
    Circle_dict = reduce_coords(bamfile, mapQ, reg, start, end)
    # intermediate dict, but keys should be values of orig dict
    Inter_dict = defaultdict(list)
    #If the number of values are larger than the unique set of values
    if len(Circle_dict.values()) > len(set(tuple(row) for row in Circle_dict.values())):
        #so here the values are the keys
        for k,v in sorted(Circle_dict.items()):
            Inter_dict[tuple(v)].append(k)

        final_dict = {tuple(v): list(k) for k, v in Inter_dict.items() if len(v) > 1}
        max_key, max_value = max(final_dict.items(), key=lambda x: len(set(x[1])))
        chr = [reg]
        new_val =  chr + max_value
        return dict([(max_key,new_val)])

    else:
        if len(Circle_dict.keys()) == 1:
            Inter_dict = defaultdict(list)
            for k, v in sorted(Circle_dict.items()):
                Inter_dict[tuple(v)].append(k)

            final_dict = {tuple(v): list(k) for k, v in Inter_dict.items()}
            Key, Value = max(final_dict.items(), key=lambda x: len(set(x[1])))
            chr = [reg]
            # min and max in case of a single read with several positions
            new_val = chr + [min(Value)] + [max(Value)]
            return dict([(Key,new_val)])
        else:
            return Circle_dict

def density_plot(bamfile,chr,start,end,savename,pos):
    """ plots the distribution of split-reads across a chromosome region """
    coord_Dict = reduce_coords(bamfile, 60, chr, start, end)
    single_pos = Single_coord(bamfile, 60, chr, start, end)
    print("Coord_dict",coord_Dict)
    print("single",single_pos)
    print("list",list(single_pos.values()))
    if len(list(single_pos.values())) != 1:
        coordinates = [item for sublist in list(single_pos.values()) for item in sublist]
        coordinates.sort()
    else:
        coordinates = list(single_pos.values())[0][1:]
    print("coo",coordinates)
    x1 = []
    x2 = []
    #append the positions to list
    for k,v in coord_Dict.items():
        x1.append(v[0])
        x2.append(v[-1])
    print(x1)
    length = coordinates[-1]-coordinates[0]
    colors=["dodgerblue","darkorange"]
    labels = ["Primary alignments","Supplementary alignments"]
    plt.figure()
    plt.hist([x1,x2],bins=20,density=True,color=colors,label=labels)
    plt.axvline(coordinates[0], color='k', linestyle='dashed', linewidth=2,label="Start pos")
    plt.axvline(coordinates[-1], color='r', linestyle='dashed', linewidth=2,label="end pos")
    plt.title("Position distribution across the circular region")
    plt.xlabel("coordinates")
    plt.ylabel("density")
    plt.legend(loc=pos)
    plt.xticks(np.arange(coordinates[0], coordinates[-1]+100,length/6))
    plt.savefig(savename)
    plt.close()

def Softclipped_NoSA_count(bamfile,tag_str):
    """Returns list with length of soft_clipped reads with and without SA"""
    SA_list = []
    No_list = []
    supp_count = 0
    s_count = 0
    SA_count = 0
    NO_SA_count = 0
    for read in bamfile.fetch():
        # the left most softclip read
        if read.is_supplementary == True:
            supp_count += 1
        if read.cigar[0][0] == 4:
            s_count += 1
            if read.has_tag("SA"):
                SA_count += 1
                SA_list.append(read.cigar[0][1])
            if not read.has_tag("SA"):
                NO_SA_count +=1
                No_list.append(read.cigar[0][1])
    print(supp_count,s_count,SA_count,NO_SA_count)
    return SA_list, No_list

def Read_dist(Read_list,bins,label,figlabel,thres1,thres2,xstart,xend):
    plt2.figure()
    plt2.hist(Read_list, bins=bins, alpha=0.5, label=label, color="blue")
    plt2.hist([], color='white',
              label= "{:.2f}".format(sum(i > thres1 for i in Read_list) / len(Read_list) * 100)
                     + "%" + " of soft-clipped part w. length above %s bp" % str(thres1))
    plt2.hist([], color = 'white',
              label="{:.2f}".format(sum(i < thres2 for i in Read_list) / len(Read_list) * 100)
                    + "%" + " of soft-clipped part w. length below %s bp" % str(thres2))
    plt2.title(str(len(Read_list)) + " Soft-clipped reads " + label)
    plt2.axvline(np.mean(Read_list),color='r', linestyle='dashed', linewidth=2,label=
    "Mean length %s" % format(np.mean(Read_list), '.2f'))
    plt2.xlim(xstart,xend)
    plt2.xlabel("Soft-clipped length")
    plt2.ylabel("Read count")
    plt2.legend()
    plt2.savefig(figlabel+".jpg")
    plt2.close()


bamfile = ps.AlignmentFile("BC05.aln_hg19.bam","rb")

density_plot(bamfile,"chr1", 243928620, 243938331,"Test_Dataset_pos_dist.png",'upper center')
#density_plot(bamfile,"chr9", 124876200, 124879026,"Test_Dataset_pos_dist.png",'upper center')
#density_plot(bamfile,"chr2", 161093561, 161141853,"Test_Dataset_pos_dist.png",'upper right')
#density_plot(bamfile,"chrX", 75680381, 75682701,"Test_Dataset_pos_dist.png",'upper center')
#"chr1", 243928620, 243938331

Supp_list, No_supp = Softclipped_NoSA_count(bamfile,"SA")
print(np.mean(Supp_list))
print(np.mean((No_supp)))
#(Read_list,bins,label,figlabel,thres1,thres2,xstart,xend)
#Read_dist(Supp_list,100,"with supplementary alignment","lol",1000,100,0,None)
Read_dist(Supp_list,100,"with supplementary alignment","Supp_dist",1000,100,0,5000)
#Read_dist(No_supp,50,"with supplementary alignment","lol",1000,100,0,None)
Read_dist(No_supp,200,"with supplementary alignment","None_supp_dist",1000,100,0,1000)