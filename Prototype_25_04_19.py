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
        new_val = chr + max_value
        return dict([(max_key,new_val)])

    else:
        print("Circle_dict")
        return Circle_dict

def complex_reads(bamfile,mapQ,reg, start, end,reads_IDs):
    """ extract those reads which is found using the single coordinate dictionary"""
    Read_list = []
    for read in bamfile.fetch(reg, start, end, multiple_iterators=True):
        if read.query_name in reads_IDs:
            if IS_SA(read, mapQ) == True:
                Read_list.append(read)
    return Read_list

def chr_coord_sort(chrlist,coordlist):
    """ Sort a list of chromosomes and their coordinates using the index of the numerically sorted coordinates"""
    coord_idx = np.argsort(coordlist)
    Coord_sort = [coordlist[i] for i in coord_idx]
    chr_sort = [chrlist[i] for i in coord_idx]
    return Coord_sort,chr_sort

def Grouping_chr(Chr_sort,Group_size):
    Grouped_chr = []
    Step = 0
    for i in Group_size:
        Grouped_chr.append(Chr_sort[Step:Step + i])
        Step += i
    return Grouped_chr

def Grouping(Coord_list,overlap_bp):
    """ Groups a given list, for which elements within the overlap range is grouped together"""
    First = [[Coord_list[0]]]
    for i in Coord_list[1:]:
        # sets the previous and current elem as the first elem
        prev = First[-1]
        current = prev[-1]
        dist = abs(i-current)

        # if dist between is below overlapping bp, then it it belong to the same group
        if dist < overlap_bp:
            prev.append(i)
        else:
            First.append([i])
    return First

def Complex(bamfile,mapQ,reg, start, end):
    Sing_dict = Single_coord(bamfile, mapQ, reg, start, end)
    Complex_dict = {}
    Total_coord = []
    Total_chr = []
    Total_overlap = []
    print(Sing_dict)
    if len(Sing_dict) == 1:
        reads_IDs = list(Sing_dict.keys())[0]
        Read_list = complex_reads(bamfile, mapQ, reg, start, end, reads_IDs)
        print("READ_LIST",Read_list)
    else:
        reads_IDs = list(Sing_dict.keys())
        Read_list = complex_reads(bamfile,mapQ,reg, start, end,reads_IDs)


    for read in Read_list:
        if len(Sing_dict) == 1:
            coord1 = list(Sing_dict.values())[0][1]
            coord2 = list(Sing_dict.values())[0][-1]
            Total_chr.extend((reg,reg))
            Total_coord.extend((coord1,coord2))
        else:
            coord1 = Sing_dict[read.query_name][0]
            coord2 = Sing_dict[read.query_name][-1]
            Total_chr.extend((reg,reg))
            Total_coord.extend((coord1, coord2))
        #print("READ",read.query_name)
        #print("length",CIGAR_len(read.cigarstring))
        Tag = read.get_tag("SA").split(';')[:-1]
        #print("tag",Tag)
        Coord_list = []
        chroms = []
        cigar_len = []
        for Tagelem in Tag:
            print("tagelem",Tagelem)
            #splitting the individual aligment info up
            Column_list = Tagelem.split(',')
            chrom = Column_list[0]
            length = CIGAR_len(Column_list[3])
            cigar_len.append(length)
            pos_start = int(Column_list[1]) - 1  #1-based
            pos_end = int(Column_list[1]) + length - 1
            overlap = sum(cigar_len)
            Total_overlap.append(overlap)

            # if the supp align is in between the circular input region then we already know the breakpoint from Sing_dict
            if chrom == reg and start - overlap <= pos_start <= end + overlap and start - overlap <= pos_end <= end + overlap:
                continue

            elif int(Column_list[4]) >= mapQ:
                #creates a coordinate list
                print("name",read.query_name)
                Coord_list.append(pos_start)
                Coord_list.append(pos_end)
                #append chr twice to ensure same length as coord_list
                chroms.append(chrom)
                chroms.append(chrom)
                #print("chomrs",chroms)
                Total_chr.extend((chrom, chrom))
                Total_coord.extend((pos_start, pos_end))

            #print("COORD",Coord_list)
            #print("CHR",chroms)
            #print("GROUPED",Grouping(Coord_list, overlap))
            if Coord_list != []:
                #sorts the chr and coordinates
                print("chroms",chroms)
                Coord_sort, chr_sort = chr_coord_sort(chroms, Coord_list)

                # first entry with the input region
                Complex_dict[read.query_name] = [reg, coord1, coord2]

                if len(Coord_sort) == 2:
                    print("lenght")
                    #one coordinate pair supporting another region
                    Complex_dict[read.query_name].extend(
                        [chr_sort[0], min(Coord_sort), max(Coord_sort)])
                else:
                    Grouped_coord = Grouping(Coord_sort, overlap)
                    Group_size = [len(i) for i in Grouped_coord]
                    Grouped_chr = Grouping_chr(chr_sort,Group_size)
                    for i in range(len(Grouped_coord)):
                        Complex_dict[read.query_name].extend(
                            [Grouped_chr[i][0], min(Grouped_coord[i]), max(Grouped_coord[i])])
        print("----------------------new_read-----------------------")

    if Complex_dict == {}:
        print("The given input region is not forming a complex circle")
        return Single_coord(bamfile, mapQ, reg, start,end)
    else:
        Coord_sort, chr_sort = chr_coord_sort(Total_chr, Total_coord)
        Grouped_coord = Grouping(Coord_sort, max(Total_overlap))
        Group_size = [len(i) for i in Grouped_coord]
        Grouped_chr = Grouping_chr(chr_sort, Group_size)
        print("The given input region forms a complex circle")
        return Complex_dict,Grouped_coord,Grouped_chr

bamfile = ps.AlignmentFile("BC01.aln_hg19.bam","rb")
complex_circ,Total_coord,Total_chr = Complex(bamfile,60,"chr2", 46846937, 46847710)

print("complex",complex_circ)
#Read_Complex_BED(complex_circ,"Complex_test.bed")
print(Total_coord)
print(len(Total_chr))
print(len(complex_circ))

def Complex_circ_BED(dict,coord_full,chr_full,Circ_no,savename):
    d = {}
    tot_len = 0
    for i in range(len(coord_full)):
        if len(d) == 0:
            d["test"] = [chr_full[i][0], min(coord_full[i]), max(coord_full[i]),
                         max(coord_full[i]) - min(coord_full[i])]
            tot_len += (max(coord_full[i]) - min(coord_full[i]))
        else:
            d["test"].extend([chr_full[i][0], min(coord_full[i]), max(coord_full[i]),
                              max(coord_full[i]) - min(coord_full[i])])
            tot_len += (max(coord_full[i]) - min(coord_full[i]))
    print(tot_len)
    df_col = ["Chr", "Start", "End", "Length"]
    rep = len(coord_full)
    first_col = 'Circle No.'
    Coord_Col = [j + "_no._" + str(i) for i in range(1, rep + 1) for j in df_col]

    complex_df = pd.DataFrame.from_dict(d, orient='index', columns=Coord_Col)
    complex_df.insert(loc=0, column=first_col, value=Circ_no)

    print(list(dict.keys()))
    add_col = ['Total_len','Read_No', 'Read_IDs']
    add_val = [tot_len,len(list(dict.keys())),[list(dict.keys())]]

    for i in range(len(add_col)):
        complex_df.insert(loc=len(complex_df.columns), column=add_col[i], value=add_val[i])

    bedtest = bt.BedTool.from_dataframe(complex_df)
    bedtest.saveas(savename)

Complex_circ_BED(complex_circ,Total_coord,Total_chr,1,"complex_full.bed")

def Read_Complex_BED(dict,savename):
    #finds the longest value and the corresponding key to create number of col in .BED
    max_key, max_value = max(dict.items(), key = lambda x: len(set(x[1])))
    df_col = ["Chr","Start","End"]
    rep = int(len(max_value)/len(df_col))
    first_col = 'Read_ID'
    Coord_Col = [j + "_no._"+str(i) for i in range(1, rep+1) for j in df_col]
    complex_df = pd.DataFrame.from_dict(dict, orient='index',columns = Coord_Col)

    complex_df.insert(loc=0, column=first_col, value=list(dict.keys()))
    All_col = [first_col]+Coord_Col

    complex_df = complex_df.sort_values(by=All_col)
    bedtest = bt.BedTool.from_dataframe(complex_df)
    bedtest.saveas(savename)



def Simple_circ_BED(dict,Circ_no,savename):
    df_col = ["Chr","Start","End"]
    simple_df = pd.DataFrame.from_dict(dict, orient='index',columns = df_col)
    simple_df = simple_df.sort_values(by=df_col)
    simple_df['Length'] = simple_df.apply(lambda x: x['End'] - x['Start'], axis=1)

    #create first column
    simple_df.insert(loc=0, column="Circle No.", value=Circ_no)

    #add these columns
    add_col = ['Read_No','Read_IDs']
    add_val = [len(list(dict.keys())[0]),[list(dict.keys())[0]]]

    for i in range(len(add_col)):
        simple_df.insert(loc=len(simple_df.columns), column=add_col[i], value=add_val[i])

    bedtest = bt.BedTool.from_dataframe(simple_df)
    bedtest.saveas(savename)

"""
bamfile = ps.AlignmentFile("BC05.aln_hg19.bam","rb")
Simple_circ = Single_coord(bamfile,60,"chr1", 243928620, 243938331)
Simple_circ_BED(Simple_circ,1,"Simple_test.bed")
"""




"""
coord=[82083498, 82085394, 24648701, 24650212, 24650308, 24651024, 82083498, 82086411, 24649837, 24651075, 82083498, 82086770, 24648862, 24651372, 24650143, 24651111, 82083710, 82084835, 82083760, 82085941, 82083769, 82086089, 82084010, 82087080, 46846937, 46847687, 53577491, 53578063, 82083498, 82086991, 24650457, 24651098, 82083787, 82085248, 24650826, 24651111, 82084355, 82085814, 82085262, 82087080, 53577491, 53578059, 46846937, 46847403, 82083795, 82086672, 82086008, 82086723, 82085836, 82087080, 46846960, 46847710, 53577493, 53578029, 53577491, 53577997, 46846937, 46847110]
chr=['chr2', 'chr2', 'chr5', 'chr5', 'chr5', 'chr5', 'chr2', 'chr2', 'chr5', 'chr5', 'chr2', 'chr2', 'chr5', 'chr5', 'chr5', 'chr5', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr16', 'chr16', 'chr2', 'chr2', 'chr5', 'chr5', 'chr2', 'chr2', 'chr5', 'chr5', 'chr2', 'chr2', 'chr2', 'chr2', 'chr16', 'chr16', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr16', 'chr16', 'chr16', 'chr16', 'chr2', 'chr2']

Coord_sort, chr_sort = chr_coord_sort(chr, coord)
Grouped_coord = Grouping(Coord_sort, 5000)

Group_size = [len(i) for i in Grouped_coord]
Group_no = len(Group_size)
Grouped_chr = []
j = [0]
d = {}
print(len(d))
for i in Group_size:
    First = sum(j[:])
    Grouped_chr.append(chr_sort[First:First+i])
    j.append(i)
print(Grouped_chr)
for i in range(len(Grouped_coord)):
    if len(d) == 0:
        d["test"] = [Grouped_chr[i][0], min(Grouped_coord[i]), max(Grouped_coord[i]),max(Grouped_coord[i])-min(Grouped_coord[i])]
    else:
        d["test"].extend([Grouped_chr[i][0], min(Grouped_coord[i]), max(Grouped_coord[i]),max(Grouped_coord[i])-min(Grouped_coord[i])])

df_col = ["Chr", "Start", "End","Length"]
rep = 4
first_col = 'Circle No.'
Coord_Col = [j + "_no._" + str(i) for i in range(1, rep + 1) for j in df_col]
print(Coord_Col)

complex_df = pd.DataFrame.from_dict(d, orient='index', columns=Coord_Col)

complex_df.insert(loc=0, column=first_col, value=1)
bedtest = bt.BedTool.from_dataframe(complex_df)
bedtest.saveas("test.bed")
"""

"""
### OTHER COMPLEX REGIONS FOR BC05###
# Complex(bamfile,60,"chr10", 66435392, 66437394)
#Complex(bamfile,60,"chr1", 16626700, 16627600)
#Complex(bamfile,60,"chr1", 243928620, 243938331)
"""


"""
print(Complex(bamfile,60,"chr1", 16626700, 16627600))
print(Complex(bamfile,60,"chr8", 51869000, 51869500))

print(Complex(bamfile,60,"chr7", 75713114, 75717415))
print("---------")
"""

"""
bamfile = ps.AlignmentFile("BC01.aln_hg19.bam","rb")
#print(Complex(bamfile,60,"chr7", 75713114, 75717415))
#print(Complex(bamfile,60,"chr2", 82083496, 82087081))
#print(Complex(bamfile,60,"chr2", 46846937, 46847710))


#def Circle_BED(dict,savename,Chr_list,Coord_list,Read_IDs):
 
chr10   18577508        18579630        13.64481132
chr10   37375321        37381579        41.09809264
chr10   66435392        66437394        20.97920892
chr10   126777139       126791124       10.4195464
chr10   133993072       133993674       16.65946844
chr11   18452021        18456198        16.27013423
chr11   84539733        84541979        45.27070347

chr1    14682287        14687166        21.45105328
chr1    16626766        16627585        25.63736264
chr1    106597836       106598346       17.30784314
chr1    114792782       114795930       7.672073791
chr1    161001223       161012454       111.6549728
chr1    174736312       174736815       80.06560636
chr1    243928620       243938331       236.1817527
"""