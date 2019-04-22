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

def Single_coord(bamfile,mapQ,reg,start,end,coord_val):
    """ returns merged coordinates for circles w. several reads having exact same coordinates. If not
    a dictionary with all reads and potential coordinates are returned"""
    Circle_dict = reduce_coords(bamfile, mapQ, reg, start, end)
    # intermediate dict, but keys should be values of orig dict
    Inter_dict = defaultdict(list)

    # WRITE SOME PROPER NOTE FOR THIS LINE
    #bascially CHECK this but if the length is larger than unique values
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

def chr_coord_sort(chrlist,coordlist):
    """ Sort a list of chromosomes and their coordinates using the index of the numerically sorted coordinates"""
    coord_idx = np.argsort(coordlist)
    Coord_sort = [coordlist[i] for i in coord_idx]
    chr_sort = [chrlist[i] for i in coord_idx]
    return Coord_sort,chr_sort


def Complex(bamfile,mapQ,reg, start, end):
    Sing_dict = Single_coord(bamfile, mapQ, reg, start, end,0)
    Complex_dict = {}
    Total_coord = []
    Total_chr = []
    if len(Sing_dict) == 1:
        reads_IDs = list(Sing_dict.keys())[0]
        Read_list = complex_reads(bamfile, mapQ, reg, start, end, reads_IDs)
    else:
        reads_IDs = list(Sing_dict.keys())
        Read_list = complex_reads(bamfile,mapQ,reg, start, end,reads_IDs)

    for read in Read_list[:-1]:
        if len(Sing_dict) == 1:
            coord1 = list(Sing_dict.values())[0][0]
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
            #print("tagelem",Tagelem)
            #splitting the individual aligment info up
            Column_list = Tagelem.split(',')
            chrom = Column_list[0]
            length = CIGAR_len(Column_list[3])
            cigar_len.append(length)
            pos_start = int(Column_list[1]) - 1  #1-based
            pos_end = int(Column_list[1]) + length - 1
            overlap = sum(cigar_len)

            # if the supp align is in between the circular input region then we already know the breakpoint from Sing_dict
            if chrom == reg and start - overlap <= pos_start <= end + overlap and start - overlap <= pos_end <= end + overlap:
                continue

            elif int(Column_list[4]) >= mapQ:
                #creates a coordinate list
                Coord_list.append(pos_start)
                Coord_list.append(pos_end)
                #append chr twice to ensure same length as coord_list
                chroms.append(chrom)
                chroms.append(chrom)
                Total_chr.extend((chrom, chrom))
                Total_coord.extend((pos_start, pos_end))

            #print("COORD",Coord_list)
            #print("CHR",chroms)
            #print("GROUPED",Grouping(Coord_list, overlap))
            if Coord_list != []:
                #sorts the chr and coordinates
                Coord_sort, chr_sort = chr_coord_sort(chroms, Coord_list)

                # first entry with the input region
                Complex_dict[read.query_name] = [reg, coord1, coord2]

                if len(Coord_sort) == 2:
                    #one coordinate pair supporting another region
                    Complex_dict[read.query_name].extend(
                        [chr_sort[0], min(Coord_sort), max(Coord_sort)])
                else:
                    Grouped_coord = Grouping(Coord_sort, overlap)
                    Group_size = [len(i) for i in Grouped_coord]
                    Grouped_chr = []
                    j = [0]
                    for i in Group_size:
                        First = j[-1]
                        Grouped_chr.append(chr_sort[First:First+i])
                        #print("teststst", chr_sort[First:First+i])
                        j.append(i)
                    for i in range(len(Grouped_coord)):
                        Complex_dict[read.query_name].extend(
                            [Grouped_chr[i][0], min(Grouped_coord[i]), max(Grouped_coord[i])])
    print("TOTALT LISTE",Total_coord)
    print("TOTAL CHR",Total_chr)
    if Complex_dict == {}:
        print("The given input region is not forming a complex circle")
        return Single_coord(bamfile, mapQ, reg, start,end,1)
    else:
        print("The given input region forms a complex circle")
        return Complex_dict

bamfile = ps.AlignmentFile("BC04.aln_hg19.bam","rb")
#File_identification("BC05.ge_mean5.bdg",60)


"""
print(Complex(bamfile,60,"chr1", 16626700, 16627600))
print(Complex(bamfile,60,"chr8", 51869000, 51869500))

print(Complex(bamfile,60,"chr7", 75713114, 75717415))
print("---------")
"""

"""
coord=[82083498, 82085394, 24648701, 24650212, 24650308, 24651024, 82083498, 82086411, 24649837, 24651075, 82083498, 82086770, 24648862, 24651372, 24650143, 24651111, 82083710, 82084835, 82083760, 82085941, 82083769, 82086089, 82084010, 82087080, 46846937, 46847687, 53577491, 53578063, 82083498, 82086991, 24650457, 24651098, 82083787, 82085248, 24650826, 24651111, 82084355, 82085814, 82085262, 82087080, 53577491, 53578059, 46846937, 46847403, 82083795, 82086672, 82086008, 82086723, 82085836, 82087080, 46846960, 46847710, 53577493, 53578029, 53577491, 53577997, 46846937, 46847110]
chr=['chr2', 'chr2', 'chr5', 'chr5', 'chr5', 'chr5', 'chr2', 'chr2', 'chr5', 'chr5', 'chr2', 'chr2', 'chr5', 'chr5', 'chr5', 'chr5', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr16', 'chr16', 'chr2', 'chr2', 'chr5', 'chr5', 'chr2', 'chr2', 'chr5', 'chr5', 'chr2', 'chr2', 'chr2', 'chr2', 'chr16', 'chr16', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr16', 'chr16', 'chr16', 'chr16', 'chr2', 'chr2']

Coord_sort, chr_sort = chr_coord_sort(chr, coord)
print(chr_sort)
Grouped_coord = Grouping(Coord_sort, 5000)
print(Grouped_coord)

Group_size = [len(i) for i in Grouped_coord]
print(Group_size)
print(len(chr_sort))
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
        d["test"] = [Grouped_chr[i][0], min(Grouped_coord[i]), max(Grouped_coord[i])]
    else:
        d["test"].extend([Grouped_chr[i][0], min(Grouped_coord[i]), max(Grouped_coord[i])])

print(d)
"""

bamfile = ps.AlignmentFile("BC01.aln_hg19.bam","rb")
#print(Complex(bamfile,60,"chr2", 82083496, 82087081))

"""
a = Complex(bamfile,60,"chr2", 82083496, 82087081)
print(a)

print("------------")
for k,v in a.items():
    print("reads",k)
    print("val",v)
print("------------")
test = sum(a.values(), [])
print(test)
chr = [item for item in test[::3] for i in range(2)]
print(chr)
del test[0::3]
coord = test
a,b = chr_coord_sort(chr,coord)
print(a)
print(b)

"""

def panda_test(dict):
    #finds the longest value and the corresponding key
    max_key, max_value = max(dict.items(), key = lambda x: len(set(x[1])))
    df_col = ["Chr","Start","End"]
    rep = int(len(max_value)/len(df_col))
    print(rep)
    All_Col = [j + "_no._"+str(i) for i in range(1, rep+1) for j in df_col]
    complex_df = pd.DataFrame.from_dict(dict, orient='index',columns = All_Col)
    print(complex_df.shape)
    complex_df = complex_df.sort_values(by=All_Col)
    print("min",complex_df.min(axis=0))
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(complex_df)

#a = Complex(bamfile,60,"chr2", 82083496, 82087081)
#panda_test(d)

"""
def File_identification(file_name,mapQ):
    count = 0
    with open(file_name) as f:
        for line in f:
            line_value = line.strip().split()
            coord = line_value[0]
            start = line_value[1]
            end = line_value[2]
            cov = line_value[3]
            coord_dict = Complex(bamfile,mapQ, str(coord), int(start), int(end))
            if bool(coord_dict) == True:
                print(coord+":"+start+"-"+end)
                print("dictionary",coord_dict)
                count += 1
            else:
                #those dict empy, possible due to complex circle formation
                print(coord + ":" + start + "-" + end)
                print("empy dictionary")
                continue
    print("the number of circles is",count)
    """

"""
def Genome_Cov_identification(bamfile,bam_name,overlap,reverse_str):
    Bedfile = pybedtools.example_bedtool(str(os.getcwd()) +"/"+bam_name)
    Cov = Bedfile.genome_coverage(bg=True)
    Merged = Cov.merge(d=overlap)

    for region in Merged:
        Chr = region[0]
        start = region[1]
        end = region[2]
        if len(Single_coord(bamfile, reverse_str, str(Chr), int(start), int(end))) == 0:
            continue
        else:
            print(Chr + ":" + start + "-" + end)
            print(Single_coord(bamfile, reverse_str, str(Chr), int(start), int(end)))

Single_coord(bamfile,reverse_str,mapQ,reg,start,end)

bamfile = ps.AlignmentFile("BC05.aln_hg19.bam","rb")

File_identification("BC05.ge_mean5.bdg",True,60)
print("----------------")
Genome_Cov_identification(bamfile,"BC05.aln_hg19.bam",500,True)
"""
def Count_reads(bamfile,reg,start,end):
    "Count the number of Soft-clipped reads with a supplementary alignment"
    s_count = 0
    SA_count = 0
    NO_SA_count = 0
    supp_count = 0
    for read in bamfile.fetch(reg,start,end,multiple_iterators=True):
        if read.is_supplementary == True:
            supp_count += 1
        #OVERVEJ LIGE OM MAN IKKE FAKTISK MISTER MANGE VED IKKE AT BRUGE cigar[-1][0]
        if read.cigar[0][0] == 4:
            s_count += 1
            if read.has_tag("SA"):
                SA_count += 1
            if not read.has_tag("SA"):
                NO_SA_count +=1
    print("total number of reads in the file",bamfile.count(reg,start,end,read_callback='nofilter'))
    print("total number of soft-clipped",s_count)
    print("Number of soft-clipped reads with SA",SA_count)
    print("Number of soft-clipped reads wo. SA",NO_SA_count)
    print("supp",supp_count)

Count_reads(bamfile,None,None,None)

def Read_length(bamfile,reg,start,end,filename,read_type):
    read_list = []
    for read in bamfile.fetch(reg, start, end, multiple_iterators=True):
        if read_type == "Soft":
            if read.cigar[0][0] == 4:
                read_list.append(read.query_length)
        elif read_type == "SA":
            if read.cigar[0][0] == 4 and read.has_tag("SA"):
                Tag = read.get_tag("SA").split(';')[:-1]
                for Tagelem in Tag:
                    Column_list = Tagelem.split(',')
                    read_list.append(CIGAR_len(Column_list[3]))
        elif read_type == "supp":
            if read.is_supplementary == True:
                read_list.append(read.query_length)
    with open(filename, 'w') as f:
        for item in read_list:
            f.write("%s\n" % item)

Read_length(bamfile,None,None,None,"Soft.txt","Soft")
Read_length(bamfile,None,None,None,"SA.txt","SA")
Read_length(bamfile,None,None,None,"Supp.txt","supp")
