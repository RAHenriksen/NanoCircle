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
        #print("REDUCTERET")
        #new_val =chr + max_value
        new_val =  chr + max_value
        return dict([(max_key,new_val)])

    else:
        if len(Circle_dict.keys()) == 1:
            Inter_dict = defaultdict(list)
            for k, v in sorted(Circle_dict.items()):
                Inter_dict[tuple(v)].append(k)
                #print("k",k)
            final_dict = {tuple(v): list(k) for k, v in Inter_dict.items()}
            Key, Value = max(final_dict.items(), key=lambda x: len(set(x[1])))
            chr = [reg]
            new_val = chr + Value
            return dict([(Key,new_val)])
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

def chr_coord_sort(chrlist,coordlist):
    """ Sort a list of chromosomes and their coordinates using the index of the numerically sorted coordinates"""
    coord_idx = np.argsort(coordlist)
    Coord_sort = [coordlist[i] for i in coord_idx]
    chr_sort = [chrlist[i] for i in coord_idx]
    return Coord_sort,chr_sort

def Grouping_chr(Chr_sort,Group_size):
    """ Groups the chromosomes in to match the grouped coordinates """
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
    if len(Sing_dict) == 1:
        reads_IDs = list(Sing_dict.keys())[0]
        Read_list = complex_reads(bamfile, mapQ, reg, start, end, reads_IDs)
    else:
        reads_IDs = list(Sing_dict.keys())
        Read_list = complex_reads(bamfile,mapQ,reg, start, end,reads_IDs)

    for read in Read_list:
        #print("Name",read.query_name)
        if len(Sing_dict) == 1:
            coord1 = list(Sing_dict.values())[0][-2]
            coord2 = list(Sing_dict.values())[0][-1]
            Total_chr.extend((reg,reg))
            Total_coord.extend((coord1,coord2))
        else:
            coord1 = Sing_dict[read.query_name][-2]
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
            #splitting the individual aligment info up
            Column_list = Tagelem.split(',')
            chrom = Column_list[0]
            length = CIGAR_len(Column_list[3])
            #print("lenght",length)
            cigar_len.append(length)
            pos_start = int(Column_list[1]) - 1  #1-based
            pos_end = int(Column_list[1]) + length - 1
            #print("pos",pos_start,pos_end)
            overlap = sum(cigar_len)*2
            Total_overlap.append(overlap)

            # if the supp align is in between the circular input region then we already know the breakpoint from Sing_dict
            if chrom == reg and start - overlap <= pos_start <= end + overlap and start - overlap <= pos_end <= end + overlap:
                continue

            elif int(Column_list[4]) >= mapQ:
                #creates a coordinate list
                print(read.query_name)
                print(Tagelem)
                #print("POS",pos_start,pos_end)
                Coord_list.append(pos_start)
                Coord_list.append(pos_end)
                #append chr twice to ensure same length as coord_list
                chroms.append(chrom)
                chroms.append(chrom)
                #print("chomrs",chroms)
                Total_chr.extend((chrom, chrom))
                #print("NEW POS",Coord_list)
                Total_coord.extend((pos_start, pos_end))
            if Coord_list != []:
                #sorts the chr and coordinates
                Coord_sort, chr_sort = chr_coord_sort(chroms, Coord_list)

                # first entry with the input region
                Complex_dict[read.query_name] = [reg, coord1, coord2]

                if len(Coord_sort) == 2:
                    #one coordinate pair supporting another region
                    Complex_dict[read.query_name].extend(
                        [chr_sort[0], int(min(Coord_sort)), int(max(Coord_sort))])
                else:
                    Grouped_coord = Grouping(Coord_sort, max(Total_overlap)*2)
                    Group_size = [len(i) for i in Grouped_coord]
                    Grouped_chr = Grouping_chr(chr_sort,Group_size)
                    for i in range(len(Grouped_coord)):
                        Complex_dict[read.query_name].extend(
                            [Grouped_chr[i][0], int(min(Grouped_coord[i])), int(max(Grouped_coord[i]))])

    if Complex_dict == {}:
        if len(Sing_dict) >= 1:
            print("The given input region is not forming a complex circle")
            final_dict = {}
            chr = [reg]
            #1 set of coord
            #print(list(Sing_dict.keys()))
            if len(list(Sing_dict.keys())) == 1:
                Output_type = 1
                #print("SING_DICT OUT")
                return Sing_dict, Output_type, None, None
            #several sets of coord
            elif len(list(Sing_dict.keys())) >= 2:
                #simple circ with several coordinate sets
                Coord_val = [v for k,v in Sing_dict.items()]
                #print("values",Coord_val)
                keys = [k for k,v in Sing_dict.items()]
                interval = [val for sublist in Coord_val for val in sublist]
                val = chr + [min(interval),max(interval)]
                Output_type = 2
                final_dict[tuple(keys)] = val
                #print("FINAL DICT")
                return final_dict, Output_type, Coord_val, None
        else:
            #empty dict
            print("The given input region is not forming a circle")
            Output_type = 0
            return Sing_dict, Output_type, None, None
    else:
        #Sorting again to create the coordinate and chromosome list
        #print("Total_coord",Total_coord)
        Coord_sort, chr_sort = chr_coord_sort(Total_chr, Total_coord)
        #print("Coord_sort",Coord_sort)
        Grouped_coord = Grouping(Coord_sort, max(Total_overlap)*2)
        Group_size = [len(i) for i in Grouped_coord]
        Grouped_chr = Grouping_chr(chr_sort, Group_size)
        print("The given input region forms a complex circle")
        Output_type = 3
        return Complex_dict,Output_type,Grouped_coord,Grouped_chr

####
# OKAY PROBLEMET LIGGER I AT NOGLE AF DE KOMPLEKSE CIRKLER ER IKKE SINGLE_COORD, SÅ DER BLIVER RETURNERET CIRCLE_DICT
# OG DERFOR IKKE TILFØJET KOORDINAT MED DET SAMME, DERFOR VED SLICING OPPE I KOMPLEKS SÅ BLIVER DE FORKERTE VÆRDIER TAGET UD
####

def Complex_full_length(file_name,bamname,mapQ):
    bamfile = ps.AlignmentFile(bamname, "rb")
    total_col_len = 0
    with open(file_name) as f:
        for line in f:
            line_value = line.strip().split()
            coord = line_value[0]
            start = line_value[1]
            end = line_value[2]
            complex_dict, circ_type, circ_coord, circ_chr = Complex(bamfile,mapQ, str(coord), int(start), int(end))
            if circ_type == 3:
                if len(circ_coord) > total_col_len:
                    total_col_len = len(circ_coord)
    return total_col_len

def Complex_circ_BED(dict,coord_full,chr_full,length,Circ_no,circ_type,savename,Filename):
    d = {}
    tot_len = 0
    SampleID = Filename.split('.')[0]
    for i in range(len(coord_full)):
        if len(d) == 0:
            d["test"] = [chr_full[i][0], int(min(coord_full[i])), int(max(coord_full[i])),
                         int(max(coord_full[i])) - int(min(coord_full[i]))]
            tot_len += (max(coord_full[i]) - min(coord_full[i]))
        else:
            d["test"].extend([chr_full[i][0], int(min(coord_full[i])), int(max(coord_full[i])),
                         int(max(coord_full[i])) - int(min(coord_full[i]))])
            tot_len += (int(max(coord_full[i])) - int(min(coord_full[i])))
    df_col = ["Chr", "Start", "End", "Length"]
    rep = len(coord_full)
    tot = length
    Coord_Col = [j + "_no._" + str(i) for i in range(1, tot + 1) for j in df_col]
    if tot > rep:
        Nan_list = ["nan"]
        extend_list = Nan_list*(4*int(tot-rep))
        d["test"].extend(extend_list)
    print("d2",d)
    complex_df = pd.DataFrame.from_dict(d, orient='index', columns=Coord_Col)

    add_col = ['Total_len','Read_No', 'Read_IDs', 'Circle_type', 'Circle_ID']

    Read_ID_str = str(list(list(dict.keys()))).replace(" ", "")
    add_val = [int(tot_len),len(list(dict.keys())), Read_ID_str, "Circle_type_%s" % circ_type,
           "%s_Complex_circ_%d" % (SampleID, Circ_no)]

    for i in range(len(add_col)):
        complex_df.insert(loc=len(complex_df.columns), column=add_col[i], value=add_val[i])


    bedtest = bt.BedTool.from_dataframe(complex_df)
    bedtest.saveas(savename)
    #return complex_df

bamfile = ps.AlignmentFile("BC01.aln_hg19.bam","rb")
#"chr2", 82083496, 82087081
#"chr2", 46846900, 46849100

Complex_dict,Output_type,Grouped_coord,Grouped_chr = Complex(bamfile,40,"chr5", 24649000,24651100)
print(len(Complex_dict))
#DICT LÆNGDE 9 passer ikke med bed dad et kun er 8
print("-----------------------------")
#print(Complex(bamfile,60,"chr2", 46846900, 46849100))
#LEN = 2 passer med bed
print("---------------------------")

Complex_dict,Output_type,Grouped_coord,Grouped_chr = Complex(bamfile,40,"chr5", 24649000,24651100)
print(Complex_dict)
print(len(Complex_dict))
print(Complex_dict.keys())
print(Complex_circ_BED(Complex_dict,Grouped_coord,Grouped_chr,len(Complex_dict),1,1,"lol.bed","lol.bed"))
#længde er 5, passer heller ikke med bed som er 4