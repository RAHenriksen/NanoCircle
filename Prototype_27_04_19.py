
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
        #chr = [reg]
        #new_val = chr + max_value
        return dict([(max_key,max_value)])

    else:
        print("choice 2")
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
    coord_val = [[Coord_list[0]]]
    for i in Coord_list[1:]:
        # sets the previous and current elem as the first elem
        prev = coord_val[-1]
        current = prev[-1]
        dist = abs(i-current)

        # if dist between is below overlapping bp, then it it belong to the same group
        if dist < overlap_bp:
            prev.append(i)
        else:
            coord_val.append([i])
    return coord_val

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
            pos_start = int(Column_list[1]) - 1  #1-based
            pos_end = int(Column_list[1]) + length - 1
            overlap = sum(cigar_len)
            Total_overlap.append(overlap)

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
                #print("chomrs",chroms)
                Total_chr.extend((chrom, chrom))
                Total_coord.extend((pos_start, pos_end))

            if Coord_list != []:
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
                    Grouped_chr = Grouping_chr(chr_sort,Group_size)
                    for i in range(len(Grouped_coord)):
                        Complex_dict[read.query_name].extend(
                            [Grouped_chr[i][0], min(Grouped_coord[i]), max(Grouped_coord[i])])

    if Complex_dict == {}:
        print("The given input region is not forming a complex circle")
        if len(Sing_dict) >= 1:
            final_dict = {}
            chr = [reg]
            #HER INDSÆTTER JEG CHROMOSOME TALLET SÅ OVERVEJ LIGE OM MAN KUNNE GØRE DET FØR
            for k, v in Sing_dict.items():
                val = chr + v
                final_dict[k] = val
            Output_type = 1
            return final_dict, Output_type, None, None
        else:
            Output_type = 0
            return Sing_dict, Output_type, None, None
    else:
        Coord_sort, chr_sort = chr_coord_sort(Total_chr, Total_coord)
        Grouped_coord = Grouping(Coord_sort, max(Total_overlap))
        Group_size = [len(i) for i in Grouped_coord]
        Grouped_chr = Grouping_chr(chr_sort, Group_size)
        print("The given input region forms a complex circle")
        Output_type = 2
        return Complex_dict,Output_type,Grouped_coord,Grouped_chr

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
    #print(tot_len)
    df_col = ["Chr", "Start", "End", "Length"]
    rep = len(coord_full)
    first_col = 'Circle No.'
    Coord_Col = [j + "_no._" + str(i) for i in range(1, rep + 1) for j in df_col]

    complex_df = pd.DataFrame.from_dict(d, orient='index', columns=Coord_Col)
    complex_df.insert(loc=0, column=first_col, value=Circ_no)

    #print(list(dict.keys()))
    add_col = ['Total_len','Read_No', 'Read_IDs']
    add_val = [tot_len,len(list(dict.keys())),[list(dict.keys())]]

    for i in range(len(add_col)):
        complex_df.insert(loc=len(complex_df.columns), column=add_col[i], value=add_val[i])

    return complex_df
    #bedtest = bt.BedTool.from_dataframe(complex_df)
    #bedtest.saveas(savename)
    # return complex_df

def Read_bed(dict,savename,circ_no):
    #finds the longest value and the corresponding key to create number of col in .BED
    max_key, max_value = max(dict.items(), key = lambda x: len(set(x[1])))
    df_col = ["Chr","Start","End"]
    rep = int(len(max_value)/len(df_col))
    first_col = 'Read_ID'
    Coord_Col = [j + "_no._"+str(i) for i in range(1, rep+1) for j in df_col]
    complex_df = pd.DataFrame.from_dict(dict, orient='index',columns = Coord_Col)

    complex_df.insert(loc=0, column=first_col, value=list(dict.keys()))
    complex_df.insert(loc=0, column='Circ_no', value=circ_no)
    print(complex_df)
    return complex_df
    #bedtest = bt.BedTool.from_dataframe(complex_df)
    #bedtest.saveas(savename)

def Simple_circ_BED(beddict,Circ_no,savename):
    if len(beddict) > 1:
        coord = []
        read = []
        reg = str()
        for k, v in beddict.items():
            coord.extend(v[1:])
            read.append(k)
            reg = str(v[0])
        value = [str(reg), min(coord), max(coord)]
        beddict = dict([(tuple(read), value)])

    df_col = ["Chr","Start","End"]
    simple_df = pd.DataFrame.from_dict(beddict, orient='index',columns = df_col)
    simple_df = simple_df.sort_values(by=df_col)
    simple_df['Length'] = simple_df.apply(lambda x: x['End'] - x['Start'], axis=1)

    #create first column
    simple_df.insert(loc=0, column="Circle No.", value=Circ_no)

    #add these columns
    add_col = ['Read_No','Read_IDs']
    add_val = [len(list(beddict.keys())[0]),[list(beddict.keys())[0]]]

    for i in range(len(add_col)):
        simple_df.insert(loc=len(simple_df.columns), column=add_col[i], value=add_val[i])

    return simple_df
    #bedtest = bt.BedTool.from_dataframe(simple_df)
    #bedtest.saveas(savename)
    #return simple_df

#Simple_circ_BED(beddict,Circ_no,savename)

def Simple_read_BED(dict,savename,circ_no):
    max_key, max_value = max(dict.items(), key=lambda x: len(set(x[1])))
    df_col = ["Chr", "Start"]
    end_col = ["End"]
    rep = int(len(max_value) / 2)
    Coord_Col = [j + "_no._" + str(i) for i in range(1, rep + 1) for j in end_col]
    ful_col = df_col + Coord_Col
    simple_df = pd.DataFrame.from_dict(dict, orient='index', columns=ful_col)
    first_col = 'Read_ID'
    simple_df.insert(loc=0, column=first_col, value=list(dict.keys()))
    simple_df.insert(loc=0, column='Circ_no', value=circ_no)
    return simple_df

def BED_file_creation(file_name,mapQ):
    Simple_count = 1
    Simple_circ = pd.DataFrame()
    Simple_read = pd.DataFrame()
    Complex_count = 1
    Complex_circ = pd.DataFrame()
    Complex_read = pd.DataFrame()
    with open(file_name) as f:
        for line in f:
            line_value = line.strip().split()
            coord = line_value[0]
            start = line_value[1]
            end = line_value[2]
            cov = line_value[3]
            print("coord",coord,start,end)
            circle_dict, circ_type, circ_coord, circ_chr = Complex(bamfile,mapQ, str(coord), int(start), int(end))
            if circ_type == 1:
                circ_bed = Simple_circ_BED(circle_dict,Simple_count,"lol")
                rows = pd.concat([Simple_circ,circ_bed])
                Simple_circ = rows
                #circ_reads = Simple_read_BED(circle_dict, "lol", Simple_count)
                #reads = pd.concat([Simple_read,circ_reads],sort=True)
                #Simple_read = reads
                Simple_count += 1
            elif circ_type == 2:
                complex_bed = Complex_circ_BED(circle_dict,circ_coord,circ_chr,Complex_count,"lol")
                rows = pd.concat([Complex_circ,complex_bed],sort=False)
                Complex_circ = rows
                Complex_count +=1
            else:
                continue
    Simple_bed = bt.BedTool.from_dataframe(Simple_circ)
    Simple_bed.saveas("Simple_circles.bed")
    Complex_bed = bt.BedTool.from_dataframe(Complex_circ)
    Complex_bed.saveas("Complex_circles.bed")
    #Simple_read_bed = bt.BedTool.from_dataframe(Simple_read)
    #Simple_read_bed.saveas("Simple_circ_read.bed")


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
            print("Coord_dict",coord_dict)
            print("length",len(coord_dict))
            if bool(coord_dict) == True:
                print(coord+":"+start+"-"+end)
                #print("dictionary",coord_dict)
                count += 1
            else:
                #those dict empy, possible due to complex circle formation
                print(coord + ":" + start + "-" + end)
                print("empty")
                print("-------------")
                continue
            print("-------------")
    print("the number of circles is",count)

bamfile = ps.AlignmentFile("BC05.aln_hg19.bam","rb")
File_identification("BC05.ge_mean5.bdg",60)


"""
bamfile = ps.AlignmentFile("BC05.aln_hg19.bam","rb")

#BED_file_creation("BC05.ge_mean5.bdg",60)


#("BC05.ge_mean5.bdg",60)
#chr16 57495224 57495389
#Type 1: - simple one coord
# 18577508 18579630
#type 0
Simple_circ1,type1,coord1,chr1 = Complex(bamfile,60,"chr10", 18577508, 18579630)
print(Simple_circ1)
#Type 1: - simple several coord
Simple_circ2,type2,coord2,chr2 = Complex(bamfile,60,"chr5", 121967700, 121972207)
print(Simple_circ2)
#Type 2: - Complex coordinates


complex_circ,type3,Total_coord,Total_chr = Complex(bamfile,60,"chr1", 16626700, 16627600)
print(complex_circ)
print(coord1,coord2,Total_coord)
print(type1,type2,type3)
print(Complex_circ_BED(complex_circ,Total_coord,Total_chr,1,"lol"))


"""
"""
df1 = Simple_circ_BED(Simple_circ1,1,"test1.bed")
df2 = Simple_circ_BED(Simple_circ2,1,"test2.bed")

final = pd.concat([df1,df2])
bedtest = bt.BedTool.from_dataframe(final)
bedtest.saveas("final.bed")
# bedtest = bt.BedTool.from_dataframe(simple_df)
# return bedtest
# bedtest.saveas(savename)
"""

"""
Read_Complex_BED(sim_dict,"sim_read.bed",1)
Complex(bamfile,60,"chr1", 16626700, 16627600)
bamfile = ps.AlignmentFile("BC01.aln_hg19.bam","rb")
#complex_circ,Total_coord,Total_chr = Complex(bamfile,60,"chr5", 24649225, 24651095)
#Complex_circ_BED(complex_circ,Total_coord,Total_chr,3,"complex_full3.bed")
#Read_Complex_BED(complex_circ,"Complex_test3.bed",3)
"""
"""
def BED_file_creation(file_name,mapQ):
    count = 1
    Simple_BED = pd.DataFrame({'col' : []})
    Simple_read_BED = pd.DataFrame({'col' : []})
    Complex_BED = pd.DataFrame({'col' : []})
    Complex_read_BED = pd.DataFrame({'col' : []})
    with open(file_name) as f:
        for line in f:
            line_value = line.strip().split()
            coord = line_value[0]
            start = line_value[1]
            end = line_value[2]
            cov = line_value[3]
            coord_dict, circ_type, circ_coord, circ_chr = Complex(bamfile,mapQ, str(coord), int(start), int(end))
            print("type",circ_type)
            print("cord",coord,start,end)
            if bool(coord_dict) == True:
                print(coord_dict)
                count +=1
            else:
                continue
            print("-------------")
    print("the number of circles is",count)
    #bedtest = bt.BedTool.from_dataframe(Simple_BED)
    #bedtest.saveas("BEDCREATION_test.bed")"""