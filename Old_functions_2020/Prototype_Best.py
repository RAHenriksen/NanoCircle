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

                #prim_pos[0] = prim_pos[0]+1
                #supp_pos[0] = supp_pos[0]+1
                #prim_pos[-1] = prim_pos[-1] + 1
                #supp_pos[-1] = supp_pos[-1] + 1

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
    """ Check for potential complex circles by reads aligning across the genome. Afterwards it returns the simple circles
     with 1 set of coordinates (type I) several set of coordinates (type II) and the complex circles (type III)"""
    Sing_dict = Single_coord(bamfile, mapQ, reg, start, end)

    Complex_dict = {}
    print("Complex_dictbeg",Complex_dict)
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
            coord1 = list(Sing_dict.values())[0][-2]
            coord2 = list(Sing_dict.values())[0][-1]
            Total_chr.extend((reg,reg))
            Total_coord.extend((coord1,coord2))
        else:
            coord1 = Sing_dict[read.query_name][-2]
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
            # 1 - based positions
            pos_start = int(Column_list[1]) - 1
            pos_end = int(Column_list[1]) + length - 1

            # the overlaps between coordinates for grouping
            overlap = sum(cigar_len)*4
            print("overlap",overlap)
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

                Total_chr.extend((chrom, chrom))

                Total_coord.extend((pos_start, pos_end))

            if Coord_list != []:

                #sorts the chr and coordinates
                Coord_sort, chr_sort = chr_coord_sort(chroms, Coord_list)
                print("coord",Coord_sort,chr_sort)
                #first entry with the input region
                Complex_dict[read.query_name] = [reg, coord1, coord2]

                if len(Coord_sort) == 2:
                    print("TRUE",Coord_sort)
                    #one coordinate pair supporting another region
                    Complex_dict[read.query_name].extend(
                        [chr_sort[0], int(min(Coord_sort)), int(max(Coord_sort))])
                else:
                    # max(Total_overlap) * 2 er fjernet da jeg ændrede overlap = sum(cigar_len)*2 til *4
                    Grouped_coord = Grouping(Coord_sort, max(Total_overlap))
                    Group_size = [len(i) for i in Grouped_coord]
                    Grouped_chr = Grouping_chr(chr_sort,Group_size)
                    for i in range(len(Grouped_coord)):
                        Complex_dict[read.query_name].extend(
                            [Grouped_chr[i][0], int(min(Grouped_coord[i])), int(max(Grouped_coord[i]))])

    print("COMPLEX",Complex_dict)

    if Complex_dict == {}:
        if len(Sing_dict) >= 1:
            final_dict = {}
            chr = [reg]

            # circle type 1, simple circle with 1 set of coordinate supported by 1 or several reads
            if len(list(Sing_dict.keys())) == 1:
                Output_type = 1
                return Sing_dict, Output_type, None, None

            # circle type 2, simple circle with several sets of coordinates from several reads
            elif len(list(Sing_dict.keys())) >= 2:
                Coord_val = [v for k,v in Sing_dict.items()]
                keys = [k for k,v in Sing_dict.items()]
                interval = [val for sublist in Coord_val for val in sublist]
                val = chr + [min(interval),max(interval)]
                Output_type = 2
                final_dict[tuple(keys)] = val
                return final_dict, Output_type, Coord_val, None
        else:
            #empty dict
            Output_type = 0
            return Sing_dict, Output_type, None, None
    else:
        #Sorting again to create the coordinate and chromosome list
        Coord_sort, chr_sort = chr_coord_sort(Total_chr, Total_coord)
        Grouped_coord = Grouping(Coord_sort, max(Total_overlap))
        Group_size = [len(i) for i in Grouped_coord]
        Grouped_chr = Grouping_chr(chr_sort, Group_size)

        # circle type 3, complex circles comprised of several chromosomal regions across the genome.
        Output_type = 3
        return Complex_dict,Output_type,Grouped_coord,Grouped_chr
    print("output", Output_type)


def Simple_circ_BED(beddict,Circ_no,circ_type,Filename):
    """ returns a dataframe with circular information for the simple circles, type 1 and 2 """
    df_col = ["Chr","Start","End"]
    simple_df = pd.DataFrame.from_dict(beddict, orient='index',columns = df_col)
    simple_df = simple_df.sort_values(by=df_col)
    simple_df['Length'] = simple_df.apply(lambda x: x['End'] - x['Start'], axis=1)
    SampleID = Filename.split('.')[0]
    add_col = ['Read_No','Read_IDs','Circle_type','Circle_ID']

    #convert list of ID's to 1 long string as to insert it as a single column in the df
    Read_ID_str = str(list(list(beddict.keys())[0])).replace(" ","")
    add_val = [len(list(beddict.keys())[0]),Read_ID_str,"Circle_type_%s" % circ_type,"%s_simple_circ_%d" % (SampleID,Circ_no)]

    for i in range(len(add_col)):
        simple_df.insert(loc=len(simple_df.columns), column=add_col[i], value=add_val[i])
    return simple_df

def Simple_reads(dict,Circ_no,coord,circ_type,Filename):
    """ returns a dataframe with specific read information for the simple circles type 2 """

    df_col = ["Chr","Start","End",'Read_ID','Circle_type','Circle_ID']
    SampleID = Filename.split('.')[0]
    Ref_df = pd.DataFrame(columns=df_col)
    keys = list(list(dict.keys())[0])
    chr = list(dict.values())[0][0]
    for i in range(len(keys)):
        Ref_df.loc[len(Ref_df)] = [str(chr), coord[i][0], coord[i][1],keys[i],"Circle_type_%s" % circ_type,"%s_simple_circ_%d" % (SampleID,Circ_no)]
    return Ref_df

def Complex_full_length(file_name,bamname,mapQ):
    """ returns the total number of columns needed for the most complex circle"""
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

def Complex_circ_BED(dict,coord_full,chr_full,length,Circ_no,circ_type,bamname):
    """ returns a dataframe with circular information for the complex circles, type 3 """

    #temporary dictionary
    d = {}
    tot_len = 0
    SampleID = bamname.split('.')[0]
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
    #number of times df_col needed to be repeated
    rep = len(coord_full)
    tot = length
    Coord_Col = [j + "_no._" + str(i) for i in range(1, tot + 1) for j in df_col]
    #insert nan for those complex circles comprised of fewer regions that the most complex
    if tot > rep:
        Nan_list = ["nan"]
        extend_list = Nan_list*(4*int(tot-rep))
        d["test"].extend(extend_list)
    complex_df = pd.DataFrame.from_dict(d, orient='index', columns=Coord_Col)

    add_col = ['Total_len','Read_No', 'Read_IDs', 'Circle_type', 'Circle_ID']

    Read_ID_str = str(list(list(dict.keys()))).replace(" ", "")
    add_val = [int(tot_len),len(list(dict.keys())), Read_ID_str, "Circle_type_%s" % circ_type,
           "%s_Complex_circ_%d" % (SampleID, Circ_no)]

    for i in range(len(add_col)):
        complex_df.insert(loc=len(complex_df.columns), column=add_col[i], value=add_val[i])

    return complex_df

def Complex_reads(dict,length,circ_no,circ_type,bamname):
    """ returns the total number of columns needed for the most complex circle"""

    #finds the longest value and the corresponding key to create number of col in .BED
    SampleID = bamname.split('.')[0]
    max_key, max_value = max(dict.items(), key = lambda x: len(set(x[1])))
    df_col = ["Chr","Start","End"]
    rep = int(len(max_value)/len(df_col))
    tot = length
    new_dict = {}
    for k,v in dict.items():
        if tot > rep:
            Nan_list = ["nan"]
            extend_list = Nan_list * (3 * int(tot - rep - 1))
            new_dict[k] = v + extend_list

    Coord_Col = [j + "_no._"+str(i) for i in range(1, tot) for j in df_col]
    complex_df = pd.DataFrame.from_dict(new_dict, orient='index',columns = Coord_Col)

    add_col = ['Read_IDs', 'Circle_type', 'Circle_ID']
    add_val = [list(new_dict.keys()),"Circle_type_%s" % circ_type,"%s_Complex_circ_%d" % (SampleID, circ_no)]
    for i in range(len(add_col)):
        complex_df.insert(loc=len(complex_df.columns), column=add_col[i], value=add_val[i])
    return complex_df

def FULL(file_name,bamname,mapQ):
    bamfile = ps.AlignmentFile(bamname, "rb")
    Simple_count = 1
    Simple_circ = pd.DataFrame()
    Simple_read = pd.DataFrame()

    Complex_count = 1
    Complex_df = pd.DataFrame()
    read_df_full = pd.DataFrame()
    Complex_col_no = Complex_full_length(file_name, bamname, mapQ)

    Sample_ID = bamname.split('.')[0]
    with open(file_name) as f:
        for line in f:
            line_value = line.strip().split()
            coord = line_value[0]
            start = line_value[1]
            end = line_value[2]
            circle_dict, circ_type, circ_coord, circ_chr = Complex(bamfile,mapQ, str(coord), int(start), int(end))

            if circ_type == 1:
                print("region", coord, start, end, "circ_type", circ_type)
                circ_bed = Simple_circ_BED(circle_dict,Simple_count,circ_type,bamname)
                rows = pd.concat([Simple_circ,circ_bed])
                Simple_circ = rows
                Simple_count += 1

            #if the circle has more than 1 potential coordinate set
            elif circ_type == 2:
                print("region", coord, start, end, "circ_type", circ_type)
                circ_bed = Simple_circ_BED(circle_dict, Simple_count,circ_type,bamname)
                circ_read = Simple_reads(circle_dict,Simple_count,circ_coord,circ_type,bamname)
                row_red = pd.concat([Simple_read,circ_read])
                rows = pd.concat([Simple_circ, circ_bed])
                Simple_circ = rows
                Simple_read = row_red
                Simple_count += 1

            elif circ_type == 3:
                print("region", coord, start, end, "circ_type", circ_type)
                df = Complex_circ_BED(circle_dict, circ_coord, circ_chr, Complex_col_no, Complex_count, circ_type, file_name)
                Complex_df = Complex_df.append(df, sort=False).fillna("nan")
                read_df = Complex_reads(circle_dict,Complex_col_no,Complex_count,circ_type,file_name)
                read_df_full = read_df_full.append(read_df,sort=False).fillna("nan")
                Complex_count += 1
            else:
                continue
    Simple_bed = bt.BedTool.from_dataframe(Simple_circ)
    Simple_bed.saveas("Simple_circles_%s_v2.bed" % Sample_ID)
    Simple_read_bed = bt.BedTool.from_dataframe(Simple_read)
    Simple_read_bed.saveas("Simple_reads_%s_v2.bed" % Sample_ID)

    bedtest = bt.BedTool.from_dataframe(Complex_df)
    bedtest.saveas("Complex_circles_%s_v2.bed" % Sample_ID)
    bedtest = bt.BedTool.from_dataframe(read_df_full)
    bedtest.saveas("Complex_reads_%s_v2.bed" % Sample_ID)

def BED_Coverage(bamname,overlap,mapQ):
    bamfile = ps.AlignmentFile(bamname, "rb")
    Sample_ID = bamname.split('.')[0]

    Bedfile = pybedtools.example_bedtool(str(os.getcwd()) +"/"+bamname)
    Cov = Bedfile.genome_coverage(bg=True)
    Merged = Cov.merge(d=overlap)

    Simple_count = 1
    Simple_circ = pd.DataFrame()
    Simple_read = pd.DataFrame()
    Complex_count = 1
    Complex_df = pd.DataFrame()
    read_df_full = pd.DataFrame()

    total_col_len = 0

    for region in Merged:
        coord = region[0]
        start = region[1]
        end = region[2]
        circle_dict, circ_type, circ_coord, circ_chr = Complex(bamfile,mapQ, str(coord), int(start), int(end))
        print("Circle type",circ_type)
        if circ_type == 1:
            circ_bed = Simple_circ_BED(circle_dict,Simple_count,circ_type,bamname)
            rows = pd.concat([Simple_circ,circ_bed])
            Simple_circ = rows
            Simple_count += 1
        #if the circle has more than 1 potential coordinate set
        elif circ_type == 2:
            circ_bed = Simple_circ_BED(circle_dict, Simple_count,circ_type,bamname)
            circ_read = Simple_reads(circle_dict,Simple_count,circ_coord,circ_type,bamname)
            row_red = pd.concat([Simple_read,circ_read])
            rows = pd.concat([Simple_circ, circ_bed])
            Simple_circ = rows
            Simple_read = row_red
            Simple_count += 1
        elif circ_type == 3:
            #det skal lige ændres med 6.
            df = Complex_circ_BED(circle_dict, circ_coord, circ_chr, 6, Complex_count, circ_type, bamname)
            Complex_df = Complex_df.append(df, sort=False).fillna("nan")
            read_df = Complex_reads(circle_dict,6,Complex_count,circ_type,bamname)
            read_df_full = read_df_full.append(read_df,sort=False).fillna("nan")
            Complex_count += 1
        else:
            continue
    Simple_bed = bt.BedTool.from_dataframe(Simple_circ)
    Simple_bed.saveas("Simple_circles_%s.bed" % Sample_ID)
    Simple_read_bed = bt.BedTool.from_dataframe(Simple_read)
    Simple_read_bed.saveas("Simple_reads_%s.bed" % Sample_ID)
    bedtest = bt.BedTool.from_dataframe(Complex_df)
    bedtest.saveas("Complex_circles_%s.bed" % Sample_ID)
    bedtest = bt.BedTool.from_dataframe(read_df_full)
    bedtest.saveas("Complex_reads_%s.bed" % Sample_ID)

bamfile = ps.AlignmentFile("BC05.aln_hg19.bam","rb")

#Complex_dict,Output_type,Grouped_coord,Grouped_chr = Complex(bamfile,60,"chr5",24908937,24914000)
#print(Complex_dict)
#print(Output_type)
#print(Grouped_coord)

#FULL("BC05.ge_mean5.bdg","BC05.aln_hg19.bam",60)
#BED_Coverage("BC03.aln_hg19.bam",1000,60)

#Complex_dict,Output_type,Grouped_coord,Grouped_chr = Complex(bamfile,60,"chr1",184459828,184462249)
#Complex_dict,Output_type,Grouped_coord,Grouped_chr = Complex(bamfile,60,"chr1",184455265,184458079)
#region chr1 184455265 184458079 circ_type 3
#region chr1 184459828 184462249 circ_type 1

#Supplement_dict(bamfile,60,"chr1",243928620,243938331)

print(len(list(Single_coord(bamfile,60,"chr1",243928620,243938331).keys())[0]))
print(Single_coord(bamfile,60,"chrM",0,16572))
print(Single_coord(bamfile,60,"chr1",243928620,243938331))