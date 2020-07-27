import argparse
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
import operator

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

def df_float_to_int(df):
    "Ensures the coordinate values in the dataframe are integers"
    row_no = df.shape[0]
    for i in df.columns:
        for j in range(row_no):
            if isinstance(df[str(i)][j], str):
                pass
            elif isinstance(df[str(i)][j], float):
                df.at[j, str(i)] = int(df[str(i)][j])
    return df

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

def reduce_end_coords(bamfile,mapQ,reg,start,end):
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

def most_frequent(List,pos):
    count1 = Counter(List)
    dic_val = list(count1.values())
    max_val = max(count1.values())
    occurence = dic_val.count(max_val)

    #single most occuring coordinate
    if occurence == 1:
        #print("COMMMON VALUE", count1.most_common(1)[0][0])
        occ = 1
        return occ, count1.most_common(1)[0][0]

    #creates the interval immediately
    else:
        most_common_no = list(map(operator.itemgetter(1), count1.most_common()))
        most_common_idx = [i for i, x in enumerate(most_common_no) if x == max(most_common_no)]


        #if there are several equal common values
        if max(most_common_no) > 1:
            #print(max(most_common_no))
            #print(most_common_idx)
            #creating the most maximum value
            common_val = [count1.most_common()[i][0] for i in most_common_idx]
            occ = 1
            #print("COMMMON VALUE",common_val)
            return occ, max(common_val)

        #if there isnt a most occuring coordinate, the interval is created
        else:
            occ = 2
            #start
            if pos == 0:
                return occ, min(List)
            #end
            if pos == 1:
                return occ, max(List)

def Single_coord(bamfile,mapQ,reg,start,end):
    """ returns the most frequent start and end coordinate for a circle in a dictionary with the key
    being all the reads. If not having a coordinate more frequent than others, a dictionary with all reads
    and the interval for potential circle """
    Circle_dict = reduce_end_coords(bamfile, mapQ, reg, start, end)
    start = [i[0] for i in Circle_dict.values()]
    #taking all possible end coordinates and merge into one list
    end = sum([i[1:] for i in Circle_dict.values()],[])
    #extract the most common
    if start or end != []:
        occ1,start_freq = most_frequent(start,0)
        #print("occurence,",occ1,"start",start_freq)
        occ2,end_freq = most_frequent(end,1)
        #print("occurence,", occ2, "start", end_freq)
        if occ1 and occ2 == 1:
            reads = []
            chr = [reg]
            #print("START COORDS", start)
            #print("END COORD",end)
            if start_freq < end_freq:
                new_val = chr + [start_freq,end_freq]
            else:
                #print("lol")
                #print(end_freq,start_freq)
                new_val = chr + [end_freq,start_freq]
            #print("Start and end",new_val)
            for k,v in Circle_dict.items():
                if any(i == start_freq or i == end_freq for i in v):
                    reads.append(k)
            final_dict = {tuple(sorted(reads)): new_val}

            #Multiple reads
            if len(list(final_dict.keys())[0]) != 1:
                type = 1
                return type,final_dict

            #Single read
            else:
                type = 2
                return type, final_dict

        #not a more supported read
        else:
            type = 2
            #mangler stadig at lave det der overlaps
            chr = [reg]
            new_val = chr + [start_freq, end_freq]
            final_dict = {tuple(sorted(Circle_dict.keys())): new_val}
            #man kunne eventuelt her bare returnere Single coord
            return type,final_dict
    else:
        type = 0
        #these dicts are empty, and serves to continue as some regions might not create a circle
        return type,Circle_dict

# saa der skal aendres lidt mere saaledes at outputtet for type 2 er det samme. dette er stadig et problem.

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
    Type,Sing_dict = Single_coord(bamfile,mapQ,reg,start,end)

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
        coord1 = list(Sing_dict.values())[0][-2]
        coord2 = list(Sing_dict.values())[0][-1]
        #print("READ QUERY",read.query_name)
        #print("SING DICT COMPLEX", Sing_dict.values())
        #print("CORD",coord1,coord2)
        Total_chr.extend((reg,reg))
        Total_coord.extend((coord1,coord2))
        #print("TOAL CHR",Total_chr,Total_coord)
        Tag = read.get_tag("SA").split(';')[:-1]
        Coord_list = []
        chroms = []
        cigar_len = []

        #examining for complexity
        for Tagelem in Tag:
            #splitting the individual aligment info up
            Column_list = Tagelem.split(',')
            #print("Column",Column_list)
            chrom = Column_list[0]
            length = CIGAR_len(Column_list[3])
            cigar_len.append(length)
            # 1 - based positions
            pos_start = int(Column_list[1]) - 1
            pos_end = int(Column_list[1]) + length - 1
            # the overlaps between coordinates for grouping
            overlap = sum(cigar_len)*4
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
                #first entry with the input region
                Complex_dict[read.query_name] = [reg, coord1, coord2]

                if len(Coord_sort) == 2:
                    #one coordinate pair supporting another region
                    Complex_dict[read.query_name].extend(
                        [chr_sort[0], int(min(Coord_sort)), int(max(Coord_sort))])
                else:
                    Grouped_coord = Grouping(Coord_sort, max(Total_overlap))
                    Group_size = [len(i) for i in Grouped_coord]
                    Grouped_chr = Grouping_chr(chr_sort,Group_size)
                    for i in range(len(Grouped_coord)):
                        Complex_dict[read.query_name].extend(
                            [Grouped_chr[i][0], int(min(Grouped_coord[i])), int(max(Grouped_coord[i]))])

    #the simple circles
    if Complex_dict == {}:
        if Type == 0 or 1 or 2:
            return Sing_dict, Type, None, None
    else:
        #Sorting again to create the coordinate and chromosome list
        Coord_sort, chr_sort = chr_coord_sort(Total_chr, Total_coord)
        #print("SORT",Coord_sort,chr_sort)
        Grouped_coord = Grouping(Coord_sort, max(Total_overlap))
        #print("GROUPED",Grouped_coord)
        Group_size = [len(i) for i in Grouped_coord]
        Grouped_chr = Grouping_chr(chr_sort, Group_size)
        # circle type 3, complex circles comprised of several chromosomal regions across the genome.
        Type = 3
        #print("SINFG", Sing_dict)
        print("COMPLEX", Complex_dict,Type,Grouped_coord,Grouped_chr)
        return Complex_dict,Type,Grouped_coord,Grouped_chr

def Simple_circ_df(beddict,Circ_no,circ_type,Filename):
    """ returns a dataframe with circular information for the simple circles, type 1 and 2 """
    df_col = ["Chr","Start","End"]
    simple_df = pd.DataFrame.from_dict(beddict, orient='index',columns = df_col)
    simple_df = simple_df.sort_values(by=df_col)
    simple_df['Length'] = simple_df.apply(lambda x: x['End'] - x['Start'], axis=1)
    SampleID = Filename.split('.')[0]
    #add_col = ['Read_No','Read_IDs','Circle_type','Circle_ID']
    add_col = ['Read_No', 'Circle_type', 'Circle_ID']
    #convert list of ID's to 1 long string as to insert it as a single column in the df

    #Read_ID_str = str(list(list(beddict.keys())[0])).replace(" ","")
    #add_val = [len(list(beddict.keys())[0]),Read_ID_str,circ_type,"%s_simple_circ_%d" % (SampleID,Circ_no)]

    add_val = [len(list(beddict.keys())[0]),circ_type,"%s_simple_circ_%d" % (SampleID,Circ_no)]

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

    if len(keys) == 1:
        Ref_df.loc[len(Ref_df)] = [str(chr), coord[0], coord[1], keys[0], circ_type,
                                   "%s_simple_circ_%d" % (SampleID, Circ_no)]
        return Ref_df
    elif len(keys) > 1:
        for i in range(len(keys)):

            Ref_df.loc[len(Ref_df)] = [str(chr), coord[i][0], coord[i][1],keys[i],circ_type,"%s_simple_circ_%d" % (SampleID,Circ_no)]
        return Ref_df

def Complex_full_length(file_name,bamfile,bamname,mapQ):
    """ returns the total number of columns needed for the most complex circle"""
    total_col_len = 0
    with open(file_name) as f:
        for line in f:
            line_value = line.strip().split()
            coord = line_value[0]
            start = line_value[1]
            end = line_value[2]
            complex_dict, circ_type, circ_coord, circ_chr = Complex(bamfile,mapQ, str(coord), int(start), int(end))
            print("COMPLEX", complex_dict)
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
            d["temp"] = [chr_full[i][0], int(min(coord_full[i])), int(max(coord_full[i])),
                         int(max(coord_full[i])) - int(min(coord_full[i]))]
            tot_len += (max(coord_full[i]) - min(coord_full[i]))
        else:
            d["temp"].extend([chr_full[i][0], int(min(coord_full[i])), int(max(coord_full[i])),
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
        d["temp"].extend(extend_list)
    complex_df = pd.DataFrame.from_dict(d, orient='index', columns=Coord_Col)

    add_col = ['Total_len','Read_No', 'Read_IDs', 'Circle_type', 'Circle_ID']

    Read_ID_str = str(list(list(dict.keys()))).replace(" ", "")
    add_val = [int(tot_len),len(list(dict.keys())), Read_ID_str, circ_type,
           "%s_Complex_circ_%d" % (SampleID, Circ_no)]

    for i in range(len(add_col)):
        complex_df.insert(loc=len(complex_df.columns), column=add_col[i], value=add_val[i])

    return complex_df


def Complex_read_info(dict,circ_no,circ_type,bamname):
    """ returns two dataframe with information for each read"""

    #finds the longest value and the corresponding key to create number of col in .BED
    SampleID = bamname.split('.')[0]
    df_col = ["Chr", "Start", "End"]
    max_key, max_value = max(dict.items(), key=lambda x: len(set(x[1])))
    rep = int(len(max_value)/len(df_col))
    new_dict = {}
    for k,v in dict.items():
        if len(max_value) == len(v):
            new_dict[k] = v
        if len(max_value) > len(v):
            #add nan for those dict keys to ensure same length of the dicts values
            Nan_list = ["nan"]
            Nan_ext = (len(max_value)-len(v))
            Extend_list = Nan_list * int(Nan_ext)
            new_dict[k] = v + Extend_list

    #number of times df_col needed to be repeated
    Coord_Col = [j + "_no._" + str(i) for i in range(1, rep + 1) for j in df_col]
    read_df = pd.DataFrame.from_dict(new_dict, orient='index',columns=Coord_Col)
    read_df = df_float_to_int(read_df)

    add_col = ['Read_IDs', 'Circle_type', 'Circle_ID']
    #add_val = [list(new_dict.keys()),"Circle_type_%s" % circ_type,"%s_Complex_circ_%d" % (SampleID, circ_no)]

    df_final3 = pd.DataFrame(columns=add_col)
    for i in range(len(new_dict.keys())):
        chr_col = {add_col[0]:list(new_dict.keys())[i],add_col[1]:"%s" % circ_type,add_col[2]:"%s_Complex_circ_%d" % (SampleID, circ_no)}
        df_final3 = df_final3.append(chr_col, ignore_index=True)

    return read_df, df_final3

def Complex_circ_df(dict,coord_full,chr_full,Circ_no,circ_type,bamname):
    """ returns a dataframe with circular information for the complex circles, type 3 """
    d = {}
    tot_len = 0
    SampleID = bamname.split('.')[0]
    for i in range(len(coord_full)):
        if len(d) == 0:
            d["Circle_temp"] = [chr_full[i][0], int(min(coord_full[i])), int(max(coord_full[i])),
                         int(max(coord_full[i])) - int(min(coord_full[i]))]
            tot_len += (max(coord_full[i]) - min(coord_full[i]))
        else:
            d["Circle_temp"].extend([chr_full[i][0], int(min(coord_full[i])), int(max(coord_full[i])),
                              int(max(coord_full[i])) - int(min(coord_full[i]))])
            tot_len += (int(max(coord_full[i])) - int(min(coord_full[i])))
    df_col = ["Chr", "Start", "End", "Length"]

    #number of times df_col needed to be repeated
    rep = len(coord_full)

    #creates a dataframe with chromosome information
    Coord_Col = [j + "_no._" + str(i) for i in range(1, rep + 1) for j in df_col]
    complex_df = pd.DataFrame.from_dict(d, orient='index',columns=Coord_Col)

    #creates a dataframe with additional information
    add_col = ['Total_len', 'Read_No', 'Read_IDs', 'Circle_type', 'Circle_ID']
    Read_ID_str = str(list(list(dict.keys()))).replace(" ", "")
    add_val = [int(tot_len), len(list(dict.keys())), Read_ID_str, "%s" % circ_type,
               "%s_Complex_circ_%d" % (SampleID, Circ_no)]
    chr_col = {add_col[0]: add_val[0], add_col[1]: add_val[1],add_col[2]: add_val[2],
         add_col[3]: add_val[3],add_col[4]: add_val[4]}
    df_final5 = pd.DataFrame(columns=add_col)
    df_final5 = df_final5.append(chr_col,ignore_index=True)

    return complex_df, df_final5

def BED_Coverage(dir,bamfile,bamname,overlap,savedir,mapQ):
    """ Creates several bed file containing information regarding, both simple and complex circles"""
    Bedfile = pybedtools.example_bedtool(dir+bamname)
    Cov = Bedfile.genome_coverage(bg=True)
    Merged = Cov.merge(d=overlap)

    Simple_count = 1
    Simple_circ = pd.DataFrame()
    Complex_count = 1
    Complex_df = pd.DataFrame()
    read_df_full = pd.DataFrame()
    Final_5 = pd.DataFrame()
    Read_Final3 = pd.DataFrame()

    for region in Merged:
        coord = region[0]
        start = region[1]
        end = region[2]
        #print("THE CHROMOSOMAL COORDINATES",coord,start,end)
        circle_dict, circ_type, circ_coord, circ_chr = Complex(bamfile,mapQ, str(coord), int(start), int(end))
        if circ_type == 1:
            circ_bed = Simple_circ_df(circle_dict,Simple_count,"high_conf",bamname)
            rows = pd.concat([Simple_circ,circ_bed])
            Simple_circ = rows
            Simple_count += 1

        elif circ_type == 2:
            if len(list(circle_dict.keys())[0]) > 1:
                circ_bed = Simple_circ_df(circle_dict, Simple_count, "conf", bamname)
                rows = pd.concat([Simple_circ, circ_bed])
                Simple_circ = rows
                Simple_count += 1
            else:
                circ_bed = Simple_circ_df(circle_dict, Simple_count, "low_conf", bamname)
                rows = pd.concat([Simple_circ, circ_bed])
                Simple_circ = rows
                Simple_count += 1
        elif circ_type == 3:
            cird_df, Last5 = Complex_circ_df(circle_dict, circ_coord, circ_chr, Complex_count, "complex", bamname)
            read_df, Last_3 = Complex_read_info(circle_dict, Complex_count, "complex", bamname)

            Complex_df = Complex_df.append(cird_df, sort=False).fillna("nan")
            Final_5 = Final_5.append(Last5)

            read_df_full = read_df_full.append(read_df, sort=False).fillna("nan")
            Read_Final3 = Read_Final3.append(Last_3)

            Complex_count += 1
        else:
            continue

    Simple_bed = bt.BedTool.from_dataframe(Simple_circ)
    Simple_bed.saveas("{0}Simple_circles_cov3.bed".format(str(savedir)))

    resulting_df = pd.concat([Complex_df.reset_index(drop=True), Final_5.reset_index(drop=True)], axis=1)
    Bed_test = bt.BedTool.from_dataframe(resulting_df)
    Bed_test.saveas("{0}Complex_circles.bed".format(str(savedir)))

    resulting_read = pd.concat([read_df_full.reset_index(drop=True), Read_Final3.reset_index(drop=True)], axis=1)
    resulting_read = df_float_to_int(resulting_read)

    Bed_reads = bt.BedTool.from_dataframe(resulting_read)
    Bed_reads.saveas("{0}Complex_reads.bed".format(str(savedir)))

def FULL(file_name,bamfile,bamname,savedir,mapQ):
    Simple_count = 1
    Simple_circ = pd.DataFrame()
    Complex_count = 1
    Complex_df = pd.DataFrame()
    read_df_full = pd.DataFrame()
    Final_5 = pd.DataFrame()
    Read_Final3 = pd.DataFrame()

    with open(file_name) as f:
        for line in f:
            region = line.strip().split()
            coord = region[0]
            start = region[1]
            end = region[2]
            #print("THE CHROMOSOMAL COORDINATES", coord, start, end)
            circle_dict, circ_type, circ_coord, circ_chr = Complex(bamfile, mapQ, str(coord), int(start), int(end))
            print("DICT",circle_dict)
            if circ_type == 1:
                circ_bed = Simple_circ_df(circle_dict, Simple_count, "high_conf", bamname)
                rows = pd.concat([Simple_circ, circ_bed])
                Simple_circ = rows
                Simple_count += 1

            elif circ_type == 2:
                if len(list(circle_dict.keys())[0]) > 1:
                    circ_bed = Simple_circ_df(circle_dict, Simple_count, "conf", bamname)
                    rows = pd.concat([Simple_circ, circ_bed])
                    Simple_circ = rows
                    Simple_count += 1
                else:
                    circ_bed = Simple_circ_df(circle_dict, Simple_count, "low_conf", bamname)
                    rows = pd.concat([Simple_circ, circ_bed])
                    Simple_circ = rows
                    Simple_count += 1
            elif circ_type == 3:
                print("CIRCLE_DICT", circle_dict,circ_coord,circ_chr,Complex_count)
                cird_df, Last5 = Complex_circ_df(circle_dict, circ_coord, circ_chr, Complex_count, "complex", bamname)
                read_df, Last_3 = Complex_read_info(circle_dict, Complex_count, "complex", bamname)

                Complex_df = Complex_df.append(cird_df, sort=False).fillna("nan")
                Final_5 = Final_5.append(Last5)

                read_df_full = read_df_full.append(read_df, sort=False).fillna("nan")
                Read_Final3 = Read_Final3.append(Last_3)

                Complex_count += 1
            else:
                continue

        Simple_bed = bt.BedTool.from_dataframe(Simple_circ)
        Simple_bed.saveas("{0}Simple_circles_1000.bed".format(str(savedir)))

        resulting_df = pd.concat([Complex_df.reset_index(drop=True), Final_5.reset_index(drop=True)], axis=1)
        Bed_test = bt.BedTool.from_dataframe(resulting_df)
        Bed_test.saveas("{0}Complex_circles_1000.bed".format(str(savedir)))

        resulting_read = pd.concat([read_df_full.reset_index(drop=True), Read_Final3.reset_index(drop=True)], axis=1)
        resulting_read = df_float_to_int(resulting_read)

        Bed_reads = bt.BedTool.from_dataframe(resulting_read)
        Bed_reads.saveas("{0}Complex_reads_1000.bed".format(str(savedir)))

if __name__ == '__main__':
     bamfile = ps.AlignmentFile("/isdata/common/wql443/NanoCircle/BC01.bam","rb")
     FULL("/isdata/common/wql443/NanoCircle/BC01_1000_cov.bed",bamfile,"BC01bam","/isdata/common/wql443/NanoCircle/",60)

"""
#aln_hg19_Qiagen, Qiagen_1000_cov.bed
ID = "HS45"
samples = "5.0samples"
dir = "/isdata/common/wql443/Projects/Sperm_2019/Analysis/Samples/{0}/{1}/".format(samples,ID)
bamfile = ps.AlignmentFile(dir+"{0}.aln_hg19.bam".format(ID),"rb")
FULL("/isdata/common/wql443/Projects/Sperm_2019/Analysis/Samples/{0}/{1}/{1}_1000_cov.bed".format(samples,ID),bamfile,
     "{0}.aln_hg19.bam".format(ID),dir,60)

#(file_name,bamfile,bamname,savedir,mapQ)
#FULL("/isdata/common/wql443/Projects/Sperm_2019/Analysis/Samples/5.0samples/HS45/HS45_1000_cov.bed",bamfile,
#     "{0}.aln_hg19.bam".format(ID),"/isdata/common/wql443/{0}_".format(ID),60)

bamfile = ps.AlignmentFile("/isdata/common/wql443/NanoCircle/ArgParse_test/Real_test_data/chr1_243928620_243938331_region.bam","rb")
FULL("/isdata/common/wql443/NanoCircle/ArgParse_test/Real_test_data/Cov.bed",bamfile,
     "chr1_243928620_243938331_region.bam","/isdata/common/wql443/NanoCircle/ArgParse_test/Real_test_data/",60)"""