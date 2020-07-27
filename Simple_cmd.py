#!/usr/bin/env python3

import Utils
import pysam as ps
import numpy as np
from collections import Counter
from itertools import groupby
import pandas as pd
import numpy as np
import pybedtools as bt
import operator

class Simple_circ:
    #Command 1
    def __init__(self,regions,bamfile,output,MapQ):
        self.regions = regions
        self.bamfile = bamfile
        self.output = output
        self.MapQ = MapQ

    def Prim_dict(self,bamfile,reg, start, end):
        """Creates a dictionary of the primary alignments for the soft-clipped reads"""
        print("PRIMARY WORK")

        Prim_dict = {}
        for read in bamfile.fetch(reg, start, end, multiple_iterators=True):
            # checks if soft-clipped
            if Utils.IS_SA(read, self.MapQ) == True:
                pos = ps.AlignedSegment.get_reference_positions(read)

                # creates dict with key being read_id and value being read ref positions
                Prim_dict[read.query_name] = [pos[0], pos[-1]]
        return Prim_dict

    def Supplement_dict(self,bamfile,reg, start, end):
        """Creates a dictionary of the primary soft-clipped reads and
        their corresponding supplementary alignments"""

        Primary_dict = self.Prim_dict(bamfile,reg, start, end)
        SA_dict = {}

        for read in bamfile.fetch(reg, start, end, multiple_iterators=True):

            # only extract those supp in primary
            if read.query_name in Primary_dict.keys():
                if Utils.IS_supp(read, self.MapQ) == True:
                    print("LOL")
                    # extracts positions
                    supp_pos = ps.AlignedSegment.get_reference_positions(read)
                    prim_pos = Primary_dict[read.query_name]

                    # Adding to positions based on direction to the empty dict
                    if read.query_name not in SA_dict.keys():

                        # From left to right
                        if Utils.Right_dir(prim_pos, supp_pos) == True:
                            SA_dict[read.query_name] = [prim_pos[0], supp_pos[-1]]

                        # From right to left
                        if Utils.Left_dir(prim_pos, supp_pos) == True:
                            SA_dict[read.query_name] = [supp_pos[0], prim_pos[-1]]

                        # From left to right once
                        if Utils.Right_circ_dir(prim_pos, supp_pos) == True:
                            SA_dict[read.query_name] = [prim_pos[0], prim_pos[-1]]

                        # From right to left once
                        if Utils.Left_circ_dir(prim_pos, supp_pos) == True:
                            SA_dict[read.query_name] = [prim_pos[0], prim_pos[-1]]

                    # Appends for reads with several supplementary alignments their position
                    elif read.query_name in SA_dict.keys():

                        if Utils.Right_dir(prim_pos, supp_pos) == True:
                            if supp_pos[-1] not in SA_dict[read.query_name]:
                                SA_dict[read.query_name].append(supp_pos[-1])
                        # From right to left
                        if Utils.Left_dir(prim_pos, supp_pos) == True:
                            if prim_pos[-1] not in SA_dict[read.query_name]:
                                SA_dict[read.query_name].append(prim_pos[-1])

                        # From left to right once
                        if Utils.Right_circ_dir(prim_pos, supp_pos) == True:
                            if prim_pos[-1] not in SA_dict[read.query_name]:
                                SA_dict[read.query_name].append(prim_pos[-1])

                        # From right to left once
                        if Utils.Left_circ_dir(prim_pos, supp_pos) == True:
                            if prim_pos[-1] not in SA_dict[read.query_name]:
                                SA_dict[read.query_name].append(prim_pos[-1])

        return SA_dict

    def reduce_end_coords(self,bamfile,reg, start, end):
        """ Reduce the number of end coordinates for those reads with several
        supplementary alignments by comparing the position to the other end positions"""

        Coord = self.Supplement_dict(bamfile,reg, start, end)
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

    def most_frequent(self,List, pos):
        count1 = Counter(List)
        dic_val = list(count1.values())
        max_val = max(count1.values())
        occurence = dic_val.count(max_val)

        # single most occuring coordinate
        if occurence == 1:
            # print("COMMMON VALUE", count1.most_common(1)[0][0])
            occ = 1
            return occ, count1.most_common(1)[0][0]

        # creates the interval immediately
        else:
            most_common_no = list(map(operator.itemgetter(1), count1.most_common()))
            most_common_idx = [i for i, x in enumerate(most_common_no) if x == max(most_common_no)]

            # if there are several equal common values
            if max(most_common_no) > 1:
                # print(max(most_common_no))
                # print(most_common_idx)
                # creating the most maximum value
                common_val = [count1.most_common()[i][0] for i in most_common_idx]
                occ = 1
                # print("COMMMON VALUE",common_val)
                return occ, max(common_val)

            # if there isnt a most occuring coordinate, the interval is created
            else:
                occ = 2
                # start
                if pos == 0:
                    return occ, min(List)
                # end
                if pos == 1:
                    return occ, max(List)

    def Reads(self,bamfile, reg, start, end, reads_IDs):
        """ extract those reads which is found using the single coordinate dictionary"""
        Read_list = []
        for read in bamfile.fetch(reg, start, end, multiple_iterators=True):
            if read.query_name in reads_IDs:
                if Utils.IS_SA(read, self.MapQ) == True:
                    Read_list.append(read)
        return Read_list

    def chr_coord_sort(self,chrlist, coordlist):
        """ Sort a list of chromosomes and their coordinates using the index of the numerically sorted coordinates"""
        coord_idx = np.argsort(coordlist)
        Coord_sort = [coordlist[i] for i in coord_idx]
        chr_sort = [chrlist[i] for i in coord_idx]
        return Coord_sort, chr_sort

    def Single_coord(self, bamfile,reg, start, end):
        """ returns the most frequent start and end coordinate for a circle in a dictionary with the key
        being all the reads. If not having a coordinate more frequent than others, a dictionary with all reads
        and the interval for potential circle """
        Circle_dict = self.reduce_end_coords(bamfile,reg, start, end)
        start = [i[0] for i in Circle_dict.values()]
        # taking all possible end coordinates and merge into one list
        end = sum([i[1:] for i in Circle_dict.values()], [])
        # extract the most common
        if start or end != []:
            occ1, start_freq = self.most_frequent(start, 0)
            # print("occurence,",occ1,"start",start_freq)
            occ2, end_freq = self.most_frequent(end, 1)
            # print("occurence,", occ2, "start", end_freq)
            if occ1 and occ2 == 1:
                reads = []
                chr = [reg]
                # print("START COORDS", start)
                # print("END COORD",end)
                if start_freq < end_freq:
                    new_val = chr + [start_freq, end_freq]
                else:
                    # print("lol")
                    # print(end_freq,start_freq)
                    new_val = chr + [end_freq, start_freq]
                # print("Start and end",new_val)
                for k, v in Circle_dict.items():
                    if any(i == start_freq or i == end_freq for i in v):
                        reads.append(k)
                final_dict = {tuple(sorted(reads)): new_val}

                # Multiple reads
                if len(list(final_dict.keys())[0]) != 1:
                    type = 1
                    return type, final_dict

                # Single read
                else:
                    type = 2
                    return type, final_dict

            # not a more supported read
            else:
                type = 2
                # mangler stadig at lave det der overlaps
                chr = [reg]
                new_val = chr + [start_freq, end_freq]
                final_dict = {tuple(sorted(Circle_dict.keys())): new_val}
                # man kunne eventuelt her bare returnere Single coord
                return type, final_dict
        else:
            type = 0
            # these dicts are empty, and serves to continue as some regions might not create a circle
            return type, Circle_dict

    def Is_Simple(self, bamfile,reg, start, end):
        """ Check for potential complex circles by reads aligning across the genome. Afterwards it returns the simple circles
         with 1 set of coordinates (type I) several set of coordinates (type II) and the complex circles (type III)"""
        Type, Sing_dict = self.Single_coord(bamfile,reg, start, end)

        Complex_dict = {}
        Total_coord = []
        Total_chr = []
        Total_overlap = []

        if len(Sing_dict) == 1:
            reads_IDs = list(Sing_dict.keys())[0]
            Read_list = self.Reads(bamfile, reg, start, end, reads_IDs)
        else:
            reads_IDs = list(Sing_dict.keys())
            Read_list = self.Reads(bamfile, reg, start, end, reads_IDs)

        for read in Read_list:
            coord1 = list(Sing_dict.values())[0][-2]
            coord2 = list(Sing_dict.values())[0][-1]
            Total_chr.extend((reg, reg))
            Total_coord.extend((coord1, coord2))
            # print("TOAL CHR",Total_chr,Total_coord)
            Tag = read.get_tag("SA").split(';')[:-1]
            Coord_list = []
            chroms = []
            cigar_len = []

            # examining for complexity
            for Tagelem in Tag:
                # splitting the individual aligment info up
                Column_list = Tagelem.split(',')
                # print("Column",Column_list)
                chrom = Column_list[0]
                length = Utils.CIGAR_len(Column_list[3])
                cigar_len.append(length)
                # 1 - based positions
                pos_start = int(Column_list[1]) - 1
                pos_end = int(Column_list[1]) + length - 1
                # the overlaps between coordinates for grouping
                overlap = sum(cigar_len) * 4
                Total_overlap.append(overlap)
                # if the supp align is in between the circular input region then we already know the breakpoint from Sing_dict
                if chrom == reg and start - overlap <= pos_start <= end + overlap and start - overlap <= pos_end <= end + overlap:
                    continue

                elif int(Column_list[4]) >= self.MapQ:
                    # creates a coordinate list
                    Coord_list.append(pos_start)
                    Coord_list.append(pos_end)

                    # append chr twice to ensure same length as coord_list
                    chroms.append(chrom)
                    chroms.append(chrom)

                    Total_chr.extend((chrom, chrom))

                    Total_coord.extend((pos_start, pos_end))

                if Coord_list != []:
                    continue

        # the simple circles
        if Complex_dict == {}:
            if Type == 0 or 1 or 2:
                return Sing_dict, Type, None, None
        else:
            Type = 3
            return None, Type, None, None

    def Simple_circ_df(self,beddict, Circ_no, circ_type):
        """ returns a dataframe with circular information for the simple circles, type 1 and 2 """
        df_col = ["Chr", "Start", "End"]
        simple_df = pd.DataFrame.from_dict(beddict, orient='index', columns=df_col)
        simple_df = simple_df.sort_values(by=df_col)
        simple_df['Length'] = simple_df.apply(lambda x: x['End'] - x['Start'], axis=1)
        # add_col = ['Read_No','Read_IDs','Circle_type','Circle_ID']
        add_col = ['Read_No', 'Circle_type', 'Circle_ID']
        # convert list of ID's to 1 long string as to insert it as a single column in the df

        # Read_ID_str = str(list(list(beddict.keys())[0])).replace(" ","")
        # add_val = [len(list(beddict.keys())[0]),Read_ID_str,circ_type,"%s_simple_circ_%d" % (SampleID,Circ_no)]

        add_val = [len(list(beddict.keys())[0]), circ_type, "simple_circ_%d" % Circ_no]

        for i in range(len(add_col)):
            simple_df.insert(loc=len(simple_df.columns), column=add_col[i], value=add_val[i])
        return simple_df

    def Circle_output(self):
        Simple_count = 1
        Simple_circ = pd.DataFrame()

        Read_File = Utils.SamBamParser(self.bamfile)

        with open(self.regions) as f:
            for line in f:
                region = line.strip().split()
                coord = region[0]
                start = region[1]
                end = region[2]
                # print("THE CHROMOSOMAL COORDINATES", coord, start, end)
                circle_dict, circ_type, circ_coord, circ_chr = self.Is_Simple(Read_File,str(coord), int(start), int(end))
                print("DICT", circle_dict)
                if circ_type == 1:
                    circ_bed = self.Simple_circ_df(circle_dict, Simple_count, "high_conf")
                    rows = pd.concat([Simple_circ, circ_bed])
                    Simple_circ = rows
                    Simple_count += 1

                elif circ_type == 2:
                    if len(list(circle_dict.keys())[0]) > 1:
                        circ_bed = self.Simple_circ_df(circle_dict, Simple_count, "conf")
                        rows = pd.concat([Simple_circ, circ_bed])
                        Simple_circ = rows
                        Simple_count += 1
                    else:
                        circ_bed = self.Simple_circ_df(circle_dict, Simple_count, "low_conf")
                        rows = pd.concat([Simple_circ, circ_bed])
                        Simple_circ = rows
                        Simple_count += 1
                elif circ_type == 3:
                    continue
                else:
                    continue

            Simple_bed = bt.BedTool.from_dataframe(Simple_circ)
            Simple_bed.saveas("{0}.bed".format(str(self.output)))

if __name__ == '__main__':
    print("main")

