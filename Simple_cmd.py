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
    def __init__(self,region,bamfile,output,MapQ):
        self.region = region
        self.bamfile = bamfile
        self.output = output
        self.MapQ = MapQ

    def Prim_dict(self,reg, start, end):
        """Creates a dictionary of the primary alignments for the soft-clipped reads"""
        Prim_dict = {}
        for read in self.bamfile.fetch(reg, start, end, multiple_iterators=True):
            # checks if soft-clipped
            if IS_SA(read, self.MapQ) == True:
                pos = ps.AlignedSegment.get_reference_positions(read)

                # creates dict with key being read_id and value being read ref positions
                Prim_dict[read.query_name] = [pos[0], pos[-1]]
        return Prim_dict

    def Supplement_dict(self, reg, start, end):
        """Creates a dictionary of the primary soft-clipped reads and
        their corresponding supplementary alignments"""

        Primary_dict = Prim_dict(self, reg, start, end)
        SA_dict = {}

        for read in self.bamfile.fetch(reg, start, end, multiple_iterators=True):

            # only extract those supp in primary
            if read.query_name in Primary_dict.keys():
                if IS_supp(read, self.MapQ) == True:

                    # extracts positions
                    supp_pos = ps.AlignedSegment.get_reference_positions(read)
                    prim_pos = Primary_dict[read.query_name]

                    # Adding to positions based on direction to the empty dict
                    if read.query_name not in SA_dict.keys():

                        # From left to right
                        if Right_dir(prim_pos, supp_pos) == True:
                            SA_dict[read.query_name] = [prim_pos[0], supp_pos[-1]]

                        # From right to left
                        if Left_dir(prim_pos, supp_pos) == True:
                            SA_dict[read.query_name] = [supp_pos[0], prim_pos[-1]]

                        # From left to right once
                        if Right_circ_dir(prim_pos, supp_pos) == True:
                            SA_dict[read.query_name] = [prim_pos[0], prim_pos[-1]]

                        # From right to left once
                        if Left_circ_dir(prim_pos, supp_pos) == True:
                            SA_dict[read.query_name] = [prim_pos[0], prim_pos[-1]]

                    # Appends for reads with several supplementary alignments their position
                    elif read.query_name in SA_dict.keys():

                        if Right_dir(prim_pos, supp_pos) == True:
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



    def Read_position(self,reg,start,end):
        """Creates a dictionary of the primary alignments for the soft-clipped reads"""
        Prim_dict = {}
        Supp_dict = {}

        Read_File = Utils.SamBamParser(self.bamfile)

        for read in Read_File.fetch(reg, start, end, multiple_iterators=True):
            # checks if soft-clipped and has a supplementary tag and a mapping quality equal to threshold
            if read.cigar[0][0] == 4 and read.has_tag("SA") and read.mapping_quality >= self.MapQ:

                #print("READ_NAME",read.query_name, int(read.reference_start), int(read.reference_end)-1)
                #print(ps.AlignedSegment.get_reference_positions(read)[0],ps.AlignedSegment.get_reference_positions(read)[-1])

                # Creating a dictionary with primary alignment coordinates
                #read.reference_end moves one past the last aligned residue
                prim_pos = [int(read.reference_start)+1, int(read.reference_end)]
                Prim_dict[read.query_name] = prim_pos

                # supplementary alignment
                Tag = read.get_tag("SA").split(',')

                # Creating a dictionary with primary alignment and supplementary alignment coordinates

                # a single supplementary alignment
                if len(Tag) == 6:
                    #print("Single supp")
                    if int(Tag[4]) >= self.MapQ:
                        supp_pos = [int(Tag[1]), int(Tag[1]) + Utils.CIGAR_len(Tag[3])]
                        #print(supp_pos)

                        # From left to right across breakpoint
                        if Utils.Right_dir(prim_pos, supp_pos) == True:
                            Supp_dict[read.query_name] = {'start':[prim_pos[0]],'end':[supp_pos[-1]]}

                        # From right to left
                        if Utils.Left_dir(prim_pos, supp_pos) == True:
                            Supp_dict[read.query_name] = {'start':[supp_pos[0]],'end':[prim_pos[-1]]}

                        # primary read covers entire circle, and supp is between prim, but use both anyways due to
                        # perhaps having repeat alignment. Which is helpful for most common
                        if Utils.Circ_once(prim_pos, supp_pos) == True:
                            Supp_dict[read.query_name] = {'start':[prim_pos[0],supp_pos[0]],
                                                          'end':[prim_pos[-1],supp_pos[-1]]}
                    else:
                        continue

                # with multiple supplementary we need to keep the possible start or ends.
                elif len(Tag) > 6:
                    #print("Multiple supp")
                    MapQ_val = Tag[4::5]
                    #print("MAPPIGN QUAL",MapQ_val)
                    Supp_start = Tag[1::5]

                    mult_coord = []

                    for i in range(len(MapQ_val)):
                        if int(MapQ_val[i]) >= self.MapQ:

                            #creates a list of list, with all possible supp coordinats
                            mult_coord.append([int(Supp_start[i]),
                                               int(Supp_start[i]) + Utils.CIGAR_len(Tag[3::5][i])])

                            # takes the first set of supp coordinates within the multiple coordinates
                            # and creates and entry in the Supp dictionary
                            if read.query_name not in Supp_dict.keys():
                                #coordinates based on the read orientation

                                if Utils.Right_dir(prim_pos, mult_coord[0]) == True:
                                    Supp_dict[read.query_name] = {'start':[prim_pos[0]], 'end':[mult_coord[0][-1]]}

                                if Utils.Left_dir(prim_pos, mult_coord[0]) == True:
                                    Supp_dict[read.query_name] = {'start':[mult_coord[0][0]], 'end':[prim_pos[-1]]}

                                #read passes one over the circle, use both supp and prim align
                                if Utils.Circ_once(prim_pos, mult_coord[0]) == True:
                                    Supp_dict[read.query_name] = {'start':[prim_pos[0]],'end':[prim_pos[-1]]}
                                    #print("Primary and multiple",prim_pos,mult_coord)

                            #considers the additional supplementary coordinates
                            if read.query_name in Supp_dict.keys():
                                for i in mult_coord[1:]:
                                    #print("Next coordinate set",i)

                                    if Utils.Right_dir(prim_pos, i) == True:
                                        Supp_dict[read.query_name]['end'].append(i[-1])
                                        #print("right v3")
                                    if Utils.Left_dir(prim_pos, i) == True:
                                        Supp_dict[read.query_name]['start'].append(i[0])
                                        #print("left v3")
                                    if Utils.Circ_once(prim_pos, i) == True:
                                        #print("full ONCE v3")
                                        Supp_dict[read.query_name] = {'start': [prim_pos[0]], 'end': [prim_pos[-1]]}
                                #print("Final_Supp",Supp_dict)
                        else:
                            continue
        return Supp_dict

    def Reduce_coord(self,pos_dict,pos):

        pos_point = []
        for k, v in pos_dict.items():
            pos_point.extend(v[pos])
        # find the most common end point

        most_common = Counter(pos_point).most_common()

        # the frequency value for each individual coordinates
        freq_val = [item[-1] for item in most_common]

        if any(x > 1 for x in freq_val):
            # Adjust the end points if the most common end point is found
            for k, v in pos_dict.items():
                if most_common[0][0] in v[pos]:
                    pos_dict[k][pos] = [most_common[0][0]]
            return pos_dict

        else:
            return pos_dict

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

    def Single_coord(self,reg,start,end):

        pos_dict = self.Read_position(reg, start, end)
        #print("lenght",len(pos_dict.keys()),pos_dict)
        # reduce start and end coordinates to those which are most common
        temp_dict = self.Reduce_coord(pos_dict, 'start')
        modify_dict = self.Reduce_coord(temp_dict, 'end')
        #print("lenght reduce start", len(temp_dict.keys()), temp_dict)
        #print("lenght reduce end", len(modify_dict.keys()), modify_dict)

        all_soft =  len(pos_dict.keys())
        # creates a single list with all s tart and end coordinates
        start_list = []
        end_list = []
        for v in modify_dict.values():
            start_list.extend(v['start'])
            end_list.extend(v['end'])

        if start_list and end_list != []:
            occ1, start_freq = self.most_frequent(start_list, 0)
            occ2, end_freq = self.most_frequent(end_list, 1)
            #print("start,end",start_freq,end_freq)
            if occ1 and occ2 == 1:
                reads = []
                chr = [reg]
                if start_freq < end_freq:
                    new_val = chr + [start_freq, end_freq]
                else:
                    #print("LLLOLLLLLLLLLLLLLL")
                    new_val = chr + [end_freq, start_freq]

                for k, v in modify_dict.items():
                    if any(i == start_freq or i == end_freq for i in v['start']):
                        reads.append(k)
                    elif any(i == start_freq or i == end_freq for i in v['end']):
                        reads.append(k)

                final_dict = {tuple(sorted(reads)): new_val}

                if len(list(final_dict.keys())[0]) != 1:
                    type = 1
                    return all_soft,type, final_dict

                else:
                    type = 2
                    return all_soft,type, final_dict

            else:
                type = 2
                chr = [reg]
                new_val = chr + [start_freq, end_freq]
                final_dict = {tuple(sorted(modify_dict.keys())): new_val}
                return all_soft,type, final_dict

        else:
            type = 0
            # these dicts are empty, and serves to continue as some regions might not create a circle
            return all_soft,type, modify_dict

    def Simple_circ_df(self,Circ_dict, circ_type,soft_no):
        """ returns a dataframe with circular information for the simple circles"""
        df_col = ["Chr", "Start", "End"]
        simple_df = pd.DataFrame.from_dict(Circ_dict, orient='index', columns=df_col)
        simple_df = simple_df.sort_values(by=df_col)
        simple_df['Length'] = simple_df.apply(lambda x: x['End'] - x['Start'], axis=1)
        # add_col = ['Read_No','Read_IDs','Circle_type','Circle_ID']
        add_col = ['Read_No', 'Circle_type',"Soft_clip"]
        # convert list of ID's to 1 long string as to insert it as a single column in the df

        add_val = [len(list(Circ_dict.keys())[0]), circ_type,soft_no]

        for i in range(len(add_col)):
            simple_df.insert(loc=len(simple_df.columns), column=add_col[i], value=add_val[i])
        return simple_df


    def Region(self):
        Simple_circ = pd.DataFrame()
        with open(self.region) as f:
            for line in f:
                region = line.strip().split()
                chr = region[0]
                start = region[1]
                end = region[2]
                #print("-------------------------------")
                #print("REGION",str(chr),int(start), int(end))
                soft_no,circ_type,circle_dict = self.Single_coord(str(chr),int(start), int(end))
                #print("modify",circle_dict)
                if circ_type == 1:
                    circ_bed = self.Simple_circ_df(circle_dict, "high_conf",soft_no)
                    rows = pd.concat([Simple_circ, circ_bed])
                    Simple_circ = rows

                elif circ_type == 2:
                    if len(list(circle_dict.keys())[0]) > 1:
                        circ_bed = self.Simple_circ_df(circle_dict, "conf",soft_no)
                        rows = pd.concat([Simple_circ, circ_bed])
                        Simple_circ = rows

                    else:
                        circ_bed = self.Simple_circ_df(circle_dict, "low_conf",soft_no)
                        rows = pd.concat([Simple_circ, circ_bed])
                        Simple_circ = rows
                else:
                    continue
        Simple_bed = bt.BedTool.from_dataframe(Simple_circ)
        Simple_bed.saveas(self.output)


if __name__ == '__main__':
    print("main")

