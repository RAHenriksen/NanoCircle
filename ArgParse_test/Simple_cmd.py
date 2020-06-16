#!/usr/bin/env python3

import Utils
import pysam as ps
import numpy as np
from collections import Counter
from itertools import groupby

class Simple_circ:
    #Command 1
    def __init__(self,region,bamfile,output,MapQ):
        self.region = region
        self.bamfile = bamfile
        self.output = output
        self.MapQ = MapQ

    def Primary(self,bamfile, mapQ, reg, start, end):
        """Creates a dictionary of the primary alignments for the soft-clipped reads"""
        Prim_dict = {}

        Read_File = Utils.SamBamParser(bamfile)

        for read in Read_File.fetch(reg, start, end, multiple_iterators=True):
            # checks if soft-clipped
            if Utils.IS_SA(read, mapQ) == True:
                prim_pos = [int(read.reference_start), int(read.reference_end) - 1]

                # creates dict with key being read_id and value being read ref positions
                Prim_dict[read.query_name] = prim_pos
        return Prim_dict

    def Supplementary(self,bamfile,MapQ,reg,start,end):
        """Creates a dictionary of the primary alignments for the soft-clipped reads"""
        Prim_dict = {}
        Supp_dict = {}

        Read_File = Utils.SamBamParser(bamfile)

        for read in Read_File.fetch(reg, start, end, multiple_iterators=True):
            # checks if soft-clipped and has a supplementary tag and a mapping quality equal to threshold
            if read.cigar[0][0] == 4 and read.has_tag("SA") and read.mapping_quality >= MapQ:
                print("-------------------------------")
                print("READ_NAME",read.query_name)
                # Creating a dictionary with primary alignment coordinates
                #read.reference_end moves one past the last aligned residue
                prim_pos = [int(read.reference_start), int(read.reference_end)-1]
                Prim_dict[read.query_name] = prim_pos

                # supplementary alignment
                Tag = read.get_tag("SA").split(',')

                # Creating a dictionary with primary alignment and supplementary alignment coordinates

                # a single supplementary alignment
                if len(Tag) == 6:
                    print("Single supp")
                    if int(Tag[4]) >= MapQ:
                        supp_pos = [int(Tag[1]), int(Tag[1]) + Utils.CIGAR_len(Tag[3])]

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

                # with multiple supplementary we need to keep the possible start or ends.
                elif len(Tag) > 6:
                    print("Multiple supp")
                    MapQ_val = Tag[4::5]
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
                            print("Primary and multiple",prim_pos,mult_coord)

                        print("First_supp dict", Supp_dict)

                    #considers the additional supplementary coordinates
                    if read.query_name in Supp_dict.keys():
                        for i in mult_coord[1:]:
                            print("Next coordinate set",i)

                            if Utils.Right_dir(prim_pos, i) == True:
                                Supp_dict[read.query_name]['end'].append(i[-1])
                                print("right v3")
                            if Utils.Left_dir(prim_pos, i) == True:
                                Supp_dict[read.query_name]['start'].append(i[0])
                                print("left v3")
                            if Utils.Circ_once(prim_pos, i) == True:
                                print("full ONCE v3")
                                Supp_dict[read.query_name] = {'start': [prim_pos[0]], 'end': [prim_pos[-1]]}
                        print("Final_Supp",Supp_dict)
        return Supp_dict

    def Region(self):
        with open(self.region) as f:
            for line in f:
                region = line.strip().split()
                chr = region[0]
                start = region[1]
                end = region[2]
                pos_dict = self.Supplementary(self.bamfile,self.MapQ,str(chr),int(start), int(end))


if __name__ == '__main__':
    print("lol")

    print("multiple", read.query_name)
    print("TAG", Tag)
    print("prim", prim_pos)
    print("all supp", mult_coord)
    # for each possible supp align check the coordinates
    print(mult_coord[0])

    if read.query_name in Supp_dict.keys():
        for i in mult_coord[1:]:
            print("Next coordinate set", i)

            # the minimum and maximum value for the created start and end coordinates
            # using both the primary and supp coordinates from above
            supp_coord = [min(list(Supp_dict[read.query_name].values())[0]),
                          max(list(Supp_dict[read.query_name].values())[1])]
            print("Coordinate interval", supp_coord)

            if Utils.Right_dir(supp_coord, i) == True:
                Supp_dict[read.query_name]['end'].append(i[-1])
                print("right v3")
            if Utils.Left_dir(supp_coord, i) == True:
                Supp_dict[read.query_name]['start'].append(i[0])
                print("left v3")
            if Utils.Circ_once(supp_coord, i) == True:
                print("full ONCE v3")
                Supp_dict[read.query_name]['start'].append(i[0])
                Supp_dict[read.query_name]['end'].append(i[-1])
        print("Final_Supp", Supp_dict)