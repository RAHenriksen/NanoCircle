#!/usr/bin/env python3

import Utils
import pysam as ps
import numpy as np
from collections import Counter
from itertools import groupby

print("Imports the simple script")

class Simple_circ:
    #Command 1
    def __init__(self,region,bamfile,output,MapQ):
        self.region = region
        self.bamfile = bamfile
        self.output = output
        self.MapQ = MapQ

    def Read_pos(self,bamfile,MapQ,reg,start,end):
        """Creates a dictionary of the primary alignments for the soft-clipped reads"""
        Pos_dict = {}

        Read_File = Utils.SamBamParser(bamfile)

        for read in Read_File.fetch(reg, start, end, multiple_iterators=True):
            # checks if soft-clipped
            if Utils.IS_SA(read,MapQ) == True:
                # primary alignment
                prim_chr = reg
                #read.reference_end moves one past the last aligned residue
                prim_pos = [int(read.reference_start), int(read.reference_end)-1]

                # supplementary alignment
                Tag = read.get_tag("SA").split(',')

                # a single supplementary alignment
                if len(Tag) == 6:
                    #print("Single", read.query_name)
                    #print("reg", reg, start, end)
                    #print("TAG",Tag)
                    #print("Prim",prim_pos)
                    if int(Tag[4]) >= self.MapQ:
                        supp_pos = [int(Tag[1]), int(Tag[1]) + Utils.CIGAR_len(Tag[3])]

                        if Utils.Right_dir(prim_pos, supp_pos) == True:
                            Pos_dict[read.query_name] = [prim_pos[0], supp_pos[-1]]

                            # From right to left
                        if Utils.Left_dir(prim_pos, supp_pos) == True:
                            Pos_dict[read.query_name] = [supp_pos[0], prim_pos[-1]]

                            # From left to right once
                        if Utils.Right_circ_dir(prim_pos, supp_pos) == True:
                            Pos_dict[read.query_name] = [prim_pos[0], prim_pos[-1]]

                            # From right to left once
                        if Utils.Left_circ_dir(prim_pos, supp_pos) == True:
                            Pos_dict[read.query_name] = [prim_pos[0], prim_pos[-1]]

                else:
                    MapQ_val = Tag[4::5]
                    Supp_start = Tag[1::5]


                    mult_coord = []

                    for i in range(len(MapQ_val)):
                        if int(MapQ_val[i]) >= self.MapQ:

                            #creates a list of list, with all possible supp coordinats
                            mult_coord.append([int(Supp_start[i]),
                                               int(Supp_start[i]) + Utils.CIGAR_len(Tag[3::5][i])])

                    #takes the first set of supp coordinates
                    supp_pos = mult_coord[0]
                    print("reg", reg, start, end)
                    print("multiple", read.query_name)
                    print("TAG", Tag)
                    print("prim", prim_pos)
                    print("all supp",mult_coord)
                    print("First supp",supp_pos)
                    # Creates a possible coordinate set comparing the primary and the first supp alignment
                    if Utils.Right_dir(prim_pos, supp_pos) == True:
                        Pos_dict[read.query_name] = [prim_pos[0], supp_pos[-1]]
                        print("left1.1")
                    if Utils.Left_dir(prim_pos, supp_pos) == True:
                        Pos_dict[read.query_name] = [supp_pos[0], prim_pos[-1]]
                        print("right.1")
                    if Utils.Right_circ_dir(prim_pos, supp_pos) == True:
                        print("full1.1")
                        Pos_dict[read.query_name] = [prim_pos[0], prim_pos[-1]]
                    if Utils.Left_circ_dir(prim_pos, supp_pos) == True:
                        print("full1.2")
                        Pos_dict[read.query_name] = [prim_pos[0], prim_pos[-1]]


                    print("TEST",Pos_dict[read.query_name])


                    # Creates the final coordinate set by comparing the coordinate set from above using primary and
                    # first supp, and appends the additional supplementary.

                    new_coord = Pos_dict[read.query_name]
                    for supp_i in mult_coord[1:]:
                        if Utils.Right_dir(new_coord, supp_i) == True:
                            if supp_i[-1] not in Pos_dict[read.query_name]:
                                Pos_dict[read.query_name].append(supp_i[-1])
                                print("left2.1")
                            # From right to left
                        if Utils.Left_dir(new_coord, supp_i) == True:
                            if new_coord[-1] not in Pos_dict[read.query_name]:
                                Pos_dict[read.query_name].append(new_coord[-1])

                            # From left to right once
                        if Utils.Right_circ_dir(new_coord, supp_i) == True:
                            if new_coord[-1] not in Pos_dict[read.query_name]:
                                Pos_dict[read.query_name].append(new_coord[-1])
                            print("full2.1")
                            # From right to left once
                        if Utils.Left_circ_dir(new_coord, supp_i) == True:
                            if new_coord[-1] not in Pos_dict[read.query_name]:
                                Pos_dict[read.query_name].append(new_coord[-1])
                        print("TEST2", Pos_dict[read.query_name])

                    print("-----------------")
        return Pos_dict

    def Region(self):
        with open(self.region) as f:
            for line in f:
                region = line.strip().split()
                chr = region[0]
                start = region[1]
                end = region[2]
                pos_dict = self.Read_pos(self.bamfile,self.MapQ,str(chr),int(start), int(end))
                #print(pos_dict)
                #print("-------------------------------")

if __name__ == '__main__':
    print("lol")