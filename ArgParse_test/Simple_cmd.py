import Utils
import pysam as ps
import numpy as np
from collections import Counter

print("Imports the simple script")

class Simple_circ:
    #Command 1
    def __init__(self, bamfile,output,MapQ):
        self.bamfile = bamfile
        self.output = output
        self.MapQ = MapQ

    def Read_pos(self):
        read_file = Utils.SamBamParser(self.bamfile)
        Pos_dict = {}
        counter = 0

        for read in read_file.fetch():
            if counter == 10:
                return Pos_dict
                break

            if Utils.IS_SA(read, self.MapQ) == True:
                #primary pos
                prim_chr = read.reference_name
                prim_pos = [int(read.reference_start), int(read.reference_end)]

                #supplementary alignment
                Tag = read.get_tag("SA").split(',')

                # single supp alignment
                if len(Tag) == 6:
                    if int(Tag[4]) >= self.MapQ:
                        counter += 1
                        supp_pos = [int(Tag[1]),int(Tag[1])+Utils.CIGAR_len(Tag[3])]
                        Pos_dict[read.query_name] = [prim_chr,
                                                     prim_pos,
                                                     supp_pos]
                # Multiple supp alignments
                else:
                    counter += 1

                    MapQ_val = Tag[4::5] # possible mapq
                    Start_possible = Tag[1::5] # possible start coordinates

                    mult_coord = []

                    #for i in range for same index of the mapQ, start pos, cigarstring
                    for i in range(len(MapQ_val)):
                        if int(MapQ_val[i]) >= self.MapQ: # if the different mapQ == 60

                            #creates a list of list, with each sublist being the coordinate interval
                            mult_coord.append([int(Start_possible[i]),
                                               int(Start_possible[i])+Utils.CIGAR_len(Tag[3::5][i])])

                    #flatten the entire mult_coord list of list, into a single list
                    supp_pos =  list(np.array(mult_coord).flat)

                    Pos_dict[read.query_name] = [Tag[0], prim_pos, supp_pos]

        return Pos_dict

    def flatten(self,irreg_list):
        for i in irreg_list:
            if isinstance(i, (list, tuple)):
                for j in self.flatten(i):
                    yield j
            else:
                yield i

    def Circ_possible(self):
        Pos_dict = self.Read_pos()
        Circ_pos = {}

        #determine if the prim pos is upstream or downstream of supp pos
        for read_ID,coord_pos in Pos_dict.items():
            # primary read upstream of supp align
            if min(coord_pos[1]) <= min(coord_pos[2]) and max(coord_pos[1]) <= max(coord_pos[2]):
                dic_val = [coord_pos[0],
                           int(min(coord_pos[1])), # minimum prim coord (start coord)
                           coord_pos[2][1::2]]     # every second supp coord from index 1 (possible end coord)

                # flattens the irrecular list (dic_val) with structure [chr,coord,[coord1,coord2]]
                Circ_pos[read_ID] = list(self.flatten(dic_val))

            # primary read spanning entire circle, but weird supp is inside interval
            if min(coord_pos[1]) <= min(coord_pos[2]) and max(coord_pos[1]) >= max(coord_pos[2]):
                dic_val = [coord_pos[0],
                           int(min(coord_pos[1])),
                           int(max(coord_pos[1]))]

                Circ_pos[read_ID] = list(self.flatten(dic_val))

            # supplementary align upstream of prim align
            # i could use the strategy of keeping the multiple supp alignments for start coord
            if min(coord_pos[2]) <= min(coord_pos[1]) and max(coord_pos[2]) <= max(coord_pos[1]):
                dic_val = [coord_pos[0],
                           coord_pos[2][0::2], #every second supp coord from index 0 (possible start coord)
                           int(max(coord_pos[1]))]
                Circ_pos[read_ID] = list(self.flatten(dic_val))

            # supplementary read spanning entire circle, but weird prim is inside interval
            if min(coord_pos[2]) <= min(coord_pos[1]) and max(coord_pos[2]) >= max(coord_pos[1]):
                dic_val = [coord_pos[0],
                           int(min(coord_pos[2])),
                           int(max(coord_pos[2]))]
                Circ_pos[read_ID] = list(self.flatten(dic_val))

        return Circ_pos

    def reduce_end_coords(self):
        """ Reduce the number of end coordinates for those reads with several
        supplementary alignments by comparing the position to the other end positions"""

        Coord = self.Circ_possible()
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

if __name__ == '__main__':

    Circ_class=Simple_circ("/isdata/common/wql443/NanoCircle/BC02.aln_hg19.bam","lol",60)
    dict = Circ_class.Circ_possible()

    list = []
    print([{k: v} for (k, v) in dict.items()])

    for k,v in dict.items():
        if
        print(k,v)




"""
    dict_test = {'key1': [19010495, 19014658, 19014064],
                 'key2': [19010495, 19014658],
                 'key3': [19010495, 19014658],
                 'key4': [19010502, 19014641]}

    Circ_class.reduce_end_coords(dict_test)"""

