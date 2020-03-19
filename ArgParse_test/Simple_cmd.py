import Utils
import pysam as ps
print("Imports the simple script")

class Simple_circ:
    #Command 1
    def __init__(self, bamfile,output,MapQ):
        self.bamfile = bamfile
        self.output = output
        self.MapQ = MapQ

    def Reads(self):
        read_file = Utils.SamBamParser(self.bamfile)
        read_dict = {}
        counter=0
        for read in read_file.fetch():
            if counter == 10:
                return read_dict
                break

            #print("READ NAME",read.query_name)
            # checks if soft-clipped
            if Utils.IS_SA(read, self.MapQ) == True:

                pos = ps.AlignedSegment.get_reference_positions(read)

                Tag = read.get_tag("SA").split(',')
                if len(Tag) == 6:
                    if int(Tag[4]) >= self.MapQ:
                        counter += 1
                        print("primary", read.query_name, read.reference_name, pos[0], pos[-1])
                        print("supp", read.get_tag("SA").split(','),Utils.CIGAR_len(Tag[3]))
                        read_dict[read.query_name] = [read.reference_name, pos[0], pos[-1]]
                else:
                    counter += 1
                    MapQ_val = Tag[4::5] # possible mapq
                    Start_possible = Tag[1::5] # possible start coordinates

                    start_coord = []
                    end_coord = []
                    #print("secondayr", read.query_name, read.reference_name,Tag,pos[0], pos[-1])
                    #for i in range for same index
                    for i in range(len(MapQ_val)):
                        if int(MapQ_val[i]) >= self.MapQ: # if the different mapQ == 60
                            start_coord.append(Start_possible[i]) # takes the coresponding val
                            end_coord.append(int(Start_possible[i])+Utils.CIGAR_len(Tag[3::5][i])) # add the cigar str

                    read_dict[read.query_name] = [read.reference_name, start_coord, end_coord]

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
                    #for i in range for same index
                    for i in range(len(MapQ_val)):
                        if int(MapQ_val[i]) >= self.MapQ: # if the different mapQ == 60
                            mult_coord.append([int(Start_possible[i]),
                                               int(Start_possible[i])+Utils.CIGAR_len(Tag[3::5][i])])

                    Pos_dict[read.query_name] = [Tag[0], prim_pos, mult_coord]

        return Pos_dict

print("---------------")
if __name__ == '__main__':

    Circ_class=Simple_circ("/isdata/common/wql443/NanoCircle/BC02.aln_hg19.bam","lol",60)
    dict = Circ_class.Read_pos()
    print(dict)
    #dict = Circ_class.Reads()

#    for k,v in dict.items():
#        print(k,v)

"""
59860e12-1576-4e0c-b190-482547539045 ['chr10', ['19013902', '19013547'], [19014590, 19014064]] #multiple start and end
b4f70162-08e4-40dc-aa10-73eec0961ef4 ['chr10', 19010495, 19011733] # single start and single end
22ec4a10-9c10-40e3-9ca1-7abb521dbb59 ['chr10', 19010502, 19014290]
4691e064-2290-4067-8082-3c7b285dcad1 ['chr10', ['19013847', '19012994', '19010507'], [19014666, 19013518, 19010927]]
8acd93e4-566e-4f1f-86ff-2f9ea2b7bb39 ['chr10', 19010805, 19013008]
c0ca6bbd-8d5a-4a5d-87fc-52c46854347b ['chr10', ['19010499', '19013814'], [19011604, 19014672]]
ed79dc24-b59a-44d0-b1d1-9760ed3c2c1c ['chr10', 19013260, 19014018]
f8d8742c-0b62-43c6-b33a-d223b9e595b3 ['chr10', ['37375766'], [37377130]] #multiple start and end, but only one with q above 60
96f5cb5f-5cb8-4cb1-8625-a723b774e91c ['chr10', ['76314667', '76311112', '76312785'], [76315990, 76312224, 76312963]]
8908600e-54a0-409f-9dcc-0c99d9869fc5 ['chr11', 14806147, 14807426]
"""