import Utils
import pysam as ps
import numpy as np
from collections import Counter
from itertools import groupby
import pandas as pd
import numpy as np
import pybedtools as bt
import operator


print("Imports the chimeric script")
class Chimeric_circ:
    #Command 2
    def __init__(self,region,bamfile,output,MapQ):
        self.region = region
        self.bamfile = bamfile
        self.output = output
        self.MapQ = MapQ

    def Read_position(self,bamfile,MapQ,reg,start,end):
        """Creates a dictionary of the primary alignments for the soft-clipped reads"""
        Pos_dict = {}

        Read_File = Utils.SamBamParser(bamfile)

        Total_coord = []
        Total_chr = []
        Total_overlap = []

        for read in Read_File.fetch(reg, start, end, multiple_iterators=True):

            prim_pos = [int(read.reference_start) + 1, int(read.reference_end)]

            #creates a dict of dict, with key = read and its value the regions to which it aligns also found with
            #keys
            Pos_dict[read.query_name] = {'region1':[reg,[prim_pos[0]],[prim_pos[-1]]]}

            Total_chr.extend((reg, reg))
            Total_coord.extend((prim_pos[0], prim_pos[1]))

            #splits all alignment into individual alignments, opposite to the simple_cmd since they
                # align to different chrom
            Tag = read.get_tag("SA").split(';')[:-1]
            cigar_len = []

            i=1
            for Tagelem in Tag:
                # splitting the individual aligment info up
                Align_info = Tagelem.split(',')
                if int(Align_info[4]) >= MapQ:
                    Chrom = Align_info[0]
                    Length = Utils.CIGAR_len(Align_info[3])
                    cigar_len.append(Length)

                    # 1 - based positions
                    pos_start = int(Align_info[1]) - 1
                    pos_end = int(Align_info[1]) + Length - 1
                    supp_pos = [pos_start,pos_end]
                    #print(Chrom, pos_start, pos_end)
                    # the overlaps between coordinates for grouping
                    overlap = sum(cigar_len) * 4
                    Total_overlap.append(overlap)

                    # if the supp align is in between the circular input region then it aligns back to the primary region
                    if Chrom == reg and start - overlap <= pos_start <= end + overlap and start - overlap <= pos_end <= end + overlap:
                        #print("YOYLO",prim_pos, pos_start, pos_end, overlap, Total_overlap)

                        # From left to right across breakpoint
                        if Utils.Right_dir(prim_pos, supp_pos) == True:
                            #appends to the end coordinate
                            Pos_dict[read.query_name]['region1'][2].append(supp_pos[-1])

                        # From right to left
                        elif Utils.Left_dir(prim_pos, supp_pos) == True:
                            #start coordinate
                            Pos_dict[read.query_name]['region1'][1].append(supp_pos[0])

                        # primary read covers entire circle, and supp is between prim, but use both anyways due to
                        # perhaps having repeat alignment. Which is helpful for most common
                        elif Utils.Circ_once(prim_pos, supp_pos) == True:
                            Pos_dict[read.query_name]['region1'][1].append(supp_pos[0])
                            Pos_dict[read.query_name]['region1'][2].append(supp_pos[-1])

                    # add supp which arent in the primary align chromosome region
                    else:
                        # if the supp chromosome is the same as already in a region, check the coordinates
                        reg_chr = Pos_dict[read.query_name]['region{0}'.format(i)][0]
                        reg_start = Pos_dict[read.query_name]['region{0}'.format(i)][1]
                        reg_end = Pos_dict[read.query_name]['region{0}'.format(i)][2]
                        # checking the coordinates, if an overlap exist, the coordinates are appended to the existing coordinates
                            # depending on the direction
                        if reg_chr == Chrom and min(reg_start) - overlap <= pos_start <= max(reg_end) + overlap and \
                                min(reg_start) - overlap <= pos_end <= max(reg_end) + overlap:
                            if Utils.Right_dir(prim_pos, supp_pos) == True:
                                # appends to the end coordinate
                                Pos_dict[read.query_name]['region{0}'.format(i)][2].append(supp_pos[-1])

                            # From right to left
                            elif Utils.Left_dir(prim_pos, supp_pos) == True:
                                # start coordinate
                                Pos_dict[read.query_name]['region{0}'.format(i)][1].append(supp_pos[0])

                            # primary read covers entire circle, and supp is between prim, but use both anyways due to
                            # perhaps having repeat alignment. Which is helpful for most common
                            if Utils.Circ_once(prim_pos, supp_pos) == True:
                                Pos_dict[read.query_name]['region{0}'.format(i)][1].append(supp_pos[0])
                                Pos_dict[read.query_name]['region{0}'.format(i)][2].append(supp_pos[-1])

                        # Otherwise the supp aligns to another region
                        else:
                            i += 1
                            Pos_dict[read.query_name]['region{0}'.format(i)] = [Chrom, [pos_start], [pos_end]]

        return Pos_dict

    def chr_coord_sort(self, chrlist, coordlist):
        """ Sort a list of chromosomes and their coordinates using the index of the numerically sorted coordinates"""
        coord_idx = np.argsort(coordlist)
        Coord_sort = [coordlist[i] for i in coord_idx]
        chr_sort = [chrlist[i] for i in coord_idx]

        return Coord_sort, chr_sort

    def Grouping_chr(self, Chr_sort, Group_size):
        """ Groups the chromosomes in to match the grouped coordinates """
        Grouped_chr = []
        Step = 0

        for i in Group_size:
            Grouped_chr.append(Chr_sort[Step:Step + i])
            Step += i

        return Grouped_chr

    def Grouping(self,Coord_list, overlap_bp):
        """ Groups a given list, for which elements within the overlap range is grouped together"""

        First = [[Coord_list[0]]]

        for i in Coord_list[1:]:

            # sets the previous and current elem as the first elem
            prev = First[-1]
            current = prev[-1]
            dist = abs(i - current)

            # if dist between is below overlapping bp, then it it belong to the same group
            if dist < overlap_bp:
                prev.append(i)
            else:
                First.append([i])

        return First

    def Circle_Fragment(self, dict):
        """Creates list with the full region of start and end coordinates (all_start) and list of
        equal length with the corresponding chromosomes"""
        chr_start = []
        all_start = []

        chr_end = []
        all_end = []

        for k1, v1 in dict.items():
            for v2 in dict[k1].values():
                #Creates a single list with the full regions with start coordinates and end coordinates
                chr_start.extend(len(v2[1]) * [v2[0]])
                all_start.extend(v2[1])

                chr_end.extend(len(v2[2]) * [v2[0]])
                all_end.extend(v2[2])

        return chr_start,all_start,chr_end,all_end

    def Common_coord(self,chr_list,coord_list,overlap,pos):
        """ Identifies the most common coordinates for a given fragment for all reads,
        identified by grouping potential coordinates"""

        #Still some issues if the cases with multiple different chromosomes which have coordinates overlapping are
            # clustered together

        Sort_coord, Sort_chr = self.chr_coord_sort(chr_list, coord_list)
        Grouped_coord = self.Grouping(Sort_coord, overlap)
        Group_size = [len(i) for i in Grouped_coord]
        Grouped_start = self.Grouping_chr(Sort_chr, Group_size)


        final_reg = []
        for i in range(len(Group_size)):
            #if there is only one coordinate
            if Group_size[i] == 1:
                final_reg.extend([Grouped_start[i][0],Grouped_coord[i][0]])

            # multiple coordinates
            elif Group_size[i] >= 1:
                occurence_count = Counter(Grouped_coord[i])

                # select the most frequent coordinate
                if list(occurence_count.values())[0] > 1:
                    final_reg.extend([Grouped_start[i][0], occurence_count.most_common(1)[0][0]])

                else:
                    #if the pos is start then we take the minimum
                    if pos == 0:
                        final_reg.extend([Grouped_start[i][0], min(Grouped_coord[i])])

                    # if the pos is max then we take the maximum
                    if pos == 1:
                        final_reg.extend([Grouped_start[i][0], max(Grouped_coord[i])])


        return final_reg

    def Circle_dict(self,dict,Start_chr,Start_coord,End_chr,End_coord):
        """ Creates dictionary for each circule with all reads, and the most common coordinates"""
        start = self.Common_coord(Start_chr, Start_coord, 10000, 0)
        end = self.Common_coord(End_chr, End_coord, 10000, 1)

        d = {}
        d[tuple(dict.keys())] = []

        try:
            assert len(start) == len(end)
        except AssertionError as e:
            print("Clustering in Circle_dict didn't make equal clusters. Increasing clustering distance to 100000bp.")
            start = self.Common_coord(Start_chr, Start_coord, 100000, 0)
            end = self.Common_coord(End_chr, End_coord, 100000, 1)
             

        for i in range(len(start[::2])):

            d[tuple(dict.keys())].extend([start[0::2][i], start[1::2][i], end[1::2][i],end[1::2][i]-start[1::2][i]])


        return d

    def Circle_dataframe(self,circle_dict,ID):
        "Creates two dataframes, the first is the fragment information and the second general information"

        df_frag = pd.DataFrame.from_dict(circle_dict, orient='index')

        #The initial columns for the final dataframe, but the fragments can vary for individual chimeric eccDNA
        df_col = ["Chr", "Start", "End", "Length"]
        rep = int(len(list(circle_dict.values())[0]) / 4)
        df_frag.columns = [j + "_" + str(i) for i in range(1, rep + 1) for j in df_col]

        #Always the last four columns of the last dataframe
        df_general = pd.DataFrame(columns=['Total_len', 'No_fragments', 'No_reads', 'Config_ID'])

        # The sum of all the Length columns giving total length
        df_general['Total_len'] = df_frag.filter(regex=("Length.*")).sum(axis=1)
        df_general['No_fragments'] = rep
        df_general['No_reads'] = len(list(circle_dict.keys())[0])
        df_general['Config_ID'] = "Chimeric_%d" % (ID)

        return df_frag, df_general

    def Region(self):
        Complex_df = pd.DataFrame()
        Final_4_df = pd.DataFrame()
        i = 1
        with open(self.region) as f:
            for line in f:
                region = line.strip().split()
                chr = region[0]
                start = region[1]
                end = region[2]

                pos_dict = self.Read_position(self.bamfile,self.MapQ,str(chr),int(start), int(end))

                if pos_dict == {}:
                    pass
                else:
                    #print(pos_dict)
                    chr_start, all_start, chr_end, all_end = self.Circle_Fragment(pos_dict)
                    Final_dict = self.Circle_dict(pos_dict,chr_start, all_start, chr_end, all_end)
                    #print(Final_dict)
                    df_circ, last4_df = self.Circle_dataframe(Final_dict,i)

                    #without .fillna("") the columns with have deafault NaN and type == float64. instead of object
                    Complex_df = Complex_df.append(df_circ, sort=False)

                    Final_4_df = Final_4_df.append(last4_df)
                    i+=1

        # concatenates the last 4 columns to the regions
        resulting_df = pd.concat([Complex_df.reset_index(drop=True), Final_4_df.reset_index(drop=True)], axis=1)
        pd.set_option("display.max_rows", None, "display.max_columns", None)

        cols = list(resulting_df)[0:-4]
        del cols[0::4]

        #ensuring coordinates are int
        resulting_df[cols] = resulting_df[cols].astype(float).astype('Int64')


        Bed_File = bt.BedTool.from_dataframe(resulting_df)
        Bed_File.saveas(self.output)

if __name__ == '__main__':
    Read_class=Chimeric_circ("BC10_1000_cov.bed","/isdata/common/wql443/NanoCircle/ArgParse_test/temp_reads/Chimeric_reads.bam","bedtest.bed",60)
    Read_class.Region()
