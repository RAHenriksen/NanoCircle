#!/usr/bin/env python
import edlib as ed
import pandas as pd
import seaborn as sns
import numpy as np
from Bio import SeqIO
import re
import os
import pysam as ps
import matplotlib.pyplot as plt
import pylab as plt2
import time
from collections import defaultdict
from matplotlib.patches import Rectangle
import pybedtools
#from pybedtools import BedTool

def read_count(bamfile):
    """ Counts all the reads in a bamfile"""
    return bamfile.count()

def Full_CIGAR(cigar_str):
    """Returns the actual length of the sequence based on the CIGAR int values"""
    count = 0
    for i in re.findall(r'(\d+)([A-Z]{1})', cigar_str):
        if i[-1] == "M" or i[-1] == "I" or i[-1] == "H" or i[-1] =="S":
            count += int(i[0])
        elif i[-1] == "D":
            continue
        #    count -= int(i[0])
    return count

def CIGAR_primary(cigar_str):
    """Returns the actual length of the sequence based on the CIGAR int values"""
    count = 0
    for i in re.findall(r'(\d+)([A-Z]{1})', cigar_str):
        if i[-1] == "M" or i[-1] == "I" or i[-1] == "H":
            count += int(i[0])
        elif i[-1] == "D" or i[-1] =="S":
            continue
        #    count -= int(i[0])
    return count

def CIGAR_pos(cig_str,char_str):
    """Returns the position in the Cigar string for the given character string"""
    return [pos for pos, char in enumerate(cig_str) if char == char_str]

def IS_SA(read,mapQ):
    """returns true if the read in the bamfile is soft-clipped"""
    if read.cigar[0][0] == 4:
        if read.has_tag("SA"):
            if read.mapping_quality == mapQ:
                return True

def IS_supp(read,mapQ):
    """returns true if the read in the bamfile is supplementary"""
    if read.is_supplementary == True:
        if read.mapping_quality == mapQ:
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

def Count_suppl(bamfile):
    "Count the number of supplementary alignments on either strand using boolean values for strand_str"
    count = 0
    for read in bamfile:
        if read.is_supplementary == True:
            if read.mapping_quality == 60:
                count += 1
    return count

def Count_Soft(bamfile):
    "Count the number of Soft-clipped reads with a supplementary alignment"
    count = 0
    for read in bamfile:
        if read.cigar[0][0] == 4:
            if read.has_tag("SA"):
                if read.mapping_quality == 60:
                    count += 1
    return count

def Coord_sort(input_bam,input_name):
    "checks if the input bame is sorted by coordinate, if not it sort it by coordinate"
    #PERHAPS ALSO ADD A INDEX FILE
    if 'HD' in input_bam.header: #HD first line if present
        if input_bam.header['HD']['SO'] != 'coordinate': #if the sorting order is not by coordinates, sort by coordinates
            print("not sorted")
            input_bam.close() #i guess its closed as to not save it in the memory

            # it makes sense to sort them after coordinates
            ps.sort("-o","coordsort_%s" % input_name, input_name)
            sorted_bam = ps.AlignmentFile("coordsort_%s.bam" % input_name[:-4])
        else:
            print("sorted")
            #input_bam.close()
            return input_bam
    return sorted_bam

def Support_circ(Dict):
    read_count = 0
    for k, v in sorted(Dict.items()):
        read_count += 1
        print(k, v)
    print("Supporting_reads",read_count)
    #return [print(k,v) for k,v in sorted(Dict.items())]

def Count_reads(bamfile,reg,start,end):
    "Count the number of Soft-clipped reads with a supplementary alignment"
    s_count = 0
    SA_count = 0
    NO_SA_count = 0
    supp_count = 0
    for read in bamfile.fetch(reg,start,end,multiple_iterators=True):
        if read.is_supplementary == True:
            supp_count += 1
        #OVERVEJ LIGE OM MAN IKKE FAKTISK MISTER MANGE VED IKKE AT BRUGE cigar[-1][0]
        if read.cigar[0][0] == 4:
            s_count += 1
            if read.has_tag("SA"):
                SA_count += 1
            if not read.has_tag("SA"):
                NO_SA_count +=1
    print("total number of reads in the file",bamfile.count(reg,start,end,read_callback='nofilter'))
    print("total number of soft-clipped",s_count)
    print("Number of soft-clipped reads with SA",SA_count)
    print("Number of soft-clipped reads wo. SA",NO_SA_count)
    print("supp",supp_count)

#Count_reads(bamfile,None,None,None)

def Read_length(bamfile,reg,start,end,filename,read_type):
    read_list = []
    for read in bamfile.fetch(reg, start, end, multiple_iterators=True):
        if read_type == "Soft":
            #full length of the entire soft-clipped
            if read.cigar[0][0] == 4:
                #print("length",read.query_length)
                #print("CIGARLen",Full_CIGAR(str(read.cigarstring)))
                #print("str",read.cigarstring)
                read_list.append(read.query_length)
        elif read_type == "SA":
            #length without the soft clipped part, so primary part
            if read.cigar[0][0] == 4 and read.has_tag("SA"):
                Tag = read.get_tag("SA").split(';')[:-1]
                for Tagelem in Tag:
                    Column_list = Tagelem.split(',')
                    print(Column_list)
                    print("cigar",CIGAR_primary(Column_list[3]))
                    read_list.append(CIGAR_primary(Column_list[3]))
        elif read_type == "supp":
            if read.is_supplementary == True:
                print("cigar",read.cigarstring)
                print("len",read.query_length)
                print("test",CIGAR_primary(str(read.cigarstring)))
                read_list.append(read.query_length)
    with open(filename, 'w') as f:
        for item in read_list:
            f.write("%s\n" % item)

sample = "BC08"
bamfile = ps.AlignmentFile("/isdata/common/wql443/Thesis_2019/Analysis/2.8samples/%s/%s.aln_hg19.bam" % (sample,sample),"rb")
### DER ER NOGLE PROBLEMER MED CIGAR LENGTH
Read_length(bamfile,None,None,None,"Read_length/2.8/%s_Soft.txt" % sample,"Soft")
#Read_length(bamfile,None,None,None,"Read_length/%s_SA.txt" % sample,"SA")
#Read_length(bamfile,None,None,None,"Read_length/%s_Supp.txt" % sample,"supp")

#print(CIGAR_len("496S345M4D"))
#print(CIGAR_len("159S119M3I2S"))
#print(CIGAR_len("1346S68M"))
#print(CIGAR_len("466S250M11D671S"))