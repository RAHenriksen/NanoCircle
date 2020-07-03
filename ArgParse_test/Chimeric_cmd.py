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
        Prim_dict = {}
        Supp_dict = {}

        Read_File = Utils.SamBamParser(bamfile)

        for read in Read_File.fetch(reg, start, end, multiple_iterators=True):
            print(read.query_name)

    def Region(self):
        Simple_circ = pd.DataFrame()
        with open(self.region) as f:
            for line in f:
                region = line.strip().split()
                chr = region[0]
                start = region[1]
                end = region[2]
                self.Read_position(self.bamfile,self.MapQ,str(chr),int(start), int(end))