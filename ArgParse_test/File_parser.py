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

class Input_Parser:
    """
       Reads and open the tab seperated file, defined by the -i/--input parameter
    """
    def __init__(self,input):
        self.input = input

    def parse_data(self):
        # Parse the input filename
        data = []
        with open(self.input) as f:
            for line in f:
                # Reads each line to create the circular regions
                line_value = line.strip().split()
                data.append(line_value)
        return data

class Output_Parser:
    """
    Saves the data from the methods, in a format defined by the -o/--output parameter
    """
    def __init__(self,input,name):
        self.input = input
        self.name = name

    def dataframe(self):
        df = pd.DataFrame(self.input, columns=['chr', 'start', 'end'])
        bed_df = bt.BedTool.from_dataframe(df)
        bed_df.saveas(self.name)

