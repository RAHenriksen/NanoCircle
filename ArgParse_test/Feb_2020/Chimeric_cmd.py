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

print("Imports the chimeric script")
class Chimeric_circ:
    #Command 2
    def __init__(self, input):
        self.input = input

    def multiply(self,val=1):
        return self.input*val

