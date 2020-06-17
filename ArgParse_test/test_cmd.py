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
import operator
from collections import Counter

start_dict1 = {
    'Key1': {'start': [243928620], 'end': [243938319]},
    'Key2': {'start': [243928620,243928625], 'end': [243938319]},
    'Key3': {'start': [243928620], 'end': [243931757, 243938319,243938354]},
    'Key4': {'start': [243928620], 'end': [243938319, 243938323]},
    'Key5': {'start': [243928620,243928634], 'end': [243938316]},
    'Key6': {'start': [243928620], 'end': [243938319]},
    'Key7': {'start': [243928634], 'end': [243938317]},
    'Key8': {'start': [243928620], 'end': [243938329,243938387]}}

start_dict2 = {
    'Key1': {'start': [243928620], 'end': [243938319]}}

start_dict3 = {
    'Key1': {'start': [243928620], 'end': [243938319,243938313]}}

start_dict4 = {'3c0df367-ce16-4ea2-b8fd-b3d40c1a630f': {'start': [2397, 2765, 3228], 'end': [4391, 3735, 3966]},
               '6d3e80e0-aa26-4ea3-9920-d2a5d7c8f30b': {'start': [10381], 'end': [12601]},
               '3a0bf9dd-d51b-4eae-b2b2-9a0509fa5d2e': {'start': [6401, 10513], 'end': [14622]},
               '4d917672-6e42-4d20-b3eb-1eac25cb8e4c': {'start': [11643, 11929, 3179], 'end': [14935, 13419]}}

a = list(start_dict3['Key1'].values())

for k,v in start_dict1.items():
    if any(i == 243928620 for i in v['start']):
        print("TETET",k,v)

exit()


def Reduce_coord(dict,pos):
    pos_point = []
    for k, v in dict.items():
        pos_point.extend(v[pos])
    #find the most common end point

    most_common = Counter(pos_point).most_common()

    #the frequency value for each individual coordinates
    freq_val = [item[-1] for item in most_common]

    if any(x > 1 for x in freq_val):
        #Adjust the end points if the most common end point is found
        for k, v in dict.items():
            if most_common[0][0] in v[pos]:
                dict[k][pos] = [most_common[0][0]]
        return dict

    else:
        return dict

def most_frequent(List,pos):
    count1 = Counter(List)
    dic_val = list(count1.values())
    max_val = max(count1.values())
    occurence = dic_val.count(max_val)

    #single most occuring coordinate
    if occurence == 1:
        #print("COMMMON VALUE", count1.most_common(1)[0][0])
        occ = 1
        return occ, count1.most_common(1)[0][0]

    #creates the interval immediately
    else:
        most_common_no = list(map(operator.itemgetter(1), count1.most_common()))
        most_common_idx = [i for i, x in enumerate(most_common_no) if x == max(most_common_no)]


        #if there are several equal common values
        if max(most_common_no) > 1:
            #print(max(most_common_no))
            #print(most_common_idx)
            #creating the most maximum value
            common_val = [count1.most_common()[i][0] for i in most_common_idx]
            occ = 1
            #print("COMMMON VALUE",common_val)
            return occ, max(common_val)

        #if there isnt a most occuring coordinate, the interval is created
        else:
            occ = 2
            #start
            if pos == 0:
                return occ, min(List)
            #end
            if pos == 1:
                return occ, max(List)



def Single_coord(start_dict):

    #reduce start and end coordinates to those which are most common
    temp_dict = Reduce_coord(start_dict, 'start')
    modify_dict = Reduce_coord(temp_dict, 'end')

    #creates a single list with all start and end coordinates
    start_list = []
    end_list = []
    for v in modify_dict.values():
        start_list.extend(v['start'])
        end_list.extend(v['end'])


    if start_list and end_list != []:
        occ1,start_freq = most_frequent(start_list,0)
        print("occurence,",occ1,"start",start_freq)
        occ2,end_freq = most_frequent(end_list,1)
        print("occurence,", occ2, "end", end_freq)

        if occ1 and occ2 == 1:
            reads = []
            chr = ['chr1']
            if start_freq < end_freq:
                new_val = chr + [start_freq,end_freq]
            else:
                print(end_freq,start_freq)
                new_val = chr + [end_freq,start_freq]

            for k,v in modify_dict.items():
                print("lol")
                if any(i == start_freq or i == end_freq for i in v['start']) or any(i == end_freq for i in v['end']):
                    reads.append(k)

            final_dict = {tuple(sorted(reads)): new_val}

            if len(list(final_dict.keys())[0]) != 1:
                type = 1
                return type, final_dict

            else:
                type = 2
                return type, final_dict

        else:
            type = 2
            chr = ['chr']
            new_val = chr + [start_freq, end_freq]
            final_dict = {tuple(sorted(modify_dict.keys())): new_val}
            return type,final_dict

type,final_dict = Single_coord(start_dict1)
print(type,final_dict)
print("_-------------")
type,final_dict = Single_coord(start_dict2)
print(type,final_dict)
print("_-------------")
type,final_dict = Single_coord(start_dict3)
print(type,final_dict)
print("_-------------")
type,final_dict = Single_coord(start_dict4)
print(type,final_dict)
print("_-------------")

a = [2,3,45]
b = [3,4,5]

for i in a and b:
    print(i)