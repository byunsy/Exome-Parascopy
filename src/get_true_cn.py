import sys
import math
import argparse
import pandas as pd
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, help="path to res.samples.bed from Parascopy")
parser.add_argument("-d", "--index", required=False, help="specific index of copy number intervals")
args = parser.parse_args()

fp = args.input
df = pd.read_csv(fp, sep="\t", skiprows=2)

for name in df['sample'].unique():
    cn_list = list(df.loc[df['sample']==name]['agCN'])
    mid = len(cn_list)//2
    if args.index:
        mid = int(args.index)
    print(name, cn_list[mid], cn_list)
