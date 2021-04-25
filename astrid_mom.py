# ASTRID median of means
import argparse
import asterid as ad
import numpy as np
import math
from astrid_base import *

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def median_of_means(seq, block):
    if block >= len(seq):
        return np.mean(seq)
    np.random.shuffle(seq)
    means = [np.mean(c) for c in chunks(seq, (len(seq) // block))]
    return np.median(means)

def r_median_of_means(seq, block = 20, tries = 20):
    return np.mean([median_of_means(seq, block) for i in range(tries)])

def parse_blocksize(seqlen, src):
    if src[0] == "c":
        return int(src[1:])
    if src[0] == "s":
        return math.floor(math.sqrt(seqlen))

parser = argparse.ArgumentParser(description='run ASTRID')
parser.add_argument("-i", "--input", type=str,
                    help="Input tree list file", required=True)
parser.add_argument("-o", "--output", type=str,
                    help="Output tree file", required=True)
parser.add_argument("-b", "--block", type=str, default="c20")
parser.add_argument("-m", "--methods", type=str, help="Methods to run", required=True)
args = parser.parse_args()

trees = open(args.input, "r").readlines()
methods = args.methods

ts = ad.get_ts(trees)
Ds = all_matrices(ts, trees)
D = matrix_elementwise(ts, Ds, lambda coll: r_median_of_means(coll, parse_blocksize(len(coll), args.block)))
t = run_iterations(ts, D, methods)
if args.output == "-":
    print(t)
else:
    with open(args.output, "w+") as fh:
        fh.write(t)