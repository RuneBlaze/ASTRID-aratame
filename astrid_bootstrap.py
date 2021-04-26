import argparse
import asterid as ad
from astrid_base import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run ASTRID')
    parser.add_argument("-i", "--input", type=str,
                        help="Input tree list file", required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Output tree file", required=True)
    parser.add_argument("-m", "--methods", type=str, help="Methods to run", required=True)
    parser.add_argument("-r", "--rate", type=float, help="deletion rate", default = 0.01)
    parser.add_argument("-b", "--bootstrap", type=int, help="bootstrap iterations", default=10)
    args = parser.parse_args()

    trees = open(args.input, "r").readlines()
    methods = args.methods

    ts = ad.get_ts(trees)
    res = []
    for i in range(args.bootstrap + 1):
        D = ad.mk_distance_matrix(ts, trees)
        if i > 0:
            random_deletion(ts, D, args.rate)
        t = run_iterations(ts, D, methods)
        res.append(t)
    t = consensus(res)
    if args.output == "-":
        print(t)
    else:
        with open(args.output, "w+") as fh:
            fh.write(t)