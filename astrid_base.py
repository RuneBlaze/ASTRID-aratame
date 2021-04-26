import asterid as ad
import asterid as astrid
import numpy as np
from random import randrange

def run_iterations(ts, D, methods):
    fns = {
        "u": lambda ts, D: astrid.upgma_star(ts, D) + ";",
        "f": lambda ts, D: astrid.fastme_balme(ts, D, 0, 0),
        "n": lambda ts, D: astrid.fastme_balme(ts, D, 1, 0),
        "s": lambda ts, D: astrid.fastme_balme(ts, D, 1, 1),
        "j": lambda ts, D: astrid.fastme_nj(ts, D, 0, 0),
        "N": lambda ts, D: astrid.fastme_nj(ts, D, 1, 0),
        "S": lambda ts, D: astrid.fastme_nj(ts, D, 1, 1),
    }
    t = None
    for m in methods:
        if m not in fns:
            raise f"{m} not found as a method!"
        f = fns[m]
        t = f(ts, D)
        D.fill_in_transient(ts, t)
    return t

def reorder(x, y):
    if x > y:
        return y, x
    return x, y

def random_deletion(ts, D, rate):
    n = len(ts)
    totalmissing = round((n * (n - 1) / 2)*rate)
    
    cnt = 0
    while cnt < totalmissing:
        i = randrange(n)
        j = randrange(n)
        i, j = reorder(i, j)
        if D.has(i, j):
            D.setmask((i, j), 0)
            cnt += 1

def all_matrices(ts, trees):
    return [ad.DistanceMatrix(ts, t) for t in trees]

def taxon_pairs(ts):
    for i in range(len(ts)):
        for j in range(i, len(ts)):
            yield i, j

def consensus(trees, minfreq=0.5):
    import dendropy
    res = dendropy.TreeList()
    for treenewick in trees:
        res.read(data=treenewick, schema="newick", rooting='force-unrooted')
    # print(trees)
    con = res.consensus(min_freq = minfreq)
    con.is_rooted = False
    return con.as_string(schema="newick")

def matrix_elementwise(ts, Ds, f):
    R = ad.DistanceMatrix(ts)
    for i, j in taxon_pairs(ts):
        if i == j:
            R[i, j] = 0
            R.setmask((i, j), len(Ds))
            continue
        l = []
        for D in Ds:
            if D.has(i, j):
                l.append(D[i, j])
        if l:
            R[i, j] = f(l)
            R.setmask((i, j), len(l))
        else:
            R[i, j] = 0
            R.setmask((i, j), 0)
    return R