import msprime
import numpy as np
import random
import json
import math

def all_ibd_segments(ts):
    n = ts.num_samples
    trees_iter = ts.trees()
    tree = next(trees_iter)
    last_mrca_m = np.zeros((n,n))
    for i in range(n):
        for j in range(i,n):
            last_mrca_m[i][j] = tree.mrca(i,j)
    last_left_m = np.zeros((n,n))
    segment_lengths_m = [[[]for x in range(n)]for y in range(n)]
    for tree in trees_iter:
        for i in range(n):
            for j in range(i,n):
                mrca = tree.mrca(i,j)
                last_mrca = last_mrca_m[i][j]
                if mrca!= last_mrca:
                    left = tree.interval[0]
                    last_left = last_left_m[i][j]
                    segment_lengths_m[i][j].append((left-last_left)/ts.sequence_length)
                    last_mrca_m[i][j] = mrca
                    last_left_m[i][j] = left
    for i in range(n):
        for j in range(i,n):
            segment_lengths_m[i][j].append((ts.sequence_length-last_left_m[i][j])/ts.sequence_length)
    return segment_lengths_m

def popn1_simple_data_morgan(N,L,m,length):
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=N)
    bb = m*length
    kk = m*100
    ts = msprime.sim_ancestry(
        samples={"A": 15}, 
        demography=demography, 
        recombination_rate = 1/length,
        sequence_length = bb
    )
    all = all_ibd_segments(ts)
    out = {'y':[],'number':[],'u':[],'v':[]}
    for i in range(len(L)-1):
        u = L[i]
        v = L[i+1]

        o = []
        for j in range(30):
            for k in range(j+1,30):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out['y'].append(sum(o)/len(o))
        out['number'].append(len(o))
        out['u'].append(u)
        out['v'].append(v)
    return out


