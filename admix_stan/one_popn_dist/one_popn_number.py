import msprime
import numpy as np
import random
import json
import seaborn
import matplotlib.pyplot as plt
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

def popn1_data(N,m,length,u):
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=N)
    bb = m*length
    kk = m*100
    ts = msprime.sim_ancestry(
        samples={"A": 20}, 
        demography=demography, 
        recombination_rate = 1/length,
        sequence_length = bb
    )
    all = all_ibd_segments(ts)
    out = {'N_obs':0,'y':[]}
    for j in range(40):
            for k in range(j+1,40):
                a = [l*kk for l in all[j][k] if u<l*kk]
                out['y'] += a
                out['N_obs'] += len(a)
    return out

def popn1_number(N,m,length,L):
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=N)
    bb = m*length
    kk = m*100
    ts = msprime.sim_ancestry(
        samples={"A": 10}, 
        demography=demography, 
        recombination_rate = 1/length,
        sequence_length = bb
    )
    all = all_ibd_segments(ts)
    out = {'N_obs':0,'y':[],'u':[],'v':[]}
    for i in range(len(L)-1):
        u = L[i]
        v = L[i+1]
        for j in range(20):
            for k in range(j+1,20):
                a = [l*kk for l in all[j][k] if u<l*kk<v]
                out['y'] += [len(a)]
                out['N_obs'] += 1
                out['u'] += [u]
                out['v'] += [v]
    return out




