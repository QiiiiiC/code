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

def popn2_simple_data_morgan(N,T,L,m,length):
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=N[0])
    demography.add_population(name="B", initial_size=N[1])
    demography.add_population(name="AB", initial_size=N[2])
    demography.add_population_split(time=T, derived=["A", "B"], ancestral="AB")
    
    bb = m*length
    kk = m*100
    ts = msprime.sim_ancestry(
        samples={"A": 5, "B": 5}, 
        demography=demography, 
        recombination_rate = 1/length,
        sequence_length = bb
    )
    all = all_ibd_segments(ts)
    out = {'N_obs':0,'y':[],'u':[],'v':[],'group':[]}

    for i in range(len(L)-1):
        u = L[i]
        v = L[i+1]
        for j in range(10):
            for k in range(j+1,10):
                a = [l*kk for l in all[j][k] if u<l*kk<v]
                out['N_obs'] += len(a)
                out['y'] += a
        
                out['u']+=([u]*len(a))
                out['v']+=([v]*len(a))
                out['group']+=([1]*len(a))

        for j in range(10):
            for k in range(j+10,20):
                a = [l*kk for l in all[j][k] if u<l*kk<v]
                out['N_obs'] += len(a)
                out['y'] += a
        
                out['u']+=([u]*len(a))
                out['v']+=([v]*len(a))
                out['group']+=([2]*len(a))      

        for j in range(10,20):
            for k in range(j+1,20):
                a = [l*kk for l in all[j][k] if u<l*kk<v]
                out['N_obs'] += len(a)
                out['y'] += a
        
                out['u']+=([u]*len(a))
                out['v']+=([v]*len(a))
                out['group']+=([3]*len(a))
        
    return out



