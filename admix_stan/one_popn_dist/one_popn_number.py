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
        samples={"A": 10}, 
        demography=demography, 
        recombination_rate = 1/length,
        sequence_length = bb
    )
    all = all_ibd_segments(ts)
    out = {'N_obs':0,'y':[]}
    for j in range(20):
            for k in range(j+1,20):
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

N_list = [500,1000,2000,3000,4000,5000]
L = [0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.25,1.5,2,3,5,10,200]
for i in range(len(N_list)):
    N = N_list[i]
    out = popn1_number(N,2,1e7,L)
    out['m'] = 200
    with open(f'one_popn_dist/one_popn_number_data/one_popn_number_{N_list[i]}.json', 'w') as json_file:
            json.dump(out, json_file)


