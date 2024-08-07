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

def popn1_simple_data_morgan(N,T,L,B,m,length):
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
    out = {'y':[],'sigma':[],'u':[],'v':[]}
    for i in range(len(L)-1):
        u = L[i]
        v = L[i+1]

        o = []
        for j in range(20):
            for k in range(j+1,20):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out['y'].append(sum(o)/len(o))
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        out['sigma'].append(math.sqrt(sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo)))
        out['u'].append(u)
        out['v'].append(v)
    return out

L = [0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.25,1.5,2,3,5,10]
N = 2500
B = 3000
output = popn1_simple_data_morgan(N,10,L,B,2,1e7)
output['N_obs'] = len(L)-1
with open('one_popn_simulated_data.json', 'w') as json_file:
    json.dump(output, json_file)

