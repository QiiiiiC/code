import msprime
import numpy as np
import random
import json

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
        o = []
        for j in range(10):
            for k in range(j+1,10):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out['y'].append(sum(o)/len(o))
        out['u'].append(u)
        out['v'].append(v)
        out['group'].append(1)
        out['N_obs'] += 1

        o = []
        for j in range(10):
            for k in range(j+10,20):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out['y'].append(sum(o)/len(o))
        out['u'].append(u)
        out['v'].append(v)
        out['group'].append(2)
        out['N_obs'] += 1


        o = []
        for j in range(10,20):
            for k in range(j+1,20):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out['y'].append(sum(o)/len(o))
        out['u'].append(u)
        out['v'].append(v)
        out['group'].append(3)
        out['N_obs'] += 1
        
    return out

N = [1400,2800,2500]
T = 20
L = [0.5,2,200]
out = popn2_simple_data_morgan(N,T,L,2,1e6)
with open('two_popn_expected_data.json', 'w') as json_file:
    json.dump(out, json_file)