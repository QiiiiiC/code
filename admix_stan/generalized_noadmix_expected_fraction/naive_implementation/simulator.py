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

def popn3_simple_data_morgan(N,T,L,m,length,n):
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=N[0])
    demography.add_population(name="B", initial_size=N[1])
    demography.add_population(name="C", initial_size=N[2])
    demography.add_population(name="BC", initial_size=N[3])
    demography.add_population(name="ABC", initial_size=N[4])
    demography.add_population_split(time=T[0], derived=["B", "C"], ancestral="BC")
    demography.add_population_split(time=T[1], derived=["A", "BC"], ancestral="ABC")
    
    bb = m*length
    kk = m*100
    ts = msprime.sim_ancestry(
        samples={"A": n, "B": n, "C": n}, 
        demography=demography, 
        recombination_rate = 1/length,
        sequence_length = bb
    )
    all = all_ibd_segments(ts)
    out = {'y':[],'u':[],'v':[],'group':[]}

    for i in range(len(L)-1):
        u = L[i]
        v = L[i+1]

        o = []
        for j in range(n*2):
            for k in range(j+1,n*2):
                a = [l*kk for l in all[j][k] if u<l*kk<v]
                out['y']+=a
                out['u']+=[u]*len(a)
                out['v']+=[v]*len(a)
                out['group']+=[1]*len(a)

        o = []
        for j in range(n*2):
            for k in range(j+n*2,n*4):
                a = [l*kk for l in all[j][k] if u<l*kk<v]
                out['y']+=a
                out['u']+=[u]*len(a)
                out['v']+=[v]*len(a)
                out['group']+=[2]*len(a)


        o = []
        for j in range(n*2):
            for k in range(j+n*4,n*6):
                a = [l*kk for l in all[j][k] if u<l*kk<v]
                out['y']+=a
                out['u']+=[u]*len(a)
                out['v']+=[v]*len(a)
                out['group']+=[3]*len(a)

        o = []
        for j in range(n*2,n*4):
            for k in range(j+1,n*4):
                a = [l*kk for l in all[j][k] if u<l*kk<v]
                out['y']+=a
                out['u']+=[u]*len(a)
                out['v']+=[v]*len(a)
                out['group']+=[4]*len(a)


        o = []
        for j in range(n*2,n*4):
            for k in range(j+n*2,n*6):
                a = [l*kk for l in all[j][k] if u<l*kk<v]
                out['y']+=a
                out['u']+=[u]*len(a)
                out['v']+=[v]*len(a)
                out['group']+=[5]*len(a)

        o = []
        for j in range(n*4,n*6):
            for k in range(j+1,n*6):
                a = [l*kk for l in all[j][k] if u<l*kk<v]
                out['y']+=a
                out['u']+=[u]*len(a)
                out['v']+=[v]*len(a)
                out['group']+=[6]*len(a)

    return out