import msprime
import numpy as np


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


# This return the avg IBD sharing fraction between populations(so 3*3) for each bin.
# N and T are list of variables specify the parameters in population topology.
# Topology is fixed under this simulation.

def popn3_simple_data(N,T,L):
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=N[0])
    demography.add_population(name="B", initial_size=N[1])
    demography.add_population(name="C", initial_size=N[2])
    demography.add_population(name="AB", initial_size=N[3])
    demography.add_population(name="ABC", initial_size=N[4])
    demography.add_population_split(time=T[0], derived=["A", "B"], ancestral="AB")
    demography.add_population_split(time=T[1], derived=["AB", "C"], ancestral="ABC")
    ts = msprime.sim_ancestry(
        samples={"A": 5, "B": 5, "C": 5}, 
        demography=demography, 
        recombination_rate = 1e-6,
        sequence_length = 1e6
    )
    all = all_ibd_segments(ts)
    out = [np.zeros(3,3) for i in range(len(L)-1)]
    for i in range(len(L)-1):
        u = L[i]
        v = L[i+1]
        o = []
        for j in range(10):
            for k in range(j+1,10):
                o.append(sum([l for l in all[j][k] if u<l*100<v]))
                out[i][0][0] = sum(o)/len(o)
        o = []
        for j in range(10):
            for k in range(j+10,20):
                o.append(sum([l for l in all[j][k] if u<l*100<v]))
                out[i][0][1] = out[i][1][0] = sum(o)/len(o)
        o = []
        for j in range(10):
            for k in range(j+20,30):
                o.append(sum([l for l in all[j][k] if u<l*100<v]))
                out[i][0][2] = out[i][2][0] = sum(o)/len(o)
        o = []
        for j in range(10,20):
            for k in range(j+1,20):
                o.append(sum([l for l in all[j][k] if u<l*100<v]))
                out[i][1][1] = sum(o)/len(o)
        o = []
        for j in range(10,20):
            for k in range(j+10,30):
                o.append(sum([l for l in all[j][k] if u<l*100<v]))
                out[i][1][2] = out[i][2][1] = sum(o)/len(o)
        o = []
        for j in range(20,30):
            for k in range(j+1,30):
                o.append(sum([l for l in all[j][k] if u<l*100<v]))
                out[i][2][2] = sum(o)/len(o)
    return out



