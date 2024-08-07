import msprime
import numpy as np
import random


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


# This returns the avg IBD sharing fraction between populations(so 3*3) for each bin, and the bootstrap estimation of the variance.
# N and T are list of variables specify the parameters in population topology.
# Topology is fixed under this simulation.

def popn3_simple_data_morgan(N,T,L,B,m,length):
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
        samples={"A": 5, "B": 5, "C": 5}, 
        demography=demography, 
        recombination_rate = 1/length,
        sequence_length = bb
    )
    all = all_ibd_segments(ts)
    out = [np.zeros((3,3)) for i in range(len(L)-1)]
    boot_var = [np.zeros((3,3)) for i in range(len(L)-1)]
    for i in range(len(L)-1):
        u = L[i]
        v = L[i+1]

        o = []
        for j in range(10):
            for k in range(j+1,10):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out[i][0][0] = sum(o)/len(o)
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        boot_var[i][0][0] = sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo)

        o = []
        for j in range(10):
            for k in range(j+10,20):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out[i][0][1] = out[i][1][0] = sum(o)/len(o)
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        boot_var[i][0][1] = boot_var[i][1][0] = sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo)


        o = []
        for j in range(10):
            for k in range(j+20,30):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out[i][0][2] = out[i][2][0] = sum(o)/len(o)
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        boot_var[i][0][2] = boot_var[i][2][0] = sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo)


        o = []
        for j in range(10,20):
            for k in range(j+1,20):
                o.append(sum([l for l in all[j][k] if u<l*kk<v]))
        out[i][1][1] = sum(o)/len(o)
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        boot_var[i][1][1] = sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo)


        o = []
        for j in range(10,20):
            for k in range(j+10,30):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out[i][1][2] = out[i][2][1] = sum(o)/len(o)
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        boot_var[i][1][2] = boot_var[i][2][1] = sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo)


        o = []
        for j in range(20,30):
            for k in range(j+1,30):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out[i][2][2] = sum(o)/len(o)
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        boot_var[i][2][2] = sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo)


    return out, boot_var




def popn3_simple_data_recombination(N,T,L,B,rho):
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=N[0])
    demography.add_population(name="B", initial_size=N[1])
    demography.add_population(name="C", initial_size=N[2])
    demography.add_population(name="BC", initial_size=N[3])
    demography.add_population(name="ABC", initial_size=N[4])
    demography.add_population_split(time=T[0], derived=["B", "C"], ancestral="BC")
    demography.add_population_split(time=T[1], derived=["A", "BC"], ancestral="ABC")
    
    m=2
    kk=200
    ts = msprime.sim_ancestry(
        samples={"A": 5, "B": 5, "C": 5}, 
        demography=demography, 
        recombination_rate = rho,
        sequence_length = round(m/rho)
    )
    all = all_ibd_segments(ts)
    out = [np.zeros((3,3)) for i in range(len(L)-1)]
    boot_var = [np.zeros((3,3)) for i in range(len(L)-1)]
    for i in range(len(L)-1):
        u = L[i]
        v = L[i+1]

        o = []
        for j in range(10):
            for k in range(j+1,10):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out[i][0][0] = sum(o)/len(o)
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        boot_var[i][0][0] = sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo)

        o = []
        for j in range(10):
            for k in range(j+10,20):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out[i][0][1] = out[i][1][0] = sum(o)/len(o)
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        boot_var[i][0][1] = boot_var[i][1][0] = sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo)


        o = []
        for j in range(10):
            for k in range(j+20,30):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out[i][0][2] = out[i][2][0] = sum(o)/len(o)
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        boot_var[i][0][2] = boot_var[i][2][0] = sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo)


        o = []
        for j in range(10,20):
            for k in range(j+1,20):
                o.append(sum([l for l in all[j][k] if u<l*kk<v]))
        out[i][1][1] = sum(o)/len(o)
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        boot_var[i][1][1] = sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo)


        o = []
        for j in range(10,20):
            for k in range(j+10,30):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out[i][1][2] = out[i][2][1] = sum(o)/len(o)
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        boot_var[i][1][2] = boot_var[i][2][1] = sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo)


        o = []
        for j in range(20,30):
            for k in range(j+1,30):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out[i][2][2] = sum(o)/len(o)
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        boot_var[i][2][2] = sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo)


    return out, boot_var