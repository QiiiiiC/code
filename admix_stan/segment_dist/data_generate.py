import numpy as np
import msprime


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

def popn1_simple_morgan(N,L,m,length,n):
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=N)
    
    bb = m*length
    kk = m*100

    ts = msprime.sim_ancestry(
        samples={"A": n}, 
        demography=demography, 
        recombination_rate = 1/length,
        sequence_length = bb
    )
    all = all_ibd_segments(ts)
    out = {'y':[],'u':[],'v':[]}

    for i in range(len(L)-1):
        u = L[i]
        v = L[i+1]
        for j in range(n*2):
            for k in range(j+1,n*2):
                a = [l*kk for l in all[j][k] if u<l*kk<v]
                out['y'] += sum(a)/len(a)
        
                out['u']+=[u]
                out['v']+=[v]
    return out
    
def popn2_simple_data_morgan(N,T,L,m,length,n):
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=N[0])
    demography.add_population(name="B", initial_size=N[1])
    demography.add_population(name="AB", initial_size=N[2])
    demography.add_population_split(time=T, derived=["A", "B"], ancestral="AB")
    
    bb = m*length
    kk = m*100

    ts = msprime.sim_ancestry(
        samples={"A": n, "B": n}, 
        demography=demography, 
        recombination_rate = 1/length,
        sequence_length = bb
    )
    all = all_ibd_segments(ts)
    out = {'y':[],'u':[],'v':[],'group':[]}

    for i in range(len(L)-1):
        u = L[i]
        v = L[i+1]
        for j in range(n*2):
            for k in range(j+1,n*2):
                a = [l*kk for l in all[j][k] if u<l*kk<v]
                out['y'] += sum(a)/len(a)
        
                out['u']+=[u]
                out['v']+=[v]
                out['group']+=[1]

        for j in range(n*2):
            for k in range(j+n*2,n*4):
                a = [l*kk for l in all[j][k] if u<l*kk<v]
                out['y'] += sum(a)/len(a)
        
                out['u']+=[u]
                out['v']+=[v]
                out['group']+=[2]     

        for j in range(n*2,n*4):
            for k in range(j+1,n*4):
                a = [l*kk for l in all[j][k] if u<l*kk<v]
                out['y'] += sum(a)/len(a)
        
                out['u']+=[u]
                out['v']+=[v]
                out['group']+=[3]
        
    return out


