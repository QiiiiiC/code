import msprime
import numpy as np
import random
import json
import math
import pandas as pd

def create_nodes_map(nodes):
    nodes_map = {}
    it = 0
    for i in nodes.keys():
        nodes_map[i] = it
        it += 1
    return nodes_map

# create a list of migration matrices, whenever there's a event(split/merge), create a migration matrix for that particular event.
# the starting element is always an identity matrix.
# each migration matrix is of dimension k*k where k is the total number of all populations(instead of starting populations).
def mig_matrix(nodes, events, nodes_map):
    n = len(nodes)
    events_num = int(len(events)-1-len([e for e in events if e[-1]==0])/2)
    out = [np.diag([1.0 for i in range(n)]).tolist() for i in range(events_num + 1)]
    
    # cnt is for loop indices.
    cnt = 1
    for i in range(1,len(events)):
        # when there's a split backward in time, two elements are created in 'events', so skip one of them.
        if events[i][1] != events[i-1][1]:
            # this is a merge backward in time.
            if events[i][-1] == 1:
                dest = events[i][0]
                change = nodes[dest]['children']
                for j in change:
                    out[cnt][nodes_map[j]][nodes_map[j]]=0
                    out[cnt][nodes_map[j]][nodes_map[dest]]=1
                cnt += 1
            # this is a split backward in time
            else:
                change = nodes[events[i][0]]['children'][0]
                dest = nodes[change]['parents']
                out[cnt][nodes_map[change]][nodes_map[change]]=0
                for j in dest:
                    out[cnt][nodes_map[change]][nodes_map[j]]=nodes[j]['frac']
                cnt += 1
    return out

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

def popn3_simple_data_morgan(N,T,L,m,length,n,seed):
    #Set the follwing model as the default model.
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=N[0])
    demography.add_population(name="C", initial_size=N[1])
    demography.add_population(name="ADMIX", initial_size=N[2])
    demography.add_population(name="ANC", initial_size=N[3])
    demography.add_admixture(time=T[0], derived='ADMIX', ancestral=["A","C"],proportions = [0.25,0.75])
    demography.add_population_split(time=T[1], derived=["A", "C"], ancestral="ANC")
    
    # bb is the number of bits simulated.
    bb = m*length

    # kk is measured in centiMorgans.
    kk = m*100
    ts = msprime.sim_ancestry(
        samples={"A": n, "ADMIX": n, "C": n}, 
        demography=demography, 
        recombination_rate = 1/length,
        sequence_length = bb,
        random_seed = seed
    )
    all = all_ibd_segments(ts)
    out = {'index':[],'number':[],'fraction':[],'u':[],'v':[],'group':[],'true_N':[],'true_T1':[],'true_T2':[]}

    for i in range(len(L)-1):
        u = L[i]
        v = L[i+1]

        #for within data, have (n*2)*(n*2-1)/2 pairs and for between data we have (n*2)*(n*2) pairs.
        #within A
        for j in range(n*2):
            for k in range(j+1,n*2):
                a = [l for l in all[j][k] if u<l*kk<v]
                out['index'] += [[j,k]]
                out['number'] += [len(a)]
                out['fraction'] += [sum(a)]
                out['u'] += [u]
                out['v'] += [v]
                out['group'] += [[0,0]]
                out['true_N'] += [N[0]]
                out['true_T1'] += [T[0]]
                out['true_T2'] += [T[1]]

       #between A and ADMIX
        for j in range(n*2):
            for k in range(n*2,n*4):
                a = [l for l in all[j][k] if u<l*kk<v]
                out['index'] += [[j,k]]
                out['number'] += [len(a)]
                out['fraction'] += [sum(a)]
                out['u'] += [u]
                out['v'] += [v]
                out['group'] += [[0,1]]
                out['true_N'] += [N[0]]
                out['true_T1'] += [T[0]]
                out['true_T2'] += [T[1]]

        #between A and C
        for j in range(n*2):
            for k in range(n*4,n*6):
                a = [l for l in all[j][k] if u<l*kk<v]
                out['index'] += [[j,k]]
                out['number'] += [len(a)]
                out['fraction'] += [sum(a)]
                out['u'] += [u]
                out['v'] += [v]
                out['group'] += [[0,2]]
                out['true_N'] += [N[0]]
                out['true_T1'] += [T[0]]
                out['true_T2'] += [T[1]]
        
        #within ADMIX and ADMIX
        for j in range(n*2,n*4):
            for k in range(j+1,n*4):
                a = [l for l in all[j][k] if u<l*kk<v]
                out['index'] += [[j,k]]
                out['number'] += [len(a)]
                out['fraction'] += [sum(a)]
                out['u'] += [u]
                out['v'] += [v]
                out['group'] += [[1,1]]
                out['true_N'] += [N[0]]
                out['true_T1'] += [T[0]]
                out['true_T2'] += [T[1]]

        #between ADMIX and C
        for j in range(n*2,n*4):
            for k in range(n*4,n*6):
                a = [l for l in all[j][k] if u<l*kk<v]
                out['index'] += [[j,k]]
                out['number'] += [len(a)]
                out['fraction'] += [sum(a)]
                out['u'] += [u]
                out['v'] += [v]
                out['group'] += [[1,2]]
                out['true_N'] += [N[0]]
                out['true_T1'] += [T[0]]
                out['true_T2'] += [T[1]]

        #within C and C
        for j in range(n*4,n*6):
            for k in range(j+1,n*6):
                a = [l for l in all[j][k] if u<l*kk<v]
                out['index'] += [[j,k]]
                out['number'] += [len(a)]
                out['fraction'] += [sum(a)]
                out['u'] += [u]
                out['v'] += [v]
                out['group'] += [[2,2]]
                out['true_N'] += [N[0]]
                out['true_T1'] += [T[0]]
                out['true_T2'] += [T[1]]

    return out



L = [0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.25,1.5,2,3,5,10,200]
N_list = [2000,5000,10000]
T1_list = [10,20,50,100]
T2_list = [110,200,500,1000]
seeds = np.arange(100,1600,100)
for seed in seeds:
    all_data = pd.DataFrame({'index':[],'number':[],'fraction':[],'u':[],'v':[],'group':[],'true_N':[],'true_T1':[],'true_T2':[]})
    for t1 in T1_list:
        for t2 in T2_list:
            for population_size in N_list:
                N = [population_size]*4
                T = [t1] + [t2]
                output = popn3_simple_data_morgan(N,T,L,2,1e6,5,seed)
                all_data = pd.concat([all_data,pd.DataFrame(output)],ignore_index=True)
                print(f'case_{t1}_{t2}_{population_size}_{seed} done')
    out = all_data.to_dict(orient = 'list')
    with open(f'admixture_{seed}.json','w') as json_file:
        json.dump(out,json_file)


