import msprime
import numpy as np
import random
import json
import math
from topology import nodes_simple, events_simple

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
    out = {'y':[],'sigma':[],'u':[],'v':[],'group':[]}

    for i in range(len(L)-1):
        u = L[i]
        v = L[i+1]

        o = []
        for j in range(10):
            for k in range(j+1,10):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out['y'].append(sum(o)/len(o))
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        out['sigma'].append(sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo))
        out['u'].append(u)
        out['v'].append(v)
        out['group'].append([1,1])

        o = []
        for j in range(10):
            for k in range(j+10,20):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out['y'].append(sum(o)/len(o))
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        out['sigma'].append(sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo))
        out['u'].append(u)
        out['v'].append(v)
        out['group'].append([1,2])


        o = []
        for j in range(10):
            for k in range(j+20,30):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out['y'].append(sum(o)/len(o))
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        out['sigma'].append(sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo))
        out['u'].append(u)
        out['v'].append(v)
        out['group'].append([1,3])

        o = []
        for j in range(10,20):
            for k in range(j+1,20):
                o.append(sum([l for l in all[j][k] if u<l*kk<v]))
        out['y'].append(sum(o)/len(o))
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        out['sigma'].append(sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo))
        out['u'].append(u)
        out['v'].append(v)
        out['group'].append([2,2])


        o = []
        for j in range(10,20):
            for k in range(j+10,30):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out['y'].append(sum(o)/len(o))
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        out['sigma'].append(sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo))
        out['u'].append(u)
        out['v'].append(v)
        out['group'].append([2,3])

        o = []
        for j in range(20,30):
            for k in range(j+1,30):
                a = [l for l in all[j][k] if u<l*kk<v]
                o.append(sum(a))
        out['y'].append(sum(o)/len(o))
        oo = []
        for b in range(B):
            p = random.choices(o, k=len(o))
            oo.append(sum(p)/len(p))
        out['sigma'].append(sum((x - sum(oo)/len(oo)) ** 2 for x in oo)/len(oo))
        out['u'].append(u)
        out['v'].append(v)
        out['group'].append([3,3])

    return out


L = [0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.25,1.5,2,3,5,10]
N = [2500,2500,2500,2500,2500]
T = [15,50]
B = 3000
output = popn3_simple_data_morgan(N,T,L,B,2,1e7)

#prepare the membership matrix A
A = mig_matrix(nodes_simple, events_simple, create_nodes_map(nodes_simple))
output['N_obs'] = 6*(len(L)-1)
output['N_popn'] = 5
output['N_events'] = 2
output['A'] = A

with open('generalized_popn_simulated_data.json', 'w') as json_file:
    json.dump(output, json_file)
