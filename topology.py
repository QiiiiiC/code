#Give two topologies that we may be interested in, one with admixture and one without.
#Note that the times and the fractions of admixture are assigned randomly, we only care about the topology.


nodes_simple = {
    1:{'parents':[5], 'children':[],'time':0,'frac':1},
    2:{'parents':[4], 'children':[],'time':0,'frac':1},
    3:{'parents':[4], 'children':[],'time':0,'frac':1},
    4:{'parents':[5], 'children':[2,3],'time':1,'frac':1},
    5:{'parents':[], 'children':[1,4],'time':2,'frac':1}
}
events_simple = [
    [0,0,[1,2,3],1],
    [4,1,[1,4],1],
    [5,2,[5],1]
]


nodes_admix = {
    1:{'parents':[4,5], 'children':[],'time':0,'frac':1},
    2:{'parents':[6], 'children':[],'time':0,'frac':1},
    3:{'parents':[7], 'children':[],'time':0,'frac':1},
    4:{'parents':[8], 'children':[1],'time':1,'frac':0.5},
    5:{'parents':[6], 'children':[1],'time':1,'frac':0.5},
    6:{'parents':[7], 'children':[5,2],'time':2,'frac':1},
    7:{'parents':[8], 'children':[6,3],'time':3,'frac':1},
    8:{'parents':[], 'children':[4,7],'time':4,'frac':1}
}
events_admix = [
    [0,0,[1,2,3],1],
    [4,1,[4,5,2,3],0],
    [5,1,[4,5,2,3],0],
    [6,2,[4,6,3],1],
    [7,3,[4,7],1],
    [8,4,[8],1]
]