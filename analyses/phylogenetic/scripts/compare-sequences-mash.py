import sys
import numpy as np

dist_file = sys.argv[1]
match_file = sys.argv[2]

dist = {}
with open(dist_file) as f:
    for l in f:
        x = l.strip().split('\t')
        # (query,sample)
        dist[(x[1], x[0])] = float(x[2])
        
match = {}
g1, g2 = [], []
lineage_samples = {}
with open(match_file) as f:
    for l in f:
        x = l.strip().split(' ')
        i1, i2 = x[2], x[3]
        g1.append(i1)
        g2.append(i2)
        if i1 not in match: match[i1] = []
        match[i1].append(i2)
        if x[1] not in lineage_samples: lineage_samples[x[1]] = []
        lineage_samples[x[1]].append(x[3])

for k in match.keys():
    match[k] = set(match[k])

for k in lineage_samples.keys():
    lineage_samples[k] = list(set(lineage_samples[k]))

def distance(s1, s2):
    return sum([a!=b for a,b in zip(s1, s2)])

g2 = list(set(g2))
for s in g2:
    print(s, end=' ')
    print(s, end=' ', file=sys.stderr)
#for lin, samples in lineage_samples.items():
#    print(lin, end=' ')
print()
for q in match.keys():
    print(q, end=' ')
    row_values = []
    #for lin, samples in lineage_samples.items():
    #    min_d = 10
    #    for s in samples:
    #        if s in match[q]: print('*', end='')
    #        d = dist[(q, s)] if (q,s) in dist else 9
    #        min_d = min(d, min_d)
    #    print(min_d, end=' ')
    #    row_values.append(min_d)
    for s in g2:
        if s in match[q]: print('*', end='')
        d = dist[(q, s)] if (q,s) in dist else 9
        print(d, end=' ')
        row_values.append(d)
        
    if sorted(row_values)[1] == 0.0:
        #print(f'B {q} {sum(np.array(row_values)==0.0)}', file=sys.stderr)
        for s in g2:
            if s in match[q]: print('*', end='', file=sys.stderr)
            d = dist[(q, s)] if (q,s) in dist else 9
            print(d, end=' ', file=sys.stderr)
        print(file=sys.stderr)
    print()

