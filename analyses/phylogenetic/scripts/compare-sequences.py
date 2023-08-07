import sys
import numpy as np

seq_file = sys.argv[1]
match_file = sys.argv[2]

seqs = {}
with open(seq_file) as f:
    for l in f:
        seqs[l.strip()[1:]] = f.readline().strip()

match = {}
g1, g2 = [], []
with open(match_file) as f:
    for l in f:
        x = l.strip().split(' ')
        g1.append(x[0])
        g2.append(x[1])
        if x[0] not in match: match[x[0]] = []
        match[x[0]].append(x[1])

for k in match.keys():
    match[k] = set(match[k])

def distance(s1, s2):
    return sum([a!=b for a,b in zip(s1, s2)])

g2 = list(set(g2))
for s in g2:
    print(s, end=' ')
print()
for q in match.keys():
    print(q, end=' ')
    for s in g2:
        if s in match[q]: print('*', end='')
        d = distance(seqs[q], seqs[s])
        if s in match[q] and d != 0:
            bads = [a!=b for a,b in zip(seqs[q], seqs[s])]
            bidx = np.arange(len(seqs[q]))[bads]
            print(f'B {q} {s} {d} {bidx}', file=sys.stderr)
        print(d, end=' ')
    print()

