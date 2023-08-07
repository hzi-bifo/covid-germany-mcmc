import sys
import numpy as np

threshold = float(sys.argv[1])
dist_file = sys.argv[2]
seq_file = sys.argv[3]

seqs = {}
with open(seq_file) as f:
    for l in f:
        seqs[l.strip()[1:]] = f.readline().strip()

def distance(a, b):
    return sum([x != y and x != 'n' and y != 'n' for x,y in zip(a,b)]) / len(a)
        
dist = {}
cnt, acc, rej = 0, 0, 0
with open(dist_file) as f:
    for l in f:
        if cnt % 10000 == 0:
            print(f'cnt={cnt} {acc}/{rej}', file=sys.stderr)
        x = l.strip().split('\t')
        # (query,sample)
        cnt += 1
        if float(x[2]) <= threshold:
            x[2] = str(distance(seqs[x[0]], seqs[x[1]]))
            print('\t'.join(x))
            acc += 1
        else:
            rej += 1
