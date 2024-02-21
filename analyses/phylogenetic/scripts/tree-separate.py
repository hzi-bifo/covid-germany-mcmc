import sys

cnt=0
pre=''
with open(sys.argv[1]) as f:
    for l in f:
        if l.startswith('tree '):
            x = l.split(' ')
            name = x[1].replace('STATE_', '')
            #if x[2].startswith('['): l = ' '.join(x[:2]+x[3:])
            with open(sys.argv[2].replace('NAME', name), 'w') as fout:
                print(pre, l, 'End;', sep='', file=fout)
            #print(name, file=sys.stderr)
            cnt+=1
        else:
            pre += l

print('saved ', cnt, 'files', file=sys.stderr)
