import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Compare clusterings of different subsampling')
parser.add_argument('-c', nargs='+', action='append',
                    help='cluster sample files, at least two needed')
parser.add_argument('-o',
                    help='out')
args = parser.parse_args()


sz = 50

def lineage_size(d1, d2):
  large_lineages1 = {c:i for i, c in enumerate(d1.groupby('cluster')['cluster'].count().sort_values(ascending=False).reset_index(name='count')[0:sz]['cluster'].to_list())}
  large_lineages2 = {c:i for i, c in enumerate(d2.groupby('cluster')['cluster'].count().sort_values(ascending=False).reset_index(name='count')[0:sz]['cluster'].to_list())}
  sizes1 = {k:v for k,v in d1.groupby('cluster')['cluster'].count().sort_values(ascending=False).reset_index(name='count')[0:sz].to_numpy()}
  sizes2 = {k:v for k,v in d2.groupby('cluster')['cluster'].count().sort_values(ascending=False).reset_index(name='count')[0:sz].to_numpy()}
  cnt = np.zeros(shape=(sz, sz), dtype=np.int32)
  for c in d1['Accession.ID']:
    #print(d1[d1['Accession.ID'] == c]['cluster'].values[0])
    if len(d2[d2['Accession.ID'] == c]['cluster'].values) == 0:
        continue
    c1 = d1[d1['Accession.ID'] == c]['cluster'].values[0]
    c2 = d2[d2['Accession.ID'] == c]['cluster'].values[0]
    if c1 in large_lineages1 and c2 in large_lineages2:
      cnt[large_lineages1[c1]][large_lineages2[c2]] += 1
  return cnt, sizes1, sizes2

print(f'opening files {args.c[0][0]} and {args.c[1][0]}')
d1 = pd.read_csv(args.c[0][0], delimiter='\t')
d2 = pd.read_csv(args.c[1][0], delimiter='\t')
cnt, large_lineages1, large_lineages2 = lineage_size(d1, d2)
with open(args.o, 'w') as f:
	line = '\t'.join([''] + list([f'{k}({v})' for k, v in large_lineages2.items()]))
	print(line, file=f)
	for cnt_row, c in zip(cnt, list(large_lineages1.keys())):
		line = '\t'.join([f'{c}({large_lineages1[c]})'] + list([str(v) for v in cnt_row]))
		print(line, file=f)
	f.close()


