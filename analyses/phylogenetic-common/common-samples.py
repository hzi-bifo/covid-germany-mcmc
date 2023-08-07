import sys
import csv
import pandas as pd
from alive_progress import alive_bar
from collections import Counter
import networkx as nx
import plotly.graph_objs as go
import plotly as plt
import pygraphviz as pgv
import json


def find_intersections_fast(file_names, largest_count):
    def get_cluster_samples(d1, c1):
        return d1[d1['cluster'] == c1]['Accession.ID'].to_numpy()
        
    data = {}
    for n, fn in file_names.items():
        data[n] = pd.read_csv(fn, delimiter='\t')

    for n1, d1 in data.items():
        sample_cluster_1 = {a[0]: a[1] for a in d1[['Accession.ID', 'cluster']].to_numpy()}
        cluster_size_1 = Counter([a[1] for a in d1[['Accession.ID', 'cluster']].to_numpy()])
        for n2, d2 in data.items():
            if n1 == n2: continue
            sample_cluster_2 = d2[['Accession.ID', 'cluster']].to_numpy()
            cluster_size_2 = Counter([a[1] for a in d2[['Accession.ID', 'cluster']].to_numpy()])
            common, not_in_1 = {}, 0
            for k2, c2 in sample_cluster_2:
                if k2 in sample_cluster_1:
                    kp = (sample_cluster_1[k2], c2)
                    if kp not in common: common[kp] = 0
                    common[kp] += 1
                else:
                    not_in_1 += 1
            print(n1, n2, len(d1), len(d2), len(set(sample_cluster_1.keys()).difference(set([a[0] for a in sample_cluster_2]))), not_in_1)
            #for kp, v in sorted(common.items(), key=lambda x: x[1], reverse = True):
            #    print(n1, kp[0], n2, kp[1], cluster_size_1[kp[0]], cluster_size_2[kp[1]], v)
            common_o, common_r = {}, {}
            for kp, v in common.items():
                if kp[0] not in common_o: common_o[kp[0]] = {}
                common_o[kp[0]][kp[1]] = v
                if kp[1] not in common_r: common_r[kp[1]] = {}
                common_r[kp[1]][kp[0]] = v
            G = pgv.AGraph(ranksep="0.5", oneblock=True)

            def lab(n, x):
                n = n.replace('-noid', '')
                if n == 'caseratio': n = 'cr'
                elif n == '5050': n = '5'
                else: n = '2'
                return n + '-' + x.replace('LIN-Germany-', '').replace('-20210602_DTA_MCC_', '-') 
            for i, x in enumerate(common_o.keys()):
                G.add_node(lab(n1, x), pos=(0,i))
                n = G.get_node(lab(n1, x))
                n.attr['pos'] = f'0,{i}!'
                n.attr['label'] = ''
                #n.attr['xlabel'] = lab(n1, x)
                n.attr['width'] = '0.1'
                n.attr['height'] = '0.1'
                n.attr['fillcolor'] = 'red'
                n.attr['color'] = 'red'
                n.attr['style'] = 'filled'
            for i, y in enumerate(common_r.keys()):
                G.add_node(lab(n2, y), pos=(1,i))
                n = G.get_node(lab(n2, y))
                n.attr['pos'] = f'1,{i}!'
                n.attr['label'] = ''
                #n.attr['xlabel'] = lab(n2, y)
                n.attr['width'] = '0.1'
                n.attr['height'] = '0.1'
                n.attr['fillcolor'] = 'green'
                n.attr['color'] = 'green'
                n.attr['style'] = 'filled'
            for kp, v in common.items():
                G.add_edge(lab(n1, kp[0]), lab(n2, kp[1]))
            for x, yv in sorted(common_o.items(), key=lambda x: len(x[1]), reverse = True):
                print(x, cluster_size_1[x])
                for y, v in yv.items():
                    print(' ', y, cluster_size_2[y], v, '(' + str(len(common_r[y])) + ')')
            #pos = nx.drawing.nx_agraph.graphviz_layout(G, prog='dot', args='-Grankdir=LR')
            #G.layout('neato', '-n2')
            G.layout('dot')
            G.draw('results/' + n1 + '-' + n2 + '.pdf')

            #max_x = 0
            #for i, x in enumerate(common_o.keys()):
            #    n = G.get_node(lab(n1, x))
            #    max_x = max(max_x, float(n.attr['pos'].split(',')[0]))
            #for i, y in enumerate(common_r.keys()):
            #    n = G.get_node(lab(n2, y))
            #    max_x = max(max_x, float(n.attr['pos'].split(',')[0]))
            #
            #for i, x in enumerate(common_o.keys()):
            #    n = G.get_node(lab(n1, x))
            #    cc = [float(c) for c in n.attr['pos'].split(',')]
            #    n.attr['pos'] = f'{cc[0]*2},{cc[1]}!'
            #
            #for i, y in enumerate(common_r.keys()):
            #    n = G.get_node(lab(n2, y))
            #    cc = [float(c) for c in n.attr['pos'].split(',')]
            #
            #G.layout('neato', '-n2')
            #G.draw('results/' + n1 + '-' + n2 + '.pdf')
            #json_string = G.pipe('json').decode()
            #json_dict = json.loads(json_string)
            #for obj in json_dict['objects']:
            #    print(obj['name'], '\t', obj['pos'])

            #pos = nx.get_node_attributes(G,'pos')
            #print(pos, file=sys.stderr)
            #nx.draw(G, pos=pos, with_labels=True, font_weight='bold')
            #pygraphviz.write_dot(G, 'results/' + n1 + '-' + n2 + '.dot')
            #plt.show()

find_intersections_fast({
            'caseratio-noid': '../phylogenetic-test-snake-2-nounsampled/results/beast/run/lin-ius-3/clusterSamples_DTA_MCC_0.5.tsv', 
            '5050-noid': '../phylogenetic-test-subsampling-3-nounsampled/results/beast/run/lin-ius-3/clusterSamples_DTA_MCC_0.5.tsv', 
            '10025-noid': '../phylogenetic-test-subsampling-5-nounsampled/results/beast/run/lin-ius-3/clusterSamples_DTA_MCC_0.5.tsv', 
            }, 200)

#find_intersections({
#            'caseratio-rid': '../phylogenetic-test-snake-2-ridentical/results/beast/run/lin-ius-3/clusterSamples_DTA_MCC_0.5.tsv', 
#            '5050-rid': '../phylogenetic-test-subsampling-3-ridentical/results/beast/run/lin-ius-3/clusterSamples_DTA_MCC_0.5.tsv', 
#            '10025-rid': '../phylogenetic-test-subsampling-5-ridentical/results/beast/run/lin-ius-3/clusterSamples_DTA_MCC_0.5.tsv', 
#            }, 30)
#
#find_intersections({
#            'caseratio-mash': '../phylogenetic-test-snake-2/results/beast/run/lin-ius-3/clusterSamples_DTA_MCC_0.5.tsv', 
#            '5050-mash': '../phylogenetic-test-subsampling-3/results/beast/run/lin-ius-3/clusterSamples_DTA_MCC_0.5.tsv', 
#            '10025-mash': '../phylogenetic-test-subsampling-5/results/beast/run/lin-ius-3/clusterSamples_DTA_MCC_0.5.tsv', 
#            }, 30)
