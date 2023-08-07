#!/usr/bin/env python
import argparse, sys, re
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument("--template", help="template")
parser.add_argument("--name", help="name")
parser.add_argument("--date", help="date")
#parser.add_argument("--grid_points", type=int, help="number of points to be used for grid")
#parser.add_argument("--cutoff", help="distance from beggining to the end of the analysis")
parser.add_argument("--clock_rate", help="initial clock rate")
parser.add_argument("--sample", help="samples file name")
parser.add_argument("--tree", help="tree file name")
parser.add_argument("--tree_init", help="tree file name")
parser.add_argument("--metadata", help="metadata file name")
parser.add_argument("--loc", help="state name")
parser.add_argument("--out", help="output xml file")
parser.add_argument("--chain", type=int, default=100000000, help="chain length for MCMC")
parser.add_argument("--loc_depth", help="instead of in/out print state (third element in location of metadata file)", action='store_true')
args = parser.parse_args()

def load_file(fn):
	with open(fn,'r') as f:
		return f.read()

template = load_file(args.template)
tree = load_file(args.tree).strip().replace(';','') if args.tree is not None else None

with open(args.sample, 'r') as f:
	samples = [x.strip() for x in f.readlines()]

with open(args.metadata, 'r') as f:
	accession_id = -1
	metadata = {}
	for l in f:
		x = l[:-1].split('\t')
		if accession_id == -1: accession_id = x.index('Accession.ID')
		else:
			metadata[x[accession_id]] = x
def d2y(d):
	y = int(d[:d.find('-')])
	d1 = datetime.strptime(str(y) + '-01-01', "%Y-%m-%d")
	d2 = datetime.strptime(d, "%Y-%m-%d")
	d3 = datetime.strptime(str(y+1) + '-01-01', "%Y-%m-%d")
	r = y + (d2 - d1).days / (d3-d1).days
	#print('Y: {} -> {}'.format(d, r), file=sys.stderr)
	return r

def loc_inout(s):
	x = [x.strip() for x in s.split('/')]
	l = args.loc if x[1] == args.loc else 'non' + args.loc
	#print('L: {} -> {}'.format(s, l), file=sys.stderr)
	return l

def loc_state(s):
	x = [x.strip() for x in s.split('/')]
	l = x[2] if len(x) >= 3 else 'NA'
	#print('L: {} -> {}'.format(s, l), file=sys.stderr)
	return l

loc = loc_inout
if args.loc_depth == True: 
    loc = loc_state

taxa = ''
date_max, date_min = 0, 3000
for s in samples:
	m = metadata[s]
	date_max, date_min = max(date_max, d2y(m[3])), min(date_min, d2y(m[3]))
	taxa += '\t\t\t<taxon id="{}"><date value="{}" direction="forwards" units="years" uncertainty="0.0"/><attr name="location">{}</attr></taxon>\n'.format(m[2], d2y(m[3]), loc(m[4]))

def replace_eval(template, tag, value):
	template = template.replace('<__'+tag+'__>', str(value))
	for p, f in [(re.compile('\<__'+tag+':(\d+)__\>'), lambda v, d: v // d), (re.compile('\<__'+tag+'\+(\d+)__\>'), lambda v, d: v + d)]:
		while True:
			m = p.search(template)
			if m is None: break
			#v = value // int(m.groups()[0])
			v = f(value, int(m.groups()[0]))
			template = template[:m.span()[0]] + str(v) + template[m.span()[1]:]
	return template

args.cutoff = date_max - 2020.5
args.grid_points = int(args.cutoff * 52)

#print(taxa)
template = replace_eval(template, 'TAXASIZE', len(samples))
template = replace_eval(template, 'GRIDPOINTS', args.grid_points)
template = template.replace('<__CHAIN__>', str(args.chain))
template = template.replace('<__CUTOFF__>', str(args.cutoff))
if args.clock_rate is not None:
	template = template.replace('<__CLOCK_RATE__>', args.clock_rate)
template = template.replace('<__TAXA__>', taxa)
if args.tree is not None:
	template = template.replace('<__DATA_TREE__>', load_file(args.tree))
if args.tree_init is not None:
	template = template.replace('<__STARTING_TREE__>', load_file(args.tree_init))
template = template.replace('__NAME__', args.name)
template = template.replace('__DATE__', args.date)

with open(args.out, 'w') as f:
	f.write(template)

