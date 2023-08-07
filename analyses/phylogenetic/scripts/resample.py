import sys, random
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--tree')
parser.add_argument('--out')
parser.add_argument('--burnin', type=int)
parser.add_argument('--rate', type=int)
parser.add_argument('--count', type=int)

args = parser.parse_args()

def load_file(f, tree_op, tree_op_default):
	region = 0
	for l in f:
		#print('L {}'.format(l), end='')
		if region == 3:
			if l.startswith('tree'):
				if l[0:11] != 'tree STATE_': raise RuntimeError('invalid tree string {}'.format(l[:50]))
				sp = l.find(' ', 11)
				#last_state, state = state, int(l[11:sp])
				state = int(l[11:sp])
				#if inc == -1 and last_state != -1: inc = state - last_state
				tree_op(state, l)
				continue
			else:
				pass
		elif l.startswith('Begin trees;'): 
			region = 1
		elif region == 1 and l.strip().startswith('Translate'): 
			region = 2
		elif region == 2 and l.strip() == ';': 
			region = 3
		#print(l, file=fw, end='')
		tree_op_default(l)

print('opening file {} to write to {}'.format(args.tree, args.out))
with open(args.tree, 'r') as f, open(args.out, 'w') as fw:

	def tree_print(l):
		print(l, file=fw, end='')


	
	if args.rate is not None:
		sampled_count = 0
		last_range_begin = 0

		def check_tree_print(state, l):
			global sampled_count, last_range_begin
			if state > args.burnin and state >= args.rate + last_range_begin:
				last_range_begin = state // args.rate * args.rate
				sampled_count += 1
				print(l, file=fw, end='')
				#print('state for sampling {} {}'.format(state, last_range_begin))
		load_file(f, check_tree_print, tree_print)
		print('sampled {} samples'.format(sampled_count))

	if args.count is not None:
		index = 0
		tree_index = []
		def count_tree(state, l):
			global index, tree_index
			if state > args.burnin:
				tree_index.append(index)
			index += 1
		load_file(f, count_tree, lambda l : 1)
		f = open(args.tree, 'r')
		sample_index = set(random.sample(tree_index, args.count))
		index, sampled_count = 0, 0
		def write_sampled(state, l):
			global index, sample_index, sampled_count
			if index in sample_index:
				print(l, file=fw, end='')
				sampled_count += 1
			index += 1
			
		load_file(f, write_sampled, tree_print)
		print('sampled {} samples'.format(sampled_count))

