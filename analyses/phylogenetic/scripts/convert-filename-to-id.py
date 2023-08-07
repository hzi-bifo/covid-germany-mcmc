import sys

filename_id = {}
def get_filename_id(fn):
	if fn not in filename_id:
		with open(fn) as f:
			l = f.readline()
			assert(l[-1] == '\n')
			l = l[:-1].split('|')[0]
			assert(len(l) > 0 and l[0] == '>')
			filename_id[fn] = l[1:]
	return filename_id[fn]

for l in sys.stdin:
	x = l[:-1].split()
	for i in range(2):
		x[i] = get_filename_id(x[i])
	print('\t'.join(x))

