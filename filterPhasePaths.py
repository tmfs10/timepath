
import sys
import math
import random

if len(sys.argv) < 5:
	print >> sys.stderr, "Usage:", sys.argv[0], "<sorted path file> <phase size> <phase genes file> <num paths>"
	sys.exit(1)

path_filepath = sys.argv[1]
phase_size = int(sys.argv[2])
phase_gene_file = sys.argv[3]
num_paths_to_select = int(sys.argv[4])

num_phases = 0
with open(phase_gene_file) as f:
	temp = next(f)
	if not '\t' in temp:
		print >> sys.stderr, "WARNING: No tab character found in file that has targets for each time point. Targets for different time points MUST be TAB separated"
	num_phases = len(temp.strip().split('\t'))
if num_phases == 0:
	print >> sys.stderr, "ERROR: The file containing targets for each time point has a blank line at the top"
	sys.exit(1)

phase_nodes = [set() for i in xrange(num_phases)]
for i in xrange(num_phases):
	genes_read = 0
	with open(phase_gene_file) as f:
		for line in f:
			line = [k.strip() for k in line.strip().split("\t")]
			assert len(line) == num_phases, str(line) + "\t" + str(num_phases)
			g = line[i]
			genes_read += 1

			# This is for not duplicating genes across phases
			use_this = True
			for j in xrange(i):
				if g in phase_nodes[j]:
					use_this = False
					break
			if use_this:
				phase_nodes[i].add(g)
			if genes_read >= phase_size:
				break
# This is what decides the phase-specific node rankings
# via phase_source_weights, phase_node_pp_weights, and phase_tf_weights
# and thus which nodes should appear in the figure
print >> sys.stderr, "reading paths"
num_phase_done = [{} for i in xrange(num_phases)]
with open(path_filepath) as f:
	num_done = 0
	for orig_line in f:
		line = [k.strip() for k in orig_line.strip().split("\t")]
		num_done += 1
		if num_done%1000000 == 0:
			print >> sys.stderr, num_done/1000000, "million paths done"
		score = float(line[0])
		path = line[2:]
		targ = path[0]
		path_phase = -1
		for i in xrange(num_phases):
			if targ in phase_nodes[i]:
				path_phase = i
				break
		if path_phase < 0:
			continue
		if targ not in num_phase_done[path_phase]:
			num_phase_done[path_phase][targ] = 0
		if num_paths_to_select > 0 and num_phase_done[path_phase][targ] >= num_paths_to_select:
			continue
		print str(path_phase) + "\t" + orig_line.strip()
		num_phase_done[path_phase][targ] += 1
