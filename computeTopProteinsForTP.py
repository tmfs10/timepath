
import sys
import math
import random

if len(sys.argv) < 2:
	print >> sys.stderr, "Usage:", sys.argv[0], "<config file>"
	sys.exit(1)

def arg_value(params, k, reqd, default_value=None):
	k_lower = k.lower()
	if k_lower not in params:
		if reqd:
			print >> sys.stderr, "Error: Parameter", k, "required in the config file!"
			sys.exit(1)
		return default_value
	else:
		return params[k_lower]

config_file = sys.argv[1]
params = {}
with open(config_file) as f:
	for line in f:
		temp = line.strip()
		if temp == "" or temp.startswith("#"):
			continue
		line = line.split('#')
		if len(line) >= 2 and len(line[0].strip()) == 0:
			continue
		if len(line[0].strip()) == 0:
			continue
		line = [k.strip().lower() for k in line[0].strip().split('\t')]
		if len(line) == 0:
			continue
		if len(line) == 1:
			print >> sys.stderr, "ERROR: parameter " + line[0] + " in config file " + config_file + " does not have a value!"
			print >> sys.stderr, "NOTE: ALL lines but be of the form key<tab>value. Comments after # character are ignored"
			sys.exit(1)

		k, v = line
		if k in params:
			print >> sys.stderr, "Error: Duplicate parameter", k
			sys.exit(1)
		params[k] = v

path_filepath = arg_value(params, "pathFile", True)
source_file = arg_value(params, "sourcesFilepath", True)
phase_size = int(arg_value(params, "numTPtargets", True))
tf_filepath = arg_value(params, "tfGeneFile", True)
rna_hits_filepath = arg_value(params, "rnaHitsFilepath", False)
num_genes_to_select = int(arg_value(params, "numProteinsPerPhaseToSelect", False, 50))
upregulated = True
min_ranking_fold_change = float(arg_value(params, "minRankingFoldChange", False, 1.5))
phase_gene_file = arg_value(params, "pathFile", True) + ".phasegenes.txt"

all_tf_set = set()
with open(tf_filepath) as f:
	next(f)
	for line in f:
		line = [k.strip() for k in line.strip().split("\t")]
		tf = line[0]
		all_tf_set.add(tf)


rna_hits_set = set()
if rna_hits_filepath is not None:
	with open(rna_hits_filepath) as f:
		for line in f:
			line = [k.strip() for k in line.strip().split("\t")]
			rna_hits_set.add(line[0])


source_nodes = set()
with open(source_file) as f:
	for line in f:
		line = [k.strip() for k in line.strip().split("\t")]
		source_nodes.add(line[0])

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

phase_source_weights = [{} for i in xrange(num_phases)]
phase_node_pp_weights = [{} for i in xrange(num_phases)]
phase_tf_weights = [{} for i in xrange(num_phases)]
joint_weights = [{} for i in xrange(num_phases)]

phase_tf_weights_to_targ = [{} for i in xrange(num_phases)]
joint_weights_to_targ = [{} for i in xrange(num_phases)]

total_phase_weight = [0.]*num_phases

# This is what decides the phase-specific node rankings
# via phase_source_weights, phase_node_pp_weights, and phase_tf_weights
# and thus which nodes should appear in the figure
print >> sys.stderr, "reading paths"
with open(path_filepath) as f:
	num_done = 0
	for line in f:
		line = [k.strip() for k in line.strip().split("\t")]
		num_done += 1
		if num_done%1000000 == 0:
			print >> sys.stderr, num_done/1000000, "million paths done"
		score = float(line[0])
		path = line[2:]
		targ = path[0]
		tf = path[1]
		path_phase = -1
		for i in xrange(num_phases):
			if path[0] in phase_nodes[i]:
				path_phase = i
				break
		if path_phase < 0:
			continue

		total_phase_weight[path_phase] += score

		# last_phase_gene is the current phase target if this is the 
		# first phase else the previous phase target
		last_phase_gene = len(path)-1
		if path_phase > 0:
			for i in reversed(xrange(1, len(path)-1)):
				if path[i] in phase_nodes[path_phase-1]:
					last_phase_gene = i
					break
		assert last_phase_gene > 0
		assert path_phase == 0 or last_phase_gene < len(path)-1, str(path) + "\t" + str(path[last_phase_gene]) + "\t" + str(path_phase)

		# intermediate nodes -- note that path is in reverse order so we're summing
		# path weights for nodes starting after the last phase gene
		for i in xrange(1, last_phase_gene+1):
			g = path[i]
			if g not in joint_weights[path_phase]:
				joint_weights[path_phase][g] = 0
				joint_weights_to_targ[path_phase][g] = {targ : 0. for targ in phase_nodes[path_phase]} 
			if g not in phase_node_pp_weights[path_phase]:
				phase_node_pp_weights[path_phase][g] = 0
			phase_node_pp_weights[path_phase][g] += score
			joint_weights[path_phase][g] += score
			joint_weights_to_targ[path_phase][g][targ] += score

		# source -- it's the last phase gene unless we're at phase 0 in which
		# case it's the actual source
		source = path[last_phase_gene]
		if source not in joint_weights[path_phase]:
			joint_weights[path_phase][source] = 0
			joint_weights_to_targ[path_phase][source] = {targ : 0. for targ in phase_nodes[path_phase]} 
		if source not in phase_source_weights[path_phase]:
			phase_source_weights[path_phase][source] = 0
		phase_source_weights[path_phase][source] += score
		joint_weights[path_phase][g] += score
		joint_weights_to_targ[path_phase][g][targ] += score

		# tf for current phase
		if tf not in joint_weights[path_phase]:
			joint_weights[path_phase][tf] = 0
			joint_weights_to_targ[path_phase][tf] = {targ : 0. for targ in phase_nodes[path_phase]} 
		if path[1] not in phase_tf_weights[path_phase]:
			phase_tf_weights[path_phase][path[1]] = 0
			phase_tf_weights_to_targ[path_phase][path[1]] = {targ : 0. for targ in phase_nodes[path_phase]} 
		phase_tf_weights[path_phase][path[1]] += score
		phase_tf_weights_to_targ[path_phase][path[1]][targ] += score
		joint_weights[path_phase][g] += score
		joint_weights_to_targ[path_phase][g][targ] += score

phase_source_ranking = [sorted([g for g in phase_source_weights[i]], key=lambda k: phase_source_weights[i][k], reverse=True) for i in xrange(num_phases)]

phase_node_pp_ranking = [sorted([g for g in phase_node_pp_weights[i]], key=lambda k: phase_node_pp_weights[i][k], reverse=True) for i in xrange(num_phases)]

phase_tf_ranking = [sorted([g for g in phase_tf_weights[i]], key=lambda k: phase_tf_weights[i][k], reverse=True) for i in xrange(num_phases)]

joint_ranking_list = [sorted([g for g in joint_weights[i]], key=lambda k: joint_weights[i][k], reverse=True) for i in xrange(num_phases)]
joint_ranking = [{t : i for t,i in zip(joint_ranking_list[j], xrange(1, len(joint_ranking_list[j])+1))} for j in xrange(num_phases)]

# Seems like the best bet might turn out to be genes that are differentially 
# expressed in their own time phase. We also want to prioritize factors, the
# removal of which will remove several targets.

with open(path_filepath + ".topproteins.txt", 'w') as f:
	if upregulated:
		for i in xrange(1, num_phases):
			for j in xrange(i):
				for t in joint_ranking[i]:
					if t not in joint_ranking[j]:
						joint_ranking[j][t] = 100000

		ranking_change = [{t : min((float(joint_ranking[k][t])/float(joint_ranking[i][t]) for k in xrange(i))) for t in joint_ranking[i]} for i in xrange(1, num_phases)]
		genes_selected = [{} for i in xrange(num_phases)]
		for i in xrange(num_phases):
			j = 0
			for g in joint_ranking_list[i]:
				if j >= num_genes_to_select:
					break
				use = True if i == 0 else ranking_change[i-1][g] >= min_ranking_fold_change
				in_which_phases = []
				gene_type = [1 if g in phase_source_ranking[i] else 0,
						1 if g in phase_node_pp_ranking[i] else 0,
						1 if g in phase_tf_ranking[i] else 0]
				for k in xrange(i):
					if g in genes_selected[k]:
						use = False
						break
				if not use:
					continue
				j += 1
				genes_selected[i][g] = [1. if i == 0 else ranking_change[i-1][g], joint_weights[i][g]/total_phase_weight[i], gene_type]

		for i in xrange(num_phases):
			print >> f, "Phase", i
			print >> f, "========"
			for g in genes_selected[i]:
				if joint_ranking[i][g] > 10000:
					continue
				if i == 0:
					print >> f, g, " ".join(["NP" if g not in joint_ranking[j] or joint_ranking[j][g] > 10000 else str(joint_ranking[j][g]) for j in xrange(num_phases)]), genes_selected[i][g][1], 1 if g in rna_hits_set else 0, " ".join((str(k) for k in genes_selected[i][g][2])) # gene, ranking in every phase, fraction weight for phase, if rna hit, is source for phase, intermediate for phase, tf for phase
				else:
					print >> f, g, genes_selected[i][g][0], " ".join(["NP" if g not in joint_ranking[j] or joint_ranking[j][g] > 10000 else str(joint_ranking[j][g]) for j in xrange(num_phases)]), genes_selected[i][g][1], 1 if g in rna_hits_set else 0, 1 if g in source_nodes else 0, " ".join((str(k) for k in genes_selected[i][g][2]))
			print >> f, "\n\n"
	else:
		for i in xrange(num_phases):
			for j in xrange(num_phases):
				if i == j:
					continue
				for t in joint_ranking[i]:
					if t not in joint_ranking[j]:
						joint_ranking[j][t] = 100000

		ranking_change = [{t : min((float(joint_ranking[k][t]) for k in xrange(i+1, num_phases)))/max((float(joint_ranking[k][t]) for k in xrange(i+1))) for t in joint_ranking[i]} for i in xrange(num_phases-1)]
		genes_selected = [{} for i in xrange(num_phases)]
		for i in xrange(num_phases-1):
			j = 0
			for g in joint_ranking_list[i]:
				if j >= num_genes_to_select:
					break
				use = ranking_change[i][g] >= min_ranking_fold_change
				in_which_phases = []
				gene_type = [1 if g in phase_source_ranking[i] else 0,
						1 if g in phase_node_pp_ranking[i] else 0,
						1 if g in phase_tf_ranking[i] else 0]
				for k in xrange(i):
					if g in genes_selected[k]:
						use = False
						break
				if not use:
					continue
				j += 1
				genes_selected[i][g] = [ranking_change[i][g], joint_weights[i][g]/total_phase_weight[i], gene_type]

		for i in xrange(num_phases-1):
			print >> f, "Phase", i
			print >> f, "========"
			for g in genes_selected[i]:
				if joint_ranking[i][g] > 10000:
					continue
				if i == 0:
					print >> f, g, genes_selected[i][g][0], " ".join(["NP" if joint_ranking[j][g] > 10000 else str(joint_ranking[j][g]) for j in xrange(num_phases)]), genes_selected[i][g][1], " ".join((str(k) for k in genes_selected[i][g][2])), 1 if g in rna_hits_set else 0
				else:
					print >> f, g, genes_selected[i][g][0], " ".join(["NP" if joint_ranking[j][g] > 10000 else str(joint_ranking[j][g]) for j in xrange(num_phases)]), genes_selected[i][g][1], " ".join((str(k) for k in genes_selected[i][g][2])), 1 if g in rna_hits_set else 0, 1 if g in source_nodes else 0 # gene, change in ranking vs previous phase, ranking in every phase, fraction weight for phase, is source for phase, intermediate for phase, tf for phase, rna hit?, main source?
			print >> f, "\n\n"


	print >> sys.stderr, [len(k) for k in joint_ranking_list]
