
import sys
import math
import random

def arg_value(params, k, reqd, default_value=None):
	k_lower = k.lower()
	if k_lower not in params:
		if reqd:
			print >> sys.stderr, "Error: Parameter", k, "required in the config file!"
			sys.exit(1)
		return default_value
	else:
		return params[k_lower]

if len(sys.argv) < 2:
	print >> sys.stderr, "Usage:", sys.argv[0], "<config file>"
	sys.exit(1)

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

		k, v = line[:2]
		if k in params:
			print >> sys.stderr, "Error: Duplicate parameter", k
			sys.exit(1)
		params[k] = v

path_filepath = arg_value(params, "pathFile", True)
source_file = arg_value(params, "sourcesFilepath", True)
phase_size = int(arg_value(params, "numTPtargets", True))
gene_files = arg_value(params, "pathFile", True) + ".phasegenes.txt"
num_nodes = int(arg_value(params, "numGenesToDisplay", False, 15))
node_width = int(arg_value(params, "graphNodeWidth", False, 60))
node_height = int(arg_value(params, "graphNodeHeight", False, 40))
pathway_score_threshold = float(arg_value(params, "pathwayScoreThreshold", False, 0))
expr_filepath = arg_value(params, "timeseriesFilepath", True)
tf_filepath = arg_value(params, "tfGeneFile", True)
rna_hits_filepath = arg_value(params, "rnaHitsFilepath", False)
important_genes_filepath = path_filepath + ".topproteins.txt"
expr_node_size = [float(k) for k in arg_value(params, "exprNodeSize", False, "10," + str(int(node_width/2))).split(",")]
tf_gene_interaction_file = arg_value(params, "tfGeneFile", True) 
max_y = 1080
num_clusters = 5

random.seed(0)

gene_to_id = {}
index_to_gene = []
gene_to_logfoldchange = []

tf_gene_edges = {}
with open(tf_gene_interaction_file) as f:
	for line in f:
		break
	for line in f:
		line = [k.strip() for k in line.strip().split()]
		t = line[0]
		g = line[1]
		if t not in tf_gene_edges:
			tf_gene_edges[t] = set()
		tf_gene_edges[t].add(g)

all_phase_nodes = set()

num_phases = 0
with open(gene_files) as f:
	for line in f:
		if not '\t' in line:
			continue
		num_phases = len(line.strip().split('\t'))
if num_phases == 0:
	print >> sys.stderr, "ERROR: The file containing targets for each time point has a blank line at the top"
	sys.exit(1)

phase_nodes = [set() for i in xrange(num_phases)]
for i in xrange(num_phases):
	genes_read = 0
	with open(gene_files) as f:
		for line in f:
			line = [k.strip() for k in line.strip().split("\t")]
			g = line[i]
			if g == "XXXX":
				continue
			genes_read += 1

			# This is for not duplicating genes across phases
			use_this = True
			for j in xrange(i):
				if g in phase_nodes[j]:
					use_this = False
					break
			if use_this:
				phase_nodes[i].add(g)
				all_phase_nodes.add(g)
			if genes_read >= phase_size:
				break

with open(expr_filepath) as f:
	for line in f:
		break
	for line in f:
		line = [k.strip() for k in line.strip().split()]
		if line[0].startswith('#') or len(line) == 0:
			continue
		if len(line) != num_phases+1:
			print >> sys.stderr, "Expression file MUST contain the log fold change expression for each time point (relative to the previous time point or control condition). If the time points have been MERGED into phases, then the log fold change expression must be for the COMBINED phase."
			sys.exit(1)
			
		gene = line[0]
		if gene not in gene_to_id:
			gene_to_id[gene] = len(gene_to_logfoldchange)
			index_to_gene += [gene]
			gene_to_logfoldchange += [[]]
		gene_index = gene_to_id[gene]

		gene_to_logfoldchange[gene_index] = [float(k) for k in line[1:]]

assert len(gene_to_logfoldchange) > 0
assert len(gene_to_logfoldchange) == len(gene_to_id)
assert len(gene_to_logfoldchange) == len(index_to_gene)

phase_clusters_asgn = []
phase_clusters = []
sizes = [[] for i in xrange(num_phases)]
for i in xrange(num_phases):
	temp = sorted([(j, gene_to_logfoldchange[j][i]) for j in xrange(len(gene_to_id)) if index_to_gene[j] in all_phase_nodes], key=lambda k : k[1])

	split_points = [int(float(len(temp))/num_clusters*c) for c in xrange(num_clusters)] + [len(temp)]

	Y = [None]*len(gene_to_id)
	for c in xrange(num_clusters):
		for j in xrange(split_points[c], split_points[c+1]):
			Y[temp[j][0]] = c

	avg = [0]*num_clusters
	n = [0]*num_clusters
	for j in xrange(len(Y)):
		if Y[j] is None:
			continue
		avg[Y[j]] += gene_to_logfoldchange[j][i]
		n[Y[j]] += 1
	phase_clusters_asgn += [Y]
	phase_clusters += [[float("inf") if n[j] == 0 else avg[j]/n[j] for j in xrange(num_clusters) if n[j] > 0]]
	sizes[i] = n

all_tf_set = set()
with open(tf_filepath) as f:
	next(f)
	for line in f:
		line = [k.strip() for k in line.strip().split("\t")]
		all_tf_set.add(line[0])


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

pp_edges_in_paths = {}
pd_edges_in_paths = {}

phase_pathways = [{} for i in xrange(num_phases)]
pd_phase_edges = [{} for i in xrange(num_phases)]

phase_source_weights = [{} for i in xrange(num_phases)]
phase_node_pp_weights = [{} for i in xrange(num_phases)]
phase_tf_weights = [{} for i in xrange(num_phases)]

all_non_targ_nodes = {}
all_pd_nodes = {}

# This is what decides the phase-specific node rankings
# via phase_source_weights, phase_node_pp_weights, and phase_tf_weights
# and thus which nodes should appear in the figure
print >> sys.stderr, "reading paths"
with open(path_filepath) as f:
	num_paths = 0
	for line in f:
		num_paths += 1
		if num_paths%1000000 == 0:
			print >> sys.stderr, num_paths/1000000, "million paths done"
		line = [k.strip() for k in line.strip().split("\t")]
		score = float(line[0])
		path = line[2:]
		path_phase = -1
		for i in xrange(num_phases):
			if path[0] in phase_nodes[i]:
				path_phase = i
				break
		if __debug__:
			for i in xrange(num_phases):
				if path_phase != i:
					assert path[0] not in phase_nodes[i]
		if path_phase < 0:
			continue
		last_phase_gene = len(path)-1
		if path_phase > 0:
			for i in reversed(xrange(1, len(path)-1)):
				if path[i] in phase_nodes[path_phase-1]:
					last_phase_gene = i
					break
		assert last_phase_gene > 0
		assert path_phase == 0 or last_phase_gene < len(path)-1, str(path)

		for g in path[1:]:
			if g not in all_non_targ_nodes:
				all_non_targ_nodes[g] = len(all_non_targ_nodes)
		if path[0] not in all_pd_nodes:
			all_pd_nodes[path[0]] = len(all_pd_nodes)+100000

		# intermediate nodes
		for i in xrange(1, last_phase_gene+1):
			g = path[i]
			if g not in phase_node_pp_weights[path_phase]:
				phase_node_pp_weights[path_phase][g] = 0
			phase_node_pp_weights[path_phase][g] += score

		# source
		source = path[last_phase_gene]
		if source not in phase_source_weights[path_phase]:
			phase_source_weights[path_phase][source] = 0
		phase_source_weights[path_phase][source] += score

		# tf
		if path[1] not in phase_tf_weights[path_phase]:
			phase_tf_weights[path_phase][path[1]] = 0
		phase_tf_weights[path_phase][path[1]] += score

		# record edges used in paths and also pairwise 
		# path flow from one gene to another
		for i in reversed(xrange(2, len(path))):
			s = path[i]
			t = path[i-1]
			if s not in pp_edges_in_paths:
				pp_edges_in_paths[s] = set()
			pp_edges_in_paths[s].add(t)

		for i in reversed(xrange(2, last_phase_gene+1)):
			s = path[i]
			for j in reversed(xrange(1, i)):
				t = path[j]
				if s not in phase_pathways[path_phase]:
					phase_pathways[path_phase][s] = {}
				if t not in phase_pathways[path_phase][s]:
					phase_pathways[path_phase][s][t] = 0
				phase_pathways[path_phase][s][t] += score

		# record overall and phase specific last pd pathways and edges
		t = path[0]
		for i in xrange(1, last_phase_gene+1):
			s = path[i]
			if s not in pd_edges_in_paths:
				pd_edges_in_paths[s] = set()
			pd_edges_in_paths[s].add(t)
			if s not in pd_phase_edges[path_phase]:
				pd_phase_edges[path_phase][s] = {}
			if t not in pd_phase_edges[path_phase][s]:
				pd_phase_edges[path_phase][s][t] = 0
			pd_phase_edges[path_phase][s][t] += score

imp_genes = [{} for i in xrange(num_phases)]
with open(important_genes_filepath) as f:
	for line in f:
		cur_phase = -1
		for line in f:
			line_split = [k.strip() for k in line.strip().split()]
			if "Phase" in line or len(line_split) == 0:
				continue
			if "==" in line:
				assert len(line_split) == 1
				cur_phase += 1
				continue
			if cur_phase == 0:
				imp_genes[cur_phase][line_split[0]] = [int(k) for k in line_split[-4:]]
			else:
				imp_genes[cur_phase][line_split[0]] = [int(k) for k in line_split[-5:]]

phase_source_ranking = [sorted([g for g in phase_source_weights[i]], key=lambda k: phase_source_weights[i][k], reverse=True) for i in xrange(num_phases)]

phase_node_pp_ranking = [sorted([g for g in phase_node_pp_weights[i]], key=lambda k: phase_node_pp_weights[i][k], reverse=True) for i in xrange(num_phases)]

phase_tf_ranking = [sorted([g for g in phase_tf_weights[i]], key=lambda k: phase_tf_weights[i][k], reverse=True) for i in xrange(num_phases)]

phase_source_nodes = set()
for g in phase_source_ranking[0]:
	if len(phase_source_nodes) >= num_nodes:
		break
	if g not in imp_genes[0]:
		continue
	if g not in all_tf_set:
		phase_source_nodes.add(g)

print >> sys.stderr, "writing graph"

#out_file = open(out_fileprefix, 'w')
#out_file.write("graph ["+"\n")
#out_file.write("directed 1"+"\n")
#out_file.write("defaultnodesize " + str(node_width)+"\n")

### Write graph file ###
edge_id = 0
cur_x = node_width*2
cur_y = node_height*2

gml_nodes = [{} for i in xrange(num_phases)]
gml_edges = [{} for i in xrange(num_phases)]

nodes_sel = {}
fedges_sel = {}
bedges_sel = {}

group_id = 0
for phase_id in xrange(num_phases):
	intermediate_nodes = set()
	tf_nodes = set()

	# populate intermediate nodes
	last_intermediate_added = None
	for g in phase_node_pp_ranking[phase_id]:
		if len(intermediate_nodes) >= num_nodes:
			break
		if g in all_tf_set or g in phase_source_nodes or g not in imp_genes[phase_id]:
			continue
		assert g not in all_tf_set
		intermediate_nodes.add(g)
		last_intermediate_added = g
	assert len(intermediate_nodes) <= num_nodes
	assert len(intermediate_nodes.intersection(all_tf_set)) == 0

	# populate tf nodes
	phase_nodes_covered = set()
	for g in phase_tf_ranking[phase_id]:
		assert g in all_tf_set
		if len(tf_nodes) == num_nodes:
			break

		if g not in imp_genes[phase_id]:
			continue

		if g not in pd_phase_edges[phase_id]:
			continue

		new_expr_node = False
		for t in gene_to_id:
			if t in pd_phase_edges[phase_id][g] and pd_phase_edges[phase_id][g][t] <= pathway_score_threshold:
				continue
			expr_node = phase_clusters_asgn[phase_id][gene_to_id[t]]
			if expr_node not in phase_nodes_covered:
				phase_nodes_covered.add(expr_node)
				new_expr_node = True

		if len(tf_nodes) <= (num_nodes-num_phases+len(phase_nodes_covered)) or new_expr_node:
			tf_nodes.add(g)

	assert len(tf_nodes) <= num_nodes
	assert len(tf_nodes.intersection(all_tf_set)) == len(tf_nodes)

	#if phase_id == 2:
	#	tf_nodes.add("STAT6")
	#	if len(intermediate_nodes) >= num_nodes:
	#		intermediate_nodes.remove(last_intermediate_added)
	#		intermediate_nodes.add("LCK")

	
	nodes_to_process = phase_source_nodes.union(intermediate_nodes.union(tf_nodes))

	sending_nodes_pp = set()
	recv_nodes_pp = set()
	for s in nodes_to_process:
		for t in nodes_to_process:
			if s == t:
				continue
			if s not in phase_pathways[phase_id] or t not in phase_pathways[phase_id][s]:
				continue
			if phase_pathways[phase_id][s][t] <= pathway_score_threshold:
				continue
			sending_nodes_pp.add(s)
			recv_nodes_pp.add(t)

	sending_nodes_pd = set()
	recv_nodes_pd = set()
	for s in pd_phase_edges[phase_id]:
		for t in pd_phase_edges[phase_id][s]:
			if not(s in tf_nodes):
				continue
			if pd_phase_edges[phase_id][s][t] <= pathway_score_threshold:
				continue
			sending_nodes_pd.add(s)
			recv_nodes_pd.add(t)

	assert len(tf_nodes.intersection(intermediate_nodes)) == 0
	assert len(tf_nodes.intersection(phase_source_nodes)) == 0
	assert len(intermediate_nodes.intersection(phase_source_nodes)) == 0
	#assert len(tf_nodes) <= num_nodes
	#assert len(intermediate_nodes) <= num_nodes

	node_filters = [
					lambda node : node in phase_source_nodes,
					lambda node : node in intermediate_nodes,
					lambda node : node in tf_nodes]

	colors = ["red", "cyan", "green"]

	print >> sys.stderr, i, len(phase_source_nodes), len(intermediate_nodes), len(tf_nodes)

	nodes_processed_pp = set()
	for f in xrange(len(node_filters)):
		node_filter = node_filters[f]
		for node in nodes_to_process:
			if node in nodes_processed_pp:
				continue
			if not node_filter(node):
				continue

			nodes_processed_pp.add(node)

			node_id = phase_id*len(all_non_targ_nodes)+all_non_targ_nodes[node]

			gml_nodes[phase_id][node_id] = {
					"label": node,
					"borderWidth": 0,
					"group": group_id,
					"graphics": {},
					"LabelGraphics": {},
					"clusterNode" : 0,
					"sourceNode" : 1 if f == 0 else 0,
					}

			gml_nodes[phase_id][node_id]["graphics"] = {
					"x": cur_x + random.randint(0,int(node_width)),
					"y": cur_y,
					"w": node_width,
					"h": node_height,
					"fill": colors[f]
					}
			assert node_id not in nodes_sel
			nodes_sel[node_id] = phase_id

			if node in rna_hits_set:
				gml_nodes[phase_id][node_id]["graphics"]["type"] = "\"diamond\""
			elif node in phase_nodes[phase_id]:
				gml_nodes[phase_id][node_id]["graphics"]["type"] = "\"rectangle\""
			else:
				gml_nodes[phase_id][node_id]["graphics"]["type"] = "\"ellipse\""

			gml_nodes[phase_id][node_id]["LabelGraphics"] = {
					"type": "\"text\"",
					"fontSize": 3,
					}

			cur_y += node_height*2 + random.randint(0,int(node_height))
			if cur_y > max_y:
				cur_x += node_width*2
				cur_y = node_height*2 + random.randint(0,int(node_height))

		cur_x += node_width*4
		group_id += 1
		cur_y = 0

	#cur_y = (max_y-(node_height*2))/2.0

	# create gene_to_logfoldchange cluster nodes
	cluster_order = range(num_clusters)
	cluster_order.sort(key=lambda i : phase_clusters[phase_id][i], reverse=True)
	min_cluster_expr = min(phase_clusters[phase_id])
	expr_gap = max(phase_clusters[phase_id])-min_cluster_expr

	# create gene cluster nodes
	for j in cluster_order:
		node_id = 200000+len(all_pd_nodes)+phase_id*6+j

		if (phase_clusters[phase_id][j] == float("inf")):
			continue
		gml_nodes[phase_id][node_id] = {
				"borderWidth": 0,
				"group": group_id,
				"label": str(round(phase_clusters[phase_id][j], 2)),
				"graphics": {},
				"LabelGraphics": {},
				"clusterNode": 1,
				"sourceNode" : 0,
				}
		assert node_id not in nodes_sel
		nodes_sel[node_id] = phase_id

		size = expr_node_size[0]+expr_node_size[1]*(phase_clusters[phase_id][j]-min_cluster_expr)
		gml_nodes[phase_id][node_id]["graphics"] = {
				"type": "\"ellipse\"",
				"x": cur_x,
				"y": cur_y,
				"w": size,
				"h": size,
				"label_y": cur_y + size,
				"fill": "orange",
				}

		gml_nodes[phase_id][node_id]["LabelGraphics"] = {
				"type": "\"text\"",
				"fontSize": 3
				}

		cur_y += node_height*2 + random.randint(0,int(node_height))
		if cur_y > max_y:
			cur_x += node_width*2
			group_id += 1
			cur_y = node_height*2 + random.randint(0,int(node_height))

	num_pp_edges = 0
	num_pd_edges = 0
	for s in nodes_processed_pp:
		for t in nodes_processed_pp:
			if s == t:
				continue
			if s not in phase_pathways[phase_id] or t not in phase_pathways[phase_id][s]:
				continue
			if phase_pathways[phase_id][s][t] <= pathway_score_threshold:
				continue
			edge_id += 1
			num_pp_edges += 1

			source_id = phase_id*len(all_non_targ_nodes)+all_non_targ_nodes[s]
			targ_id = phase_id*len(all_non_targ_nodes)+all_non_targ_nodes[t]
			if source_id not in fedges_sel:
				fedges_sel[source_id] = {}
			if targ_id not in bedges_sel:
				bedges_sel[targ_id] = {}
			assert targ_id not in fedges_sel[source_id]
			assert source_id not in bedges_sel[targ_id]
			fedges_sel[source_id][targ_id] = [edge_id, phase_id]
			bedges_sel[targ_id][source_id] = [edge_id, phase_id]

			gml_edges[phase_id][edge_id] = {
					"source": source_id,
					"target": targ_id,
					"weight": phase_pathways[phase_id][s][t],
					"graphics": {},
					}

			gml_edges[phase_id][edge_id]["graphics"] = {
					"targetArrow": "\"standard\""
					}

			if s in pp_edges_in_paths and t in pp_edges_in_paths[s]:
				gml_edges[phase_id][edge_id]["graphics"]["edgeLineType"] = "\"EQUAL_DASH\""

	for s in nodes_processed_pp:
		not_done = range(num_clusters)
		#for t in pd_phase_edges[phase_id][s]:
		for t in gene_to_id:
			#if s not in all_tf_set:
			#	continue
			if phase_clusters_asgn[phase_id][gene_to_id[t]] not in not_done:
				continue
			#if s not in pd_phase_edges[phase_id] or t not in pd_phase_edges[phase_id][s]:
			#	continue
			if s not in pd_phase_edges[phase_id]:
				continue
			#if t not in tf_gene_edges[s]:
			#	continue
			if t not in pd_phase_edges[phase_id][s] or pd_phase_edges[phase_id][s][t] <= pathway_score_threshold:
				continue
			edge_id += 1
			num_pd_edges += 1

			source_id = phase_id*len(all_non_targ_nodes)+all_non_targ_nodes[s]
			targ_id =  200000+len(all_pd_nodes)+phase_id*6+phase_clusters_asgn[phase_id][gene_to_id[t]]
			if source_id not in fedges_sel:
				fedges_sel[source_id] = {}
			if targ_id not in bedges_sel:
				bedges_sel[targ_id] = {}
			assert targ_id not in fedges_sel[source_id]
			assert source_id not in bedges_sel[targ_id]
			fedges_sel[source_id][targ_id] = [edge_id, phase_id]
			bedges_sel[targ_id][source_id] = [edge_id, phase_id]
			edge_type = "\"EQUAL_DASH\"" if s not in all_tf_set else "\"standard\""
			gml_edges[phase_id][edge_id] = {
					"source": source_id,
					"target": targ_id,
					"weight": pd_phase_edges[phase_id][s][t] if t in pd_phase_edges[phase_id][s] else 0.9,
					"graphics": {"targetArrow": edge_type},
					}
			not_done.remove(phase_clusters_asgn[phase_id][gene_to_id[t]])

	if phase_id > 0:
		for t in nodes_processed_pp:
			not_done = range(num_clusters)
			for s in phase_source_ranking[phase_id]:
				if phase_clusters_asgn[phase_id-1][gene_to_id[s]] not in not_done:
					continue
				if s not in phase_pathways[phase_id] or t not in phase_pathways[phase_id][s]:
					continue
				if phase_pathways[phase_id][s][t] <= pathway_score_threshold:
					continue
				not_done.remove(phase_clusters_asgn[phase_id-1][gene_to_id[s]])

				edge_id += 1
				num_pd_edges += 1

				source_id = 200000+len(all_pd_nodes)+(phase_id-1)*6+phase_clusters_asgn[phase_id-1][gene_to_id[s]]
				targ_id = phase_id*len(all_non_targ_nodes)+all_non_targ_nodes[t]
				if source_id not in fedges_sel:
					fedges_sel[source_id] = {}
				if targ_id not in bedges_sel:
					bedges_sel[targ_id] = {}
				assert targ_id not in fedges_sel[source_id]
				assert source_id not in bedges_sel[targ_id]
				fedges_sel[source_id][targ_id] = [edge_id, phase_id]
				bedges_sel[targ_id][source_id] = [edge_id, phase_id]
				gml_edges[phase_id][edge_id] = {
						"source": source_id,
						"target": targ_id,
						"weight": phase_pathways[phase_id][s][t],
						"graphics": {"targetArrow": "\"standard\""}
						}


	phase_source_nodes = set()

fpaths_sel = {}
genes = list(set(fedges_sel.keys()).union(set(bedges_sel.keys())))
for i in xrange(len(genes)):
	s = genes[i]
	if s not in fpaths_sel:
		break
	for j in xrange(len(genes)):
		if i == j:
			continue
		t = genes[j]
		nl = fedges_sel[s]
		done = False
		while done:
			nl2 = []
			for g in nl:
				if g == t:
					fpaths_sel.add(t)
					done = True
					break
				nl2 += list(fedges_sel[g])
			nl = nl2

nodes_to_remove = set()
for node in nodes_sel:
	phase_id = nodes_sel[node]
	if gml_nodes[phase_id][node]["sourceNode"]:
		if node not in fedges_sel or len(fedges_sel[node]) == 0:
			nodes_to_remove.add(node)
			if node in fedges_sel:
				for targ in fedges_sel[node]:
					edge_id, phase_id = fedges_sel[node][targ]
					if edge_id in gml_edges[phase_id]:
						del gml_edges[phase_id][edge_id]
	elif gml_nodes[phase_id][node]["clusterNode"] and phase_id == num_phases-1:
		if node not in bedges_sel or len(bedges_sel[node]) == 0:
			nodes_to_remove.add(node)
			if node in bedges_sel:
				for source in bedges_sel[node]:
					edge_id, phase_id = bedges_sel[node][source]
					if edge_id in gml_edges[phase_id]:
						del gml_edges[phase_id][edge_id]
	else:
		if node not in bedges_sel or len(bedges_sel[node]) == 0 or node not in fedges_sel or len(fedges_sel[node]) == 0:
			nodes_to_remove.add(node)
			if node in fedges_sel:
				for targ in fedges_sel[node]:
					edge_id, phase_id = fedges_sel[node][targ]
					if edge_id in gml_edges[phase_id]:
						del gml_edges[phase_id][edge_id]
			if node in bedges_sel:
				for source in bedges_sel[node]:
					edge_id, phase_id = bedges_sel[node][source]
					if edge_id in gml_edges[phase_id]:
						del gml_edges[phase_id][edge_id]

for node in nodes_to_remove:
	del gml_nodes[nodes_sel[node]][node]

y_range = [[0, 0] for i in xrange(group_id+1)]
total_weight = [0]*num_phases
for phase_id in xrange(num_phases):
	for node_id,node in gml_nodes[phase_id].iteritems():
		group_id = node["group"]
		y = node["graphics"]["y"]
		if y_range[group_id][0] > y:
			y_range[group_id][0] = y
		if y_range[group_id][1] < y:
			y_range[group_id][1] = y
		assert max_y >= y_range[group_id][1]-y_range[group_id][0], str(y_range[group_id][1]) + ", " + str(yrange[group_id][0])

for phase_id in xrange(num_phases):
	for edge_id, edge in gml_edges[phase_id].iteritems():
		total_weight[phase_id] += edge["weight"]

to_add = [(max_y-(y_range[i][1]-y_range[i][0]))/2. for i in xrange(len(y_range))]

with open(path_filepath + ".graph.gml", 'w') as f:
	print >> f, "graph ["
	print >> f, "directed 1"
	print >> f, "defaultnodesize " + str(node_width)
	for phase_id in xrange(num_phases):
		for node_id,node in gml_nodes[phase_id].iteritems():
			print >> f, "node ["
			print >> f, "id " + str(node_id)
			print >> f, "borderWidth " + str(node["borderWidth"])
			print >> f, "TrueLabel \"" + node["label"] + "\""

			g = node["graphics"]
			y = g["y"] + to_add[node["group"]]
			print >> f, "graphics ["
			print >> f, "type " + g["type"]
			print >> f, "x " + str(g["x"])
			print >> f, "y " + str(y)
			print >> f, "w " + str(g["w"])
			print >> f, "h " + str(g["h"])
			print >> f, "fill \"" + g["fill"] + "\""
			print >> f, "]"

			lg = node["LabelGraphics"]
			print >> f, "LabelGraphics ["
			print >> f, "type " + lg["type"]
			print >> f, "color " "\"black\""
			if "label_y" in g:
				print >> f, "y " + str(g["label_y"] + to_add[node["group"]])
			print >> f, "fontSize " + str(lg["fontSize"])
			print >> f, "]"
			print >> f, "]"

	for phase_id in xrange(num_phases):
		for edge_id, edge in gml_edges[phase_id].iteritems():
			print >> f, "edge ["
			print >> f, "id " + str(edge_id)
			print >> f, "source " + str(edge["source"])
			print >> f, "target " + str(edge["target"])
			print >> f, "weight " + str(edge["weight"]/total_weight[phase_id]*10000)
			print >> f, "graphics ["
			if "edgeLineType" in edge["graphics"]:
				print >> f, "edgeLineType " + edge["graphics"]["edgeLineType"]
			print >> f, "targetArrow " + edge["graphics"]["targetArrow"]
			print >> f, "]"
			print >> f, "]"

	print >> f, "]"


	print >> sys.stderr, "# pp edges", str(num_pp_edges)
	print >> sys.stderr, "# pd edges", str(num_pd_edges)

