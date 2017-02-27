
import sys
import gzip

if len(sys.argv) < 2:
	print >> sys.stderr, "Usage:", "<sorted path file>"
	sys.exit(1)

path_file = sys.argv[1]
num_top_paths = 0

path_file_obj = None
if path_file.endswith(".gz"):
    path_file_obj = gzip.open(path_file)
else:
    path_file_obj = open(path_file)

gene_scores = {}
i = 0
for line in path_file_obj:
    i += 1
    if i > num_top_paths and num_top_paths > 0:
        break
    if i%1000000 == 0:
        print >> sys.stderr, i/1000000, "million paths done"

    line = [k.strip() for k in line.strip().split("\t")]
    score = float(line[0])
    path_genes = line[2:]
    for gene in path_genes:
        if gene not in gene_scores:
            gene_scores[gene] = 0
        gene_scores[gene] += score

path_file_obj.close()


for gene in gene_scores:
	print gene + "\t" + str(gene_scores[gene])
