import sys,os

inp = open('Psa_all_hits.blast','r')
inl = inp.readline()
D = {}
while inl:
	if 'aG' in inl:
		tem = inl.strip().split('\t')
		if 'aG' in tem[0] and 'aG' in tem[1] and tem[0] != tem[1]:
			gene1 = tem[0]
			gene2 = tem[1]
			if gene1 not in D:
				D[gene1] = gene2
	inl = inp.readline()

set3 = ['02', '06', '08', '09', '13', '18']
out = open('Reciprocal_best_hist_set3_and_set4_genes.txt','w')
out2 = open('Reciprocal_best_hist_set3_and_set3_genes.txt','w')
out3 = open('Reciprocal_best_hist_set4_and_set4_genes.txt','w')
P = {}
P2 = {}
P3 = {}
for gene1 in D:
	if D[gene1] in D:
		gene2 = D[gene1]
		if D[gene2] == gene1:
			chr1 = gene1[4:6]
			chr2 = gene2[4:6]
			if chr1 in set3 and chr2 not in set3:
				if gene1 + '__' + gene2 not in P:
					P[gene1 + '__' + gene2] = 1
					out.write('%s\t%s\n'%(gene1, gene2))
			if chr1 not in set3 and chr2 in set3:
				if gene2 + '__' + gene1 not in P:
					P[gene2 + '__' + gene1] = 1
					out.write('%s\t%s\n'%(gene2, gene1))
			if chr1 in set3 and chr2 in set3 and chr1 != chr2:
				if min(gene1,gene2) + '__' + max(gene1,gene2) not in P2:
					P2[min(gene1,gene2) + '__' + max(gene1,gene2)] = 1
					out2.write('%s\t%s\n'%(gene1, gene2))
			if chr1 not in set3 and chr2 not in set3 and chr1 != chr2:
				if min(gene1,gene2) + '__' + max(gene1,gene2) not in P3:
					P3[min(gene1,gene2) + '__' + max(gene1,gene2)] = 1
					out3.write('%s\t%s\n'%(gene1, gene2))


out.close()
out2.close()
out3.close()