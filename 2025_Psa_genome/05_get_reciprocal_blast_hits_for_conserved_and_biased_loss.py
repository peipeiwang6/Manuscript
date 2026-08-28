import sys,os
import pandas as pd
import numpy as np
import math
import pickle

alle = pd.read_csv('Psa_allelic_genes_final_types_20250328.txt',header=0, index_col=None, sep='\t')
alle = alle.fillna('')
alle['Representative'] = ''
for i in range(0, alle.shape[0]):
	for j in range(0, 4):
		if alle.iloc[i,j] != '':
			gene = alle.iloc[i,j].split(',')[0]
			alle.iloc[i,6] = gene
			break

inp = open('Psa_all_hits.blast','r')
inl = inp.readline()
D = {}
while inl:
	tem = inl.strip().split('\t')
	gene1 = tem[0]
	gene2 = tem[1]
	if gene1 not in D:
		if gene1 != gene2:
			if len(np.where(alle['Representative'] == gene1)[0]) > 0:
				loc = np.where(alle['Representative'] == gene1)[0][0]
				if gene2 in alle['Representative'].tolist():
					D[gene1] = gene2
	inl = inp.readline()

f = open("Psa_reciprocal_best_match_for_representative.pkl","wb")
pickle.dump(D,f)
f.close()

# load dictionary
with open('Psa_reciprocal_best_match_for_representative.pkl', 'rb') as f1:
	D = pickle.load(f1)

def type_classification(Type):
	if Type in ['0/0/1', '0/0/0/1']:
		alle_type = 'Singleton'
	elif len(list(set(Type.split('/')))) == 1:
		if Type.split('/')[0] == '1':
			alle_type = 'Conserved'
		else:
			if len(list(set(Type.split('/')))) == 1:
				alle_type = 'Conserved_with_same_tandem_duplicates'
			else:
				alle_type = 'With_biased_tandem_duplicates'
	elif '0' in Type.split('/'):
		alle_type = 'With_biased_gene_loss'
	else:
		alle_type = 'With_biased_tandem_duplicates'
	return(alle_type)

out = open('Reciprocal_best_hist_for_different_allelic_types.txt','w')
out.write('Gene1\tGene2\tType1\tType2\tPair_type\n')
R = {}
for gene1 in D:
	if D[gene1] in D:
		gene2 = D[gene1]
		if D[gene2] == gene1:
			type1 = type_classification(alle.iloc[np.where(alle['Representative'] == gene1)[0][0],5])
			type2 = type_classification(alle.iloc[np.where(alle['Representative'] == gene2)[0][0],5])
			if min(gene1, gene2) + '_' + max(gene1, gene2) not in R:
				R[min(gene1, gene2) + '_' + max(gene1, gene2)] = 1
				out.write('%s\t%s\t%s\t%s\t%s/%s\n'%(gene1, gene2,alle.iloc[np.where(alle['Representative'] == gene1)[0][0],5], alle.iloc[np.where(alle['Representative'] == gene2)[0][0],5], type1, type2))
				out.flush()


out.close()
