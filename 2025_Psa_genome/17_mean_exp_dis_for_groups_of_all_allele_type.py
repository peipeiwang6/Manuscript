import sys,os
import pandas as pd
import numpy as np
import pickle
import math

with open('Psa_keep_pairs_for_Ks_exp_dis.pkl', 'rb') as f1:
	Keep_pairs = pickle.load(f1)

Exp = pd.read_csv('Psa_normalized_mean_TPM_clean_20250217.csv',header=0, sep=',', low_memory=False)

Allele_type = pd.read_csv('Psa_syri_allelic_fillout_using_McScanX_and_inversion_20250328_combine_tandem_duplicates_six_types.txt',header=0, sep='\t', low_memory=False)

T = {'With_biased_gene_loss':'Type2','With_biased_tandem_duplicates':'Type3','Conserved_with_same_tandem_duplicates':'Type4','Conserved_but_with_domain_variation':'Type5','Conserved':'Type6'}

Dis = {}

for i in range(0, Allele_type.shape[0]):
	if Allele_type.iloc[i,5] != 'Singleton':
		gene_list = []
		for j in range(0,4):
			if type(Allele_type.iloc[i,j]) != float:
				genes = Allele_type.iloc[i,j].split(',')
				for g in genes:
					gene_list.append(g)
		dis = []
		for m in range(0, len(gene_list)-1):
			for n in range(m, len(gene_list)):
				gene1 = gene_list[m]
				gene2 = gene_list[n]
				pair = min(gene1, gene2) + '_' + max(gene1, gene2)
				if pair in Keep_pairs[gene_list[0]]:
					gene1_exp = Exp[Exp['Unnamed: 0']==gene1].iloc[0,1:Exp.shape[1]].tolist()
					gene2_exp = Exp[Exp['Unnamed: 0']==gene2].iloc[0,1:Exp.shape[1]].tolist()
					dis.append(math.dist(gene1_exp, gene2_exp))
		Dis[gene_list[0]] = [np.mean(dis), T[Allele_type.iloc[i,5]]]

df_Dis = pd.DataFrame.from_dict(Dis, orient='index')
df_Dis.to_csv('Psa_allele_mean_exp_dist_for_all_types.txt', index=True, header=True,sep="\t") #Euclidean
