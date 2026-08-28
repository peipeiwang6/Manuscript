import sys,os
import pandas as pd
import numpy as np
import pickle

Allele_type = pd.read_csv('Psa_syri_allelic_fillout_using_McScanX_and_inversion_20250328_combine_tandem_duplicates_six_types.txt',header=0, sep='\t', low_memory=False)

T = {'With_biased_gene_loss':'Type2','With_biased_tandem_duplicates':'Type3','Conserved_with_same_tandem_duplicates':'Type4','Conserved_but_with_domain_variation':'Type5','Conserved':'Type6'}

Keep_pairs = {}
Ds = {}
Dn = {}
Omega = {}
for i in range(0, Allele_type.shape[0]):
	if Allele_type.iloc[i,5] != 'Singleton':
		gene_list = []
		for j in range(0,4):
			if type(Allele_type.iloc[i,j]) != float:
				genes = Allele_type.iloc[i,j].split(',')
				for g in genes:
					gene_list.append(g)
		if gene_list[0] not in Keep_pairs:
			Keep_pairs[gene_list[0]] = []
			if len(Keep_pairs) % 1000 == 0:
				print('Have done for %s groups!'%len(Keep_pairs))
		gene_pairs = []
		for m in range(0,len(gene_list)-1):
			for n in range(m+1,len(gene_list)):
				if gene_list[m][0:7] != gene_list[n][0:7]:
					gene_pairs.append([gene_list[m],gene_list[n]])
		Chr_pair = {}
		for pair in gene_pairs:
			gene1 = pair[0]
			gene2 = pair[1]
			genepair = min(pair) + '_' + max(pair)
			chrpair = min(gene1[6], gene2[6]) + '_' + max(gene1[6], gene2[6]) 
			if chrpair not in Chr_pair:
				Chr_pair[chrpair] = {}
			path = 'KaKs_all_allele_types/%s/%s/%s'%(T[Allele_type.iloc[i,5]], gene_list[0], genepair)
			file = open('%s/%s_raml.out'%(path, genepair),'r').readlines()
			for k in range(0, len(file)):
				if file[k].startswith('seq. seq.'):
					k = k + 2
					tem = file[k].split()
					omega = float(tem[6])
					dn = float(tem[7])
					ds = float(tem[10])
					Chr_pair[chrpair][ds] = [genepair, dn, omega]
		ds = []
		dn = []
		omega = []
		for chrpair in Chr_pair:
			minds = min(Chr_pair[chrpair].keys())
			Keep_pairs[gene_list[0]].append(Chr_pair[chrpair][minds][0])
			ds.append(minds)
			dn.append(Chr_pair[chrpair][minds][1])
			omega.append(Chr_pair[chrpair][minds][2])
		Ds[gene_list[0]] = [np.mean(ds), T[Allele_type.iloc[i,5]]]
		Dn[gene_list[0]] = [np.mean(dn), T[Allele_type.iloc[i,5]]]
		Omega[gene_list[0]] = [np.mean(omega), T[Allele_type.iloc[i,5]]]


f = open("Psa_keep_pairs_for_Ks_exp_dis.pkl","wb")
pickle.dump(Keep_pairs,f)
f.close()

df_Ds = pd.DataFrame.from_dict(Ds, orient='index')
df_Ds.to_csv('Psa_allele_mean_Ds_for_all_types.txt', index=True, header=True,sep="\t")

df_Dn = pd.DataFrame.from_dict(Dn, orient='index')
df_Dn.to_csv('Psa_allele_mean_Dn_for_all_types.txt', index=True, header=True,sep="\t")

df_Omega = pd.DataFrame.from_dict(Omega, orient='index')
df_Omega.to_csv('Psa_allele_mean_Omega_for_all_types.txt', index=True, header=True,sep="\t")


