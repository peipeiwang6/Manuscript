import sys,os
import pandas as pd
import numpy as np
import pickle


# Syri information
Syri_AB3 = pd.read_csv('Set3_A_Bsyri.out', header=None, sep = '\t', low_memory=False)
subSyri_AB3 = Syri_AB3[Syri_AB3[10].isin(['SYNAL', 'INVAL',  'TRANSAL', 'INVTRAL'])]
Syri_AB4 = pd.read_csv('Set4_A_Bsyri.out', header=None, sep = '\t', low_memory=False)
subSyri_AB4 = Syri_AB4[Syri_AB4[10].isin(['SYNAL', 'INVAL',  'TRANSAL', 'INVTRAL'])]
Syri_AB = pd.concat([subSyri_AB3, subSyri_AB4], axis=0)
Syri_AB.to_csv('Psa_Syri_AB.txt', index=False, header=True, sep='\t')

Syri_BC3 = pd.read_csv('Set3_B_Csyri.out', header=None, sep = '\t', low_memory=False)
subSyri_BC3 = Syri_BC3[Syri_BC3[10].isin(['SYNAL', 'INVAL',  'TRANSAL', 'INVTRAL'])]
Syri_BC4 = pd.read_csv('Set4_B_Csyri.out', header=None, sep = '\t', low_memory=False)
subSyri_BC4 = Syri_BC4[Syri_BC4[10].isin(['SYNAL', 'INVAL',  'TRANSAL', 'INVTRAL'])]
Syri_BC = pd.concat([subSyri_BC3, subSyri_BC4], aPsa_genes_in_syntenic_regionsxis=0)
Syri_BC.to_csv('Psa_Syri_BC.txt', index=False, header=True, sep='\t')

Syri_CD = pd.read_csv('Set4_C_Dsyri.out', header=None, sep = '\t', low_memory=False)
Syri_CD = Syri_CD[Syri_CD[10].isin(['SYNAL', 'INVAL',  'TRANSAL', 'INVTRAL'])]
Syri_CD.to_csv('Psa_Syri_CD.txt', index=False, header=True, sep='\t')

#####################
Syri_AB3 = pd.read_csv('Set3_A_Bsyri.out', header=None, sep = '\t', low_memory=False)
subSyri_AB3 = Syri_AB3[Syri_AB3[10].isin(['SYNAL', 'INVAL',  'TRANSAL', 'INVTRAL','DUPAL','INVDPAL'])]
Syri_AB4 = pd.read_csv('Set4_A_Bsyri.out', header=None, sep = '\t', low_memory=False)
subSyri_AB4 = Syri_AB4[Syri_AB4[10].isin(['SYNAL', 'INVAL',  'TRANSAL', 'INVTRAL','DUPAL','INVDPAL'])]
Syri_AB = pd.concat([subSyri_AB3, subSyri_AB4], axis=0)
Syri_AB.to_csv('Psa_Syri_AB_20250317.txt', index=False, header=True, sep='\t')

Syri_BC3 = pd.read_csv('Set3_B_Csyri.out', header=None, sep = '\t', low_memory=False)
subSyri_BC3 = Syri_BC3[Syri_BC3[10].isin(['SYNAL', 'INVAL',  'TRANSAL', 'INVTRAL','DUPAL','INVDPAL'])]
Syri_BC4 = pd.read_csv('Set4_B_Csyri.out', header=None, sep = '\t', low_memory=False)
subSyri_BC4 = Syri_BC4[Syri_BC4[10].isin(['SYNAL', 'INVAL',  'TRANSAL', 'INVTRAL','DUPAL','INVDPAL'])]
Syri_BC = pd.concat([subSyri_BC3, subSyri_BC4], axis=0)
Syri_BC.to_csv('Psa_Syri_BC_20250317.txt', index=False, header=True, sep='\t')

Syri_CD = pd.read_csv('Set4_C_Dsyri.out', header=None, sep = '\t', low_memory=False)
Syri_CD = Syri_CD[Syri_CD[10].isin(['SYNAL', 'INVAL',  'TRANSAL', 'INVTRAL','DUPAL','INVDPAL'])]
Syri_CD.to_csv('Psa_Syri_CD_20250317.txt', index=False, header=True, sep='\t')

GFF = pd.read_csv('Psa_all.gff',sep='\t',header=None)

Syri_AB = pd.read_csv('Psa_Syri_AB_20250317.txt', header=0, sep='\t')
Syri_BC = pd.read_csv('Psa_Syri_BC_20250317.txt', header=0, sep='\t')
Syri_CD = pd.read_csv('Psa_Syri_CD_20250317.txt', header=0, sep='\t')

# get overlapping Syri information for gene pair
def Get_overlaping_Syri(subSyri, Left_A, Right_A, Left_B, Right_B):
	res_A = pd.DataFrame()
	for i in range(0, subSyri.shape[0]):
		if (int(subSyri.iloc[i, 2])-Left_A)*(int(subSyri.iloc[i,1])-Right_A) <= 0:
			res_A = pd.concat([res_A,subSyri.iloc[i,:]], axis = 1)
	res_B = pd.DataFrame()
	for i in range(0, subSyri.shape[0]):
		if (int(subSyri.iloc[i, 7])-Left_B)*(int(subSyri.iloc[i,6])-Right_B) <= 0:
			res_B = pd.concat([res_B,subSyri.iloc[i,:]], axis = 1)
	if res_A.shape[0] != 0 and res_B.shape[0] != 0:
		if len(list(set(res_A.columns) & set(res_B.columns) )) > 0: # geneA and geneB are in the same syntenic block
			return(res_A.iloc[10,0])
		else:
			return('')
	else:
		return('')

D_chr = {}
for n in range(1, 20):
	for letter in ['a','b','c','d']:
		D_chr['%s_%s'%(n, letter)] = 'P%s%s'%(letter,"{:02d}".format(n))
		#D_chr[ 'P%s%s'%(letter,"{:02d}".format(n))] = '%s_%s'%(n, letter)

# for genes overlapping with syntenic regions
out = open('Psa_genes_in_syntenic_regions.txt','w')
inv_AB = Syri_AB[Syri_AB['10'] == 'SYNAL']
for i in range(0, inv_AB.shape[0]):
	subgffA = GFF[GFF[0] == D_chr[inv_AB.iloc[i,0]]]
	subgffA = subgffA[subgffA[2] < inv_AB.iloc[i,2]]
	subgffA = subgffA[subgffA[3] > inv_AB.iloc[i,1]]
	subgffB = GFF[GFF[0] == D_chr[inv_AB.iloc[i,5]]]
	subgffB = subgffB[subgffB[2] < inv_AB.iloc[i,7]]
	subgffB = subgffB[subgffB[3] > inv_AB.iloc[i,6]]
	if subgffA.shape[0] > 0 or subgffB.shape[0] > 0:
		tem = inv_AB.iloc[i,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
		out.write(tem.str.cat(sep='\t') + '\n')
		for j in range(0, subgffA.shape[0]):
			temA = subgffA.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temA.str.cat(sep='\t') + '\n')
		for j in range(0, subgffB.shape[0]):
			temB = subgffB.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temB.str.cat(sep='\t') + '\n')

inv_AB = Syri_BC[Syri_BC['10'] == 'SYNAL']
for i in range(0, inv_AB.shape[0]):
	subgffA = GFF[GFF[0] == D_chr[inv_AB.iloc[i,0]]]
	subgffA = subgffA[subgffA[2] < inv_AB.iloc[i,2]]
	subgffA = subgffA[subgffA[3] > inv_AB.iloc[i,1]]
	subgffB = GFF[GFF[0] == D_chr[inv_AB.iloc[i,5]]]
	subgffB = subgffB[subgffB[2] < inv_AB.iloc[i,7]]
	subgffB = subgffB[subgffB[3] > inv_AB.iloc[i,6]]
	if subgffA.shape[0] > 0 or subgffB.shape[0] > 0:
		tem = inv_AB.iloc[i,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
		out.write(tem.str.cat(sep='\t') + '\n')
		for j in range(0, subgffA.shape[0]):
			temA = subgffA.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temA.str.cat(sep='\t') + '\n')
		for j in range(0, subgffB.shape[0]):
			temB = subgffB.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temB.str.cat(sep='\t') + '\n')

inv_AB = Syri_CD[Syri_CD['10'] == 'SYNAL']
for i in range(0, inv_AB.shape[0]):
	subgffA = GFF[GFF[0] == D_chr[inv_AB.iloc[i,0]]]
	subgffA = subgffA[subgffA[2] < inv_AB.iloc[i,2]]
	subgffA = subgffA[subgffA[3] > inv_AB.iloc[i,1]]
	subgffB = GFF[GFF[0] == D_chr[inv_AB.iloc[i,5]]]
	subgffB = subgffB[subgffB[2] < inv_AB.iloc[i,7]]
	subgffB = subgffB[subgffB[3] > inv_AB.iloc[i,6]]
	if subgffA.shape[0] > 0 or subgffB.shape[0] > 0:
		tem = inv_AB.iloc[i,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
		out.write(tem.str.cat(sep='\t') + '\n')
		for j in range(0, subgffA.shape[0]):
			temA = subgffA.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temA.str.cat(sep='\t') + '\n')
		for j in range(0, subgffB.shape[0]):
			temB = subgffB.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temB.str.cat(sep='\t') + '\n')

out.close()

# for genes overlapping with inversion regions
out = open('Psa_genes_in_inversion_regions.txt','w')
inv_AB = Syri_AB[Syri_AB['10'] == 'INVAL']
for i in range(0, inv_AB.shape[0]):
	subgffA = GFF[GFF[0] == D_chr[inv_AB.iloc[i,0]]]
	subgffA = subgffA[subgffA[2] < inv_AB.iloc[i,2]]
	subgffA = subgffA[subgffA[3] > inv_AB.iloc[i,1]]
	subgffB = GFF[GFF[0] == D_chr[inv_AB.iloc[i,5]]]
	subgffB = subgffB[subgffB[2] < inv_AB.iloc[i,7]]
	subgffB = subgffB[subgffB[3] > inv_AB.iloc[i,6]]
	if subgffA.shape[0] > 0 or subgffB.shape[0] > 0:
		tem = inv_AB.iloc[i,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
		out.write(tem.str.cat(sep='\t') + '\n')
		for j in range(0, subgffA.shape[0]):
			temA = subgffA.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temA.str.cat(sep='\t') + '\n')
		for j in range(0, subgffB.shape[0]):
			temB = subgffB.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temB.str.cat(sep='\t') + '\n')

inv_AB = Syri_BC[Syri_BC['10'] == 'INVAL']
for i in range(0, inv_AB.shape[0]):
	subgffA = GFF[GFF[0] == D_chr[inv_AB.iloc[i,0]]]
	subgffA = subgffA[subgffA[2] < inv_AB.iloc[i,2]]
	subgffA = subgffA[subgffA[3] > inv_AB.iloc[i,1]]
	subgffB = GFF[GFF[0] == D_chr[inv_AB.iloc[i,5]]]
	subgffB = subgffB[subgffB[2] < inv_AB.iloc[i,7]]
	subgffB = subgffB[subgffB[3] > inv_AB.iloc[i,6]]
	if subgffA.shape[0] > 0 or subgffB.shape[0] > 0:
		tem = inv_AB.iloc[i,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
		out.write(tem.str.cat(sep='\t') + '\n')
		for j in range(0, subgffA.shape[0]):
			temA = subgffA.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temA.str.cat(sep='\t') + '\n')
		for j in range(0, subgffB.shape[0]):
			temB = subgffB.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temB.str.cat(sep='\t') + '\n')

inv_AB = Syri_CD[Syri_CD['10'] == 'INVAL']
for i in range(0, inv_AB.shape[0]):
	subgffA = GFF[GFF[0] == D_chr[inv_AB.iloc[i,0]]]
	subgffA = subgffA[subgffA[2] < inv_AB.iloc[i,2]]
	subgffA = subgffA[subgffA[3] > inv_AB.iloc[i,1]]
	subgffB = GFF[GFF[0] == D_chr[inv_AB.iloc[i,5]]]
	subgffB = subgffB[subgffB[2] < inv_AB.iloc[i,7]]
	subgffB = subgffB[subgffB[3] > inv_AB.iloc[i,6]]
	if subgffA.shape[0] > 0 or subgffB.shape[0] > 0:
		tem = inv_AB.iloc[i,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
		out.write(tem.str.cat(sep='\t') + '\n')
		for j in range(0, subgffA.shape[0]):
			temA = subgffA.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temA.str.cat(sep='\t') + '\n')
		for j in range(0, subgffB.shape[0]):
			temB = subgffB.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temB.str.cat(sep='\t') + '\n')

out.close()


# for genes overlapping with translocation regions
out = open('Psa_genes_in_transclocation_regions.txt','w')
inv_AB = Syri_AB[Syri_AB['10'].isin(['TRANSAL','INVTRAL'])]
for i in range(0, inv_AB.shape[0]):
	subgffA = GFF[GFF[0] == D_chr[inv_AB.iloc[i,0]]]
	subgffA = subgffA[subgffA[2] < inv_AB.iloc[i,2]]
	subgffA = subgffA[subgffA[3] > inv_AB.iloc[i,1]]
	subgffB = GFF[GFF[0] == D_chr[inv_AB.iloc[i,5]]]
	subgffB = subgffB[subgffB[2] < inv_AB.iloc[i,7]]
	subgffB = subgffB[subgffB[3] > inv_AB.iloc[i,6]]
	if subgffA.shape[0] > 0 or subgffB.shape[0] > 0:
		tem = inv_AB.iloc[i,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
		out.write(tem.str.cat(sep='\t') + '\n')
		for j in range(0, subgffA.shape[0]):
			temA = subgffA.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temA.str.cat(sep='\t') + '\n')
		for j in range(0, subgffB.shape[0]):
			temB = subgffB.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temB.str.cat(sep='\t') + '\n')

inv_AB = Syri_BC[Syri_BC['10'].isin(['TRANSAL','INVTRAL'])]
for i in range(0, inv_AB.shape[0]):
	subgffA = GFF[GFF[0] == D_chr[inv_AB.iloc[i,0]]]
	subgffA = subgffA[subgffA[2] < inv_AB.iloc[i,2]]
	subgffA = subgffA[subgffA[3] > inv_AB.iloc[i,1]]
	subgffB = GFF[GFF[0] == D_chr[inv_AB.iloc[i,5]]]
	subgffB = subgffB[subgffB[2] < inv_AB.iloc[i,7]]
	subgffB = subgffB[subgffB[3] > inv_AB.iloc[i,6]]
	if subgffA.shape[0] > 0 or subgffB.shape[0] > 0:
		tem = inv_AB.iloc[i,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
		out.write(tem.str.cat(sep='\t') + '\n')
		for j in range(0, subgffA.shape[0]):
			temA = subgffA.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temA.str.cat(sep='\t') + '\n')
		for j in range(0, subgffB.shape[0]):
			temB = subgffB.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temB.str.cat(sep='\t') + '\n')

inv_AB = Syri_CD[Syri_CD['10'].isin(['TRANSAL','INVTRAL'])]
for i in range(0, inv_AB.shape[0]):
	subgffA = GFF[GFF[0] == D_chr[inv_AB.iloc[i,0]]]
	subgffA = subgffA[subgffA[2] < inv_AB.iloc[i,2]]
	subgffA = subgffA[subgffA[3] > inv_AB.iloc[i,1]]
	subgffB = GFF[GFF[0] == D_chr[inv_AB.iloc[i,5]]]
	subgffB = subgffB[subgffB[2] < inv_AB.iloc[i,7]]
	subgffB = subgffB[subgffB[3] > inv_AB.iloc[i,6]]
	if subgffA.shape[0] > 0 or subgffB.shape[0] > 0:
		tem = inv_AB.iloc[i,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
		out.write(tem.str.cat(sep='\t') + '\n')
		for j in range(0, subgffA.shape[0]):
			temA = subgffA.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temA.str.cat(sep='\t') + '\n')
		for j in range(0, subgffB.shape[0]):
			temB = subgffB.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temB.str.cat(sep='\t') + '\n')

out.close()

# for genes overlapping with duplication regions
out = open('Psa_genes_in_duplication_regions.txt','w')
inv_AB = Syri_AB[Syri_AB['10'].isin(['DUPAL','INVDPAL'])]
for i in range(0, inv_AB.shape[0]):
	subgffA = GFF[GFF[0] == D_chr[inv_AB.iloc[i,0]]]
	subgffA = subgffA[subgffA[2] < inv_AB.iloc[i,2]]
	subgffA = subgffA[subgffA[3] > inv_AB.iloc[i,1]]
	subgffB = GFF[GFF[0] == D_chr[inv_AB.iloc[i,5]]]
	subgffB = subgffB[subgffB[2] < inv_AB.iloc[i,7]]
	subgffB = subgffB[subgffB[3] > inv_AB.iloc[i,6]]
	if subgffA.shape[0] > 0 or subgffB.shape[0] > 0:
		tem = inv_AB.iloc[i,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
		out.write(tem.str.cat(sep='\t') + '\n')
		for j in range(0, subgffA.shape[0]):
			temA = subgffA.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temA.str.cat(sep='\t') + '\n')
		for j in range(0, subgffB.shape[0]):
			temB = subgffB.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temB.str.cat(sep='\t') + '\n')

inv_AB = Syri_BC[Syri_BC['10'].isin(['DUPAL','INVDPAL'])]
for i in range(0, inv_AB.shape[0]):
	subgffA = GFF[GFF[0] == D_chr[inv_AB.iloc[i,0]]]
	subgffA = subgffA[subgffA[2] < inv_AB.iloc[i,2]]
	subgffA = subgffA[subgffA[3] > inv_AB.iloc[i,1]]
	subgffB = GFF[GFF[0] == D_chr[inv_AB.iloc[i,5]]]
	subgffB = subgffB[subgffB[2] < inv_AB.iloc[i,7]]
	subgffB = subgffB[subgffB[3] > inv_AB.iloc[i,6]]
	if subgffA.shape[0] > 0 or subgffB.shape[0] > 0:
		tem = inv_AB.iloc[i,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
		out.write(tem.str.cat(sep='\t') + '\n')
		for j in range(0, subgffA.shape[0]):
			temA = subgffA.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temA.str.cat(sep='\t') + '\n')
		for j in range(0, subgffB.shape[0]):
			temB = subgffB.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temB.str.cat(sep='\t') + '\n')

inv_AB = Syri_CD[Syri_CD['10'].isin(['DUPAL','INVDPAL'])]
for i in range(0, inv_AB.shape[0]):
	subgffA = GFF[GFF[0] == D_chr[inv_AB.iloc[i,0]]]
	subgffA = subgffA[subgffA[2] < inv_AB.iloc[i,2]]
	subgffA = subgffA[subgffA[3] > inv_AB.iloc[i,1]]
	subgffB = GFF[GFF[0] == D_chr[inv_AB.iloc[i,5]]]
	subgffB = subgffB[subgffB[2] < inv_AB.iloc[i,7]]
	subgffB = subgffB[subgffB[3] > inv_AB.iloc[i,6]]
	if subgffA.shape[0] > 0 or subgffB.shape[0] > 0:
		tem = inv_AB.iloc[i,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
		out.write(tem.str.cat(sep='\t') + '\n')
		for j in range(0, subgffA.shape[0]):
			temA = subgffA.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temA.str.cat(sep='\t') + '\n')
		for j in range(0, subgffB.shape[0]):
			temB = subgffB.iloc[j,].apply(lambda elt: str(int(elt)) if isinstance(elt, float) else str(elt))
			out.write(temB.str.cat(sep='\t') + '\n')

out.close()

import sys,os
import pandas as pd
# syntenic
inp = open('Psa_genes_in_syntenic_regions.txt','r').readlines()
i = 0
G = []
R = {}
while i < len(inp):
	if inp[i].startswith('P'):
		gene = inp[i].split('\t')[1]
		if gene not in G:
			G.append(gene)
		i += 1
		while inp[i].startswith('P'):
			gene = inp[i].split('\t')[1]
			if gene not in G:
				G.append(gene)
			i += 1
		tem = '__'.join(sorted(G))
		if tem not in R:
			R[tem] = 1
		i = i - 1
		G = []
	i += 1

blast = pd.read_csv('Psa_Vs_Psa_between_alle_chr.blast.txt',header=None, sep='\t')
gff = pd.read_csv('../06_collinear_analysis/Psa_gff_for_only_genes.txt',header=0, sep='\t')
out = open('Psa_syntenic_genes.txt','w')
for pair in R:
	genes = pair.split('__')
	for i in range(0,len(genes)-1):
		tem1 = blast[blast[0] == genes[i]]
		for j in range(i + 1, len(genes)):
			if genes[j] in tem1[1].tolist():
				out.write('%s\t%s\n'%(genes[i], genes[j]))
				out.flush()

out.close()


# inversion
inp = open('Psa_genes_in_inversion_regions.txt','r').readlines()
i = 0
G = []
R = {}
while i < len(inp):
	if inp[i].startswith('P'):
		gene = inp[i].split('\t')[1]
		if gene not in G:
			G.append(gene)
		i += 1
		while inp[i].startswith('P'):
			gene = inp[i].split('\t')[1]
			if gene not in G:
				G.append(gene)
			i += 1
		tem = '__'.join(sorted(G))
		if tem not in R:
			R[tem] = 1
		i = i - 1
		G = []
	i += 1

blast = pd.read_csv('../04_allelic_genes_based_on_updated_gff/Psa_Vs_Psa_between_alle_chr.blast.txt',header=None, sep='\t')
gff = pd.read_csv('Psa_gff_for_only_genes.txt',header=0, sep='\t')
out = open('Psa_inversed_genes.txt','w')
for pair in R:
	genes = pair.split('__')
	for i in range(0,len(genes)-1):
		tem1 = blast[blast[0] == genes[i]]
		for j in range(i + 1, len(genes)):
			if genes[j] in tem1[1].tolist():
				dir1 = gff[gff['GeneID'] == genes[i]].iloc[0,3]
				dir2 = gff[gff['GeneID'] == genes[j]].iloc[0,3]
				if dir1 != dir2:
					out.write('%s\t%s\n'%(genes[i], genes[j]))

out.close()

# translocation
inp = open('Psa_genes_in_transclocation_regions.txt','r').readlines()
i = 0
G = []
R = {}
while i < len(inp):
	if inp[i].startswith('P'):
		gene = inp[i].split('\t')[1]
		if gene not in G:
			G.append(gene)
		i += 1
		while inp[i].startswith('P'):
			gene = inp[i].split('\t')[1]
			if gene not in G:
				G.append(gene)
			i += 1
		tem = '__'.join(sorted(G))
		if tem not in R:
			R[tem] = 1
		i = i - 1
		G = []
	i += 1

blast = pd.read_csv('../04_allelic_genes_based_on_updated_gff/Psa_Vs_Psa_between_alle_chr.blast.txt',header=None, sep='\t')
gff = pd.read_csv('Psa_gff_for_only_genes.txt',header=0, sep='\t')
out = open('Psa_translocated_genes.txt','w')
for pair in R:
	genes = pair.split('__')
	for i in range(0,len(genes)-1):
		tem1 = blast[blast[0] == genes[i]]
		for j in range(i + 1, len(genes)):
			if genes[j] in tem1[1].tolist():
				dir1 = gff[gff['GeneID'] == genes[i]].iloc[0,3]
				dir2 = gff[gff['GeneID'] == genes[j]].iloc[0,3]
				if dir1 != dir2:
					out.write('%s\t%s\tinversed_translocated\n'%(genes[i], genes[j]))
				else:
					out.write('%s\t%s\ttranslocated\n'%(genes[i], genes[j]))

out.close()


# duplication
inp = open('Psa_genes_in_duplication_regions.txt','r').readlines()
i = 0
G = []
R = {}
while i < len(inp):
	if inp[i].startswith('P'):
		gene = inp[i].split('\t')[1]
		if gene not in G:
			G.append(gene)
		i += 1
		while inp[i].startswith('P'):
			gene = inp[i].split('\t')[1]
			if gene not in G:
				G.append(gene)
			i += 1
		tem = '__'.join(sorted(G))
		if tem not in R:
			R[tem] = 1
		i = i - 1
		G = []
	i += 1

blast = pd.read_csv('../04_allelic_genes_based_on_updated_gff/Psa_Vs_Psa_between_alle_chr.blast.txt',header=None, sep='\t')
gff = pd.read_csv('Psa_gff_for_only_genes.txt',header=0, sep='\t')
out = open('Psa_duplication_genes.txt','w')
for pair in R:
	genes = pair.split('__')
	for i in range(0,len(genes)-1):
		tem1 = blast[blast[0] == genes[i]]
		for j in range(i + 1, len(genes)):
			if genes[j] in tem1[1].tolist():
				dir1 = gff[gff['GeneID'] == genes[i]].iloc[0,3]
				dir2 = gff[gff['GeneID'] == genes[j]].iloc[0,3]
				if dir1 != dir2:
					out.write('%s\t%s\tinversed_duplication\n'%(genes[i], genes[j]))
				else:
					out.write('%s\t%s\tduplication\n'%(genes[i], genes[j]))

out.close()
