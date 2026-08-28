import sys,os
import pandas as pd
import numpy as np
import pickle

# distinguish collinear blocks between allelic chromosomes or different chromosomes or within the same chromosome
inp = open('Psa_all.collinearity','r').readlines()
Alle = {}
Diff = {}
Same = {}
i = 0
while i < len(inp):
	if not inp[i].startswith('#'):
		tem = inp[i].split()
		align = int(tem[0].split('-')[0])
		gene1 = inp[i].split('\t')[1]
		gene2 = inp[i].split('\t')[2]
		chr1 = gene1.split('G')[0]
		chr2 = gene2.split('G')[0]
		if chr1[:-1] == chr2[:-1] and chr1 != chr2:  # belong to different allelic chromsomes but not the same chromosome
			if align not in Alle:
				Alle[align] = []
			Alle[align].append([gene1,gene2])
		if chr1 == chr2 : # within the same chromosome
			if align not in Same:
				Same[align] = []
			Same[align].append([gene1,gene2])
		if chr1[:-1] != chr2[:-1]  : # belong to different chromosomes
			if align not in Diff:
				Diff[align] = []
			Diff[align].append([gene1,gene2])
	i = i + 1

D = {}  # save the number of syntenic gene pairs, and proportion of syntenic gene pairs over the number of all genes in the syntenic block
GFF = pd.read_csv('Psa_all.gff',sep='\t',header=None)
for align in Alle:
	df = pd.DataFrame(Alle[align])
	index1_1 = np.where(GFF[1] == df.iloc[0, 0])[0][0]
	index1_2 = np.where(GFF[1] == df.iloc[-1, 0])[0][0]
	subgff1 = GFF.iloc[index1_1: (index1_2 + 1), :]
	chr1 = df.iloc[0,0][4:7]
	index2_1 = np.where(GFF[1] == df.iloc[0, 1])[0][0]
	index2_2 = np.where(GFF[1] == df.iloc[-1, 1])[0][0]
	subgff2 = GFF.iloc[index2_1: (index2_2 + 1), :]
	chr2 = df.iloc[0,1][4:7]
	chr_combination = min(chr1, chr2) + '_' + max(chr1, chr2)
	n = df.shape[0]
	prop = n/(subgff1.shape[0] + subgff2.shape[0])
	for gene in df[0]:
		if gene not in D:
			D[gene] = {}
		if chr_combination not in D[gene]:
			D[gene][chr_combination] = {}
		if align not in D[gene][chr_combination]:
			D[gene][chr_combination][align] = [n, prop]
	for gene in df[1]:
		if gene not in D:
			D[gene] = {}
		if chr_combination not in D[gene]:
			D[gene][chr_combination] = {}
		if align not in D[gene][chr_combination]:
			D[gene][chr_combination][align] = [n, prop]

S = {} # save the selected syntenic block
for gene in D:
	for chr_combination in D[gene]:
		if len(D[gene][chr_combination]) > 1:
			tem = pd.DataFrame.from_dict(D[gene][chr_combination], orient='index')
			tem[2] = tem[0]*tem[1]
			#winner = np.where(tem.iloc[0,:] == np.max(tem.iloc[0,:]))[0][0] # based on number of syntenic pairs
			winner = np.where(tem.iloc[:,2] == np.max(tem.iloc[:,2]))[0][0] # based on proportion of syntenic pairs
			winner = tem.index[winner]
			for align in D[gene][chr_combination]:
				if align not in S:
					S[align] = [0,0,0,0]
					S[align][2] = D[gene][chr_combination][align][0]
					S[align][3] = D[gene][chr_combination][align][1]
				if align == winner:
					S[align][0] += 1
				else:
					S[align][1] += 1
		else:
			for align in D[gene][chr_combination]:
				if align not in S:
					S[align] = [1,0,0,0]
					S[align][2] = D[gene][chr_combination][align][0]
					S[align][3] = D[gene][chr_combination][align][1]


res = pd.DataFrame.from_dict(S, orient='index')
align_winner = res[res[1] == 0]
winner_02 = res[res[0] > res[1]]
winner_02 = winner_02[winner_02[1] != 0]
winner_02['Combined'] = winner_02[2] * winner_02[3]
winner_02 = winner_02.sort_values('Combined', ascending = False)
winner_02 = winner_02[winner_02['Combined'] > 40]
align_winner = pd.concat([align_winner, winner_02], axis=0)
res.to_csv('Psa_Collinear_blocks_allelic_syntenic_proportion.txt', index=True, header=True, sep='\t')
align_winner.to_csv('Psa_Collinear_blocks_allelic_winner.txt', index=True, header=True, sep='\t')


# make dataframe for collinear blocks
Collinear = pd.DataFrame()
for align in Alle:
	tem_df = pd.DataFrame({'Align':align, 'Type':'allelic_chr'}, index=range(0, len(Alle[align])))
	tem_df_02 = pd.DataFrame(Alle[align])
	tem = pd.concat([tem_df, tem_df_02], axis=1, ignore_index=True)
	Collinear = pd.concat([Collinear, tem])

for align in Same:
	tem_df = pd.DataFrame({'Align':align, 'Type':'within_same'}, index=range(0, len(Same[align])))
	tem_df_02 = pd.DataFrame(Same[align])
	tem = pd.concat([tem_df, tem_df_02], axis=1, ignore_index=True)
	Collinear = pd.concat([Collinear, tem])

for align in Diff:
	tem_df = pd.DataFrame({'Align':align, 'Type':'between_diff'}, index=range(0, len(Diff[align])))
	tem_df_02 = pd.DataFrame(Diff[align])
	tem = pd.concat([tem_df, tem_df_02], axis=1, ignore_index=True)
	Collinear = pd.concat([Collinear, tem])

Collinear.columns = ['Align_ID' , 'Type', 'Gene1', 'Gene2']
Collinear.to_csv('Psa_Collinear_blocks.txt', index=False, header=True, sep='\t')
# make a dictionary for tandem duplicate groups
inp = open('Psa_all.tandem_more.txt','r').readlines()
Group = {}
n = 0
for inl in inp:
	tem = inl.strip().split(',')
	gene1 = tem[0]
	gene2 = tem[1]
	if n == 0:
		Group[n] = []
		Group[n].append(gene1)
		Group[n].append(gene2)
		n = n + 1
	else:
		exist = 0
		for group in Group:
			if gene1 in Group[group] and gene2 not in Group[group]:
				Group[group].append(gene2)
				exist = 1
			if gene1 not in Group[group] and gene2 in Group[group]:
				Group[group].append(gene1)
				exist = 1
		if exist == 0:
			Group[n] = []
			Group[n].append(gene1)
			Group[n].append(gene2)	
			n = n + 1	


# make dictionary for tandem duplicates
Tandem = {}
for n in Group:
	for gene in Group[n]:
		if gene not in Tandem:
			Tandem[gene] = []
		for other in Group[n]:
			if other != gene and other not in Tandem[gene]:
				Tandem[gene].append(other)

f = open("Psa_tandem_genes.pkl","wb")
pickle.dump(Tandem,f)
f.close()

