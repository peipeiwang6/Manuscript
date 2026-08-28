import sys,os
import pandas as pd
import numpy as np
import pickle

syntenic = pd.read_csv('Psa_syntenic_genes.txt',header=None, sep='\t')  # results from 02_get_genes_in_syntenic_inversion_translocation_regions_based_on_syri.py
syntenic[2] = 'SYNAL'
Syrigenes = syntenic
Syrigenes = Syrigenes.drop_duplicates()
Syrigenes['Chr'] = Syrigenes[0].str[4:6:]


def Order_allelic(df):
	df_ordered = pd.DataFrame()
	df_ordered = pd.concat([df_ordered, df[df['GeneA'].notna()]])
	df_ordered = df_ordered.sort_values(by='GeneA')
	df_GeneA_na = df[df['GeneA'].isna()]
	df_GeneB_notna = df_GeneA_na[df_GeneA_na['GeneB'].notna()]
	df_GeneB_notna = df_GeneB_notna.sort_values(by='GeneB')
	df_GeneB_na = df_GeneA_na[df_GeneA_na['GeneB'].isna()]
	df_GeneC_notna = df_GeneB_na[df_GeneB_na['GeneC'].notna()]
	df_GeneC_notna = df_GeneC_notna.sort_values(by='GeneC')
	df_GeneC_na = df_GeneB_na[df_GeneB_na['GeneC'].isna()]
	df_GeneD_notna = df_GeneC_na[df_GeneC_na['GeneD'].notna()]
	df_GeneD_notna = df_GeneD_notna.sort_values(by='GeneD')
	# insert B genes
	for i in range(0, df_GeneB_notna.shape[0]):
		try:
			loc = np.where(df_ordered['GeneB'] > df_GeneB_notna.iloc[i,1])[0][0]
			df_ordered = pd.concat([df_ordered.iloc[0:loc, :], df_GeneB_notna.iloc[[i], :], df_ordered.iloc[loc:, :]], axis=0).reset_index(drop=True)
		except:
			loc = df_ordered.shape[0]
			df_ordered = pd.concat([df_ordered.iloc[0:loc, :], df_GeneB_notna.iloc[[i], :]], axis=0).reset_index(drop=True)
	# insert C genes
	for i in range(0, df_GeneC_notna.shape[0]):
		try:
			loc = np.where(df_ordered['GeneC'] > df_GeneC_notna.iloc[i,2])[0][0]
			df_ordered = pd.concat([df_ordered.iloc[0:loc, :], df_GeneC_notna.iloc[[i], :], df_ordered.iloc[loc:, :]], axis=0).reset_index(drop=True)
		except:
			loc = df_ordered.shape[0]
			df_ordered = pd.concat([df_ordered.iloc[0:loc, :], df_GeneC_notna.iloc[[i], :]], axis=0).reset_index(drop=True)
	# insert D genes
	for i in range(0, df_GeneD_notna.shape[0]):
		try:
			loc = np.where(df_ordered['GeneD'] > df_GeneD_notna.iloc[i,3])[0][0]
			df_ordered = pd.concat([df_ordered.iloc[0:loc, :], df_GeneD_notna.iloc[[i], :], df_ordered.iloc[loc:, :]], axis=0).reset_index(drop=True)
		except:
			loc = df_ordered.shape[0]
			df_ordered = pd.concat([df_ordered.iloc[0:loc, :], df_GeneD_notna.iloc[[i], :]], axis=0).reset_index(drop=True)
	return(df_ordered)


def Allelic_finder_set3(chra, chrb, chrc, subsyri, i_a, i_b, i_c, increasing, remaining):
	if remaining == 'n':
		geneA = chra.iloc[i_a:(i_a + increasing), 1]
		geneB = chrb.iloc[i_b:(i_b + increasing), 1]
		geneC = chrc.iloc[i_c:(i_c + increasing), 1]
	else:
		geneA = chra.iloc[i_a:chra.shape[0], 1]
		geneB = chrb.iloc[i_b:chrb.shape[0], 1]
		geneC = chrc.iloc[i_c:chrc.shape[0], 1]
	temsyriAB = subsyri[subsyri[0].isin(geneA.tolist()) & subsyri[1].isin(geneB.tolist())]
	temsyriAB = temsyriAB.merge(pd.DataFrame(geneA), how='outer', left_on=0, right_on=1)
	temsyriBC = subsyri[subsyri[0].isin(geneB.tolist()) & subsyri[1].isin(geneC.tolist())]
	temsyriBC = temsyriBC.merge(pd.DataFrame(geneB), how='outer', left_on=0, right_on=1)
	temsyriABC = temsyriAB.loc[:,[0,'1_x',2]].merge(temsyriBC.loc[:,[0,'1_x',2]], how='outer', left_on='1_x', right_on=0)
	temsyriABC = temsyriABC.loc[:,['0_x','0_y','1_x_y','2_x','2_y']]
	temsyriABC = temsyriABC.merge(pd.DataFrame(geneC), how='outer', left_on='1_x_y', right_on=1)
	temsyriABC = temsyriABC.loc[:,['0_x','0_y','1_x_y','2_x','2_y']]
	temsyriABC['GeneD'] = np.nan
	temsyriABC['C_D'] = np.nan
	temsyriABC = temsyriABC.loc[:,['0_x','0_y','1_x_y','GeneD', '2_x','2_y','C_D']]
	temsyriABC.columns = ['GeneA', 'GeneB', 'GeneC', 'GeneD', 'A_B', 'B_C', 'C_D']
	temsyriABC = temsyriABC.drop_duplicates()
	temsyriABC = Order_allelic(temsyriABC)
	rowcount = temsyriABC.iloc[:,0:3].count(axis = 'columns')
	loc = np.where(rowcount == 3)[0][-1]
	genea = temsyriABC.iloc[loc,0]
	geneb = temsyriABC.iloc[loc,1]
	genec = temsyriABC.iloc[loc,2]
	if all(temsyriABC.iloc[(loc+1):,0].dropna() > genea) and all(temsyriABC.iloc[(loc+1):,1].dropna() > geneb) and all(temsyriABC.iloc[(loc+1):,2].dropna() > genec):
		return(temsyriABC.iloc[0:(loc+1),:], temsyriABC.iloc[loc,0], temsyriABC.iloc[loc,1], temsyriABC.iloc[loc,2])
	else:
		tosave = temsyriABC[(temsyriABC['GeneA'] <= genea) | (temsyriABC['GeneB'] <= geneb)  | (temsyriABC['GeneC'] <= genec)]
		return(tosave, temsyriABC.iloc[loc,0], temsyriABC.iloc[loc,1], temsyriABC.iloc[loc,2])

def Allelic_finder_set4(chra, chrb, chrc, chrd, subsyri, i_a, i_b, i_c, i_d, increasing, remaining):
	if remaining == 'n':
		geneA = chra.iloc[i_a:(i_a + increasing), 1]
		geneB = chrb.iloc[i_b:(i_b + increasing), 1]
		geneC = chrc.iloc[i_c:(i_c + increasing), 1]
		geneD = chrd.iloc[i_d:(i_d + increasing), 1]
	else:
		geneA = chra.iloc[i_a:chra.shape[0], 1]
		geneB = chrb.iloc[i_b:chrb.shape[0], 1]
		geneC = chrc.iloc[i_c:chrc.shape[0], 1]
		geneD = chrd.iloc[i_d:chrd.shape[0], 1]
	temsyriAB = subsyri[subsyri[0].isin(geneA.tolist()) & subsyri[1].isin(geneB.tolist())]
	temsyriAB = temsyriAB.merge(pd.DataFrame(geneA), how='outer', left_on=0, right_on=1)
	temsyriBC = subsyri[subsyri[0].isin(geneB.tolist()) & subsyri[1].isin(geneC.tolist())]
	temsyriBC = temsyriBC.merge(pd.DataFrame(geneB), how='outer', left_on=0, right_on=1)
	temsyriCD = subsyri[subsyri[0].isin(geneC.tolist()) & subsyri[1].isin(geneD.tolist())]
	temsyriCD = temsyriCD.merge(pd.DataFrame(geneC), how='outer', left_on=0, right_on=1)
	temsyriABC = temsyriAB.loc[:,[0,'1_x',2]].merge(temsyriBC.loc[:,[0,'1_x',2]], how='outer', left_on='1_x', right_on=0)
	temsyriABCD = temsyriABC.loc[:,['0_x','0_y','1_x_y', '2_x', '2_y']].merge(temsyriCD.loc[:,[0,'1_x',2]], how='outer', left_on='1_x_y', right_on=0)
	temsyriABCD = temsyriABCD.loc[:,['0_x','0_y','1_x_y','1_x','2_x','2_y',2]]
	temsyriABCD = temsyriABCD.merge(pd.DataFrame(geneD), how='outer', left_on='1_x', right_on=1)
	temsyriABCD = temsyriABCD.loc[:,['0_x','0_y','1_x_y','1_x','2_x','2_y',2]]
	temsyriABCD.columns = ['GeneA', 'GeneB', 'GeneC', 'GeneD', 'A_B', 'B_C', 'C_D']
	temsyriABCD = temsyriABCD.drop_duplicates()
	temsyriABCD = Order_allelic(temsyriABCD)
	rowcount = temsyriABCD.iloc[:,0:4].count(axis = 'columns')
	loc = np.where(rowcount == 4)[0][-1]
	genea = temsyriABCD.iloc[loc,0]
	geneb = temsyriABCD.iloc[loc,1]
	genec = temsyriABCD.iloc[loc,2]
	gened = temsyriABCD.iloc[loc,3]
	if all(temsyriABCD.iloc[(loc+1):,0].dropna() > genea) and all(temsyriABCD.iloc[(loc+1):,1].dropna() > geneb) and all(temsyriABCD.iloc[(loc+1):,2].dropna() > genec) and all(temsyriABCD.iloc[(loc+1):,3].dropna() > gened):
		return(temsyriABCD.iloc[0:(loc+1),:], temsyriABCD.iloc[loc,0], temsyriABCD.iloc[loc,1], temsyriABCD.iloc[loc,2], temsyriABCD.iloc[loc,3])
	else:
		tosave = temsyriABCD[(temsyriABCD['GeneA'] <= genea) | (temsyriABCD['GeneB'] <= geneb)  | (temsyriABCD['GeneC'] <= genec)  | (temsyriABCD['GeneD'] <= gened)]
		return(tosave, temsyriABCD.iloc[loc,0], temsyriABCD.iloc[loc,1], temsyriABCD.iloc[loc,2], temsyriABCD.iloc[loc,3])

GFF = pd.read_csv('Psa_all.gff',sep='\t',header=None)
GFF['Chr'] = GFF[0].str[-2:]

res = pd.DataFrame()
n = 0
for Chr in GFF['Chr'].unique():
	print(Chr)
	increasing = 20
	subgff = GFF[GFF['Chr'] == Chr]
	subsyri = Syrigenes[Syrigenes['Chr']==Chr]
	chra = subgff[subgff[0]=='Pa%s'%Chr]
	chrb = subgff[subgff[0]=='Pb%s'%Chr]
	chrc = subgff[subgff[0]=='Pc%s'%Chr]
	i_a = 0
	i_b = 0
	i_c = 0
	chrset = 3
	if len(subgff[0].unique())==4:
		chrd = subgff[subgff[0]=='Pd%s'%Chr]
		i_d = 0
		chrset = 4
	if chrset == 3:
		while i_a < chra.shape[0] - increasing and i_b < chrb.shape[0] - increasing and i_c < chrc.shape[0] - increasing:
			try:
				temres, genea, geneb, genec = Allelic_finder_set3(chra, chrb, chrc, subsyri, i_a, i_b, i_c, increasing, 'n')
				res = pd.concat([res, temres], axis = 0, ignore_index=True)
				i_a = np.where(chra[1] == genea)[0][0] + 1
				i_b = np.where(chrb[1] == geneb)[0][0] + 1
				i_c = np.where(chrc[1] == genec)[0][0] + 1
				n += 1
				if n % 100 == 0:
					print(res.shape[0])
				increasing = 20
			except:
				increasing += 10
				print(increasing)
		geneA = chra.iloc[i_a:chra.shape[0], 1]
		geneB = chrb.iloc[i_b:chrb.shape[0], 1]
		geneC = chrc.iloc[i_c:chrc.shape[0], 1]
		temsyriAB = subsyri[subsyri[0].isin(geneA.tolist()) & subsyri[1].isin(geneB.tolist())]
		temsyriAB = temsyriAB.merge(pd.DataFrame(geneA), how='outer', left_on=0, right_on=1)
		temsyriBC = subsyri[subsyri[0].isin(geneB.tolist()) & subsyri[1].isin(geneC.tolist())]
		temsyriBC = temsyriBC.merge(pd.DataFrame(geneB), how='outer', left_on=0, right_on=1)
		temsyriABC = temsyriAB.loc[:,[0,'1_x',2]].merge(temsyriBC.loc[:,[0,'1_x',2]], how='outer', left_on='1_x', right_on=0)
		temsyriABC = temsyriABC.loc[:,['0_x','0_y','1_x_y','2_x','2_y']]
		temsyriABC = temsyriABC.merge(pd.DataFrame(geneC), how='outer', left_on='1_x_y', right_on=1)
		temsyriABC = temsyriABC.loc[:,['0_x','0_y','1_x_y','2_x','2_y']]
		temsyriABC['GeneD'] = np.nan
		temsyriABC['C_D'] = np.nan
		temsyriABC = temsyriABC.loc[:,['0_x','0_y','1_x_y','GeneD', '2_x','2_y','C_D']]
		temsyriABC.columns = ['GeneA', 'GeneB', 'GeneC', 'GeneD', 'A_B', 'B_C', 'C_D']
		temsyriABC = temsyriABC.drop_duplicates()
		temsyriABC = Order_allelic(temsyriABC)
		res = pd.concat([res, temsyriABC], axis = 0, ignore_index=True)
	else:
		while i_a < chra.shape[0] - increasing and i_b < chrb.shape[0] - increasing and i_c < chrc.shape[0] - increasing and i_d < chrd.shape[0] - increasing:
			try:
				temres, genea, geneb, genec, gened = Allelic_finder_set4(chra, chrb, chrc, chrd, subsyri, i_a, i_b, i_c, i_d, increasing, 'n')
				res = pd.concat([res, temres], axis = 0, ignore_index=True)
				i_a = np.where(chra[1] == genea)[0][0] + 1
				i_b = np.where(chrb[1] == geneb)[0][0] + 1
				i_c = np.where(chrc[1] == genec)[0][0] + 1
				i_d = np.where(chrd[1] == gened)[0][0] + 1
				n += 1
				if n % 100 == 0:
					print(res.shape[0])
				increasing = 20
			except:
				increasing += 10
				print(increasing)
		geneA = chra.iloc[i_a:chra.shape[0], 1]
		geneB = chrb.iloc[i_b:chrb.shape[0], 1]
		geneC = chrc.iloc[i_c:chrc.shape[0], 1]
		geneD = chrd.iloc[i_d:chrd.shape[0], 1]
		temsyriAB = subsyri[subsyri[0].isin(geneA.tolist()) & subsyri[1].isin(geneB.tolist())]
		temsyriAB = temsyriAB.merge(pd.DataFrame(geneA), how='outer', left_on=0, right_on=1)
		temsyriBC = subsyri[subsyri[0].isin(geneB.tolist()) & subsyri[1].isin(geneC.tolist())]
		temsyriBC = temsyriBC.merge(pd.DataFrame(geneB), how='outer', left_on=0, right_on=1)
		temsyriCD = subsyri[subsyri[0].isin(geneC.tolist()) & subsyri[1].isin(geneD.tolist())]
		temsyriCD = temsyriCD.merge(pd.DataFrame(geneC), how='outer', left_on=0, right_on=1)
		temsyriABC = temsyriAB.loc[:,[0,'1_x',2]].merge(temsyriBC.loc[:,[0,'1_x',2]], how='outer', left_on='1_x', right_on=0)
		temsyriABCD = temsyriABC.loc[:,['0_x','0_y','1_x_y', '2_x', '2_y']].merge(temsyriCD.loc[:,[0,'1_x',2]], how='outer', left_on='1_x_y', right_on=0)
		temsyriABCD = temsyriABCD.loc[:,['0_x','0_y','1_x_y','1_x','2_x','2_y',2]]
		temsyriABCD = temsyriABCD.merge(pd.DataFrame(geneD), how='outer', left_on='1_x', right_on=1)
		temsyriABCD = temsyriABCD.loc[:,['0_x','0_y','1_x_y','1_x','2_x','2_y',2]]
		temsyriABCD.columns = ['GeneA', 'GeneB', 'GeneC', 'GeneD', 'A_B', 'B_C', 'C_D']
		temsyriABCD = temsyriABCD.drop_duplicates()
		temsyriABCD = Order_allelic(temsyriABCD)
		res = pd.concat([res, temsyriABCD], axis = 0, ignore_index=True)

res.to_csv('Psa_Syri_allelic_genes_based_on_GFF_20250324.txt',header=True, index=False, sep='\t')

# fill out missing allelic genes using McScanX

McScanX = pd.read_csv('Psa_Collinear_blocks.txt',header=0, sep='\t',low_memory=False) # results from 01_parse_McScanX_and_tandem_information.py
McScanX = McScanX[McScanX['Type'] == 'allelic_chr']
McScanX['Chr'] = McScanX['Gene1'].str[4:6:]
McScanX = McScanX.drop('Type', axis=1)

C = {'a':0, 'b':1, 'c':2, 'd':3, 'ab':4, 'bc': 5, 'cd':6, 'ac':7, 'ad':8, 'bd':9}

def Combine_two_list(list1, list2):
	res = ['','','','','','','','','','','']
	res2 = ['','','','','','','','','','','']
	for i in range(0, len(list1)):
		if list1[i] != '' and list2[i] == '':
			res[i] = list1[i]
			res2[i] = list1[i]
		if list1[i] == '' and list2[i] != '':
			res[i] = list2[i]
			res2[i] = list2[i]
		if list1[i] != '' and list2[i] != '':
			if list1[i] == list2[i]:
				res[i] = list1[i]
				res2[i] = list1[i]
			else:
				if str(list1[i]) in str(list2[i]):
					res[i] = list2[i]
					res2[i] = list2[i]
				elif str(list1[i]) in str(list2[i]):
					res[i] = list1[i]
					res2[i] = list1[i]
				else:
					if i > 3:
						if str(list1[i]).endswith('AL'):
							res[i] = list1[i]
							res2[i] = list1[i]
						elif str(list2[i]).endswith('AL'):
							res[i] = list1[i]
							res2[i] = list1[i]
						else:
							res[i] = list1[i]
							res2[i] = list1[i]
					else:
						res[i] = list1[i]
						res2[i] = list2[i]
	if res != res2:
		res = pd.concat([pd.DataFrame(res), pd.DataFrame(res2)], axis=1).transpose()
	else:
		res = pd.DataFrame(res).transpose()
	return(res)

def fill_missing_allelic(region, regionMcS):
	temres = pd.DataFrame(['','','','','','','','','','','']).transpose()
	temres.columns = region.columns
	temres.iloc[0, C[regionMcS.iloc[0,1][6:7]]] = regionMcS.iloc[0,1]
	temres.iloc[0, C[regionMcS.iloc[0,2][6:7]]] = regionMcS.iloc[0,2]
	temres.iloc[0, C['%s%s'%(min(regionMcS.iloc[0,1][6:7], regionMcS.iloc[0,2][6:7]), max(regionMcS.iloc[0,1][6:7], regionMcS.iloc[0,2][6:7]))]] = regionMcS.iloc[0,0]
	for i in range(1, regionMcS.shape[0]):
		try: # fill the missing allelic genes
			loc = np.where(temres['Gene%s'%regionMcS.iloc[i,1][6:7].upper()] == regionMcS.iloc[i,1])[0][0]
			if temres.iloc[loc, C[regionMcS.iloc[i,2][6:7]]] == '':
				temres.iloc[loc, C[regionMcS.iloc[i,2][6:7]]] = regionMcS.iloc[i,2]
				temres.iloc[loc, C['%s%s'%(min(regionMcS.iloc[i,1][6:7], regionMcS.iloc[i,2][6:7]), max(regionMcS.iloc[i,1][6:7], regionMcS.iloc[i,2][6:7]))]] = regionMcS.iloc[i,0]
		except: # add new line
			temres_02 = pd.DataFrame(['','','','','','','','','','','']).transpose()
			temres_02.columns = region.columns
			temres = pd.concat([temres, temres_02], axis=0)
			temres.iloc[-1, C[regionMcS.iloc[i,1][6:7]]] = regionMcS.iloc[i,1]
			temres.iloc[-1, C[regionMcS.iloc[i,2][6:7]]] = regionMcS.iloc[i,2]
			temres.iloc[-1, C['%s%s'%(min(regionMcS.iloc[i,1][6:7], regionMcS.iloc[i,2][6:7]), max(regionMcS.iloc[i,1][6:7], regionMcS.iloc[i,2][6:7]))]] = regionMcS.iloc[i,0]
	# compare two regions
	region = region.fillna('')
	combined = pd.concat([region, temres], axis=0, ignore_index=True)
	combined = combined.sort_values(['GeneA', 'GeneB', 'GeneC', 'GeneD'])
	combined_update = combined.copy()
	save = pd.DataFrame(['','','','','','','','','','','']).transpose()
	save.columns = region.columns
	for i in range(0, combined.shape[0] - 1):
		for j in range(i + 1, combined.shape[0]):
			shared = list(set(combined.iloc[i,0:4]) & set(combined.iloc[j,0:4]))
			if len(shared) > 1 or (len(shared) == 1 and shared[0] != ''):
				save_line = Combine_two_list(combined.iloc[i,:].tolist(), combined.iloc[j,:].tolist())
				save_line.columns = region.columns
				save = pd.concat([save, save_line], axis=0)
				try:
					combined_update = combined_update.drop(combined.index[i])
				except:
					None
				try:
					combined_update = combined_update.drop(combined.index[j])
				except:
					None
	save = pd.concat([save, combined_update],axis=0, ignore_index=True)
	return(save[1:].drop_duplicates())

res = pd.read_csv('Psa_Syri_allelic_genes_based_on_GFF_20250324.txt',header=0, sep='\t')
res['A_C'] = ''
res['A_D'] = ''
res['B_D'] = ''
res['Chr'] = ''
for i in range(0, res.shape[0]):
	for j in range(0, 4):
		try:
			if res.iloc[i,j].startswith('Psa'):
				Chr = res.iloc[i,j][4:6]
				res.iloc[i,10] = Chr
				break
		except:
			None

set3 = ['02', '06', '08', '09', '13', '18']

Final_save = pd.DataFrame(['','','','','','','','','','','']).transpose()
Final_save.columns = res.columns

for Chr in res['Chr'].unique():
	subres = res[res['Chr'] == Chr]
	subMcS = McScanX[McScanX['Chr'] == Chr]
	gene_count = subres.iloc[:,0:4].count(axis = 'columns')
	if Chr in set3:
		loci = np.where(gene_count == 3)[0]
	else:
		loci = np.where(gene_count == 4)[0]
	start = 0
	for end in loci:
		if end - start > 0:
			region = subres.iloc[start:end,:]
			genes = list(set(list(set(region['GeneA'])) + list(set(region['GeneB'])) + list(set(region['GeneC'])) + list(set(region['GeneD']))))
			regionMcS = subMcS[subMcS['Gene1'].isin(genes) & subMcS['Gene2'].isin(genes)]
			if regionMcS.shape[0] > 0:
				Final_save = pd.concat([Final_save, fill_missing_allelic(region, regionMcS)], axis=0, ignore_index=True)
			else:
				Final_save = pd.concat([Final_save, region], axis=0, ignore_index=True)
		Final_save = pd.concat([Final_save, subres.iloc[end:(end+1),:]], axis=0, ignore_index=True)
		start = end + 1
		Final_save = Final_save.drop_duplicates()
		print(Final_save.shape[0])
	if start != subres.shape[0]:
		end = subres.shape[0]
		region = subres.iloc[start:end,:]
		genes = list(set(list(set(region['GeneA'])) + list(set(region['GeneB'])) + list(set(region['GeneC'])) + list(set(region['GeneD']))))
		regionMcS = subMcS[subMcS['Gene1'].isin(genes) & subMcS['Gene2'].isin(genes)]
		if regionMcS.shape[0] > 0:
			Final_save = pd.concat([Final_save, fill_missing_allelic(region, regionMcS)], axis=0, ignore_index=True)
		else:
			Final_save = pd.concat([Final_save, region], axis=0, ignore_index=True)
		Final_save = Final_save.drop_duplicates()
		print(Final_save.shape[0])

Final_save.to_csv('Psa_syri_allelic_fillout_using_McScanX_20250325.txt',header=True, index=False, sep='\t')

# fillout allelic genes in inversion
# add allelic genes in inversion information
inverse = pd.read_csv('../06_collinear_analysis/Psa_inversed_genes.txt',header=None, sep='\t') # results from 02_get_genes_in_syntenic_inversion_translocation_regions_based_on_syri.py
inverse['Align_ID'] = 'INVAL'
inverse['Chr'] = inverse[0].str[4:6:]
inverse = inverse.loc[:,['Align_ID', 0, 1, 'Chr']]
inverse.columns = McScanX.columns

translocated = pd.read_csv('../06_collinear_analysis/Psa_translocated_genes.txt',header=None, sep='\t')
translocated['Align_ID'] = 'TRANSAL'
translocated['Chr'] = translocated[0].str[4:6:]
translocated = translocated.loc[:,['Align_ID', 0, 1, 'Chr']]
translocated.columns = McScanX.columns

inverse = pd.concat([inverse,translocated], axis = 0, ignore_index = True)

res = pd.read_csv('Psa_syri_allelic_fillout_using_McScanX_20250325.txt',header=0, sep='\t')
res = res[1:]
for i in range(0, res.shape[0]):
	for j in range(0, 4):
		try:
			if res.iloc[i,j].startswith('Psa'):
				Chr = res.iloc[i,j][4:6]
				res.iloc[i,10] = Chr
				break
		except:
			None

set3 = ['02', '06', '08', '09', '13', '18']

Final_save = pd.DataFrame(['','','','','','','','','','','']).transpose()
Final_save.columns = res.columns

for Chr in res['Chr'].unique():
	subres = res[res['Chr'] == Chr]
	subMcS = inverse[inverse['Chr'] == Chr]
	gene_count = subres.iloc[:,0:4].count(axis = 'columns')
	if Chr in set3:
		loci = np.where(gene_count == 3)[0]
	else:
		loci = np.where(gene_count == 4)[0]
	start = 0
	for end in loci:
		if end - start > 0:
			region = subres.iloc[start:end,:]
			genes = list(set(list(set(region['GeneA'])) + list(set(region['GeneB'])) + list(set(region['GeneC'])) + list(set(region['GeneD']))))
			regionMcS = subMcS[subMcS['Gene1'].isin(genes) | subMcS['Gene2'].isin(genes)]
			if regionMcS.shape[0] > 0:
				Final_save = pd.concat([Final_save, fill_missing_allelic(region, regionMcS)], axis=0, ignore_index=True)
			else:
				Final_save = pd.concat([Final_save, region], axis=0, ignore_index=True)
		Final_save = pd.concat([Final_save, subres.iloc[end:(end+1),:]], axis=0, ignore_index=True)
		start = end + 1
		Final_save = Final_save.drop_duplicates()
		print(Final_save.shape[0])
	if start != subres.shape[0]:
		end = subres.shape[0]
		region = subres.iloc[start:end,:]
		genes = list(set(list(set(region['GeneA'])) + list(set(region['GeneB'])) + list(set(region['GeneC'])) + list(set(region['GeneD']))))
		regionMcS = subMcS[subMcS['Gene1'].isin(genes) & subMcS['Gene2'].isin(genes)]
		if regionMcS.shape[0] > 0:
			Final_save = pd.concat([Final_save, fill_missing_allelic(region, regionMcS)], axis=0, ignore_index=True)
		else:
			Final_save = pd.concat([Final_save, region], axis=0, ignore_index=True)
		Final_save = Final_save.drop_duplicates()
		print(Final_save.shape[0])

Final_save.to_csv('Psa_syri_allelic_fillout_using_McScanX_and_inversion_20250327.txt',header=True, index=False, sep='\t')

df = pd.read_csv('Psa_syri_allelic_fillout_using_McScanX_and_inversion_20250327.txt',header=0, sep='\t')
df = df[1:]
allelic = df.iloc[:,0:4].drop_duplicates()
allelic.to_csv('Psa_syri_allelic_fillout_using_McScanX_and_inversion_20250327_drop_duplicates.txt',header=True, index=False, sep='\t')
allelic['Chr'] = ''
for i in range(0, allelic.shape[0]):
	for j in range(0, 4):
		try:
			if allelic.iloc[i,j].startswith('Psa'):
				Chr = allelic.iloc[i,j][4:6]
				allelic.iloc[i,4] = Chr
				break
		except:
			None

def Order_allelic(df):
	df_ordered = pd.DataFrame()
	df_ordered = pd.concat([df_ordered, df[df['GeneA'].notna()]])
	df_ordered = df_ordered.sort_values(by='GeneA')
	df_GeneA_na = df[df['GeneA'].isna()]
	df_GeneB_notna = df_GeneA_na[df_GeneA_na['GeneB'].notna()]
	df_GeneB_notna = df_GeneB_notna.sort_values(by='GeneB')
	df_GeneB_na = df_GeneA_na[df_GeneA_na['GeneB'].isna()]
	df_GeneC_notna = df_GeneB_na[df_GeneB_na['GeneC'].notna()]
	df_GeneC_notna = df_GeneC_notna.sort_values(by='GeneC')
	df_GeneC_na = df_GeneB_na[df_GeneB_na['GeneC'].isna()]
	df_GeneD_notna = df_GeneC_na[df_GeneC_na['GeneD'].notna()]
	df_GeneD_notna = df_GeneD_notna.sort_values(by='GeneD')
	# insert B genes
	for i in range(0, df_GeneB_notna.shape[0]):
		try:
			loc = np.where(df_ordered['GeneB'] > df_GeneB_notna.iloc[i,1])[0][0]
			df_ordered = pd.concat([df_ordered.iloc[0:loc, :], df_GeneB_notna.iloc[[i], :], df_ordered.iloc[loc:, :]], axis=0).reset_index(drop=True)
		except:
			loc = df_ordered.shape[0]
			df_ordered = pd.concat([df_ordered.iloc[0:loc, :], df_GeneB_notna.iloc[[i], :]], axis=0).reset_index(drop=True)
	# insert C genes
	for i in range(0, df_GeneC_notna.shape[0]):
		try:
			loc = np.where(df_ordered['GeneC'] > df_GeneC_notna.iloc[i,2])[0][0]
			df_ordered = pd.concat([df_ordered.iloc[0:loc, :], df_GeneC_notna.iloc[[i], :], df_ordered.iloc[loc:, :]], axis=0).reset_index(drop=True)
		except:
			loc = df_ordered.shape[0]
			df_ordered = pd.concat([df_ordered.iloc[0:loc, :], df_GeneC_notna.iloc[[i], :]], axis=0).reset_index(drop=True)
	# insert D genes
	for i in range(0, df_GeneD_notna.shape[0]):
		try:
			loc = np.where(df_ordered['GeneD'] > df_GeneD_notna.iloc[i,3])[0][0]
			df_ordered = pd.concat([df_ordered.iloc[0:loc, :], df_GeneD_notna.iloc[[i], :], df_ordered.iloc[loc:, :]], axis=0).reset_index(drop=True)
		except:
			loc = df_ordered.shape[0]
			df_ordered = pd.concat([df_ordered.iloc[0:loc, :], df_GeneD_notna.iloc[[i], :]], axis=0).reset_index(drop=True)
	return(df_ordered)

# order the allelic first
df_ordered = pd.DataFrame(['','','','','']).transpose()
df_ordered.columns = allelic.columns
for Chr in allelic['Chr'].unique():
	subdf = allelic[allelic['Chr'] == Chr]
	df_ordered = pd.concat([df_ordered, Order_allelic(subdf)])

df_ordered.iloc[1:,:].to_csv('Psa_syri_allelic_fillout_using_McScanX_and_inversion_20250327_ordered.txt',header=True, index=False, sep='\t')

# add tandem information from Tandem and then combine tandem duplicates
def get_all_genes(List):
	genes = []
	for l in List:
		if ',' in l:
			for g in l.split(','):
				genes.append(g)
		else:
			if l != '':
				genes.append(l)
	uniquegenes = sorted(list(set(genes)))
	return(uniquegenes)

with open('Psa_tandem_genes.pkl', 'rb') as f1:  # results from 01_parse_McScanX_and_tandem_information.py
	Tandem = pickle.load(f1)

def Get_tandem_duplicates(genes):
	td = []
	for g in genes:
		if g in Tandem:
			td.append(g)
			td = td + Tandem[g]
		else:
			td.append(g)
	td = sorted(list(set(td)))
	return(td)

df_ordered = pd.read_csv('Psa_syri_allelic_fillout_using_McScanX_and_inversion_20250327_ordered.txt',header=0, sep='\t')
df_ordered = df_ordered.dropna(axis = 0, how = 'all')
df_ordered = df_ordered.fillna('')
df_grr = df_ordered.copy()
for i in df_grr.index:
	genes = get_all_genes(df_grr.loc[i,['GeneA', 'GeneB', 'GeneC', 'GeneD']])
	if '' in genes:
		genes.remove('')
	tem = df_grr[df_grr['GeneA'].isin(genes) | df_grr['GeneB'].isin(genes) | df_grr['GeneC'].isin(genes) | df_grr['GeneD'].isin(genes)] 
	genes = list(set(get_all_genes(tem['GeneA'].tolist())) | set(get_all_genes(tem['GeneB'].tolist())) | set(get_all_genes(tem['GeneC'].tolist())) | set(get_all_genes(tem['GeneD'].tolist())))
	genes = Get_tandem_duplicates(genes)
	subdf = df_grr[df_grr['GeneA'].isin(genes) | df_grr['GeneB'].isin(genes) | df_grr['GeneC'].isin(genes) | df_grr['GeneD'].isin(genes)] 
	while subdf.shape[0] != tem.shape[0]:
		tem = subdf
		genes = list(set(get_all_genes(tem['GeneA'].tolist())) | set(get_all_genes(tem['GeneB'].tolist())) | set(get_all_genes(tem['GeneC'].tolist())) | set(get_all_genes(tem['GeneD'].tolist())))
		genes = Get_tandem_duplicates(genes)
		subdf = df_grr[df_grr['GeneA'].isin(genes) | df_grr['GeneB'].isin(genes) | df_grr['GeneC'].isin(genes) | df_grr['GeneD'].isin(genes)] 
	if subdf.shape[0] > 1:
		for k in ['GeneA','GeneB','GeneC','GeneD']:
			df_grr.loc[i,k] = ','.join(get_all_genes(subdf.loc[:,k].tolist()))
		for j in subdf.index:
			if j != i:
				df_grr.loc[j,:] = ['','','','','']


df_grr = df_grr.dropna(axis = 0, how = 'all')
df_grr = df_grr.drop_duplicates()
df_grr.to_csv('Psa_syri_allelic_fillout_using_McScanX_and_inversion_20250328_combine_tandem_duplicates.txt',header=True, index=False, sep='\t')


# categorize the allelic genes
# count number of allelic genes
def Type_allelic(types):
	tem = types.split('/')
	if len(list(set(tem))) == 1:
		return(types)
	else:
		for i in range(0, len(tem)):
			if int(tem[i]) > 1:
				tem[i] = 'n'	
		tem.sort()
		return('/'.join(tem))

df_ordered = pd.read_csv('Psa_syri_allelic_fillout_using_McScanX_and_inversion_20250328_combine_tandem_duplicates.txt',header=0, sep='\t')
df_ordered = df_ordered.iloc[:,0:4]

df_ordered = df_ordered.dropna(axis = 0, how = 'all')
df_ordered = df_ordered.fillna('')

Type = pd.DataFrame({'Type_of_allelic':'', 'Type_ab':''}, index=range(0, df_ordered.shape[0]))
df_ordered = pd.concat([df_ordered, Type], axis=1)
chr_list = [2, 6, 8, 9, 13, 18]
for i in df_ordered.index:
	chr = ''
	if str(df_ordered.loc[i, 'GeneA']).startswith('Psa'):
		count_A = len(df_ordered.loc[i, 'GeneA'].split(','))
		chr = int(df_ordered.loc[i, 'GeneA'].split('G')[0].split('Psa_')[1][0:-1])
	else:
		count_A = 0
	if str(df_ordered.loc[i, 'GeneB']).startswith('Psa'):
		count_B = len(df_ordered.loc[i, 'GeneB'].split(','))
		chr = int(df_ordered.loc[i, 'GeneB'].split('G')[0].split('Psa_')[1][0:-1])
	else:
		count_B = 0
	if str(df_ordered.loc[i, 'GeneC']).startswith('Psa'):
		count_C = len(df_ordered.loc[i, 'GeneC'].split(','))
		chr = int(df_ordered.loc[i, 'GeneC'].split('G')[0].split('Psa_')[1][0:-1])
	else:
		count_C = 0
	if str(df_ordered.loc[i, 'GeneD']).startswith('Psa'):
		count_D = len(df_ordered.loc[i, 'GeneD'].split(','))
		chr = int(df_ordered.loc[i, 'GeneD'].split('G')[0].split('Psa_')[1][0:-1])
	else:
		count_D = 0
	if chr not in chr_list:
		df_ordered.loc[i, 'Type_of_allelic'] = '%s/%s/%s/%s'%(count_A, count_B, count_C, count_D) 
		df_ordered.loc[i, 'Type_ab'] = Type_allelic(df_ordered.loc[i, 'Type_of_allelic']) 
	else:
		df_ordered.loc[i, 'Type_of_allelic'] = '%s/%s/%s'%(count_A, count_B, count_C) 
		df_ordered.loc[i, 'Type_ab'] = Type_allelic(df_ordered.loc[i, 'Type_of_allelic']) 

df_ordered.to_csv('Psa_allelic_genes_final_types_20250328.txt', index=False, header=True, sep='\t')
counts = df_ordered['Type_of_allelic'].value_counts()
counts.to_csv('Psa_allelic_genes_final_types_counts_20250328.txt', index=True, header=True, sep='\t')
counts_ab = df_ordered['Type_ab'].value_counts()
counts_ab.to_csv('Psa_allelic_genes_final_types_counts_20250328_higher_level.txt', index=True, header=True, sep='\t')


# functional annotation
Blast = pd.read_csv('Psa_Vs_Psa_top_non-self_blast_hit.txt', header=None, sep='\t')  # results from 00_parse_Psa_allbyall_blast_results.py
D_Blast = {}
for i in range(0, Blast.shape[0]):
	if Blast.iloc[i,0] not in D_Blast:
		D_Blast[Blast.iloc[i,0]] = Blast.iloc[i,1]

path = '/public/home/wangpeipei/08_Psa_genome_annotation_v2_20250114/2_19_9_total_functional_annotation_based_on_longest_pep/'
COG = pd.read_csv(path + '000_final_cog_anno_longest_transcripts.tsv', header=0, sep='\t')
COG['GID'] = COG['GID'].str[0:14]
COG = COG[COG['COG_Name'] != 'Function unknown']
D_COG = {}
for i in range(0, COG.shape[0]):
	if COG.iloc[i,0] not in D_COG:
		D_COG[COG.iloc[i,0]] = 1

InterProScan = pd.read_csv(path + '00_final_interproscan_final.func_info_gene_refer_to_longest_transcripts.txt', header=None, sep='\t')
D_InterProScan = {}
for i in range(0, InterProScan.shape[0]):
	if InterProScan.iloc[i,0] not in D_InterProScan:
		D_InterProScan[InterProScan.iloc[i,0]] = 1

GO = pd.read_csv(path + '01_final_gene_GO_information_longest_transcripts.txt', header=None, sep='\t')
GO[0] = GO[0].str[0:14]
D_GO = {}
for i in range(0, GO.shape[0]):
	if GO.iloc[i,0] not in D_GO:
		D_GO[GO.iloc[i,0]] = 1

Pfam = pd.read_csv(path + '02_final_annotation_Pfam_longest_transcripts.txt', header=None, sep='\t')
Pfam[0] = Pfam[0].str[0:14]
D_Pfam = {}
for i in range(0, Pfam.shape[0]):
	if Pfam.iloc[i,0] not in D_Pfam:
		D_Pfam[Pfam.iloc[i,0]] = 1

TAIR = pd.read_csv(path + '04_final_mapping_P.san_AT_pep_and_description_longest_transcripts.txt', header=None, sep='\t')
TAIR[0] = TAIR[0].str[0:14]
D_TAIR = {}
for i in range(0, TAIR.shape[0]):
	if TAIR.iloc[i,0] not in D_TAIR:
		D_TAIR[TAIR.iloc[i,0]] = 1

NR = pd.read_csv(path + '05_final_nr_longest_transcripts.txt', header=None, sep='\t')
NR[0] = NR[0].str[0:14]
D_NR = {}
for i in range(0, NR.shape[0]):
	if NR.iloc[i,0] not in D_NR:
		D_NR[NR.iloc[i,0]] = 1

Unipro = pd.read_csv(path + '06_final_uniprot_description_longest_transcripts.txt', header=None, sep='\t')
Unipro[0] = Unipro[0].str[0:14]
D_Unipro = {}
for i in range(0, Unipro.shape[0]):
	if Unipro.iloc[i,0] not in D_Unipro:
		D_Unipro[Unipro.iloc[i,0]] = 1

EC = pd.read_csv(path + '07_final_P.san_EC_longest_transcripts.txt', header=None, sep='\t')
EC[0] = EC[0].str[0:14]
D_EC = {}
for i in range(0, EC.shape[0]):
	if EC.iloc[i,0] not in D_EC:
		D_EC[EC.iloc[i,0]] = 1

IPR = pd.read_csv(path + '08_final_IPR_info_longest_transcripts.txt', header=None, sep='\t')
IPR[0] = IPR[0].str[0:14]
D_IPR = {}
for i in range(0, IPR.shape[0]):
	if IPR.iloc[i,0] not in D_IPR:
		D_IPR[IPR.iloc[i,0]] = 1

def Functional_annotation(gene):
	res = []
	if gene in D_COG:
		res.append('COG')
	if gene in D_InterProScan:
		res.append('InterProScan')
	if gene in D_GO:
		res.append('GO')
	if gene in D_Pfam:
		res.append('Pfam')
	if gene in D_TAIR:
		res.append('TAIR')
	if gene in D_NR:
		res.append('NR')
	if gene in D_Unipro:
		res.append('Unipro')
	if gene in D_EC:
		res.append('EC')
	if gene in D_IPR:
		res.append('IPR')
	if len(res) == 0:
		return('No functional annotation')
	else:
		return(';'.join(res))

# protein length
protein_length = pd.read_csv('/public/home/wangpeipei/08_Psa_genome_annotation_v2_20250114/Psa_longest_pep_length_20250114.txt', header=0, sep = ',')
PL = {}
for i in range(0, protein_length.shape[0]):
	PL[protein_length.iloc[i,0]] = int(protein_length.iloc[i,1])

Expression = pd.read_csv('/public/home/wangpeipei/07_analysis_after_genome/07_transcriptome_analysis/Psa_normalized_mean_TPM_clean_20250217.csv', header=0, sep=',')

df_ordered = pd.read_csv('Psa_allelic_genes_final_types_20250328.txt', header=0, sep='\t')

# classify all genes
out = open('Psa_allelic_genes_final_functional_annotations_20250408.txt','w')
out.write('Gene\tType\tBest_blast_hit\tFunctional_annotation\tExpression\tProtein_length\n')
for i in range(0, df_ordered.shape[0]):
	if df_ordered.loc[i, 'Type_ab'] in ['0/0/1', '0/0/0/1']:
		alle_type = 'Singleton'
	elif len(list(set(df_ordered.loc[i, 'Type_ab'].split('/')))) == 1:
		if df_ordered.loc[i, 'Type_ab'].split('/')[0] == '1':
			alle_type = 'Conserved'
		else:
			if len(list(set(df_ordered.loc[i, 'Type_of_allelic'].split('/')))) == 1:
				alle_type = 'Conserved_with_same_tandem_duplicates'
			else:
				alle_type = 'With_biased_tandem_duplicates'
	elif '0' in df_ordered.loc[i, 'Type_ab'].split('/'):
		alle_type = 'With_biased_gene_loss'
	else:
		alle_type = 'With_biased_tandem_duplicates'
	for j in range(0, 4):
		if str(df_ordered.iloc[i, j]).startswith('Psa'):
			for gene in df_ordered.iloc[i, j].split(','):
				res = gene + '\t' + alle_type
				# for blast
				if gene in D_Blast:
					hit = D_Blast[gene]
					chr1 = gene.split('Psa_')[1].split('G')[0]
					chr2 = hit.split('Psa_')[1].split('G')[0]
					if chr1 == chr2: # within the same chr
						res = res + '\t' + 'Best_hit_on_same_chr'
					if chr1 != chr2 and chr1[0:-1] == chr2[0:-1]:
						res = res + '\t' + 'Best_hit_on_alle_chr'
					if chr1[0:-1] != chr2[0:-1]:
						res = res + '\t' + 'Best_hit_on_diff_chr'
				else:
					res = res + '\t' +  'No_blast_hit'
				# for functional annotation
				res = res + '\t' +  Functional_annotation(gene)
				# for expression
				exp = Expression[Expression['Unnamed: 0'] == gene]
				if any(exp.iloc[0, 1:] > 2):   # 2 is the threshold to call whether a gene is expressed or not
					if any(exp.iloc[0, 1:] > 100): 
						res = res + '\t' +  'highly_Expressed'
					else:
						res = res + '\t' +  'Expressed'
				elif all(exp.iloc[0, 1:] == 0):
					res = res + '\t' +  'Not_Expressed'
				else:
					res = res + '\t' +  'Almost_not_expressed'
				# for protein length
				res = res + '\t%s'%PL[gene]
				out.write(res + '\n')
				out.flush()

out.close()

