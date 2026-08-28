import sys,os
import pandas as pd
import numpy as np
import pickle

# functional annotation
Blast = pd.read_csv('Psa_Vs_Psa_top_non-self_blast_hit.txt', header=None, sep='\t')
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

def Detect_domain_variation(gene_list):
	Domain = []
	res = 0
	for gene in gene_list:
		Domain.append(Pfam[Pfam[0]==gene][1].tolist())
	if len(Domain) > 1:
		for i in range(0, len(Domain)-1):
			for j in range(i+1, len(Domain)):
				diff = list(set(set(Domain[i]).difference(set(Domain[j]))) | set(set(Domain[j]).difference(set(Domain[i]))))
				if len(diff) > 0:
					res = 1
	return(res)


# protein length
protein_length = pd.read_csv('/public/home/wangpeipei/08_Psa_genome_annotation_v2_20250114/Psa_longest_pep_length_20250114.txt', header=0, sep = ',')
PL = {}
for i in range(0, protein_length.shape[0]):
	PL[protein_length.iloc[i,0]] = int(protein_length.iloc[i,1])

Expression = pd.read_csv('/public/home/wangpeipei/07_analysis_after_genome/07_transcriptome_analysis/Psa_normalized_mean_TPM_clean_20250217.csv', header=0, sep=',')

df_ordered = pd.read_csv('Psa_allelic_genes_final_types_20250328.txt', header=0, sep='\t')


# classify all genes
out = open('Psa_allelic_genes_final_functional_annotations_six_classes_20260109.txt','w')
out.write('Gene\tType\tBest_blast_hit\tFunctional_annotation\tExpression\tProtein_length\n')
for i in range(0, df_ordered.shape[0]):
	if df_ordered.loc[i, 'Type_ab'] in ['0/0/1', '0/0/0/1']:
		alle_type = 'Singleton'
	elif len(list(set(df_ordered.loc[i, 'Type_ab'].split('/')))) == 1:
		if df_ordered.loc[i, 'Type_ab'].split('/')[0] == '1':
			gene_list = df_ordered.iloc[i,0:4]
			if Detect_domain_variation(gene_list) == 1:
				alle_type = 'Conserved_but_with_domain_variation'
			else:
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
					if any(exp.iloc[0, 1:] > 10): 
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





























