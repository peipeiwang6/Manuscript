import sys,os
import pandas as pd
import numpy as np
import pickle
with open('Psa_representative_pep.pkl', 'rb') as f1:
	Pep = pickle.load(f1)

Pep_update = {}
for gene in Pep:
	gene_update = gene.split('.')[0]
	Pep_update[gene_update] = Pep[gene]

with open('Psa_representative_cds.pkl', 'rb') as f2:
	CDS = pickle.load(f2)

CDS_update = {}
for gene in CDS:
	gene_update = gene.split('.')[0]
	CDS_update[gene_update] = CDS[gene]

Allele = pd.read_csv('Psa_syri_allelic_fillout_using_McScanX_and_inversion_20250328_combine_tandem_duplicates.txt',header=0, sep='\t', low_memory=False)
Type = pd.read_csv('Psa_allelic_genes_final_functional_annotations_six_classes_20260109.txt',header=0, sep='\t', low_memory=False)

Allele_type = Allele.copy()
Allele_type['Type'] = ''
for i in range(0, Allele_type.shape[0]):
	for j in range(0,4):
		if type(Allele_type.iloc[i,j]) != float:
			gene = Allele_type.iloc[i,j].split(',')[0]
			allele_type = Type[Type['Gene']==gene]['Type'].iloc[0]
			break
	Allele_type.iloc[i,5] = allele_type

Allele_type.to_csv('Psa_syri_allelic_fillout_using_McScanX_and_inversion_20250328_combine_tandem_duplicates_six_types.txt', index=False, header=True,sep="\t")
Allele_type = pd.read_csv('Psa_syri_allelic_fillout_using_McScanX_and_inversion_20250328_combine_tandem_duplicates_six_types.txt',header=0, sep='\t', low_memory=False)


T = {'With_biased_gene_loss':'Type2','With_biased_tandem_duplicates':'Type3','Conserved_with_same_tandem_duplicates':'Type4','Conserved':'Type5','Conserved_but_with_domain_variation':'Type6'}
for i in range(0, Allele_type.shape[0]):
	if Allele_type.iloc[i,5] != 'Singleton':
		gene_list = []
		for j in range(0,4):
			if type(Allele_type.iloc[i,j]) != float:
				genes = Allele_type.iloc[i,j].split(',')
				for g in genes:
					gene_list.append(g)
		os.system('mkdir KaKs_all_allele_types/%s/%s'%(T[Allele_type.iloc[i,5]], gene_list[0]))
		gene_pairs = []
		for m in range(0,len(gene_list)-1):
			for n in range(m+1,len(gene_list)):
				if gene_list[m][0:7] != gene_list[n][0:7]:
					gene_pairs.append([gene_list[m],gene_list[n]])
		for pair in gene_pairs:
			os.system('mkdir KaKs_all_allele_types/%s/%s/%s'%(T[Allele_type.iloc[i,5]], gene_list[0],min(pair) + '_' + max(pair)))
			gene1 = pair[0]
			gene2 = pair[1]
			path = 'KaKs_all_allele_types/%s/%s/%s'%(T[Allele_type.iloc[i,5]], gene_list[0],min(pair) + '_' + max(pair))
			pep_seq = open('%s/%s_%s_pep.fas'%(path, min(gene1, gene2), max(gene1, gene2)),'w')
			pep_seq.write('>%s\n%s\n>%s\n%s\n'%(gene1, Pep_update[gene1], gene2, Pep_update[gene2]))
			pep_seq.close()
			cds_seq = open('%s/%s_%s_cds.fas'%(path,min(gene1, gene2), max(gene1, gene2)),'w')
			cds_seq.write('>%s\n%s\n>%s\n%s\n'%(gene1, CDS_update[gene1], gene2, CDS_update[gene2]))
			cds_seq.close()
			for file in os.listdir(path):
				if file.endswith('_pep.fas'):
					os.system('/public/home/likangli/miniforge3/envs/mafft/bin/mafft --anysymbol --maxiterate 1000 --localpair %s/%s > %s/%s_aligned.fas'%(path, file, path, file.split('.fas')[0]))
					os.system('python 13_convert_pep_alignment_to_CDS_alignment.py %s/%s_pep_aligned.fas %s/%s_cds.fas'%(path, file.split('_pep.fas')[0], path, file.split('_pep.fas')[0]))
					os.system('cp /public/home/wangpeipei/Software/PAML/paml4.9j/yn00.ctl %s'%path)
					os.system('awk \'{gsub(\"examples/abglobin.nuc\",\"%s/%s_cds_aligned.fas\",$0); print $0}\' %s/yn00.ctl > %s/yn00.ctl_tem'%(path, file.split('_pep.fas')[0], path, path))
					os.system('awk \'{gsub(\"= yn\",\"= %s/%s_raml.out\",$0); print $0}\' %s/yn00.ctl_tem > %s/yn00.ctl'%(path, file.split('_pep.fas')[0], path, path))
					os.system('rm %s/yn00.ctl_tem'%path)
					os.system('/public/home/wangpeipei/Software/PAML/paml4.9j/bin/yn00 %s/yn00.ctl'%path)

