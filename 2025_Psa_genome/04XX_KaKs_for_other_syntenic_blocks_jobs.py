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

collinear = pd.read_csv('Psa_Collinear_blocks.txt',header=0, sep='\t', low_memory=False)
align_winner = pd.read_csv('Psa_Collinear_blocks_allelic_winner.txt',sep='\t',header=0)
align_winner.columns = ['Align','Winning','Lossing','n','Proportion']

# only calculate for allelic syntenic gene pairs
collinear = collinear[-collinear['Align_ID'].isin(align_winner['Align'].tolist())]
collinear = collinear[collinear['Type'].isin(['within_same','between_diff'])]
collinear = collinear[collinear['Gene1'].str.contains('aG')]
collinear = collinear[collinear['Gene2'].str.contains('aG')] #79546

start = int(sys.argv[1])
end = int(sys.argv[2])

D = {}
for i in range(start, end):
	if collinear.iloc[i, 0] not in align_winner['Align'].tolist():
		if collinear.iloc[i, 0] not in D:
			D[collinear.iloc[i, 0]] = 1
			os.system('mkdir Ks_for_non_allelic_syntenic_blocks/%s'%(collinear.iloc[i, 0]))
		gene1 = collinear.iloc[i,2]
		gene2 = collinear.iloc[i,3]
		os.system('mkdir Ks_for_non_allelic_syntenic_blocks/%s/%s_%s'%(collinear.iloc[i, 0], min(gene1, gene2), max(gene1, gene2)))
		path = 'Ks_for_non_allelic_syntenic_blocks/%s/%s_%s'%(collinear.iloc[i, 0], min(gene1, gene2), max(gene1, gene2))
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
				# os.system('cp /public/home/wangpeipei/Software/PAML/paml4.9j/codeml.ctl %s'%path)
				# os.system('awk \'{gsub(\"stewart.aa\",\"%s/%s_cds_aligned.fas\",$0); print $0}\' %s/codeml.ctl > %s/codeml.ctl_tem'%(path, file.split('_pep.fas')[0], path, path))
				# os.system('awk \'{gsub(\"mlc\",\"%s/%s\",$0); print $0}\' %s/codeml.ctl_tem > %s/codeml.ctl'%(path, file.split('_pep.fas')[0], path, path))
				# os.system('rm %s/codeml.ctl_tem'%path)
				# os.system('/public/home/wangpeipei/Software/PAML/paml4.9j/bin/codeml %s/codeml.ctl'%path)
				os.system('cp /public/home/wangpeipei/Software/PAML/paml4.9j/yn00.ctl %s'%path)
				os.system('awk \'{gsub(\"examples/abglobin.nuc\",\"%s/%s_cds_aligned.fas\",$0); print $0}\' %s/yn00.ctl > %s/yn00.ctl_tem'%(path, file.split('_pep.fas')[0], path, path))
				os.system('awk \'{gsub(\"= yn\",\"= %s/%s_raml.out\",$0); print $0}\' %s/yn00.ctl_tem > %s/yn00.ctl'%(path, file.split('_pep.fas')[0], path, path))
				os.system('rm %s/yn00.ctl_tem'%path)
				os.system('/public/home/wangpeipei/Software/PAML/paml4.9j/bin/yn00 %s/yn00.ctl'%path)

