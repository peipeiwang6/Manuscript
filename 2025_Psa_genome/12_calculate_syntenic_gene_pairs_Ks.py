import sys,os
import pandas as pd
import numpy as np
import pickle

pep = open('/public/home/wangpeipei/08_Psa_genome_annotation_v2_20250114/Psa_longest_pep_20250114.fa','r').readlines()
Pep = {}
i = 0
while i < len(pep):
	inl = pep[i].strip()
	if inl.startswith('>'):
		Pep[inl[1:]] = pep[i + 1].strip()
		i += 1
	i += 1

cds = open('/public/home/wangpeipei/08_Psa_genome_annotation_v2_20250114/Psa_CDS_20250113.fa','r').readlines()
CDS = {}
i = 0
while i < len(cds):
	inl = cds[i].strip()
	if inl.startswith('>'):
		if inl[1:] in Pep:
			CDS[inl[1:]] = cds[i + 1].strip()
			i += 1
			i += 1
			while not cds[i].startswith('>'):
				CDS[inl[1:]] = CDS[inl[1:]] + cds[i].strip()
				i += 1
			i = i - 1
	i += 1

f = open("Psa_representative_pep.pkl","wb")
pickle.dump(Pep,f)
f.close()

f = open("Psa_representative_cds.pkl","wb")
pickle.dump(CDS,f)
f.close()

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

collinear = collinear[collinear['Align_ID'].isin(align_winner['Align'].tolist())]

D = {}
for i in range(0, collinear.shape[0]):
	if collinear.iloc[i, 0] in align_winner['Align'].tolist():
		if collinear.iloc[i, 0] not in D:
			D[collinear.iloc[i, 0]] = 1
			os.system('mkdir Ks_for_syntenic_blocks/%s'%(collinear.iloc[i, 0]))
		gene1 = collinear.iloc[i,2]
		gene2 = collinear.iloc[i,3]
		os.system('mkdir Ks_for_syntenic_blocks/%s/%s_%s'%(collinear.iloc[i, 0], min(gene1, gene2), max(gene1, gene2)))
		path = 'Ks_for_syntenic_blocks/%s/%s_%s'%(collinear.iloc[i, 0], min(gene1, gene2), max(gene1, gene2))
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
				# os.system('awk \'{gsub(\"stewart.aa\",\"%s/%s_cds_align.fas\",$0); print $0}\' %s/codeml.ctl > %s/codeml.ctl_tem'%(path, file.split('_pep.fas')[0], path, path))
				# os.system('awk \'{gsub(\"mlc\",\"%s/%s\",$0); print $0}\' %s/codeml.ctl_tem > %s/codeml.ctl'%(path, file.split('_pep.fas')[0], path, path))
				# os.system('rm %s/codeml.ctl_tem'%path)
				# os.system('/public/home/wangpeipei/Software/PAML/paml4.9j/bin/codeml %s/codeml.ctl'%path)
				os.system('cp /public/home/wangpeipei/Software/PAML/paml4.9j/yn00.ctl %s'%path)
				os.system('awk \'{gsub(\"examples/abglobin.nuc\",\"%s/%s_cds_align.fas\",$0); print $0}\' %s/yn00.ctl > %s/yn00.ctl_tem'%(path, file.split('_pep.fas')[0], path, path))
				os.system('awk \'{gsub(\"= yn\",\"= %s/%s_raml.out\",$0); print $0}\' %s/yn00.ctl_tem > %s/yn00.ctl'%(path, file.split('_pep.fas')[0], path, path))
				os.system('rm %s/yn00.ctl_tem'%path)
				os.system('/public/home/wangpeipei/Software/PAML/paml4.9j/bin/yn00 %s/yn00.ctl'%path)

align_winner = pd.read_csv('Psa_Collinear_blocks_allelic_winner.txt',sep='\t',header=0)
align_winner.columns = ['Align','Winning','Lossing','n','Proportion']

out = open('Psa_winner_allelic_gene_pairs_KaKs.txt','w')
out.write('Align\tChr_pair\tGene1\tGene2\tdN\tdS\tOmega\n')
for align in align_winner['Align'].tolist():
	for pair in os.listdir('Ks_for_syntenic_blocks/%s/'%align):
		file = open('Ks_for_syntenic_blocks/%s/%s/%s_raml.out'%(align, pair, pair),'r').readlines()
		for i in range(0, len(file)):
			if file[i].startswith('seq. seq.'):
				i = i + 2
				tem = file[i].split()
				omega = tem[6]
				dn = tem[7]
				ds = tem[10]
		out.write('%s\t%s_%s\t%s\t%s\t%s\t%s\t%s\n'%(align, pair[4:7], pair[19:22], pair.split('_P')[0], 'P'+pair.split('_P')[1], dn, ds, omega))
		out.flush()

out.close()

