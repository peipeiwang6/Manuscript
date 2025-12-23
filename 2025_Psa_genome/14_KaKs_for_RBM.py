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

RBM = pd.read_csv('Reciprocal_best_hist_set3_and_set3_genes.txt',header=None, sep='\t', low_memory=False)

for i in range(0, RBM.shape[0]):
	gene1 = RBM.iloc[i,0]
	gene2 = RBM.iloc[i,1]
	os.system('mkdir Ks_for_RBM_set3_vs_set3/%s_%s'%(min(gene1, gene2), max(gene1, gene2)))
	path = 'Ks_for_RBM_set3_vs_set3/%s_%s'%(min(gene1, gene2), max(gene1, gene2))
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

RBM = pd.read_csv('Reciprocal_best_hist_set3_and_set4_genes.txt',header=None, sep='\t', low_memory=False)

for i in range(0, RBM.shape[0]):
	gene1 = RBM.iloc[i,0]
	gene2 = RBM.iloc[i,1]
	os.system('mkdir Ks_for_RBM_set3_vs_set4/%s_%s'%(min(gene1, gene2), max(gene1, gene2)))
	path = 'Ks_for_RBM_set3_vs_set4/%s_%s'%(min(gene1, gene2), max(gene1, gene2))
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

RBM = pd.read_csv('Reciprocal_best_hist_set4_and_set4_genes.txt',header=None, sep='\t', low_memory=False)

for i in range(0, RBM.shape[0]):
	gene1 = RBM.iloc[i,0]
	gene2 = RBM.iloc[i,1]
	os.system('mkdir Ks_for_RBM_set4_vs_set4/%s_%s'%(min(gene1, gene2), max(gene1, gene2)))
	path = 'Ks_for_RBM_set4_vs_set4/%s_%s'%(min(gene1, gene2), max(gene1, gene2))
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

out = open('Psa_RBM_KaKs.txt','w')
out.write('RBM_type\tChr_pair\tGene1\tGene2\tdN\tdS\tOmega\n')
for pair in os.listdir('Ks_for_RBM_set3_vs_set3/'):
	file = open('Ks_for_RBM_set3_vs_set3/%s/%s_raml.out'%(pair, pair),'r').readlines()
	for i in range(0, len(file)):
		if file[i].startswith('seq. seq.'):
			i = i + 2
			tem = file[i].split()
			omega = tem[6]
			dn = tem[7]
			ds = tem[10]
	out.write('%s\t%s_%s\t%s\t%s\t%s\t%s\t%s\n'%('Set3_vs_set3', pair[4:7], pair[19:22], pair.split('_P')[0], 'P'+pair.split('_P')[1], dn, ds, omega))
	out.flush()

for pair in os.listdir('Ks_for_RBM_set3_vs_set4/'):
	file = open('Ks_for_RBM_set3_vs_set4/%s/%s_raml.out'%(pair, pair),'r').readlines()
	for i in range(0, len(file)):
		if file[i].startswith('seq. seq.'):
			i = i + 2
			tem = file[i].split()
			omega = tem[6]
			dn = tem[7]
			ds = tem[10]
	out.write('%s\t%s_%s\t%s\t%s\t%s\t%s\t%s\n'%('Set3_vs_set4', pair[4:7], pair[19:22], pair.split('_P')[0], 'P'+pair.split('_P')[1], dn, ds, omega))
	out.flush()

for pair in os.listdir('Ks_for_RBM_set4_vs_set4/'):
	file = open('Ks_for_RBM_set4_vs_set4/%s/%s_raml.out'%(pair, pair),'r').readlines()
	for i in range(0, len(file)):
		if file[i].startswith('seq. seq.'):
			i = i + 2
			tem = file[i].split()
			omega = tem[6]
			dn = tem[7]
			ds = tem[10]
	out.write('%s\t%s_%s\t%s\t%s\t%s\t%s\t%s\n'%('Set4_vs_set4', pair[4:7], pair[19:22], pair.split('_P')[0], 'P'+pair.split('_P')[1], dn, ds, omega))
	out.flush()

out.close()

