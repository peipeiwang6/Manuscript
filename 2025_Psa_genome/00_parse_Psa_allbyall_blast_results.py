import sys,os
import pandas as pd
import numpy as np
import pickle

distance = 10

# blast hit within the same chr
out1 = open('Psa_Vs_Psa_within_same_chr.blast.txt','w')
out2 = open('Psa_Vs_Psa_between_alle_chr.blast.txt','w')
out3 = open('Psa_Vs_Psa_between_diff_chr.blast.txt','w')
inp = open('Psa_all_hits.blast','r')
inl = inp.readline()
while inl:
	tem = inl.split('\t')
	chr1 = tem[0].split('Psa_')[1].split('G')[0]
	chr2 = tem[1].split('Psa_')[1].split('G')[0]
	if chr1 == chr2 and tem[0] != tem[1]:
		out1.write(inl)
	if chr1[0:-1] == chr2[0:-1] and chr1[-1] != chr2[-1]:
		out2.write(inl)
	if chr1[0:-1] != chr2[0:-1]:
		out3.write(inl)
	inl = inp.readline()

out1.close()
out2.close()
out3.close()

# the toppest non-self blast hit
inp = open('Psa_all_hits.blast','r')
out4 = open('Psa_Vs_Psa_top_non-self_blast_hit.txt','w')
R = {}
inl = inp.readline()
while inl:
	tem = inl.split('\t')
	if tem[0] not in R:
		if tem[0] != tem[1]:
			out4.write(inl)
			R[tem[0]] = 1
	inl = inp.readline()

out4.close()

GFF = pd.read_csv('Psa_all.gff',sep='\t',header=None)
Blast = pd.read_csv('Psa_Vs_Psa_within_same_chr.blast.txt',sep='\t',header=None)
tandem = open('Psa_all.tandem_more.txt','w')
for chr in GFF[0].unique():
	subGFF = GFF[GFF[0] == chr]
	for i in range(0, subGFF.shape[0] - 1):
		for j in range(i+1, i + distance + 2):
			if j < subGFF.shape[0]:
				tem = Blast[Blast[0] == subGFF.iloc[i, 1]]
				if tem.shape[0] > 0:
					if subGFF.iloc[j, 1] in tem[1].tolist():
						tandem.write('%s,%s\n'%(subGFF.iloc[i, 1], subGFF.iloc[j, 1]))
				else:
					tem = Blast[Blast[0] == subGFF.iloc[j, 1]]
					if tem.shape[0] > 0:
						if subGFF.iloc[i, 1] in tem[1].tolist():
							tandem.write('%s,%s\n'%(subGFF.iloc[i, 1], subGFF.iloc[j, 1]))
			tandem.flush()

tandem.close()

