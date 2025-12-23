import sys,os
file1 = sys.argv[1]
file2 = sys.argv[2]
pep = open(file1,'r').readlines()
cds = open(file2,'r').readlines()

PEP = {}
i = 0
while i < len(pep):
	inl = pep[i].strip()
	if inl.startswith('>'):
		PEP[inl[1:]] = pep[i + 1].strip()
		i += 1
		i += 1
		while i < len(pep) and not pep[i].startswith('>'):
			PEP[inl[1:]] = PEP[inl[1:]] + pep[i].strip()
			i += 1
		i = i - 1
	i += 1

CDS = {}
i = 0
while i < len(cds):
	inl = cds[i].strip()
	if inl.startswith('>'):
		CDS[inl[1:]] = cds[i + 1].strip()
		i += 1
		i += 1
		while i < len(cds) and not cds[i].startswith('>'):
			CDS[inl[1:]] = CDS[inl[1:]] + cds[i].strip()
			i += 1
		i = i - 1
	i += 1

out = open('%s_align.fas'%file2.split('.fas')[0], 'w')
for gene in PEP:
	nucle = ''
	n = 0
	for i in range(0,len(PEP[gene])):
		if PEP[gene][i] == '-':
			nucle = nucle + '---'
		else:
			nucle = nucle + CDS[gene][(3*n):(3*n+3)]
			n += 1
	out.write('>%s\n%s\n'%(gene, nucle))

out.close()