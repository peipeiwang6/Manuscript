import os,sys
inp=sys.argv[1]
top=sys.argv[2]
inp_f=open(inp,"r").readlines()
out=open(inp+"_top"+top+".txt","w")
gene={}
for line in inp_f:
    if  not line.startswith("#"):
        query=line.split("\t")[0].strip()
        subject = line.split("\t")[1].strip()
        if query not in gene.keys():
            gene[query]=[subject]
            out.write("%s\n"%line.strip())
        else:
            if subject not in gene[query]:
                gene[query].append(subject)
                if len(gene[query])<=int(top):
                    out.write("%s\n"%line.strip())
            else:
                continue

