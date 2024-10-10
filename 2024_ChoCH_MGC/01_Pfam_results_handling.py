# this script was used to parse the hmmscan results for protein sequences.
import sys,os
import pandas as pd
df = pd.read_csv('Tair_dom_02.tbl',sep='\t')
gene = list(df['Gene'].unique())
out = open("Tair_domain_information_parsed_final.txt",'w')
out.write('%s\n'%df.columns.str.cat(sep="\t"))
for g in gene:
	sub_df = df[df['Gene']==g]
	sub_df = sub_df.sort_values(by='Ali_start') # sort the table based on the Ali_start
	index_kept = [] # saving the kept index
	sub_df_original = sub_df
	sub_df_next = sub_df
	while len(set(index_kept)) != sub_df_original.shape[0]: # till there is no overlapping regions in the kept table
		sub_df_original = sub_df_next # fix the initial table, which is used to check whether the table was changed after the loop
		sub_df = sub_df_next
		index_kept = []
		while sub_df.shape[0] > 0:
			if sub_df.shape[0] == 1:
				index_kept.append(sub_df.index[0])
				break
			else:
				L1 = sub_df.iloc[0,:]
				n = 0
				sub_df_tem = sub_df
				for i in range(1,sub_df_tem.shape[0]):
					L2 = sub_df_tem.iloc[i,:]
					if (L2.iloc[6] - L1.iloc[7]) * (L2.iloc[7] - L1.iloc[6]) < 0 : # columns 6,7 are the ali coordinates
						if L1.iloc[3] < L2.iloc[3]: # 3 is the column number of E-value
							sub_df = sub_df.drop(sub_df_tem.index[i])
						if L1.iloc[3] > L2.iloc[3]:
							sub_df = sub_df.drop(sub_df_tem.index[0])
							break
						if L1.iloc[3] == L2.iloc[3]: # when the evalues are the same, keep the one with longer overlapping with hmm domain
							if abs(L1.iloc[5] - L1.iloc[4]) > abs(L2.iloc[5] - L2.iloc[4]):
								sub_df = sub_df.drop(sub_df_tem.index[i])
							if abs(L1.iloc[5] - L1.iloc[4]) < abs(L2.iloc[5] - L2.iloc[4]):
								sub_df = sub_df.drop(sub_df_tem.index[0])
								break
							if abs(L1.iloc[5] - L1.iloc[4]) == abs(L2.iloc[5] - L2.iloc[4]): # if the overlapped regions have the same length, then keep the first one
								sub_df = sub_df.drop(sub_df_tem.index[i])
						n = 1
				if n == 0: #indicate the first line was not overlapping with any other lines
					index_kept.append(sub_df.index[0])
					sub_df = sub_df.drop(sub_df.index[0]) # and then remove it from the table
		sub_df_next = df.loc[set(index_kept),:] # after each loop, update the table
	sub_df_next = sub_df_next.sort_values(by='Ali_start')
	sub_df_next = sub_df_next.astype(str)
	for i in range(0,sub_df_next.shape[0]):
		out.write('%s\n'%sub_df_next.iloc[i,:].str.cat(sep="\t")) 
		out.flush() # reflesh the results after every loop

out.close()

