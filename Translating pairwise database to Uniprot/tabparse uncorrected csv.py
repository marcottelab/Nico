from pandas import Series,DataFrame
import pandas as pd
import csv


dic = {}
listp = []

with open("arath_uniprot.tab",'rb') as ref:
	for row in csv.reader(ref, dialect="excel-tab"):
		entry,entryname, yourlist = row
		#dic[yourlist]=entry
		listp.append([yourlist,entry])
dp = DataFrame(listp,columns=['A','Uniprot A'])

dp2 = dp.set_index('A')
#print dp2
proteins = []
with open("arath_biogrid.tab","rb") as grid:
	for row in csv.reader(grid, dialect="excel-tab"):	
		systa, systb = row
		proteins.append([systa,systb])
		#if systa in dic and systb in dic:
		#	interactions.append([dic[systa],dic[systb]])
		#else: interactions.append([systa,systb])

prot = DataFrame(proteins,columns=['A','B'])
prot2 = prot[(prot.A != '-') & (prot.B != '-') & (prot.B.str[0:2] == 'AT') 
						& (prot.A.str[0:2] == 'AT')]


prot3 = prot2.set_index('A')
#print prot3
#print dp2
#print prot3.columns
#print dp2.columns

join1 = prot3.join(dp2,how='left')
#print join1
join1b = join1.reset_index()
join1c = join1b.set_index('B')


dp.columns = ['B','Uniprot B']

dic2 = dp.set_index('B')

join2 = join1c.join(dic2,how='left')

final =  join2.reset_index().ix[:,[2,3]]
print final

final.to_csv('arath_pwi.tsv', sep='\t' , header=False , index=False)


