import pandas as pd
import numpy as np
from numpy import union1d,asarray
from pandas import DataFrame,Series
from regex import search,findall,finditer


table = pd.read_csv('elution_peptides_positions_arath_dropdups.csv',sep=',').ix[:,1:]




#peptides that appear more than once in a protein
proteins = list(set(list(table.Protein)))

tmp = table[table.Protein == proteins[18]]
print tmp
exit()

def coverage(table):
	#return a dataframe with the coverage of each individual peptide in a protein
	interval = (table.apply(lambda row : range(int(row['Start']),int(row['End'])+1),axis=1))#.apply(lambda row: row[:-3])
	#if there is only one peptide, return the range between its start and end positions
	if len(table) == 1: return asarray(range(int(table['Start']),int(table['End'])+1))
	#if there are more, unite all the intervals
	if len(table) > 1:
		return reduce(union1d,(list(interval)))

print


covdir = {}
i = 0 #counter
for prot in proteins:
	i+=1
	print '>'+str(i)+'\t'+prot
	data = table[table.Protein == prot]
	covdir[prot] = coverage(data)
	print coverage(data)
	print '###################################################'
exit()
outfile = open("coverage_arath.txt", "w")

for key,value in covdir.iteritems():
	output = str(key) + "\t" + str(value) +"\n"
 
        outfile.write(output)

outfile.close()
