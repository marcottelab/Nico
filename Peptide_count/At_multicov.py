import pandas as pd
import numpy as np
from numpy import union1d,asarray
from pandas import DataFrame,Series
from regex import search,findall,finditer
from ast import literal_eval

#set visualization options for the data: do not summarize arrays with less than 4000 elements (all), fit every node into a line
np.set_printoptions(threshold=4000,linewidth=10000)

table = pd.read_csv('elution_peptides_positions_arath.csv',sep=',').ix[:,1:]


#peptides that appear more than once in a protein
proteins = list(set(list(table.Protein)))
#proteins.remove('sp|Q9LSB4|NAI2_ARATH')
#proteins.remove('tr|Q9LH98|Q9LH98_ARATH')
#proteins.remove('sp|Q9LQ31|JAL4_ARATH')
#proteins.remove('sp|P0CAP5|REM13_ARATH')
#proteins.remove('sp|Q9FKA5|Y5957_ARATH')
#proteins.remove('sp|Q9SHF2|AGO3_ARATH')
#proteins.remove('tr|F4JW79|F4JW79_ARATH')
#proteins.remove('tr|Q9LW12|Q9LW12_ARATH')
 
print type(literal_eval(table[table.Protein == proteins[1]].ix[10743,'Start'])) 

def multicov(row):
	intervals = []
	for i in range(len(literal_eval(row['Start']))):
		print row['Protein'],row['Start']
		intervals.append((range(literal_eval(row['Start'])[i],literal_eval(row['End'])[i]+1)))
	return [reduce(union1d,intervals),'foobar']


#dir = {'Start':[[1,7],[14]],
#	'End':[[5,10],[18]]}
#
#df = DataFrame(dir,columns=['Start','End'])
#print df
#print df.apply(multicov,axis=1).ix[:,0]



def coverage(table):
	#return a dataframe with the coverage of each individual peptide in a protein
	interval = table.apply(multicov,axis=1)
	interval = interval.apply(lambda row : row[0])
	#if there is only one peptide, return the range between its start and end positions
	if len(table) == 1: return asarray(list(interval))[0]
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
	print coverage(data), (coverage(data)).shape
	print '###################################################'

outfile = open("coverage_arath.tsv", "w")

for key,value in covdir.iteritems():
	output = str(key) + "\t" + str(value) +"\n"
 
        outfile.write(output)

outfile.close()
