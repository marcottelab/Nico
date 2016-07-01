from sys import argv
from pandas import Series, DataFrame
from itertools import product
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

script, name = argv
	
pd.options.display.max_rows = 999

data = pd.read_csv("Hs_helaN_1003.pep_list_FDR0010",sep='\t',names=['Protein','Fraction','Peptide','Count'],index_col=['Protein','Peptide']).sortlevel(1)

# input: last 6 digits of a protein name (ENSEMBL). Output: Spectral count for each of its constituent peptides
def prot(name):
	prot = 'ENSP00000'+str(name)
	pep = data.ix[prot]
	return pep

#Input: Fraction name ending in the format "_P1A01", being "_P2C06" the 126th fraction. Output: fraction index.
def fracnum(name):
	val = {'A':0,
		'B':12,
		'C':24,
		'D':36,
		'E':48,
		'F':60,
		'G':72,
		'H':84}
	return int(int(name[-2:])+val[name[-3]]+96*(int(name[-4])-1))



def fill_zeros(list):
	final = [0]*127
	for i in list:
		final[int(i[0])] = i[1]
	return final

#replace fraction names by fraction numbers
data['Fraction'] =  data['Fraction'].map(fracnum)

#Replace the number by the last 6 ENSEMBL digits of the protein you want to test
#from here on we work with the reduced dataset of the protein being tested
test = prot(str(name))

#create a list of all the fractions
fractions = range(1,127)

#create a list of unique peptides for the protein being tested
ind = list(set((test.index.values.tolist())))
#print ind

#produce all the possible combinations in fractions x peptides
idx = [(x, fraction)  for x, fraction in product(ind,fractions)]
#print idx

#reindex the specific protein data set with all the missing fractions, fill them with zero
print test

table = test.reset_index().set_index(['Peptide','Fraction']).reindex(idx,fill_value=0)


#create a dictionary where each entry contains the y-axis data (value) for each peptide (key)
dir = {}
for key,grp in table.groupby(level=0): #given the dataset for the specific protein, group the data by peptides
	dir[key] = (table.ix[key].values.tolist()) #for each group, append its data to a different entry and name it after its peptide


##### PLOT #####
#fig, ax = sns.plt.subplots(len(dir), 1, figsize=(7,5))
#for a,key in zip(ax,dir.keys() ):
#    y = dir[key]
#    n = len(y)
#    x = np.linspace(1,n,n)
#    a.plot(x,y)
    # add labels/titles and such here

plt.show()
