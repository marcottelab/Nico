from pandas import Series,DataFrame
import pandas as pd
import csv

dict = pd.read_csv("arath_uniprot.tab",sep='\t',names=['Uniprot A','A'],usecols=
								[0,2])
dict_a = dict.set_index('A')


pairs = pd.read_csv("arath_biogrid.tab",sep='\t',names=['A','B'],usecols=
								[0,1])

pairs_arath = pairs[(pairs.A != '-') & (pairs.B != '-') & (pairs.B.str[0:2] == 'AT') 
						& (pairs.A.str[0:2] == 'AT')]

pairs_arath_a  = pairs_arath.set_index('A')


join1 = pairs_arath_a.join(dict_a,how='left').reset_index()
join1_b = join1.set_index('B')


dict.columns = ['Uniprot B','B']
dict_b = dict.set_index('B')


join2 = join1_b.join(dict_b,how='left')
print join2


final =  join2.reset_index().ix[:,[2,3]]
print final

final.to_csv('arath-biogrid-uniprot-pairwise.tsv', sep='\t' , header=False , index=False)


