import pandas as pd
from Bio import SeqIO
from regex import search,findall,finditer
from pandas import Series,DataFrame
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
proteome = list(SeqIO.parse("Hs_complex_map.fasta","fasta"))

#import experimental results from LC/MS
data = pd.read_csv("Hs_helaN_1003.pep_list_FDR0010",sep='\t',names=['Protein','Fraction','Peptide','Count'],index_col=['Protein']).sortlevel(1)


list_seq = []

for protein in proteome:
	list_seq.append([str(protein.description[-15:]),str(protein.seq)])



sequences = DataFrame(list_seq,columns=['Protein','Sequence']).set_index('Protein')


tmp = data.join(sequences)

def pep_start(row):
	pep = row['Peptide'].replace('I','J').replace('L','J').replace('J','(I|L)')
	seq = row['Sequence']
	if len(findall(pep,seq)) > 1:
		runs = finditer(pep,seq)
		coord = []
		for match in runs:
			coord.append(match.start()+1)
		return coord
	elif len(findall(pep,seq)) == 1:
		return search(pep,seq).start()
	else: return 'Not found'
				

def pep_end(pep,seq):
	if len(findall(pep,seq)) > 1:
		runs = finditer(pep,seq)
		coord = []
		for match in runs:
			coord.append(match.end())
		return coord
	elif len(findall(pep,seq)) == 1:
		return search(pep,seq).end()
	else: return 'Not found'


#tmp2 = tmp.ix[1040:1070,:]

print

#print tmp[(tmp.index =='ENSP00000453801') & (tmp.Count == 0.5)]
tmp['Start'] =  tmp.apply(pep_start,axis=1)

tmp[tmp.Start == 'Not found'].ix[:,['Peptide','Count','Start']].to_csv("notfound.tsv",sep="\t",header=False,index=False)

