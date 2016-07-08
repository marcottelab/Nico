import pandas as pd
from Bio import SeqIO
from regex import search,findall,finditer
from pandas import Series,DataFrame
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
proteome = list(SeqIO.parse("uniprot-proteome%3AUP000006548.fasta","fasta"))

#import experimental results from LC/MS
data = pd.read_csv("elution_peptides_arath.csv",sep=' ',header=0,names=['Protein','Peptide','Fraction','Count'])#.sortlevel(1)
data =  data.dropna().set_index('Protein')


#create a list with the IDs and sequences of the proteins in the given proteome
list_seq = []
for protein in proteome:
	list_seq.append([str(protein.id),str(protein.seq)])

#turn that list into a dataframe
sequences = DataFrame(list_seq,columns=['Protein','Sequence']).set_index('Protein')

#append sequences to the original dataframe
tmp = data.join(sequences)

#return starting position of peptide with respect to the protein on the same row
def pep_start(row):
	#correct leucine-isoleucine insensibility
	pep = row['Peptide'].replace('I','J').replace('L','J').replace('J','(I|L)')
	seq = row['Sequence']
	#if it matches more than once, return a list of positions
	if len(findall(pep,seq)) > 1:
		runs = finditer(pep,seq)
		coord = []
		for match in runs:
			coord.append(match.start()+1)
		return coord
	elif len(findall(pep,seq)) == 1:
		return [search(pep,seq).start()]
	else: return 'Not found'
				
#idem last function, end position
def pep_end(row):
	pep = row['Peptide'].replace('I','J').replace('L','J').replace('J','(I|L)')
        seq = row['Sequence']
	if len(findall(pep,seq)) > 1:
		runs = finditer(pep,seq)
		coord = []
		for match in runs:
			coord.append(match.end())
		return coord
	elif len(findall(pep,seq)) == 1:
		return [search(pep,seq).end()]
	else: return 'Not found'

#add the three columns using the functions previously defined
tmp['Start'] =  tmp.apply(pep_start,axis=1)
tmp['End'] = tmp.apply(pep_end,axis=1)
tmp['Length'] = tmp['Sequence'].map(len)

final = tmp.reset_index()

print (final)
print len(final[len(final.Start) > 1])
print len(final[len(final.Start) == 1])
#save all the output
final.to_csv('elution_peptides_positions_arath.csv',sep=',')
#save all the proteins except for those that contain a peptide sequence more than once
final[final.Start < 5000].to_csv('elution_peptides_positions_arath_dropdups.csv',sep=',')
