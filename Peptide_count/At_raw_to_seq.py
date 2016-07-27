#NOTE: This code will only run with Pandas 0.18 or higher (str.split function)
#NOTE: To compute total count, groupby rows where Appearance==1 or peptides that show up more than once in a protein will overincrease the total Count value

import sys 
#importing local updated Pandas module
sys.path.append("/home/nag2378/.local/lib/python2.7/site-packages")
import pandas as pd
from Bio import SeqIO
from regex import search,findall,finditer
from pandas import Series,DataFrame
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq 

def pep_start(row):
        #correct leucine-isoleucine insensibility
        pep = row['Peptide'].replace('I','J').replace('L','J').replace('J','(I|L)')
        seq = row['Sequence']
        #if it matches more than once, return a list of positions
        if len(findall(pep,seq)) > 1:
                runs = finditer(pep,seq)
                coord = []
		string = ''
		firstart = True
                for match in runs:
			if firstart:
				string = str(match.start()+1)
				firstart = False
			else:	string+= ','+str(match.start()+1)
                return string
        elif len(findall(pep,seq)) == 1:
                return str(search(pep,seq).start())
        else: return 'Not found'

#idem last function, end position
def pep_end(row):
        #correct leucine-isoleucine insensibility
        pep = row['Peptide'].replace('I','J').replace('L','J').replace('J','(I|L)')
        seq = row['Sequence']
        #if it matches more than once, return a list of positions
        if len(findall(pep,seq)) > 1:
                runs = finditer(pep,seq)
		string = ''
		firstend = True
                for match in runs:
			if firstend:
				string = str(match.end())
				firstend = False
			else:	
				string+= ','+str(match.end())
                return string
        elif len(findall(pep,seq)) == 1:
                return str(search(pep,seq).end())
        else: return 'Not found'

proteome = list(SeqIO.parse("uniprot-proteome%3AUP000006548.fasta","fasta"))

#import experimental results from LC/MS
data = pd.read_csv("elution_peptides_arath.csv",sep=' ',header=0,names=['Protein','Peptide','Fraction','Count'])
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
#add the three columns using the functions previously defined
tmp['Start'] =  tmp.apply(pep_start,axis=1)
tmp['End'] = tmp.apply(pep_end,axis=1)
tmp['Length'] = tmp['Sequence'].map(len)
tmp = tmp.reset_index().set_index(['Protein','Peptide'])
#save a tmp
tmp.to_csv('elution_peptides_tmp.csv',sep=',')

#open tmp
tmp = pd.read_csv("elution_peptides_tmp.csv")

#split multiple start sites
multistart = tmp['Start'].str.split(',', expand=True).stack()
multistart.index = multistart.index.droplevel(-1)
multistart.name = 'Start'

#split multiple end sites
multiend = tmp['End'].str.split(',', expand=True).stack()
multiend.index = multiend.index.droplevel(-1)
multiend.name = 'End'
tmp.drop(['Start', 'End'], inplace=True, axis=1)

#add new sites to original dataframe
multiboth = pd.DataFrame({'Start':multistart, 'End':multiend}, index=multistart.index)
final = tmp.join(multiboth)
#reorder columns
final = final[['Protein','Peptide','Fraction','Count','Sequence','Length','Start','End']]

#reindex to add Appearance
final = final.set_index(['Protein','Peptide','Fraction'])
#add 'Appearance' column for peptides that appear multiple times
final['Appearance'] = final.groupby(final.index).cumcount() + 1 
#reindex to original index
final = final.reset_index().set_index(['Protein','Peptide'])

#test peptide with multiple appearances
print final.ix['sp|Q9FKA5|Y5957_ARATH'].ix['KPSYGR']

final.to_csv("elution_peptides_numerical.csv")

print "All done! Output file: elution_peptides_numerical.csv"
