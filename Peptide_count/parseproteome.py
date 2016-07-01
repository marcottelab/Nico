from pandas import Series,DataFrame
import pandas as pd
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#import experimental results from LC/MS
data = pd.read_csv("Hs_helaN_1003.pep_list_FDR0010",sep='\t',names=['Protein','Fracti    on','Peptide','Count'],index_col=['Protein','Peptide']).sortlevel(1)

#import proteome from ENSEMBL FTP database, parse it with BioPython
proteome = list(SeqIO.parse("Homo_sapiens.GRCh38.pep.all.fa","fasta"))

unique = [protein for protein in proteome if protein.id[-1]=='1']


for protein in unique:
	protein.id = protein.id[:-2]
	print protein.id

print proteome[0]
print unique[0]

SeqIO.write(unique,'unique_proteins.fasta',"fasta")

