from pandas import Series,DataFrame
import pandas as pd
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

proteome = list(SeqIO.parse("Homo_sapiens.GRCh38.pep.all.fa","fasta"))

#filter proteins whose ID ends in 1
unique = [protein for protein in proteome if protein.id[-1]=='1']

#delete the .1
for protein in unique:
	protein.id = protein.id[:-2]
	print protein.id


#write new FASTA file
SeqIO.write(unique,'unique_proteins.fasta',"fasta")

