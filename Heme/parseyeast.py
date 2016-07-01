from sys import argv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import regex as re

script, org = argv

dna_seqs = list(SeqIO.parse(org+" heme sequences.fasta","fasta"))
orfs = []
for seq in dna_seqs:
	orfs.append(SeqRecord(seq=Seq(max(re.findall(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)',str(seq.seq)), key = len)),id=seq.id))


#translate each sequence
peps = []
for i in range(len(orfs)):
	peps.append(SeqRecord(seq=orfs[i].seq.translate(),id=orfs[i].id))

#take the ids of the proteins and find the corresponding fasta files,
uniprot = []
lengths=[]

txt = open("translated_"+org+".txt","w")
for i in peps:
	try:
		uniprot.append(SeqIO.read("All/"+i.id+".fasta","fasta").seq)
		lengths.append(len(SeqIO.read("All/"+i.id+".fasta","fasta").seq))
		print "Sequence "+i.id+"\n\tExpected: "+str(len(SeqIO.read("All/"+i.id+".fasta","fasta").seq))+"\n\tObserved: "+str(len(i.seq)-1)+'\n'
		file.write
	except:
		print str(i.id)+" not found in the master folder.\n"


with open("translated_"+org+".txt","w") as file:
	for i in range(len(peps)):
		file.write(peps[i].id)
		file.write('\n')
		file.write(str(peps[i].seq)[:-1])
		file.write('\n'*2)	
print "The translated sequences of the experimental proteins have been stored in a file called 'translated_%s.txt'\n" %(org)

print "All done!"
