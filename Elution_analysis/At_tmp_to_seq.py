#NOTE: This code will only run with Pandas 0.18 or higher (str.split function)

import sys
#importing local updated Pandas
sys.path.append("/home/nag2378/.local/lib/python2.7/site-packages")
import pandas as pd
from Bio import SeqIO
from regex import search,findall,finditer
from pandas import Series,DataFrame
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


tmp = pd.read_csv("elution_peptides_tmp.csv")
multistart = tmp['Start'].str.split(',', expand=True).stack()
multistart.index = multistart.index.droplevel(-1)
multistart.name = 'Start'

multiend = tmp['End'].str.split(',', expand=True).stack()
multiend.index = multiend.index.droplevel(-1)
multiend.name = 'End'
tmp.drop(['Start', 'End'], inplace=True, axis=1)

multiboth = pd.DataFrame({'Start':multistart, 'End':multiend}, index=multistart.index)
final = tmp.join(multiboth)
final = final[['Protein','Peptide','Fraction','Count','Sequence','Length','Start','End']]





final = final.set_index(['Protein','Peptide'])
final['Appearance'] = final.groupby(final.index).cumcount() + 1

print final.ix['sp|Q9FKA5|Y5957_ARATH'].ix['KPSYGR']

final.to_csv("elution_peptides_numerical.csv")
