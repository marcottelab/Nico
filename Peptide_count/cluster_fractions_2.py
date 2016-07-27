from sys import argv
import numpy as np
import pandas as pd
import regex as re
from pandas import Series,DataFrame,read_csv
from sklearn.preprocessing import normalize
from sklearn.cluster import MeanShift,estimate_bandwidth

script, protein_name = argv

pd.set_option('display.width',120)
data = read_csv("elution_peptides_numerical.csv",sep=',')

#default protein = tr|F4KHD5|F4KHD5_ARATH 
data_prot =  data[data.Protein==protein_name].reset_index()
data_prot2 =  data_prot[['Peptide','Fraction','Count','Start']].sort_values('Start')
data_wide = data_prot2.reset_index(drop=True).pivot('Peptide','Fraction','Count').fillna(0)
data_wide = data_wide.join(data_prot2.set_index('Peptide')['Start']).sort_values('Start').drop_duplicates()

#print data_prot2.reset_index(drop=True)#.set_index(['Peptide','Fraction'])['Start']


#must fix duplicates for it to work in other proteins

#print data_wide
print data_wide.ix[:,-1]
#observation matrix
matrix= data_wide.ix[:,:-1].as_matrix()



#run MeanShift  on matrix
ms = MeanShift(cluster_all=False,bandwidth=1.1)

print estimate_bandwidth(normalize(matrix), quantile=0.3, random_state=0)

print ms.fit(normalize(matrix))

cluster_centers = ms.cluster_centers_
labels = ms.labels_

#print centers
print (Series(labels)).shape
print (data_wide).shape

print labels
exit()
########################################################
#when I get a robust clustering I will finish this script

data_wide = data_wide.reset_index()

#add Label column according to "labels"
data_wide['Label'] = Series(labels)

data_wide.index = data_wide.Peptide 
print data_wide.head(1)['Label']

data_prot.index = data_prot.Peptide
labeled_data = data_prot.join(data_wide.Label,how='left').reset_index(drop=True)

print labeled_data.columns

file_name = re.split(r'\|',protein_name)[-1]
labeled_data.to_csv('labeled_peptides_'+file_name+'.csv',sep=',')
	
print 'Labeled data stored in labeled_peptides_'+file_name+'_ms.csv'

