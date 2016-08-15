from sys import argv
import numpy as np
import pandas as pd
import regex as re
from pandas import Series,DataFrame,read_csv
from sklearn.cluster import KMeans

script, protein_name = argv

pd.set_option('display.width',120)
data = read_csv("elution_peptides_numerical.csv",sep=',')

#default protein = tr|F4KHD5|F4KHD5_ARATH 
data_prot =  data[data.Protein==protein_name].reset_index()
data_prot2 =  data_prot[['Peptide','Fraction','Count']]
data_wide = data_prot2.pivot('Peptide','Fraction','Count').fillna(0)
data_wide = data_wide

#must fix duplicates for it to work in other proteins

#observation matrix
matrix = data_wide.reset_index().ix[:,1:].as_matrix()



#run kmeans k=2 on matrix
kmeans = KMeans(n_clusters =2)
print kmeans.fit(matrix)

centroids = kmeans.cluster_centers_
labels = kmeans.labels_

#print centroids
print (Series(labels)).shape
print (data_wide).shape

data_wide = data_wide.reset_index()
data_wide['Label'] = Series(labels)
data_wide.index = data_wide.Peptide 
print data_wide.head(1)['Label']

data_prot.index = data_prot.Peptide
labeled_data = data_prot.join(data_wide.Label,how='left').reset_index(drop=True)

print labeled_data.columns

file_name = re.split(r'\|',protein_name)[-1]
labeled_data.to_csv('labeled_peptides_'+file_name+'.csv',sep=',')
	
print 'Labeled data stored in labeled_peptides_'+file_name+'.csv'

