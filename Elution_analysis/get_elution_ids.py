import pandas as pd
import sys
import argparse


def create_tables(species, level, orthology_file, elution_file, peptide_file):


    experimentname = elution_file.replace(".csv", "")
    experimentname = experimentname.split("/")[-1]

    
    peps = pd.DataFrame(pd.read_csv(peptide_file)) 
    print "Peptide file loaded"
    #Remove undistinguishable isoleucine/leucine
    peps['Peptide'] = peps['Peptide'].str.replace('I', 'J')
    peps['Peptide'] = peps['Peptide'].str.replace('L', 'J')


    peps = peps.set_index(['ProteinID'])
 
    print peps
    
    
    frac = pd.DataFrame(pd.read_csv(elution_file))
    #Remove undistinguishable isoleucine/leucine
    frac['Peptide'] = frac['Peptide'].str.replace('I', 'J')
    frac['Peptide'] = frac['Peptide'].str.replace('L', 'J')
    frac['Peptide'] = frac['Peptide'].str.replace('*', '')

    frac = frac.set_index(['Peptide'])
    
    prot = pd.DataFrame(pd.read_csv(orthology_file, sep="\t"))
    prot = prot.set_index(['ProteinID'])
    
    
    print peps
    #print prot
    
    group_pep = prot.join(peps, how = "left")
    print "group_pep", group_pep.shape
    group_pep = group_pep.reset_index()
    
    group_pep_cols = ['GroupID', 'Peptide']
    
    group_pep = group_pep[group_pep_cols]
    
    
    #Get unique Rows, as one peptide can occur multiples in one group
    uniq_group_pep = group_pep.drop_duplicates()
    print "nonredundant_group_pep", uniq_group_pep.shape
    nonred_group = "peptide_assignments/" + species + "/nonredundant_orthogroup_" + species + "_" + level + ".csv"
    uniq_group_pep.to_csv(nonred_group)
    
    #uniq_group_pep = final_group_pep.set_index(['Peptide'])
    #Get Peptides which only occur in one group
    final_group_pep = uniq_group_pep.drop_duplicates(subset=['Peptide'], keep=False) #current Docs/version have subset
    print "final group_pep", final_group_pep.shape
    
    #This is currently same between group and protein identifications.

    ident_group = "peptide_assignments/" + species + "/identifying_orthogroup_" + species + "_" + level + ".csv"
    final_group_pep.to_csv(ident_group)
    
    final_group_pep = final_group_pep.set_index(['Peptide'])
    
    #Join peptides in experiment to unique Group Identifying peptides
    frac_group = frac.join(final_group_pep, how='left')
    print "frac_group", frac_group.shape
   
    frac_group = frac_group.reset_index()
    group_identified_peptides = frac_group[['Peptide']].drop_duplicates()
    identified_group = "peptide_assignments/" + species + "/identifiedpeps_orthogroup_" + "_" + experimentname + "_" + species + "_" + level + ".csv"
    group_identified_peptides.to_csv(identified_group)
    
        
    frac_group = frac_group[['ExperimentID', 'FractionID', 'GroupID', 'PeptideCount']]
    
    grouped_frac_group = frac_group.groupby(by=['ExperimentID', 'FractionID', 'GroupID'])['PeptideCount'].sum()
    
    final_frac_group =  grouped_frac_group.reset_index(name='Total_SpecCounts')
   
    print "number of spectral counts", final_frac_group['Total_SpecCounts'].sum()

    elut_group = "identified_elutions/" + species + "/" + experimentname +"_" + species + "_" + level + ".csv"
    final_frac_group.to_csv(elut_group)


    #ungrouped analysis    
    protein_pep = prot.join(peps, how = "left")
    print "protein_pep", protein_pep.shape
   
    protein_pep = protein_pep.reset_index()
    
    
    protein_pep_cols = ['ProteinID', 'Peptide']
    
    protein_pep = protein_pep[protein_pep_cols]
    
    
    #Get unique Rows, as one peptide can occur multiples in one protein
    uniq_protein_pep = protein_pep.drop_duplicates()
    print "nonredundant protein_pep", uniq_protein_pep.shape
    nonred_prot = "peptide_assignments/" + species +"/nonredundant_protein_" + species + ".csv"
    uniq_protein_pep.to_csv(nonred_prot)

    ###This is just dropping subsequent appearances. It's not removing anything with a duplicate. 

  
    #Get Peptides which only occur in one protein
    final_protein_pep = uniq_protein_pep.drop_duplicates(subset=['Peptide'], keep=False) #current Docs/version have subset
    print "final protein_pep", final_protein_pep.shape
    ident_prot = "peptide_assignments/" + species + "/identifying_protein_" + species + ".csv"
  
    final_protein_pep.to_csv(ident_prot)
    
    final_protein_pep = final_protein_pep.set_index(['Peptide'])
    
    #Join peptides in experiment to unique Group Identifying peptides
    frac_protein = frac.join(final_protein_pep, how='left')
    print "frac_protein", frac_protein.shape
   
    frac_protein = frac_protein.reset_index()
    protein_identified_peptides = frac_protein[['Peptide']].drop_duplicates()
    identified_prot = "peptide_assignments/" + species + "/identifiedpeps_protein_" + "_" + experimentname + "_" + species + ".csv"

    protein_identified_peptides.to_csv(identified_prot)
    
        
    frac_protein = frac_protein[['ExperimentID', 'FractionID', 'ProteinID', 'PeptideCount']]
    
    grouped_frac_protein = frac_protein.groupby(by=['ExperimentID', 'FractionID', 'ProteinID'])['PeptideCount'].sum()
    
    final_frac_protein =  grouped_frac_protein.reset_index(name='Total_SpecCounts')
    print "number of spectral counts", final_frac_protein['Total_SpecCounts'].sum()

    elut_prot = "identified_elutions/" + species + "/" + experimentname +"_" + species + "_proteins.csv"
    final_frac_protein.to_csv(elut_prot)
    
 
    
    #Test identifiable peptides
    
    #No difference in arabidopsis? 
    #Check every step
    #Grouping condences 28000 identifications to 22000
    
    #claire@Oppenheimer:~/for_mySQL$ awk -F',' '{print $4}' singleIdentified.csv | sort -u | wc -l
    #10381
    #claire@Oppenheimer:~/for_mySQL$ awk -F',' '{print $4}' GroupIdentified.csv | sort -u | wc -l
    #6015
    
    #So, grouping is getting fewer groups, but average large size groups. Grouping get 20% more identifying peptides
    
    
    
    
    
    
parser = argparse.ArgumentParser(description='Interpret mass spec experiments using orthologous groups of proteins to make identifications')

parser.add_argument('species_code', action="store", type=str)
parser.add_argument('phylogenetic_level', action="store", type=str)
parser.add_argument('orthology_file', action="store", type=str)
parser.add_argument('elution_file', action="store", type=str)
parser.add_argument('peptides_file', action="store", type=str)
 

inputs = parser.parse_args()

create_tables(inputs.species_code, inputs.phylogenetic_level, inputs.orthology_file, inputs.elution_file, inputs.peptides_file)






    
