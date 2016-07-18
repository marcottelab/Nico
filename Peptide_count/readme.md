Nico/Peptide_count
=============

**addseq.py**  

        W: This script pretended to do what At_addseq.py, but
        there were too many mismatches between the complex map
        and the experimental results (in *H. sapiens*), so what 
        it actually outputs is a list of experimentally identified 
        peptides that couldn’t be found in the original sequences        
        according to the proteome (I1).

        I1: Hs_complex_map.fasta (sequence of the identified         proteins)  
        I2: Hs_helaN_1003.pep_list_FDR0010 (LC/MS2 results post-identification)
        O1: notfound.tsv 
 
**At_addseq.py**

        W: Takes the experimental results of a LC/MS2 assay of *A.thaliana*, outputs a csv file wherein the start and end of each peptide found is added, in the         format “[a,b,c]”. O2 does not include peptides that appear more than once in a single protein.

        I1:        uniprot-proteome%3AUP000006548.fasta
        I2:        elution_peptides_arath.csv
        O1: elution_peptides_positions_arath.csv
        O2: elution_peptides_positions_arath_dropdups.csv

**At_Col_0_indark_201505_elution.csv**

**At_coverage.py**

        W:
        I: elution_peptides_positions_arath_dropdups.csv
        O: coverage_arath.txt

**At_multicov.py**

        W: Given the output of At_addseq.py, create a tsv in which it is specified how much of a protein is covered by the peptides. Unlike **At_coverage.py**, this script considers
        I: elution_peptides_positions_arath.csv
        O: coverage_arath.tsv

**At_raw_to_seq.py**

        W: Writes a CSV in which, just like in **At_addseq.py**,
        the starting and ending positions of each peptide are 
        added. Unlike in that file, multiple appearances of a single
        peptide are split into separate rows and each appearance
        is ordered by an extra column, “Appearance”. Start and
        End columns are not lists but a number.
        I1: uniprot-proteome%3AUP000006548.fasta
        I2: elution_peptides_arath.csv
        O: elution_peptides_numerical.csv

**At_raw_to_tmp.py**

        W: Same as **At_addseq.py**, but it doesn’t write 
        brackets around the numbers (**At_multicov** can’t read 
        it literally). This will be fixed.
        I: uniprot-proteome%3AUP000006548.fasta
        O: elution_peptides_tmp.csv

**At_tmp_to_seq.py**

        W: It takes the output of **At_raw_to_tmp.py** and splits the lines, adding the Appearance column. It performs the “second half” of what **At_raw_to_seq.py** does.

        I: elution_peptides_tmp.csv
        O: elution_peptides_numerical.csv

**parseproteome.py**

        W: Filter proteins whose ID ends in 1 and delete “.1” in their IDs so that it can be compared to other databases.
        I: Homo_sapiens.GRCh38.pep.all.fa
        O: unique_proteins.fasta
**peptidecount2.py**

        W: Converts fractions IDs to numbers, add missing combinations of fraction-peptide. Output a .tsv with the data frame. The method proved to be inefficient and a better result was obtained by using R’s Tidyr and Diplyr.
        I: Hs_helaN_1003.pep_list_FDR0010
        O: Hs_helaN_1003_reindexed.tsv
**peptidecount.py**

        Arg: Last six digits of a protein name (ENSEMBL).
        W: Provides the spectral count for each of the peptides of         
        the protein specified
        I: Hs_helaN_1003.pep_list_FDR0010 (Experimental results of mass-spec post assignation
        O: No files. Matplotlib graph.

get_elution_ids.py

	Origin: Claire’s analogous script to add_seq / raw_to_seq. Grouping algorithm.

Homo_sapiens.GRCh38.pep.all.fa

Hs_addseq.py

coverage_arath_dropdups.txt

	Origin: Renamed version of coverage_arath.txt. Only 
        includes peptides that appear once in a protein.

coverage_arath.txt

elution_peptides_arath.csv

elution_peptides_At_Col_0_indark_201505_elution_arath.csv

elution_peptides_numerical.csv

elution_peptides_positions_arath.csv

elution_peptides_positions_arath_dropdups.csv

elution_peptides_positions_numerical_arath.csv

	Origin: 

elution_peptides_tmp.csv


Hs_helaN_1003.pep_list_FDR0010

nohup.out

notfound.tsv

	Origin: addseq.py


Protein_list.tsv

uniprot-proteome%3AUP000006548.fasta

unique_proteins.fasta

