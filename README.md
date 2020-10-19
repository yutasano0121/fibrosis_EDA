# Code for anlyzing Adams et al., (Science Advances, 2020) scRNA-seq data

## Files
### QC'ed with gene count and mitochondrial gene ratio
* counts_QCed_1000gene10MT.mtx
* annotation_QCed_1000gene10MT.csv
### Data whose sample size was reduced to 10% after QC
* counts_QCed_1000gene10MT_reduced.mtx
* annotation_QCed_1000gene10MT_reduced.csv
### List of proteins/genes involved in ECM
* ECM_matrisome_hs_masterlist.csv
<br />\* Fetched from http://matrisomeproject.mit.edu/
<br />\* Better than MatrixDB as it has correct human gene symbols
### List of interaction involving ECM molecules in human (Uniprot keywords and gene symbols)
* ECM_interaction.csv
<br />\* Fetched from http://matrixdb.univ-lyon1.fr/
<br />\* Some gene symbols are not standard, which need to be modified
