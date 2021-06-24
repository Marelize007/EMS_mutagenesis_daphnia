# EMS_mutagenesis_daphnia

The following repository containes the scripts used to identify EMS induced mutations as well as the commands used to functionally annotate the mutations using SNPEff. The scripts are as follow:

1. filter_mutations.py
This is a custom python script. The input is a VCF file containing all the tentative mutations for the treatment lines. The number of forward and reverse read support can be adjusted as well as the amount of samples that can share a mutation. The output files are a VCF file containing the final set of mutations, a text file containing the breakdown of the types of base substitutions, and a text file containing all of the mutations that failed the filtering criteria.

2. SNPEff
The command line used to functionally annotate the final list of mutations.
