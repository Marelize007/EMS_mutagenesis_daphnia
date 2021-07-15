# EMS_mutagenesis_daphnia

The following repository containes the scripts used to identify EMS induced mutations as well as the commands used to functionally annotate the mutations using SNPEff. The scripts are as follow:

1. The following BWA (Li and Durbin, 2010) command was used along with default parameters to align the fastq files of mutant lines to the D. pulex and D. pulicaria reference genomes:

bwa mem reference_genome.fa inputfile_1.fq.gz inputfile_2.fq.gz > output.sam
    
2. The mpileup and call functions of BCFtools (Li, 2011) along with default parameters were used to generate genotype likelihoods and genotype calls in a VCF file containing all EMS mutant lines derived from each natural Daphnia isolate. We added the following FORMAT and INFO tags to the VCF file: AD (allelic depth), DP (number of high-quality bases), ADF (allelic depth on forward strand) and ADR (allelic depth on reverse strand). An examle output file, output_file.vcf, can be viewed in this repository.
    
bcftools mpileup -Ou -f reference_genome.fa -a INFO/AD,FORMAT/AD,FORMAT/DP,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,FORMAT/SCR  input_file_1.bam  input_file_2.bam input_file_3.bam | bcftools call --threads 20  -mO z -o output_file.vcf.gz

3. After aditional filtering with BCFtools, a custom python script, filter_mutations.py, was used to further filter the mutations based on the forward and reverse read support. This number can be manually adjusted in the script as well as the amount of samples that can share a mutation. The output files are a VCF file containing the final set of mutations, a text file containing the breakdown of the types of base substitutions, and a text file containing all of the mutations that failed the filtering criteria.

4. SNPEff was used to functionally annotate the final list of mutations. The input file is the filtered mutations to be annotated and the output file contains the functionally annotated mutations. For more information on SNPEff lease see https://pcingola.github.io/SnpEff/.

java -Xmx4g -jar snpEff.jar -v -cancer -cancerSamples untreated_treated_samples.txt ref_genome input.vcf  > annotated_output.vcf


