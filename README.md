# EMS_mutagenesis_daphnia

The following repository containes the scripts used to identify EMS induced mutations as well as the commands used to functionally annotate the mutations using SNPEff. The scripts are as follow:

1. The following BWA (Li and Durbin, 2010) command was used along with default parameters to align the fastq files of mutant lines to the D. pulex (Ye et al. 2017) and D. pulicaria (Jackson et al. 2021) reference genomes:

bwa mem reference_genome.fa inputfile_1.fq.gz inputfile_2.fq.gz > output.sam

2. Samtools (Li et al. 2009) was used to remove reads that mapped to multiple locations and Picard tools's (http://broadinstitute.github.io/picard/) MarkDuplicates function was used to locate and tag PCR duplicates. 

samtools view -h input.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b -o output.bam
java -Xmx2g -XX:MaxMetaspaceSize=256m -jar picard.jar MarkDuplicates INPUT=input.bam OUTPUT=output.bam METRICS_FILE=file.metric
    
3. The mpileup and call functions of BCFtools (Li, 2011) along with default parameters were used to generate genotype likelihoods and genotype calls in a VCF file containing all EMS mutant lines derived from each natural Daphnia isolate. We added the following FORMAT and INFO tags to the VCF file: AD (allelic depth), DP (number of high-quality bases), ADF (allelic depth on forward strand) and ADR (allelic depth on reverse strand). An examle output file, output_file.vcf, can be viewed in this repository.
    
bcftools mpileup -Ou -f reference_genome.fa -a INFO/AD,FORMAT/AD,FORMAT/DP,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,FORMAT/SCR  input_file_1.bam  input_file_2.bam input_file_3.bam | bcftools call --threads xx  -mO z -o output_file.vcf.gz

For this study only biallelic single nucleotide polymorphisms (SNPs) that met the following filter parameters were used:
    * Quality score (QUAL) >= 20
    * Sequencing depth (DP) >= 10
    * Distance >= 50 bp from an indel
    ** Indels were not examined in this study and were also filtered out.

4. After additional filtering with BCFtools, a custom python script, filter_mutations.py, was used to further filter mutations based on a consensus method and forward and reverse read support. For each SNP site a consensus genotype call is established using a majority rule. If a VCF file has a total of N samples, N-1 of those samples need to be in agreement to establish a consensus genotype call. If an EMS sample shows a genotype call different from the consensus, a tentative mutation is called. Since EMS has been shown to induce mutations randomly into the genome, it is highly unlikely that the same mutation will be observed amongst multiple EMS treatment lines. Each genotype call also has to be supported by at least 2 forward and 2 reverse reads in order to limit false positives due to sequencing errors. Both the number of samples needed for the consensus genotype and the read support can be manually adjusted. The output files are a VCF file containing the final set of mutations, a text file containing the breakdown of the types of base substitutions, and a text file containing all of the mutations that failed the filtering criteria.

5. SNPEff was used to functionally annotate the final list of mutations. The input file is the filtered mutations to be annotated and the output file contains the functionally annotated mutations. For more information on SNPEff please see https://pcingola.github.io/SnpEff/.

java -Xmx4g -jar snpEff.jar -v -cancer -cancerSamples untreated_treated_samples.txt ref_genome input.vcf  > annotated_output.vcf


