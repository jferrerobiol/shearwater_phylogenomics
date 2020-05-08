#!/bin/sh
#$ -V
#$ -cwd

# 1) trimmomatic.sh: Script to remove the Illumina adapters #

```cd /ddn/data/sbvd77/UCE/raw ## move to folder cotaining the raw files to filter
path=/ddn/data/sbvd77/UCE
mkdir ../cleaned
TruSeq3=$path/info/TruSeq3-PE-2.fa```

## Define the invariable parts of the illumina adapters

i51=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
i52=GTGTAGATCTCGGTGGTCGCCGTATCATT
i71=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
i72=ATCTCGTATGCCGTCTTCTGCTTG

for file in *READ1.fastq.gz; do ## start a for loop to treat all the raw files in the folder

### Define some variables (file2 and name of sample) ###
file2=$(echo $file | sed 's/-READ1/-READ2/')
sample=$(echo $file | cut -d "-" -f1)

### Create the directory jerarchy as required by phyluce ###

mkdir ../cleaned/$sample
mkdir ../cleaned/$sample/raw-reads
mkdir ../cleaned/$sample/split-adapter-quality-trimmed
mkdir ../cleaned/$sample/stats

### Define input and output directories ###

input=$path/raw
output=$path/cleaned/$sample/split-adapter-quality-trimmed

### Create the adapters.fasta file ###

## Define the barcodes ##

i5barcode=$(echo $(zcat $file | grep "^@J00" | cut -d ":" -f 10 | sort | uniq -c | sort -r | head -1 | sed -E 's/[0-9]+ [A-Z]+\+([A-Z]+)/\1/'))
i7barcode=$(echo $(zcat $file | grep "^@J00" | cut -d ":" -f 10 | sort | uniq -c | sort -r | head -1 | sed -E 's/[0-9]+ ([A-Z]+)\+[A-Z]+/\1/'))

## Define the adapter sequences to filter for ##

i5=$i51$i5barcode$i52
i7=$i71$i7barcode$i72
i5revcomp=$(echo $i5 | rev | tr ATCG TAGC)
i7revcomp=$(echo $i7 | rev | tr ATCG TAGC)

## Write the adapters.fasta file ##

cat $TruSeq3 > ../cleaned/$sample/adapters.fasta
printf '\n'\>i5'\n'$i5'\n'\>i7'\n'$i7'\n'\>i5revcomp'\n'$i5revcomp'\n'\>i7revcomp'\n'$i7revcomp'\n' >> ../cleaned/$sample/adapters.fasta

### Run trimmomatic to generate the cleaned files ###

java -jar /ddn/apps/Cluster-Apps/miniconda2/4.3.27/jar/trimmomatic.jar PE -threads 8 -phred33 $input/$file $input/$file2 $output/$file $output/$sample\-READ1-single.fastq.gz $output/$file2 $output/$sample\-READ2-single.fastq.gz ILLUMINACLIP:$path/cleaned/$sample/adapters.fasta:2:30:10 LEADING:5 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:40 &> $path/cleaned/$sample/stats/$sample\-adapter-contam.txt

### Concatenate the 2 files of unpaired reads and remove them afterwards ###

cat $output/$sample\-READ1-single.fastq.gz > $output/$sample\-READ-singleton.fastq.gz
cat $output/$sample\-READ2-single.fastq.gz >> $output/$sample\-READ-singleton.fastq.gz
rm $output/$sample\-READ1-single.fastq.gz $output/$sample\-READ2-single.fastq.gz

done

#########################################
## 2) UCE contig assembly with Trinity ##
#########################################

### Array job to run Trinity assemblies ###

#!/bin/bash
#$ -V
#$ -cwd
#SBATCH -J trinity_rest
#SBATCH -p par7.q
#SBATCH -n 120
#SBATCH -N 10 --ntasks-per-node=12
#SBATCH -o /ddn/data/sbvd77/UCE/oe_files/array.log
#SBATCH -a 1-35%10

path=/ddn/data/sbvd77/UCE
FILE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $path/info/list_samples_rest)

module load dbl/phyluce/1.5.0
export PATH=`echo ${PATH} | awk -v RS=: -v ORS=: '/jdk1.8.0/ {next} {print}'`

phyluce_assembly_assemblo_trinity --config $path/info/$(echo ${FILE} | cut -d ":" -f1).conf --output $path/assemblies --subfolder split-adapter-quality-trimmed --cores 12 --clean --log-path $path/info/trinity_log &> $path/oe_files/trinity$SLURM_ARRAY_TASK_ID.log

###################################
## 3) Matching contigs to probes ##
###################################

path=/ddn/data/sbvd77/UCE
phyluce_assembly_match_contigs_to_probes \
    --contigs $path/assemblies/contigs/ \
    --probes $path/probes/uce-5k-probes.fasta \
    --output $path/contigs_aligned_to_probes &> $path/info/match_contigs_to_probes.log

###############################################################################
## 4) Build a monolithic fasta file containing all selected samples and loci ##
###############################################################################

path=/ddn/data/sbvd77/UCE

phyluce_assembly_get_match_counts \
    --locus-db $path/contigs_aligned_to_probes/probe.matches.sqlite \
    --taxon-list-config $path/info/datasets.conf \
    --taxon-group 'dataset_out' \
    --output $path/match_counts/dataset_out.conf \
    --incomplete-matrix

phyluce_assembly_get_fastas_from_match_counts \
    --contigs $path/assemblies/contigs/ \
    --locus-db $path/contigs_aligned_to_probes/probe.matches.sqlite \
    --match-count-output $path/match_counts/dataset_out.conf \
    --incomplete-matrix $path/fastas_from_match_counts/dataset_out.incomplete \
    --output $path/fastas_from_match_counts/dataset_out.fasta

##################################
## 5) Align UCE loci with mafft ##
##################################

phyluce_align_seqcap_align \
    --fasta $path/fastas_from_match_counts/dataset_out.fasta \
    --output $path/alignments/dataset_out/ \
    --output-format fasta \
    --taxa 52 \
    --aligner mafft \
    --no-trim \
    --cores 12 \
    --incomplete-matrix \
    --log-path $path/info/

#####################################
## 6) Trim alignments with GBlocks ##
#####################################

phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed \
    --alignments $path/alignments/dataset_out/ \
    --output $path/alignments/dataset_out_trimmed/ \
    --b4 8 \
    --cores 12 \
    --log-path $path/info/

###################################
## 7) Contig datasets generation ##
###################################

### Script to remove locus names ###

phyluce_align_remove_locus_name_from_nexus_lines \
    --alignments $path/alignments/dataset_out_trimmed \
    --output $path/alignments/dataset_out_trimmed_cleaned \
    --cores 12

### Scripts to build the 75 and 95% contig datasets ###

phyluce_align_get_only_loci_with_min_taxa \
    --alignments $path/alignments/dataset_out_trimmed_cleaned \
    --taxa 52 \
    --percent 0.75 \
    --output $path/alignments/dataset_out_trimmed_cleaned_selected_75 \
    --cores 12 \
    --log-path $path/info

phyluce_align_get_only_loci_with_min_taxa \
    --alignments $path/alignments/dataset_out_trimmed_cleaned \
    --taxa 52 \
    --percent 0.95 \
    --output $path/alignments/dataset_out_trimmed_cleaned_selected_95 \
    --cores 12 \
    --log-path $path/info

### Scripts to remove species with no correspondence in the ddRAD dataset and generate the phylip concatenated alignments ###

TriSeq -in $path/alignments/dataset_out_trimmed_cleaned_selected_75/uce*nexus -of phylip -rm PAAsD1 -o UCE_75
TriSeq -in $path/alignments/dataset_out_trimmed_cleaned_selected_95/uce*nexus -of phylip -rm PAAsD1 -o UCE_95

###########################################
## 8) IUPAC consensus dataset generation ##
###########################################

### Script to create fasta files with all the uce loci per sample ###

phyluce_align_explode_alignments\
    --alignments $path/alignments/dataset_out_trimmed\
    --input-format nexus \
    --output $path/sample_specific_fastas\
    --by-taxon

### Script to build the phasing.conf ###

cd $path

printf \[references\]'\n' > $path/info/phasing.conf
find `pwd` -name *.fasta | grep "sample_specific" | sed -E 's/(.*\/([a-zA-Z0-9]+)\_contigs\.fasta)/\2\:\1/' >> $path/info/phasing.conf
printf '\n'\[individuals\]'\n' >> $path/info/phasing.conf
cd cleaned
for file in $(ls -1); do
  fastq=/ddn/data/sbvd77/UCE/cleaned/$file
  printf $file\:$fastq'\n' >> $path/info/phasing.conf
done
printf '\n'\[flowcell\]'\n' >> $path/info/phasing.conf
for file in $(ls -1); do
  flowcell=$(zcat $file/split*/$file\-READ1.fastq.gz | head -1 | cut -d ":" -f3)
  printf $file\:$flowcell'\n' >> $path/info/phasing.conf
done

### Script to align the raw reads to the aligned and trimmed contigs in an individual basis ###

phyluce_snp_bwa_multiple_align \
--config $path/info/phasing.conf \
--output $path/bwamem_mapping \
--subfolder split-adapter-quality-trimmed \
--mem \
--cores 10 \
--log-path $path/info

### Script to phase the bam files and write fasta files with both haplotypes per loci ###

phyluce_snp_phase_uces \
    --config $path/info/phasing.conf \
    --bams $path/bwamem_mapping \
    --output $path/phased_uces \
    --conservative \
    --cores 10 \
    --log-path $path/info

### Script to merge both alleles to create IUPAC consensus sequences ###

cd /ddn/data/sbvd77/UCE/phased_uces/fastas

cat joined_allele_sequences_all_samples.fasta | sed -n -e '/\_0/,/>/ p' | grep -v "\_1" | sed -E 's/\_0 \|.*//g' > 0_file
cat joined_allele_sequences_all_samples.fasta | sed -n -e '/\_1/,/>/ p' | grep -v "\_0" | sed -E 's/\_1 \|.*//g' > 1_file
seqtk mergefa -m 0_file 1_file > joined_allele_sequences_all_samples_iupac.fasta
rm 0_file 1_file


### Script to align IUPAC sequences of UCEs ###

phyluce_align_seqcap_align \
    --fasta $path/phased_uces/fastas/joined_allele_sequences_all_samples_iupac.fasta \
    --output $path/alignments/dataset_out_phased_iupac_trimmed \
    --taxa 58 \
    --aligner mafft \
    --cores 12 \
    --incomplete-matrix \
    --ambiguous \
    --log-path $path/info

### Script to remove locus names ###

phyluce_align_remove_locus_name_from_nexus_lines \
    --alignments $path/alignments/dataset_out_phased_iupac_trimmed \
    --output $path/alignments/dataset_out_phased_iupac_trimmed_cleaned \
    --cores 12

### Scripts to build the 75 and 95% datasets ###

phyluce_align_get_only_loci_with_min_taxa \
    --alignments $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected \
    --taxa 52 \
    --percent 0.75 \
    --output $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_95 \
    --cores 12 \
    --log-path $path/info

phyluce_align_get_only_loci_with_min_taxa \
    --alignments $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected \
    --taxa 52 \
    --percent 0.95 \
    --output $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_75 \
    --cores 12 \
    --log-path $path/info 

### Scripts to remove species with no correspondence in the ddRAD dataset and filter out positions with more than 25% and 5% missing data ###

TriSeq -in $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_75/uce*nexus -of phylip -rm PAAsD1 --missing-filter 25 25 -o UCE_75_phased_iupac_partitions_75post
TriSeq -in $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_95/uce*nexus -of phylip -rm PAAsD1 --missing-filter 5 5 -o UCE_95_phased_iupac_partitions_95post

### Change names so as they are the same in ddRAD and UCEs ###

UCE=($(cat correspondence_ddRAD-UCE-concat.tsv | grep -v "ddRAD" | cut -f2))
concat=($(cat correspondence_ddRAD-UCE-concat.tsv | grep -v "ddRAD" | cut -f3))
for ((i = 0; i < ${#UCE[@]}; i++)); do sed -iE "s/${UCE[i]}/${concat[i]}/" UCE_75_phased_iupac_partitions_75post.phy; done
for ((i = 0; i < ${#UCE[@]}; i++)); do sed -iE "s/${UCE[i]}/${concat[i]}/" UCE_95_phased_iupac_partitions_95post.phy; done

########################################
## 9) SNPs dataset for SNAPP analyses ##
########################################

### Script to convert fasta files to vcf using snp-sites ###

for file in *.fas; do snp-sites -v -o ../dataset_out_phased_iupac_trimmed_cleaned_selected_75_vcf/$(echo $file | sed 's/.fas/.vcf/') $file; cat ../dataset_out_phased_iupac_trimmed_cleaned_selected_75_vcf/$(echo $file | sed 's/.fas/.vcf/') | sed "s/^1/$(echo $file | sed 's/.fas//')/" | grep -v "*" > ../dataset_out_phased_iupac_trimmed_cleaned_selected_75_vcf/$(echo $file | sed 's/.fas/.vcfbo/');rm ../dataset_out_phased_iupac_trimmed_cleaned_selected_75_vcf/$(echo $file | sed 's/.fas/.vcf/'); mv ../dataset_out_phased_iupac_trimmed_cleaned_selected_75_vcf/$(echo $file | sed 's/.fas/.vcfbo/') ../dataset_out_phased_iupac_trimmed_cleaned_selected_75_vcf/$(echo $file | sed 's/.fas/.vcf/'); done

### Remove non-biallelic SNPs and choose one SNP per locus ###

cd $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_75_vcf

for file in *vcf; do
  cat $file | grep "^#" > $file\_int
  cat $file | grep -v "^#" | grep -v "," >> $file\_int
  cat $file | grep -v "^#" | grep "," | grep -vE ",.*," | awk '/R/ && /A/ && /G/' >> $file\_int
  cat $file | grep -v "^#" | grep "," | grep -vE ",.*," | awk '/Y/ && /C/ && /T/' >> $file\_int
  cat $file | grep -v "^#" | grep "," | grep -vE ",.*," | awk '/M/ && /A/ && /C/' >> $file\_int
  cat $file | grep -v "^#" | grep "," | grep -vE ",.*," | awk '/K/ && /G/ && /T/' >> $file\_int
  cat $file | grep -v "^#" | grep "," | grep -vE ",.*," | awk '/S/ && /C/ && /G/' >> $file\_int
  cat $file | grep -v "^#" | grep "," | grep -vE ",.*," | awk '/W/ && /A/ && /T/' >> $file\_int
  cat $file\_int | grep "^#" > $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_75_vcf_onlybiallelic/$file
  cat $file\_int | grep -v "^#" | shuf >> $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_75_vcf_onlybiallelic/$file
  rm $file\_int
done

cd $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_75_vcf_onlybiallelic

for file in *vcf; do snps=$(cat $file | grep -vc "^#"); if [ $snps -eq 0 ]; then rm $file; fi; done; ls | wc -l

for file in *vcf; do cat $file | grep "^#" > $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_75_vcf_onlybiallelic_1SNPperlocus/$file; cat $file | grep -v "^#" | shuf | head -1 >> $path/dataset_out_phased_iupac_trimmed_cleaned_selected_75_vcf_onlybiallelic_1SNPperlocus/$file; done

### Convert to fasta ###

cd $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_75_vcf_onlybiallelic_1SNPperlocus

for file in *vcf; do
  zero=$(cat $file | grep "uce" | cut -f4)
  alt=$(cat $file | grep "uce" | cut -f5)
  two=$(echo $alt | grep ",")
  if [ -z "$two" ]; then one=$alt; else one=$(echo $alt | cut -d "," -f1); two=$(echo $alt | cut -d "," -f2); fi  
  sp=($(cat $file | grep "CHROM" | cut -f10- | tr '\t' '\n'))
  value=($(cat $file | grep "uce" | cut -f10- | tr '\t' '\n'))
  for ((i=0;i<${#sp[@]};i++)); do 
    if [ ${value[i]}  -eq 0 ]; then printf \>${sp[i]}'\n'$zero'\n'>> $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_75_fasta_onlybiallelic_1SNPperlocus/$(echo $file | sed 's/.vcf/.fa/')
    elif [ ${value[i]}  -eq 1 ]; then printf \>${sp[i]}'\n'$one'\n'>> $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_75_fasta_onlybiallelic_1SNPperlocus/$(echo $file | sed 's/.vcf/.fa/')
    else printf \>${sp[i]}'\n'$two'\n'>> $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_75_fasta_onlybiallelic_1SNPperlocus/$(echo $file | sed 's/.vcf/.fa/')
    fi
  done
done

### Concatenate the fasta files ###

cd $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_75_fasta_onlybiallelic_1SNPperlocus
TriSeq -in *.fa -of fasta -o ../UCE_phased_snapp

### Change x introduced by snp-sites by n ###

sed -i 's/x/n/g' ../UCE_phased_snapp.fas

### Convert fasta files to SNAPP nexus files ###

TriSeq -c -in $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_75_fasta_onlybiallelic_1SNPperlocus/*.fa -of snapp -o $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_75_fasta_onlybiallelic_1SNPperlocus_snapp 

### Concatenate the nexus files and make some modifications ###

python3 ./ElConcatenero.py -if nexus -of nexus -in $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_75_fasta_onlybiallelic_1SNPperlocus_snapp/*nex -o UCE_phased_snapp

cd $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_75_fasta_onlybiallelic_1SNPperlocus_snapp

cat UCE_phased_snapp.nex | sed -E 's/mixed.*\;/integerdata\ symbols\=\"012\"\ missing\=\?\ gap\=\-\;/' | head -6 > UCE_phased_snapp.nexbo

tail -n52 UCE_phased_snapp.nex | sed -E 's/n/?/g' >> UCE_phased_snapp.nexbo

rm UCE_phased_snapp.nex

mv UCE_phased_snapp.nexbo UCE_phased_snapp.nex

## I finally change e?d by end in the last line ##

#################################################################
## 10) Get summary statistics for UCE assembly for each sample ##
#################################################################

### Calculate some sum stats (num of contigs, length of contigs, etc.) ###

printf sample,contigs,total_bp,mean_length,95_CI_length,min_length,max_length,median_length,'contigs_>1kb''\n' > $path/info/get_fasta_lengths.tsv
for i in $path/assemblies/contigs/*.fasta;do
  phyluce_assembly_get_fasta_lengths --input $i --csv >> $path/info/get_fasta_lengths.tsv
done

### Compute coverage across assemblies and summarize those data by contig ###

phyluce_assembly_get_trinity_coverage --assemblies $path/assemblies/ --assemblo-config $path/info/assembly_conf/assembly.conf --cores 12 --subfolder split-adapter-quality-trimmed --clean --bwa-mem &> 
$path/info/get_trinity_coverage.tsv     

### Get a table with summary information about contigs matching probes ###

cat $path/info/match_contigs_to_probes.log | grep "%" | cut -d "-" -f6 | sed -E 's/^(.*)\:\s([0-9]+)\s.*\s([0-9]+)\s.*\s([0-9]+)\s.*\s([0-9]+)\s.*\s([0-9]+)\s.*/\1\t\2\t\3\t\4\t\5\t\6/g' | cut -f 1-3,5,6 >> info/match_contigs_to_probes.tsv 

### Get a file to plot in R ###

printf tTaxon'\t'Total_contigs'\t'Mean_length'\t'Mean_coverage'\t'Ontarget_bp_percent'\n' > UCEs_coverage.tsv
cat phyluce_assembly_get_trinity_coverage_for_uce_loci.log | grep -A5 "Processing" | grep -E 'Processing | contigs' | sed -E 's/.*\ing\s([a-zA-Z]+?[0-9]+).*/\1/' | sed -E 's/.*\ing\s([a-zA-Z]+).*/\1/' | sed -E 's/.*\INFO -//' | sed 's/ contigs, mean trimmed length = /     /' | sed 's/, mean trimmed coverage = /     /' | sed 's/x, on-target bases (uce contigs) = /    /' | sed -E 's/%.*//' | sed 's/     /        /' | tr '\n' '\t' | sed -E 's/\t([a-zA-Z])/\n\1/g' | sed 's/\t \t/\t/' >> UCEs_coverage.tsv

#######################################################
## 11) Plot summary statistics for UCE assembly in R ##
#######################################################

#!/bin/R

library(ggplot2)
library(RColorBrewer)
library(devtools)
library(ggpubr)

setwd("~/Dropbox/Tesi/Data/Durham/UCE/info/")
d = read.delim("UCEs_coverage.tsv")

### Script to make plots of the correlations between different variables ###

concov = ggplot(d, aes(x=Total_contigs, y=Mean_coverage)) +
  geom_point(color="gray30", size = 1.5) +
  labs(x = "Total UCE loci (max. 5060 loci)", y = "Mean sequencing coverage")

## Calculating the regression line and r2 ##

m = lm(Mean_coverage ~ Total_contigs, d)
eq1 = paste0("y = ", format(coef(m)[1], digits = 4), " + ", format(coef(m)[2], digits = 2), " x")
eq2 = paste0("rÂ² = ", format(summary(m)$r.squared, digits = 3))
eq3 = paste0("p-value = ", format(summary(m)$coefficients[2,4], digits = 3))
eq = c(eq1, eq2, eq3)
eq = paste(eq, collapse="\n")
concov = concov +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text", x = 3400, y = 45, label = eq, color="black", size = 2,parse = FALSE)


conlen = ggplot(d, aes(x=Total_contigs, y=Mean_length)) +
  geom_point(color="gray30", size = 1.5) +
  labs(x = "Total UCE loci (max. 5060 loci)", y = "Mean UCE loci length (bp)")

m = lm(Mean_length ~ Total_contigs, d)
eq1 = paste0("y = ", format(coef(m)[1], digits = 4), " + ", format(coef(m)[2], digits = 2), " x")
eq2 = paste0("rÂ² = ", format(summary(m)$r.squared, digits = 3))
eq3 = paste0("p-value = ", format(summary(m)$coefficients[2,4], digits = 3))
eq = c(eq1, eq2, eq3)
eq = paste(eq, collapse="\n")
conlen = conlen +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text", x = 3400, y = 800, label = eq, color="black", size = 2,parse = FALSE)


conon = ggplot(d, aes(x=Total_contigs, y=Ontarget_bp_percent)) +
  geom_point(color="gray30", size = 1.5) +
  labs(x = "Total UCE loci (max. 5060 loci)", y = "On target bp (%)")

m = lm(Ontarget_bp_percent ~ Total_contigs, d)
eq1 = paste0("y = ", format(coef(m)[1], digits = 4), " + ", format(coef(m)[2], digits = 2), " x")
eq2 = paste0("rÂ² = ", format(summary(m)$r.squared, digits = 3))
eq3 = paste0("p-value = ", format(summary(m)$coefficients[2,4], digits = 3))
eq = c(eq1, eq2, eq3)
eq = paste(eq, collapse="\n")
conon = conon +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text", x = 3400, y = 84, label = eq, color="black", size = 2,parse = FALSE)


covlen = ggplot(d, aes(x=Mean_coverage, y=Mean_length)) +
  geom_point(color="gray30", size = 1.5) +
  labs(x = "Mean sequencing coverage", y = "Mean UCE loci length (bp)")

m = lm(Mean_length ~ Mean_coverage, d)
eq1 = paste0("y = ", format(coef(m)[1], digits = 4), " + ", format(coef(m)[2], digits = 2), " x")
eq2 = paste0("rÂ² = ", format(summary(m)$r.squared, digits = 3))
eq3 = paste0("p-value = ", format(summary(m)$coefficients[2,4], digits = 3))
eq = c(eq1, eq2, eq3)
eq = paste(eq, collapse="\n")
covlen = covlen +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text", x = 15, y = 800, label = eq, color="black", size = 2,parse = FALSE)


covon = ggplot(d, aes(x=Mean_coverage, y=Ontarget_bp_percent)) +
  geom_point(color="gray30", size = 1.5) +
  labs(x = "Mean sequencing coverage", y = "On target bp (%)")

m = lm(Ontarget_bp_percent ~ Mean_coverage, d)
eq1 = paste0("y = ", format(coef(m)[1], digits = 4), " + ", format(coef(m)[2], digits = 2), " x")
eq2 = paste0("rÂ² = ", format(summary(m)$r.squared, digits = 3))
eq3 = paste0("p-value = ", format(summary(m)$coefficients[2,4], digits = 3))
eq = c(eq1, eq2, eq3)
eq = paste(eq, collapse="\n")
covon = covon +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text", x = 15, y = 84, label = eq, color="black", size = 2,parse = FALSE)


lenon = ggplot(d, aes(x=Mean_length, y=Ontarget_bp_percent)) +
  geom_point(color="gray30", size = 1.5) +
  labs(x = "Mean UCE loci length (bp)", y = "On target bp (%)")

m = lm(Ontarget_bp_percent ~ Mean_length, d)
eq1 = paste0("y = ", format(coef(m)[1], digits = 4), " + ", format(coef(m)[2], digits = 2), " x")
eq2 = paste0("rÂ² = ", format(summary(m)$r.squared, digits = 3))
eq3 = paste0("p-value = ", format(summary(m)$coefficients[2,4], digits = 3))
eq = c(eq1, eq2, eq3)
eq = paste(eq, collapse="\n")
lenon = lenon +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text", x = 470, y = 84, label = eq, color="black", size = 2,parse = FALSE)


ggarrange(concov, conlen, conon, covlen, covon, lenon, ncol = 3, nrow = 2)
