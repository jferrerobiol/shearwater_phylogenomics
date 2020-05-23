#######################################
## 1) Trimmomatic to remove adaptors ##
#######################################
```
cd $HOME/ddRAD/rawdata/raw_to_filter ## move to folder cotaining the raw files to filter

for file in *R1_001.fastq.gz ## start a for loop to treat all the raw files in the folder

do

file2="$(echo $file | sed 's/_R1_/_R2_/')"

java -jar /soft/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 24 -phred33 $HOME/ddRAD/rawdata/raw_to_filter/$file $HOME/ddRAD/
rawdata/raw_to_filter/$file2 $HOME/ddRAD/rawdata/raw_filtered_trimmomatic/$file"_paired" $HOME/ddRAD/rawdata/raw_filtered_trimmoma
tic/$file"_unpaired" $HOME/ddRAD/rawdata/raw_filtered_trimmomatic/$file2"_paired" $HOME/ddRAD/rawdata/raw_filtered_trimmomatic/$fi
le2"_unpaired" ILLUMINACLIP:/soft/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:30:10 MINLEN:150

done
```
```
rename 's/.fastq.gz_paired/_paired.fastq.gz/' *_paired
rename 's/.fastq.gz_unpaired/_unpaired.fastq.gz/' *_unpaired
```
##########################################
## 2) Remove contaminants using bowtie2 ##
##########################################
```
cd $HOME/ddRAD/rawdata/raw_filtered_trimmomatic ## move to folder cotaining the raw files to filter

for file in *_paired.fastq ## start a for loop to treat all the raw files in the folder

do

/soft/bowtie2-223/bowtie2 -x $HOME/ddRAD/bowtie2/index/contaminant_db -U $HOME/ddRAD/rawdata/raw_filtered_trimmomatic/$file --un-g
z $HOME/ddRAD/rawdata/raw_filtered_cont/$file"_nocont" --very-sensitive-local -p 24

done
```
###################################
## 3) Syncronize R1 and R2 files ##
###################################

### Sort R1 and R2 files ###
```
for file in *gz
do
zcat $file | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > ../raw_filtered_sorted/$file"_sorted.fastq"
```
### Syncronize read files after contaminant removal ###
```
prefixref=$HOME/ddRAD/rawdata/raw_to_filter/
prefixout=$HOME/ddRAD/rawdata/raw_filtered/sync/

cd $HOME/ddRAD/rawdata/raw_filtered_sorted ## move to folder containing the files to syncronize

for file in *_R1_*.fastq ## start a loop that will only consider R1 files

do

file2=$(echo $file | sed 's/R1/R2/')
fileori=$prefixref$(echo $file | sed 's/_paired_nocont//')
fileout1=$prefixout$(echo $file | sed 's/_sorted/_sorted_sync/')
fileout2=$prefixout$(echo $file | sed 's/_sorted/_sorted_sync/')

python $HOME/ddRAD/scripts/python/sync_paired_end_reads.py $fileori $file $file2 $fileout1 $fileout2

done
```
############################################
## 4) Demultiplexing with process_radtags ##
############################################
```
cd $HOME/Durham/ddRAD/cleaned

echo "Cleaning and demultiplexing the reads..."

list=$(ls -1 ../raw)

for pool in $list
do
process_radtags -p ../raw/$pool -P -o ./ -b ../info/barcodes$pool.tsv --inline_null --renz_1 ecoRI --renz_2 mspI -s 20 -r -c -q -i gzfastq &> process_radtags.$pool.oe ## The files in raw/ correspond to the syncronized files from step 3)
done
```
#####################################
## 5) Check process_radtags output ##
#####################################

### Check the per-sample coverages: extract the number of reads from the log of process_radtags ###
```
cd $HOME/Durham/ddRAD/info

echo -e '#sample\tn_reads' > n_reads_per_sample.tsv

list=$(ls -1 ../raw)

for pool in $list
do

# Retrieve the part of the log between 'Barcode...' and the next empty line,
# then discard the first and last lines and keep the 2nd and 6th columns.

sed -n '/^Barcode\tFilename\t/,/^$/ p' ../cleaned/process_radtags.$pool.log \
                | sed '1 d; $ d' \
                | cut -f2,6 \
                >> n_reads_per_sample.tsv
done
```
### Plot these numbers using 2.plot_n_reads_per_sample.R from Rochette & Catchen 2017
#### Remove bad samples (if any) or samples we will not use for the analyses from the population map
```
echo -n "\
aPMyr
aPBry01
CEdw2
" > discarded_samples

discarded_regex=$(paste -s -d'|' discarded_samples)
mv popmap.tsv complete_popmap.tsv
grep -vE "^($discarded_regex)\>" complete_popmap.tsv > popmap.tsv
```
############################################
## 6) Parameter optimization using STACKS ##
############################################

### STACKS M parameter optimization ###
####  from M=2 to M=10 with n=M+2 ##
```
denovo_map.pl -S -b 1 -m 3 -M $M -n $((M+2)) -T 10 -o $HOME/ddRAD/stacks/denovomap/stacksugM$M -O $HOME/ddRAD/stacks/popmap/popmapoutfilo --samples $HOME/Durham/ddRAD/cleaned
```
### STACKS n parameter optimization ###
#### from n=4 to n=20 with M=4 ##
```
denovo_map.pl -S -b 1 -m 3 -M 4 -n $n -T 10 -o $HOME/ddRAD/stacks/denovomap/stacksugn$n -O $HOME/ddRAD/stacks/popmap/popmapoutfilo --samples $HOME/Durham/ddRAD/cleaned
```
### STACKS m parameter optimization ###
####  from m=1 to m=10 with M=4 and n=6 ##
```
denovo_map.pl -S -b 1 -m $m -M 4 -n 6 -T 10 -o $HOME/ddRAD/stacks/denovomap/stacksugm$m -O $HOME/ddRAD/stacks/popmap/popmapoutfilo --samples $HOME/Durham/ddRAD/cleaned
```
### STACKS M parameter plots ###
 ```
#!/bin/R

## New loci plots as in Paris et al. 2017 ##

library(plotrix) # Library to make gaps at plot axes

# Plots for M iterations keeping n at M=2 (makes no sense in a phylogenetic dataset to set values of n smaller or equal to M)
setwd("/home/joan/Dropbox/Tesi/Data/ddRAD/stacks_parameter_tests")

data <- read.table("outfilo_M_parameter_tests_R.csv", header=T, sep = "\t")
data$kassembloci <- data$Num_assemb_loci/1000 # Number of assembled loci as K (the same for polymorphic loci and SNPs in the following rows)
data$kpolimloci <- data$Num_polim_loci/1000
data$kSNPs <- data$Num_SNPs/1000
data$SNPs_p_polimloci <- data$Num_SNPs/data$Num_polim_loci

d40 <- data[data$Dataset=="p40pc",] # Building the datasets to calculate the number of new assembled loci, polymorphic and SNPs for every increment of M
d40$newassemloci <- NA
d40$newpolloci <- NA
d40$newSNPs <- NA

for (i in 1:length(d40$Num_polim_loci)){
  d40$newassemloci[i] <- d40$Num_assemb_loci[i+1]-d40$Num_assemb_loci[i]
  d40$newpolloci[i] <- d40$Num_polim_loci[i+1]-d40$Num_polim_loci[i]
  d40$newSNPs[i] <- d40$Num_SNPs[i+1]-d40$Num_SNPs[i]
}
 
d40 <- d40[!is.na(d40$newSNPs),]

d60 <- data[data$Dataset=="p60pc",] # Building the datasets to calculate the number of new assembled loci, polymorphic and SNPs for every increment of M
d60$newassemloci <- NA
d60$newpolloci <- NA
d60$newSNPs <- NA

for (i in 1:length(d60$Num_polim_loci)){
  d60$newassemloci[i] <- d60$Num_assemb_loci[i+1]-d60$Num_assemb_loci[i]
  d60$newpolloci[i] <- d60$Num_polim_loci[i+1]-d60$Num_polim_loci[i]
  d60$newSNPs[i] <- d60$Num_SNPs[i+1]-d60$Num_SNPs[i]
}

d60 <- d60[!is.na(d60$newSNPs),]

d80 <- data[data$Dataset=="p80pc",] # Building the datasets to calculate the number of new assembled loci, polymorphic and SNPs for every increment of M
d80$newassemloci <- NA
d80$newpolloci <- NA
d80$newSNPs <- NA

for (i in 1:length(d80$Num_polim_loci)){
  d80$newassemloci[i] <- d80$Num_assemb_loci[i+1]-d80$Num_assemb_loci[i]
  d80$newpolloci[i] <- d80$Num_polim_loci[i+1]-d80$Num_polim_loci[i]
  d80$newSNPs[i] <- d80$Num_SNPs[i+1]-d80$Num_SNPs[i]
}

d80 <- d80[!is.na(d80$newSNPs),]

pdf(file = "M_iterations_new_loci.pdf", width=9, height=21, pointsize=20) # Plot as in Paris et al. 2017
par(mfrow=c(3,1))
plot(d40$newassemloci, type="b", col="blue", lwd=3, xaxt="n", yaxt="n", xlab="", ylab="Number of new assembled loci", ylim=c(-150,850))
axis(side=1, at=c(1:8), labels=c("M2-M3","M3-M4","M4-M5","M5-M6","M6-M7","M7-M8","M8-M9","M9-M10"))
axis(2, at=seq(-150, 850, by=50))
lines(d60$newassemloci, type="b", col="red", lwd=3)
lines(d80$newassemloci, type="b", col="green", lwd=3)
legend(x=5, y=850, y.intersp=1, bty="n", c("40% of individuals", "60% of individuals", "80% of individuals"), lty=c(1,1,1), lwd=c(3,3,3), col=c("blue","red","green"))

plot(d40$newpolloci, type="b", col="blue", lwd=3, xaxt="n", yaxt="n", xlab="", ylab="Number of new polymorphic loci", ylim=c(-150,850))
axis(side=1, at=c(1:8), labels=c("M2-M3","M3-M4","M4-M5","M5-M6","M6-M7","M7-M8","M8-M9","M9-M10"))
axis(2, at=seq(-150, 850, by=50))
lines(d60$newpolloci, type="b", col="red", lwd=3)
lines(d80$newpolloci, type="b", col="green", lwd=3)

plot(d40$newSNPs, type="b", col="blue", lwd=3, xaxt="n", yaxt="n", xlab="", ylab="Number of new variable positions", ylim=c(500,25000))
axis(side=1, at=c(1:8), labels=c("M2-M3","M3-M4","M4-M5","M5-M6","M6-M7","M7-M8","M8-M9","M9-M10"))
axis(2, at=seq(500, 25000, by=500))
lines(d60$newSNPs, type="b", col="red", lwd=3)
lines(d80$newSNPs, type="b", col="green", lwd=3)
dev.off()

## SNP plots as in Rochette & Catchen 2017 ##

# Read the M r80 file
Mr80 = read.delim('./n_snps_per_locus_M=n-2_r80.tsv')

Mr80 = subset(Mr80, M==n-2 & m==3)

# Make sure the table is ordered by number of snps.
Mr80 = Mr80[order(Mr80$n_snps),]

Mn_values80 = sort(unique(Mr80$M))

# Write the counts in a matrix.
m80 = matrix(NA, nrow=length(Mn_values80), ncol=max(Mr80$n_snps)+1)
for(i in 1:nrow(Mr80)) {
  m80[Mr80$M[i]-1,Mr80$n_snps[i]+1] = Mr80$n_loci[i] # [M-1] as row 1 is for M2, [n_snps+1] as column 1 is for loci with 0 SNPs
}

m80 = m80 / rowSums(m80, na.rm=T)

# Read the M r40 file
Mr40 = read.delim('./n_snps_per_locus_M=n-2_r40.tsv')

Mr40 = subset(Mr40, M==n-2 & m==3)

# Make sure the table is ordered by number of snps.
Mr40 = Mr40[order(Mr40$n_snps),]

Mn_values40 = sort(unique(Mr40$M))

# Write the counts in a matrix.
m40 = matrix(NA, nrow=length(Mn_values40), ncol=max(Mr40$n_snps)+1)
for(i in 1:nrow(Mr40)) {
  m40[Mr40$M[i]-1,Mr40$n_snps[i]+1] = Mr40$n_loci[i] # [M-1] as row 1 is for M2, [n_snps+1] as column 1 is for loci with 0 SNPs
}

max_n_snps = 30
m40[,max_n_snps+2] = rowSums(m40[,(max_n_snps+2):ncol(m40)], na.rm=T)
m40 = m40[,1:(max_n_snps+2)]
m40 = m40 / rowSums(m40, na.rm=T)

# Draw the barplots for M iterations
pdf('n_snps_per_locus_M=n-2.pdf', height = 15, width = 15)
par(mfrow=c(2,1))
col = rev(heat.colors(length(Mn_values80)))

barplot(m80,
        beside=T, col=col, las=1,
        names.arg=c(0:30, '>30'),
        ylim=c(0,0.08),
        xlab='Number of SNPs',
        ylab='Percentage of loci',
        main='Distributions of the number of SNPs per locus\nfor a range of M==n-2 values\nr80'
)
legend('topright', title='M==n-2', legend=Mn_values80, fill=col)

barplot(m40,
        beside=T, col=col, las=1,
        names.arg=c(0:30, '>30'),
        ylim=c(0,0.08),
        xlab='Number of SNPs',
        ylab='Percentage of loci',
        main='r40'
)
legend('topright', title='M==n-2', legend=Mn_values40, fill=col)

null=dev.off()
```

### STACKS n parameter plots ###

#### New loci plots as in Paris et al. 2017 ##
```
# Plot for n iterations at M=4
setwd("/home/joan/Dropbox/Tesi/Data/ddRAD/stacks_parameter_tests")

data <- read.table("outfilo20_n_parameter_tests_R.csv", header=T, sep = "\t")
data$kassembloci <- data$Num_assemb_loci/1000 # Number of assembled loci as K (the same for polymorphic loci and SNPs in the following rows)
data$kpolimloci <- data$Num_polim_loci/1000
data$kSNPs <- data$Num_SNPs/1000
data$SNPs_p_polimloci <- data$Num_SNPs/data$Num_polim_loci

d40 <- data[data$Dataset=="p40pc",] # Building the datasets to calculate the number of new assembled loci, polymorphic and SNPs for every increment of M
d40$newassemloci <- NA
d40$newpolloci <- NA
d40$newSNPs <- NA

for (i in 1:length(d40$Num_polim_loci)){
  d40$newassemloci[i] <- d40$Num_assemb_loci[i+1]-d40$Num_assemb_loci[i]
  d40$newpolloci[i] <- d40$Num_polim_loci[i+1]-d40$Num_polim_loci[i]
  d40$newSNPs[i] <- d40$Num_SNPs[i+1]-d40$Num_SNPs[i]
}

d40 <- d40[!is.na(d40$newSNPs),]

d60 <- data[data$Dataset=="p60pc",] # Building the datasets to calculate the number of new assembled loci, polymorphic and SNPs for every increment of M
d60$newassemloci <- NA
d60$newpolloci <- NA
d60$newSNPs <- NA

for (i in 1:length(d60$Num_polim_loci)){
  d60$newassemloci[i] <- d60$Num_assemb_loci[i+1]-d60$Num_assemb_loci[i]
  d60$newpolloci[i] <- d60$Num_polim_loci[i+1]-d60$Num_polim_loci[i]
  d60$newSNPs[i] <- d60$Num_SNPs[i+1]-d60$Num_SNPs[i]
}

d60 <- d60[!is.na(d60$newSNPs),]

d80 <- data[data$Dataset=="p80pc",] # Building the datasets to calculate the number of new assembled loci, polymorphic and SNPs for every increment of M
d80$newassemloci <- NA
d80$newpolloci <- NA
d80$newSNPs <- NA

for (i in 1:length(d80$Num_polim_loci)){
  d80$newassemloci[i] <- d80$Num_assemb_loci[i+1]-d80$Num_assemb_loci[i]
  d80$newpolloci[i] <- d80$Num_polim_loci[i+1]-d80$Num_polim_loci[i]
  d80$newSNPs[i] <- d80$Num_SNPs[i+1]-d80$Num_SNPs[i]
}

d80 <- d80[!is.na(d80$newSNPs),]

pdf(file = "n_iterations_new_loci.pdf", width=9, height=21, pointsize=20) # Plot as in Paris et al. 2017
par(mfrow=c(3,1))
plot(d40$newassemloci, type="b", col="blue", lwd=3, xaxt="n", yaxt="n", xlab="", ylab="Number of new assembled loci", ylim=c(-150,850))
axis(side=1, at=c(1:16), labels=c("n4-n5","n5-n6","n6-n7","n7-n8","n8-n9","n9-n10","n10-n11","n11-n12","n12-n13","n13-n14","n14-n15","n15-n16","n16-n17","n17-n18","n18-n19","n19-n20"))
axis(2, at=seq(-150, 850, by=50))
lines(d60$newassemloci, type="b", col="red", lwd=3)
lines(d80$newassemloci, type="b", col="green", lwd=3)
legend(x=5, y=850, y.intersp=1, bty="n", c("40% of individuals", "60% of individuals", "80% of individuals"), lty=c(1,1,1), lwd=c(3,3,3), col=c("blue","red","green"))

plot(d40$newpolloci, type="b", col="blue", lwd=3, xaxt="n", yaxt="n", xlab="", ylab="Number of new polymorphic loci", ylim=c(-150,850))
axis(side=1, at=c(1:16), labels=c("n4-n5","n5-n6","n6-n7","n7-n8","n8-n9","n9-n10","n10-n11","n11-n12","n12-n13","n13-n14","n14-n15","n15-n16","n16-n17","n17-n18","n18-n19","n19-n20"))
axis(2, at=seq(-150, 850, by=50))
lines(d60$newpolloci, type="b", col="red", lwd=3)
lines(d80$newpolloci, type="b", col="green", lwd=3)

plot(d40$newSNPs, type="b", col="blue", lwd=3, xaxt="n", yaxt="n", xlab="", ylab="Number of new variable positions", ylim=c(500,25000))
axis(side=1, at=c(1:16), labels=c("n4-n5","n5-n6","n6-n7","n7-n8","n8-n9","n9-n10","n10-n11","n11-n12","n12-n13","n13-n14","n14-n15","n15-n16","n16-n17","n17-n18","n18-n19","n19-n20"))
axis(2, at=seq(500, 25000, by=500))
lines(d60$newSNPs, type="b", col="red", lwd=3)
lines(d80$newSNPs, type="b", col="green", lwd=3)
dev.off()

pdf(file = "n_iterations_4-20.pdf", width=28, height=21, pointsize=20) # Plot as in Paris et al. 2017
par(mfrow=c(3,3))

plot(data[data$Dataset=="p40pc" & data$M==4,]$n, data[data$Dataset=="p40pc" & data$M==4,]$kassembloci, type="p", pch=19, col="blue", xlab="", ylab="Number of assembled loci (k)", axes=FALSE, main="40% of individuals")
box()
axis(side=1, at=c(4:20), labels=c("n4","n5","n6","n7","n8","n9","n10","n11","n12","n13","n14","n15","n16","n17","n18","n19","n20"))
axis(side=2)
plot(data[data$Dataset=="p60pc" & data$M==4,]$n, data[data$Dataset=="p60pc" & data$M==4,]$kassembloci, type="p", pch=19, col="red", xlab="", ylab="", axes=FALSE, main="60% of individuals")
box()
axis(side=1, at=c(4:20), labels=c("n4","n5","n6","n7","n8","n9","n10","n11","n12","n13","n14","n15","n16","n17","n18","n19","n20"))
axis(side=2)
plot(data[data$Dataset=="p80pc" & data$M==4,]$n, data[data$Dataset=="p80pc" & data$M==4,]$kassembloci, type="p", pch=19, col="green",xlab="", ylab="", axes=FALSE, main="80% of individuals")
box()
axis(side=1, at=c(4:20), labels=c("n4","n5","n6","n7","n8","n9","n10","n11","n12","n13","n14","n15","n16","n17","n18","n19","n20"))
axis(side=2)

plot(data[data$Dataset=="p40pc" & data$M==4,]$n, data[data$Dataset=="p40pc" & data$M==4,]$kpolimloci, type="p", pch=19, col="blue", xlab="", ylab="Number of polymorphic loci (k)", axes=FALSE)
box()
axis(side=1, at=c(4:20), labels=c("n4","n5","n6","n7","n8","n9","n10","n11","n12","n13","n14","n15","n16","n17","n18","n19","n20"))
axis(side=2)
plot(data[data$Dataset=="p60pc" & data$M==4,]$n, data[data$Dataset=="p60pc" & data$M==4,]$kpolimloci, type="p", pch=19, col="red", xlab="", ylab="", axes=FALSE)
box()
axis(side=1, at=c(4:20), labels=c("n4","n5","n6","n7","n8","n9","n10","n11","n12","n13","n14","n15","n16","n17","n18","n19","n20"))
axis(side=2)
plot(data[data$Dataset=="p80pc" & data$M==4,]$n, data[data$Dataset=="p80pc" & data$M==4,]$kpolimloci, type="p", pch=19, col="green",xlab="", ylab="", axes=FALSE)
box()
axis(side=1, at=c(4:20), labels=c("n4","n5","n6","n7","n8","n9","n10","n11","n12","n13","n14","n15","n16","n17","n18","n19","n20"))
axis(side=2)

plot(data[data$Dataset=="p40pc" & data$M==4,]$n, data[data$Dataset=="p40pc" & data$M==4,]$kSNPs, type="p", pch=19, col="blue", xlab="", ylab="Number of variable positions (k)", axes=FALSE)
box()
axis(side=1, at=c(4:20), labels=c("n4","n5","n6","n7","n8","n9","n10","n11","n12","n13","n14","n15","n16","n17","n18","n19","n20"))
axis(side=2)
plot(data[data$Dataset=="p60pc" & data$M==4,]$n, data[data$Dataset=="p60pc" & data$M==4,]$kSNPs, type="p", pch=19, col="red", xlab="", ylab="", axes=FALSE)
box()
axis(side=1, at=c(4:20), labels=c("n4","n5","n6","n7","n8","n9","n10","n11","n12","n13","n14","n15","n16","n17","n18","n19","n20"))
axis(side=2)
plot(data[data$Dataset=="p80pc" & data$M==4,]$n, data[data$Dataset=="p80pc" & data$M==4,]$kSNPs, type="p", pch=19, col="green",xlab="", ylab="", axes=FALSE)
box()
axis(side=1, at=c(4:20), labels=c("n4","n5","n6","n7","n8","n9","n10","n11","n12","n13","n14","n15","n16","n17","n18","n19","n20"))
axis(side=2)
dev.off()
```
#### SNP plots as in Rochette & Catchen 2017 ##
```
# Read the n r80 file
nr80 = read.delim('./n_snps_per_locus_M4_n4-20_r80.tsv')

nr80 = subset(nr80, M==4 & m==3)

# Make sure the table is ordered by number of snps.
nr80 = nr80[order(nr80$n_snps),]

n_values80 = sort(unique(nr80$n))

# Write the counts in a matrix.
m80 = matrix(NA, nrow=length(n_values80), ncol=max(nr80$n_snps)+1)
for(i in 1:nrow(nr80)) {
  m80[nr80$n[i]-3,nr80$n_snps[i]+1] = nr80$n_loci[i] # [n-3] as row 1 is for n4, [n_snps+1] as column 1 is for loci with 0 SNPs
}

max_n_snps = 30
m80[,max_n_snps+2] = rowSums(m80[,(max_n_snps+2):ncol(m80)], na.rm=T)
m80 = m80[,1:(max_n_snps+2)]
m80 = m80 / rowSums(m80, na.rm=T)

# Read the n r40 file
nr40 = read.delim('./n_snps_per_locus_M4_n4-20_r40.tsv')

nr40 = subset(nr40, M==4 & m==3)

# Make sure the table is ordered by number of snps.
nr40 = nr40[order(nr40$n_snps),]

n_values40 = sort(unique(nr40$n))

# Write the counts in a matrix.
m40 = matrix(NA, nrow=length(n_values40), ncol=max(nr40$n_snps)+1)
for(i in 1:nrow(nr40)) {
  m40[nr40$n[i]-3,nr40$n_snps[i]+1] = nr40$n_loci[i] # [n-3] as row 1 is for n4, [n_snps+1] as column 1 is for loci with 0 SNPs
}

max_n_snps = 30
m40[,max_n_snps+2] = rowSums(m40[,(max_n_snps+2):ncol(m40)], na.rm=T)
m40 = m40[,1:(max_n_snps+2)]
m40 = m40 / rowSums(m40, na.rm=T)

pdf('n_snps_per_locus_M_4_n4-20.pdf', height = 15, width = 15)
par(mfrow=c(2,1))

col = rev(heat.colors(length(n_values80)))

barplot(m80,
        beside=T, col=col, las=1,
        names.arg=c(0:30, '>30'),
        ylim=c(0,0.08),
        xlab='Number of SNPs',
        ylab='Percentage of loci',
        main='Distributions of the number of SNPs per locus\nfor a range of n values and M=4\nr80'
)
legend('topright', title='n', legend=n_values80, fill=col)

barplot(m40,
        beside=T, col=col, las=1,
        names.arg=c(0:30, '>30'),
        ylim=c(0,0.08),
        xlab='Number of SNPs',
        ylab='Percentage of loci',
        main='r40'
)
legend('topright', title='n', legend=n_values40, fill=col)

null=dev.off()
```
### STACKS m parameter plots ###

### Script to count the number and proportion of singletons in the data sets with m values from 3 to 10 ##
```
cd $HOME/Durham/ddRAD/info
printf '#m_value\tn_singletons\tn_snps\tp_singletons\n' > singleton_count.tsv # Write header on file containing number and proportion of singletons per m value

cd $HOME/ddRAD/stacks/populations/test_outfilo_m/ # go to the directory where the directories we are going to work with are

for dire in p*40 # start a loop to enter the directories for all m values
  do
  m=$(echo $dire | sed -E 's/poutfilo_phred20_(m[0-9]+)_40/\1/') # define m value to be printed afterwards

  cd $dire # enter the directory for any m value

  vcftools --vcf batch_1.vcf --min-alleles 2 --max-alleles 2 --counts --out batch_1 # count the number of alleles per biallelic SNP

  singcount=$(cat batch_1.frq.count | grep -v "CHROM" | cut -f6 | cut -d ":" -f2 | sort -n | uniq -c | head -1 | sed -E 's/^ +//' | cut -d ' ' -f1) # number of singletons

  totcount=$(cat batch_1.frq.count | grep -vc "CHROM") # total number of SNPs

  psing=$(echo "$singcount/$totcount" | bc -l) # proportion of singletons

  printf $m'\t'$singcount'\t'$totcount'\t'$psing'\n' >> $HOME/Durham/ddRAD/info/singleton_count.tsv

  cd ..

done
```

### Plots of singleton number and proportion ##
```
#!/bin/R

d = read.delim('./singleton_count.tsv')

pdf(file = "singleton_per_m_value.pdf", width=15, height=7.5)
par(mfrow=c(1,2))

plot(d$n_singletons, type="b", lwd=3, xaxt="n", xlab="", ylab="Number of singletons")
axis(side=1, at=c(1:8), labels = c("m3","m4","m5","m6","m7","m8","m9","m10"))

plot(d$p_singletons, type="b", lwd=3, xaxt="n", ylim = c(0.323,0.328), xlab="", ylab="Proportion of singletons")
axis(side=1, at=c(1:8), labels = c("m3","m4","m5","m6","m7","m8","m9","m10"))

dev.off()
```
#################
## 7) Chosen Parameterizations ##
#####################

#### a) STACKS optimal:  m3 M5 n8: results of optimizing each of the parameters as in Paris et al. 2017 and Rochette & Catchen 2017 #
#### b) STACKS default:  m3 M2 n1: Stacks default parameterization optimized for population genomics analyses                       #
#### c) STACKS higher n: m3 M5 n15: parameterization to see the effects of adding extremely variable loci but also more paralogs    #
#### d) STACKS higher m: m7 M5 n8: parameterization to avoid including sequencing errors as real alleles                            #
#### e) PyRAD clust 89: similar to STACKS higher n                                                                                  #
#### f) PyRAD clust 94: similar to STACKS optimal                                                                                   #
#### g) STACKS refmap: integrate alignments with the higher m parameterization loci mapped to Calonectris borealis genome           #

###############################           
## 8) STACKS denovo analyses ##
###############################

#### In silico extraction of outgroups radtags (Oceanites oceanicus, Thalassarche chlororhynchos and Fulmarus glacialis ###

#### Extract the radtags with Digital_RADs.py from DaCosta & Sorenson 2014 (example Oceanites oceanicus) ##
```
python3 Digital_RADs.py Oceanites_oceanicus.genomic.fa Oceoce_radtags_150-300 2 GAATTC CCGG 150 300
```
#### Modify the output file so as it can be used in STACKS ##
```
cat Oceoce_radtags_150-300 | grep -v "contig" | cut -f 1-5 | sed -E 's/^/>/' | sed -E 's/\t/_/g' | sed -E 's/(.+_.+_.+)_1_/\1_1_R1\n/' | sed -E 's/(.+_.+_.+)_-1_/\1_-1_R1\n/' | sed 's/^GAATTC/AATTC/g'| cut -c 1-145 > Oceoce_radtags_150-300_R1.fa

cat Oceoce_radtags_150-300 | grep -v "contig" | cut -f 1-5 | sed -E 's/^/>/' | sed -E 's/\t/_/g' | sed -E 's/(.+_.+_.+)_1_/\1_1_R2\n/' | sed -E 's/(.+_.+_.+)_-1_/\1_-1_R2\n/' | while read line; do first=$(echo $line | cut -c 1); if [ $(echo $first) != '>' ]; then echo $line | rev | cut -c 2-150 | tr ATCG TAGC; else echo $line; fi; done > Oceoce_radtags_150-300_R2.fa # because R2 must have 149 bp

cat Oceoce_radtags_150-300_R1.fa | gzip > OOce.1.fa.gz

cat Oceoce_radtags_150-300_R2.fa | gzip > OOce.2.fa.gz
```
#### Subsample reads from samples with coverage >60x so as they have coverage of ~50x ###

#### I select 1630478 reads which is the mean of the 5 samples with coverages around 50x
```
cd $HOME/Durham/ddRAD/cleaned_sel/

list="AGri1 CDio1 PLBo2 PPuf2 PLBa1 PMau2 CLeu1 PMau1 PLBo1 PLBa2 CLeu2"
R2path=$HOME/Durham/ddRAD/cleaned_sel_R2/

for element in $list; do
  zcat $element".1.fq.gz" | grep -E "^@[0-9]+\_" | shuf -n 1630478 | sed 's/^@//' > $element"_whitelist.lst"
  cat $element"_whitelist.lst" | sed -E 's/1$/2/' > $element"_whitelist2.lst"
  $HOME/programari/seqtk/seqtk subseq $element".1.fq.gz" $element"_whitelist.lst" | gzip > $element"bo.1.fq.gz"
  $HOME/programari/seqtk/seqtk subseq $R2path$element".2.fq.gz" $element"_whitelist2.lst" | gzip > $R2path$element"bo.2.fq.gz"
  rm $element".1.fq.gz"
  mv $element"bo.1.fq.gz" $element".1.fq.gz"
  rm $R2path$element".2.fq.gz"
  mv $R2path$element"bo.2.fq.gz" $R2path$element".2.fq.gz"
  echo $element" R1 has "$(echo $(zcat $element".1.fq.gz" | grep -c "^@7_"))" reads"
  echo $element" R2 has "$(echo $(zcat $R2path$element".2.fq.gz" | grep -c "^@7_"))" reads"
done
```
### Script to trim reads all to the same length of 145 and 149 bp for R1 and R2 respectively ###
```
mkdir $HOME/Durham/ddRAD/cleaned_sel_bo $HOME/Durham/ddRAD/cleaned_sel_R2_bo
cd $HOME/Durham/ddRAD/cleaned_sel/

for file in *fq.gz; do
  prefix=$(echo $file | cut -d "." -f1)
  if [ $(echo $(zcat $file | head -2 | grep -v "^@" | awk '{print length}')) == '145' ]
    then cp $file ../cleaned_sel_bo
         $HOME/programari/seqtk/seqtk trimfq -e 1 ../cleaned_sel_R2/$prefix.2.fq.gz | gzip > ../cleaned_sel_R2_bo/$prefix.2.fq.gz
    else $HOME/programari/seqtk/seqtk trimfq -e 1 $file | gzip > ../cleaned_sel_bo/$file
         $HOME/programari/seqtk/seqtk trimfq -e 2 ../cleaned_sel_R2/$prefix.2.fq.gz | gzip > ../cleaned_sel_R2_bo/$prefix.2.fq.gz
  fi
done
```
#### After checking that it worked, remove original directories and rename the directories with the trimmed files ##
```
rm -r $HOME/Durham/ddRAD/cleaned_sel/ $HOME/Durham/ddRAD/cleaned_sel_R2/
mv $HOME/Durham/ddRAD/cleaned_sel_bo/ $HOME/Durham/ddRAD/cleaned_sel/
mv $HOME/Durham/ddRAD/cleaned_sel_R2_bo/ $HOME/Durham/ddRAD/cleaned_sel_R2/
```
### Run the STACKS pipeline for each parameterization (example optimal) ###

### Stacks pipeline for PE ddRAD data using outgroups extracted in silico ## 
### Outputs: concatenated phylip files for 3 different completeness values: 65 (p18_in), 75 (p21_in) and 95% (p26_in) ##

### Setting the data set and cd to it
```
dataset=optimalPE
cd $HOME/Durham/ddRAD/stacks.denovo
mkdir $dataset
cd $dataset
```
### Setting path variables
```
samples=$HOME/Durham/ddRAD/cleaned_sel/
popmap=$HOME/Durham/ddRAD/info/popmap_woutgroups_tot.tsv
popmap2=$HOME/Durham/ddRAD/info/popmap_woutgroups.tsv
```
### Setting parameter variables
```
m=3
M=5
n=8
```
### Ustacks: outgroups: m 1 as it comes from an insilico extraction, M as the rest.
```
outgroups=(OOce FGra TChl OTet FGla) #List of the outgroups
id=(1 2 3 4 5) #List of the id to give to each outgroup
for ((i = 0; i < 5; i++)); do #Start for loop to treat each outgroup with ustacks. It is better to start with 0 if we use the length of the array as the stop because we do not have to do any operations inside. This kind of loop has de advantage that we can make several variables (array type) to iterate inside the loop
  ustacks -f $samples${outgroups[i]}.1.fa.gz -i ${id[i]} --name ${outgroups[i]} -o ./ -m 1 -M $M -p 40 &> ustacks_${outgroups[i]}.oe #Using path variables defined before and array variables with suffixes that refer to the elements of the predefined lists
done
```
### Ustacks: samples with high coverage.
##### Some of these samples have already been shuffled to reduce the number of reads. Disable secondary reads for these samples as they have already enough reads and like this we make them more comparable with lower coverage samples: -N 0
```
highcov=(CLeu1 CLeu2 CDio1 CBor1 AGri1 ATen1 ATen2 PPuf2 PYel2 PMau1 PMau2 PLBa1 PLBa2 PLBa3 PLBo1 PLBo2)
id=(6 7 8 10 18 20 21 30 33 34 35 54 55 56 57 58)
for ((i = 0; i < 16; i++)); do
  ustacks -f $samples${highcov[i]}.1.fq.gz -i ${id[i]} --name ${highcov[i]} -o ./ -m $m -M $M -N 0 -p 40 &> ustacks_${highcov[i]}.oe
done
```
### Ustacks: samples with coverage 10-50x. Here we allow secondary reads as they do not have so high coverages
```
norm=(CDio2 CBor2 CEdw1 CEdw2 CEdw17 APac1 APacLMG ABul1 AGri2 ACre1 ACre8 ACar3 ACar5 AGra1 AGra2 PNat1 PNat2 XPP3 PYel1 POpi2 PNNe1 PNNe2 PGav76 PGav96 PHut8 XPLL1 PLLh1 PLLo1 PBNi45 PBNi01 PBBa54b PBBa65 PBDi1 PAHa2 PEle2 PEle3 aPEle)
id=(9 11 12 13 14 15 16 17 19 22 23 24 25 26 27 28 29 31 32 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53)
for ((i = 0; i < 37; i++)); do
  ustacks -f $samples${norm[i]}.1.fq.gz -i ${id[i]} --name ${norm[i]} -o ./ -m $m -M $M -p 40 &> ustacks_${norm[i]}.oe
done
```
### Outgroups: replace U by O to snps, tags and models files as it comes from an haploid genome and thus everything is homozygotic by definition
```
for ((i = 0; i < 5; i++)); do
  zcat ${outgroups[i]}.tags.tsv.gz | sed 's/U/O/g' | gzip > ${outgroups[i]}bo.tags.tsv.gz #Replace U by O and save in a file with the same name but with sufix bo
  zcat ${outgroups[i]}.models.tsv.gz | sed 's/U/O/g' | gzip > ${outgroups[i]}bo.models.tsv.gz
  zcat ${outgroups[i]}.snps.tsv.gz | sed 's/U/O/g' | gzip > ${outgroups[i]}bo.snps.tsv.gz
  rm ${outgroups[i]}.tags.tsv.gz #Remove the original files
  rm ${outgroups[i]}.models.tsv.gz
  rm ${outgroups[i]}.snps.tsv.gz
  mv ${outgroups[i]}bo.tags.tsv.gz ${outgroups[i]}.tags.tsv.gz #Rename the modified files as the original files
  mv ${outgroups[i]}bo.models.tsv.gz ${outgroups[i]}.models.tsv.gz
  mv ${outgroups[i]}bo.snps.tsv.gz ${outgroups[i]}.snps.tsv.gz
done
```
### Cstacks: creating the catalog
```
cstacks -P ./ -n $n -p 40 -M $popmap &> cstacks.oe #Using the predefined value for n and the complete popmap
```
### Sstacks: matching samples to the catalog to call the alleles and know which of the catalog samples are present in every sample
```
sstacks -P ./ -p 40 -M $popmap &> sstacks.oe #Using the complete popmap
```
### tsv2bam: making a bam file with the information of the matches to have the information in form of an alignment
```
tsv2bam -P ./ -M $popmap -R ../../cleaned_sel_R2/ -t 40 &> tsv2bam.oe
```
### gstacks: incorporating PE data, merging the R2 and using SNP calling method Maruki to refine the SNP calls
```
gstacks -P ./ -M $popmap -t 40 &> gstacks.oe
```
### populations: selecting the samples and the loci we want to build the phylogenies and making the final alignments
```
popmap=$HOME/Durham/ddRAD/info/popmap_ingroup.tsv # We will work with no outgroups

mkdir p18_in p21_in p26_in

samples=($(echo $(cat $popmap | cut -f1)))
lensamples=$(echo ${#samples[@]})

pops=(p18_in p21_in p26_in)
vpops=(18 21 26)
lenpops=$(echo ${#pops[@]})

samplesori=($(echo $(cat $HOME/Durham/ddRAD/info/samples_map_ingroup_sorted.tsv | cut -f1)))
samplesrep=($(echo $(cat $HOME/Durham/ddRAD/info/samples_map_ingroup_sorted.tsv | cut -f2)))

for ((i = 0; i < $lenpops; i++)); do
  populations -P ./ -O ./${pops[i]}/ -M $popmap -t 40 -p ${vpops[i]} -r 0.4 --min_maf 0.026 --fasta_samples &> populations_${pops[
i]}.oe
  cat ./${pops[i]}/populations.samples.fa | grep -v "Stacks version" > ./${pops[i]}/populations.samplesbo.fa
  rm ./${pops[i]}/populations.samples.fa
  mv ./${pops[i]}/populations.samplesbo.fa ./${pops[i]}/populations.samples.fa

  cd ./${pops[i]}
  python $HOME/ddRAD/scripts/python/alignments_from_stacks_fasta.py populations.samples.fa 2 --fasta

  cd ./out_fasta/
  for file in *.fa; do
    cat $file | grep -A1 -E "a$" | grep -v "^--" | sed -E 's/a$//'> a_$file
    cat $file | grep -A1 -E "b$" | grep -v "^--" | sed -E 's/b$//' > b_$file
    rm $file
    seqtk mergefa a_$file b_$file > $file
  done
  rm a_* b_*

  TriSeq -in *.fa -of phylip -o ../${pops[i]}_concat --missing-filter 98 98

  for ((j = 0; j < $lensamples; j++)); do
   sed -i -e "s/${samplesori[j]}/${samplesrep[j]}/g" ../${pops[i]}_concat.phy
  done
  cd ../../
done
```
#######################
## 9) PyRAD analyses ##
#######################

### Merge read pairs with PEAR to be able to run merged reads protocol in PyRAD ###
```
cd $HOME/Durham/ddRAD/pyrad/raw

for file in *R1.fq.gz; do
  /users/jferrer/programari/pear/bin/pear -f $file -r ${file/_R1.fq.gz/_R2.fq.gz} -o ./merged/${file/_R1.fq.gz/} -m 300 -n 144 -v 10 -j 10 >> pear.log 2>&1
  rm merged/*unassembled* merged/*discarded*
  gzip merged/${file/_R1.fq.gz/}*
done
```
### PyRAD analyses for merged reads with the 2 different parameterizations aand 3 different completeness (example clust 94 and 75% completeness) ###
```
cd $HOME/Durham/ddRAD/pyrad
/users/jferrer/programari/pyrad-3.0.66/pyrad/pyRAD.py -p params_clust94_percent75.txt
```
####################################################################
## 10) STACKS analyses with Calonectris borealis reference genome ##
####################################################################

### Script to map the higher m parameterization catalog to the Calonectris borealis genome ###
```
cd $HOME/Durham/ddRAD/stacks.denovo/highermPE
bwa_db=$HOME/Durham/ddRAD/genomes/index/Cbor

bwa mem -M -L 7 -t 10 $bwa_db catalog.fa.gz | /users/DB_shared/samtools-1.3.1/samtools view -bS -o catalog.bam

stacks-integrate-alignments -P ./ -B catalog.bam -O ./integrate_alignments
```

#### Then, run populations as shown before for the optimal dataset ###

################################################################################
## 11) Calculate parsimony informative sites (PIS) per locus for each dataset ##
################################################################################

### Phyluce script to calculate the number of PIS per locus ###
```
phyluce_align_get_informative_sites --alignments out_fasta/ --output PIS.tsv --input-format fasta --log-path ./ --cores 7
```
### Create a file called PIS_all.tsv concatenating all the PIS.tsv files and adding the dataset and the completeness in the end ###

##########################################################################################
## 12) Phylogenetic analyses using RAxML with only the ingroup samples for each dataset ##
##########################################################################################
```
cd $HOME/Durham/ddRAD/
find `pwd` -name *.phy  > list_phy_files

mkdir /users/jferrer/Durham/phylo_analyses/RAxML/ingroup/
output=/users/jferrer/Durham/phylo_analyses/RAxML/ingroup/
```
```
#!/bin/bash
#$ -V
#$ -cwd
#$ -t 1:21

FILE=$(sed "${SGE_TASK_ID}q;d" ./list_phy_files)
raxmlHPC-PTHREADS-SSE3 -f a -d -N 1000 -s ${FILE} -m GTRGAMMA -T 30 -x $RANDOM -p $RANDOM -n $(echo ${FILE} | tr . _ | cut -d "_" -f2-3)_N1000 -w $output
```
################################################################################################################################
## 13) Plot of sum of BS (for the nodes present in all the trees) for each tree (dataset) and violin plots of PIS per dataset ##
################################################################################################################################
```
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)

setwd("~/Dropbox/Tesi/Data/Durham/ddRAD/sum_stats/")

d <- read.delim("PIS_all.tsv")
d$dataset_completeness <- as.factor(paste(d$dataset, d$completeness, sep="_"))
d$pPIS <- d$informative_sites/d$length

g <- read.delim("datasets_phylostats.txt")

completeness <- rep(c(65,75,95),7)

BSsum <- ggplot(g, aes(x=dataset_completeness, y=BS_sum)) +
  theme_minimal() +
  geom_point(aes(color=dataset), size=5) +
  scale_color_viridis_d(name = "Dataset", labels = c("stacks default", "stacks higher m", "stacks higher n", "stacks optimal", "pyrad clust 89", "pyrad clust 94", "stacks refmap")) +
  scale_x_discrete(labels=completeness, name="") +
  ylab("Node bootstrap support sum")

pPIS <- ggplot(d, aes(x=dataset_completeness, y=pPIS)) +
  theme_minimal() +
  geom_violin(aes(fill=dataset), draw_quantiles = c(0.25,0.5,0.75), color="#757575") +
  scale_fill_viridis_d(name = "Dataset", labels = c("stacks default", "stacks higher m", "stacks higher n", "stacks optimal", "pyrad clust 89", "pyrad clust 94", "stacks refmap")) +
  scale_x_discrete(labels=completeness, name="Dataset completeness (%)") +
  ylab("Proportion of PIS")

gA=ggplot_gtable(ggplot_build(BSsum))
gB=ggplot_gtable(ggplot_build(pPIS))

maxWidth = grid::unit.pmax(gA$widths, gB$widths)
gA$widths <- as.list(maxWidth)
gB$widths <- as.list(maxWidth)

grid.newpage()
BSsum_pPIS <-grid.arrange(arrangeGrob(gA,gB,nrow=2,heights=c(1,1.5)))

ggsave("propPIS_BSsum_dataset.pdf", BSsum_pPIS, device="pdf", units="cm", width=30, height=21)
```
########################################################################
## 14) Build 75 and 95% completeness datasets for comparison with UCE ##
########################################################################

### Remove species with no correspondence in the UCE dataset and filter out positions with more than 25% and 5% missing data ###
```
TriSeq -in highermPE/p20_fgla/out_fasta/*.fa -of phylip -rm PAHa2 --missing-filter 25 25 -o ddRAD_75
TriSeq -in highermPE/p25_fgla/out_fasta/*.fa -of phylip -rm PAHa2 --missing-filter 5 5 -o ddRAD_95
```
### Change names so as they are the same in ddRAD and UCEs ###
```
ddRAD=($(cat correspondence_ddRAD-UCE-concat.tsv | grep -v "ddRAD" | cut -f1))
concat=($(cat correspondence_ddRAD-UCE-concat.tsv | grep -v "ddRAD" | cut -f3))
for ((i = 0; i < ${#ddRAD[@]}; i++)); do sed -iE "s/${ddRAD[i]}/${concat[i]}/" ddRAD_75.phy; done
for ((i = 0; i < ${#ddRAD[@]}; i++)); do sed -iE "s/${ddRAD[i]}/${concat[i]}/" ddRAD_95.phy; done
```
#########################################
## 15) SNPs dataset for SNAPP analyses ##
#########################################

### Convert the vcf file to SNAPP input file using TriSeq and the partitions file ###
```
python vcf2snapp.py highermPE/p20_fgla/populations.single_snp.vcf
```
###############
## End of script ##
###############
