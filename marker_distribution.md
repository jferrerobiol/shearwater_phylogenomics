# Marker distribution

#################################################################################
## 1) Map ddRAD fragments to the Puffinus mauretanicus and to the Calypte anna genome with blast ##
#################################################################################

### Build databases for blast ###
```
$HOME/programari/ncbi-blast-2.8.1+/bin/makeblastdb -in $HOME/ref_sequences/GCA_003957555.1_bCalAnn1_v1.p_genomic.fna -dbtype nucl -out $HOME/Durham/mapping/index/blast/Calann_genome
$HOME/programari/ncbi-blast-2.8.1+/bin/makeblastdb -in $HOME/ref_sequences/final.genome.scf.fasta.masked -dbtype nucl -out $HOME/Durham/mapping/index/blast/Pufmau_genome
```
### Create and save a file with ddRAD loci for the Puffinus mauretanicus sample PMau2 in the pyRAD clust 89 75% data set ###
```
cat Shearwaters_clust89_percent75.loci | grep "PMau2" | sed 's/-//g' | sed 's/.assembled//g'| sed -E 's/      /\n/g' > /users/jferrer/Durham/mapping/ddRAD_PMau2.fab

i=1; while read p; do if [ $p = ">PMau2" ]; then printf $p\_$i'\n'; i=$(( ++i )); else printf $p'\n'; fi; done < ddRAD_PMau2.fab > ddRAD_PMau2.fa

rm ddRAD_PMau2.fab
```
### Run blast ###
```
$HOME/programari/ncbi-blast-2.8.1+/bin/blastn -query $HOME/Durham/mapping/input/ddRAD_PMau2.fa -db $HOME/Durham/mapping/index/blast/Calann_genome -evalue 1e-3 -max_target_seqs 1 -out $HOME/Durham/mapping/blastn/ddRAD_Calann -outfmt 6 -num_threads 10
$HOME/programari/ncbi-blast-2.8.1+/bin/blastn -query $HOME/Durham/mapping/input/ddRAD_PMau2.fa -db $HOME/Durham/mapping/index/blast/Pufmau_genome -evalue 1e-3 -max_target_seqs 1 -out $HOME/Durham/mapping/blastn/ddRAD_Pufmau -outfmt 6 -num_threads 10
```
###################################################################################
## 2) Map UCE fragments to the Puffinus mauretanicus and to the Calypte anna genome with blast ##
###################################################################################

### Create and save a file with ddRAD loci for the Puffinus mauretanicus sample PMau2 from the 75% dataset ###
```
for file in *fa; do cat $file | grep -A1 "PMau2" >> ../UCE_PMau2.fab; done

i=1; while read p; do if [ $p = ">PMau2" ]; then printf $p\_$i'\n'; i=$(( ++i )); else printf $p'\n'; fi; done < UCE_PMau2.fab > UCE_PMau2.fa

rm UCE_PMau2.fab
```
### Run blast ###
```
$HOME/programari/ncbi-blast-2.8.1+/bin/blastn -query $HOME/Durham/mapping/input/UCE_PMau2.fa -db $HOME/Durham/mapping/index/blast/Calann_genome -evalue 1e-3 -max_target_seqs 1 -out $HOME/Durham/mapping/blastn/UCE_Calann -outfmt 6 -num_threads 10
$HOME/programari/ncbi-blast-2.8.1+/bin/blastn -query $HOME/Durham/mapping/input/UCE_PMau2.fa -db $HOME/Durham/mapping/index/blast/Pufmau_genome -evalue 1e-3 -max_target_seqs 1 -out $HOME/Durham/mapping/blastn/UCE_Pufmau -outfmt 6 -num_threads 10
```
##############################################################################  
## 3) Map Calypte anna CDS to Puffinus mauretanicus and Calypte anna genome ##
##############################################################################

### Build databases for gmap ###
```
/soft/gmap-2017-06-20/bin/gmap_build -d Durham/mapping/index/gmap/Cal_ann ref_sequences/GCA_003957555.1_bCalAnn1_v1.p_genomic.fna # Cal_ann
/soft/gmap-2017-06-20/bin/gmap_build -d $HOME/Durham/mapping/index/gmap/Puf_mau $HOME/ref_sequences/final.genome.scf.fasta.masked # Puf_mau
```
### Run gmap ###
```
/soft/gmap-2017-06-20/bin/gmap -D $HOME/Durham/mapping/index/gmap/Cal_ann -d Cal_ann -B 5 -t 20 $HOME/Durham/mapping/input/APP-008.cds.fa
/soft/gmap-2017-06-20/bin/gmap -D $HOME/Durham/mapping/index/gmap/Puf_mau -d Puf_mau -B 5 -t 20 $HOME/Durham/mapping/input/APP-008.cds.fa
```
############################################################################################
## 4) Plot the mappings of ddRAD, UCEs and CDS on the Puffinus mauretanicus longest scaffolds as in Harvey et al. 2016 ##
############################################################################################

### Prepare a file with the lengths of each scaffold ###

## Script to create a file with 2 columns: the first one with the name of the scaffold and the second with the length ##

## Usage: python fasta_len.py <assembly.fa> <output.txt>
```
import sys
from Bio import SeqIO

FastaFile = open(sys.argv[1], 'r')
output = open(sys.argv[2], 'w')

for rec in SeqIO.parse(FastaFile, 'fasta'):
    name = rec.id
    seq = rec.seq
    seqLen = len(rec)
    output.write(name+'\t'+str(seqLen)+'\n')

FastaFile.close()
output.close
```
## Run fasta_len.py ##
```
python fasta_len.py ../Reference_genomes/final.genome.scf.fasta.masked Pufmau_genome_scaffolds.txtp

cat Pufmau_genome_scaffolds.txtp | sort -nr  -k2 > Pufmau_genome_scaffolds.txt # Sort it

rm Pufmau_genome_scaffolds.txtp # Remove the intermediate file
```
### Modify the output of gmap to read it in R ###
```
cat CDS_Pufmau.log | grep "Path 1" | sed -E 's/.*genome //g' | tr ':' '\t' | sed 's/\.\./        /g' | sed 's/,//g' | sed -E 's/\([0-9]+\ bp\)/+/g' | sed -E 's/\(\-[0-9]+\ bp\)/-/g' | sed -E 's/\ +/  /g' > CDS_Pufmau
```
### R script to make the plots with the 20 longest scaffolds ###
```
#!/bin/R

library(dplyr)
setwd("~/Dropbox/Tesi/Data/Durham/Mapping_genome/input")

## Plot longest scaffolds with Puffinus mauretanicus as reference genome ##

# Read input files

scaffolds <- read.table("Pufmau_genome_scaffolds.txt", sep="\t", header=T) # Reference info
ddRAD_hits <- read.table("ddRAD_Pufmau", sep="\t", header=F) # NCBI Blastn hit table for ddRAD data
UCE_hits <- read.table("UCE_Pufmau", sep="\t", header=F) # NCBI Blastn hit table for UCE data
CDS_hits <- read.table("CDS_Pufmau", sep="\t", header=F) # File specifying coordinates of coding sequences modified from gmap results
scaffplot <- scaffolds[1:20,]
ddRAD_hits <- select(ddRAD_hits, V2, V9, V10)
UCE_hits <- select(UCE_hits, V2, V9, V10)

# Process input, get information for plotting

sep <- (8/length(scaffplot)) # determine distance between reference fragments
scaf.lengths <- 5*(scaffplot$Size/max(scaffplot$Size)) # determine the size for each fragment
seg.pos <- vector() # vector to record Y positions for axis labels
par(mar=c(5,9,1,3)) # set up plot margins

# Plot
plot(NA, xlim=c(0,6), ylim=c(0,20), xlab="Length (Mbp)", ylab="", axes=F) # empty plot
for(i in 1:nrow(scaffplot)) {
	segments(0, i, scaf.lengths[i], i, col="gray", lwd=6) # plot fragment
	seg.pos <- c(seg.pos, i) # record fragment Y positions for axis labels
	## Add in all hits
	#ddRAD
	ddRAD_hit.rows <- ddRAD_hits[which(ddRAD_hits$V2 == as.character(scaffplot$RefSeq[i])),]
	ddRAD_hit.pos <- ((ddRAD_hit.rows[3]-((ddRAD_hit.rows[3]-ddRAD_hit.rows[2])/2))/max(scaffplot$Size))*5
	points(ddRAD_hit.pos[,1], rep(i+0.2, length(ddRAD_hit.pos[,1])), col="#03A9F4", pch=19, cex=0.6)
	#UCE
	UCE_hit.rows <- UCE_hits[which(UCE_hits$V2 == as.character(scaffplot$RefSeq[i])),]
	UCE_hit.pos <- ((UCE_hit.rows[3]-((UCE_hit.rows[3]-UCE_hit.rows[2])/2))/max(scaffplot$Size))*5
	points(UCE_hit.pos[,1], rep(i-0.2, length(UCE_hit.pos[,1])), col="#FF5252", pch=19, cex=0.6)
	#CDS
	CDS_hit.rows <- CDS_hits[which(CDS_hits$V1 == as.character(scaffplot$RefSeq[i])),]
	CDS_hit.ini <- (CDS_hit.rows[3]/max(scaffplot$Size))*5
	CDS_hit.fin <- (CDS_hit.rows[2]/max(scaffplot$Size))*5
	
	for (j in 1:length(CDS_hit.ini[,1])) {
	  segments(CDS_hit.ini[j,1], i, CDS_hit.fin[j,1], i, col="#FFC107", lwd=6)
	}
}
interval <- 5000000/max(scaffplot$Size) # interval for x axis ticks
axis(1, at=c(0, interval, interval*2, interval*3, interval*4, interval*5, interval*6, interval*7), labels=c(0, 1, 2, 3, 4, 5, 6, 7)) # x axis ticks/labels may need to be tweaked depending on reference
axis(2, at=seg.pos, labels=scaffplot$RefSeq, las=2, col="white") # y axis ticks/labels
dev.off()
```
##########################################################################################################
## 5) Liftover of ddRAD fragments mapped to the Puffinus mauretanicus genome to the Calypte anna genome ##
##########################################################################################################

### Prepare a file with the lengths of each scaffold ###

## Run fasta_len.py ##
```
python ../scripts/fasta_len.py GCA_003957555.1_bCalAnn1_v1.p_genomic.fna ../input/Calann_genome_chromosomes.txtp

cat Calann_genome_chromosomes.txtp | sort -nr  -k2 > Calann_genome_chromosomes.txt # Sort it

rm Calann_genome_chromosomes.txtp # Remove the intermediate file
```
### Align Puffinus mauretanicus and Calypte anna genomes with Satsuma2 ###
```
SatsumaSynteny2 -q $HOME/ref_sequences/final.genome.scf.fasta.masked -t $HOME/ref_sequences/GCA_003957555.1_bCalAnn1_v1.p_genomic.fna -o $HOME/Durham/mapping/satsuma2 -slaves 2 -threads 32
```
### Prepare config file for kraken ###

#### [genomes]
```
Puf_mau        /Users/joan/Dropbox/Tesi/Data/Durham/Mapping_genome/Reference_genomes/final.genome.scf.fasta.masked
Cal_ann        /Users/joan/Dropbox/Tesi/Data/Durham/Mapping_genome/Reference_genomes/GCA_003957555.1_bCalAnn1_v1.p_genomic.fna
```

#### [pairwise-maps]
```
Cal_ann Puf_mau /Users/joan/Dropbox/Tesi/Data/Durham/Mapping_genome/input/satsuma_summary_onlychromosomes.chained.out
```
### Run the liftover with kraken ###
```
RunKraken -c kraken.config -s ddRAD_Pufmau_bo.gff -S Puf_mau -T Cal_ann
```
### Convert the output of kraken to a file readable in R ###
```
cat mapped.gtf | sed -E 's/^(CM.*\.1)\_.*sequence/\1/g' > ddRAD_kraken_liftover_Pufmau2Calann.gff
```
### Modify the output of gmap to read it in R ###
```
cat CDS_Calann.log | grep "Path 1" | sed -E 's/.*genome //g' | tr ':' '\t' | sed 's/\.\./        /g' | sed 's/,//g' | sed -E 's/\([0-9]+\ bp\)/+/g' | sed -E 's/\(\-[0-9]+\ bp\)/-/g' | sed -E 's/\ +/  /g' > CDS_Calann
```
#################################################################
## 6) Plot ddRAD, UCEs and CDS on the Calypte anna chromosomes ##
#################################################################
```
#!/bin/R

### Plot Cal_ann chromosomes and UCEs, ddRAD and CDS on them ###

# Read input files

chromosomes <- read.table("Calann_genome_chromosomes.txt", sep="\t", header=T) # Reference info
ddRAD_hits <- read.table("ddRAD_kraken_liftover_Pufmau2Calann.gff", sep="\t", header=F) # gff liftover file for ddRAD data
UCE_hits <- read.table("UCE_Calann", sep="\t", header=F) # NCBI Blastn hit table for UCE data
CDS_hits <- read.table("CDS_Calann", sep="\t", header=F) # File specifying coordinates of coding sequences modified from gmap results
ddRAD_hits <- select(ddRAD_hits, V1, V4, V5)
UCE_hits <- select(UCE_hits, V2, V9, V10)

# Process input, get information for plotting

scaf.lengths <- 5*(chromosomes$Size/max(chromosomes$Size)) # determine the size for each fragment
seg.pos <- vector() # vector to record Y positions for axis labels

par(mar=c(5,9,1,3)) # set up plot margins

# Plot
plot(NA, xlim=c(0,6), ylim=c(0,66), axes=F, xlab="", ylab="") # empty plot
for(i in 1:nrow(chromosomes)) {
  segments(0, i*2, scaf.lengths[i], i*2, col="gray", lwd=10) # plot fragment
  seg.pos <- c(seg.pos, i*2) # record fragment Y positions for axis labels
  ## Add in all hits
  #ddRAD
  ddRAD_hit.rows <- ddRAD_hits[which(ddRAD_hits$V1 == as.character(chromosomes$RefSeq[i])),]
  ddRAD_hit.pos <- ((ddRAD_hit.rows[3]-((ddRAD_hit.rows[3]-ddRAD_hit.rows[2])/2))/max(chromosomes$Size))*5
  points(ddRAD_hit.pos[,1], rep(i*2+0.4, length(ddRAD_hit.pos[,1])), col="#03A9F4", pch=19, cex=1)
  #UCE
  UCE_hit.rows <- UCE_hits[which(UCE_hits$V2 == as.character(chromosomes$RefSeq[i])),]
  UCE_hit.pos <- ((UCE_hit.rows[3]-((UCE_hit.rows[3]-UCE_hit.rows[2])/2))/max(chromosomes$Size))*5
  points(UCE_hit.pos[,1], rep(i*2-0.4, length(UCE_hit.pos[,1])), col="#FF5252", pch=19, cex=1)
  #CDS
  CDS_hit.rows <- CDS_hits[which(CDS_hits$V1 == as.character(chromosomes$RefSeq[i])),]
  CDS_hit.ini <- (CDS_hit.rows[3]/max(chromosomes$Size))*5
  CDS_hit.fin <- (CDS_hit.rows[2]/max(chromosomes$Size))*5
  
  for (j in 1:length(CDS_hit.ini[,1])) {
    segments(CDS_hit.ini[j,1], i*2, CDS_hit.fin[j,1], i*2, col="#FFC107", lwd=10)
  }
}
interval <- 100000000/max(chromosomes$Size) # interval for x axis ticks
axis(1, at=seq(0, interval*10, interval), labels=seq(0,200,20), cex.axis=3, line=-1, col="white")
axis(1, at=seq(0, interval*10, interval), lwd=3, line=-2,labels=F) # x axis ticks/labels may need to be tweaked depending on reference
axis(2, at=seg.pos, labels=chromosomes$chromosome, las=2, col="white", cex.axis=3) # y axis ticks/labels
title(xlab="Length (Mbp)", ylab="", line=3, cex.lab=3)
title(xlab="", ylab="chromosome", line=7, cex.lab=3)
dev.off()
```
##################################################
## 7) Prepare input files for regioneR analyses ##
##################################################

### Modify the input files from the previous analyses so as they can be read by regioneR ###
```
cat Calann_genome_chromosomes.txt | cut -f1,3 | grep -v "RefSeq" > ../regioneR/Calann_genome_chromosomes.txt
gff2bed < ddRAD_kraken_liftover_Pufmau2Calann.gff > ../regioneR/ddRAD_kraken_liftover_Pufmau2Calann.bed
cat UCE_Calann | cut -f2,9,10 | awk '{print $1 "      " $2 "  " $3 "  " $3-$2}' | grep -v "-" | cut -f1-3 | sed 's/$/        +/g'> ../regioneR/UCE_Calann.txt; cat UCE_Calann | cut -f2,9,10 | awk '{print $1 "   " $2 "  " $3 "  " $3-$2}' | grep "-" | awk '{print $1 "        " $3 "  " $2 "  -"}' >> ../regioneR/UCE_Calann.txt
cat CDS_Calann | grep "+" > ../regioneR/CDS_Calann.txt; cat CDS_Calann | grep "-" | awk '{print $1 "        " $3 " " $2 "  " $4}' >> ../regioneR/CDS_Calann.txt
```
### Sort and remove the information for the sex chromosomes ###
```
cat UCE_Calann.txt | sort -k1,1 -k2n > UCE_Calann_sorted.txt
cat CDS_Calann.txt | sort -k1,1 -k2n > CDS_Calann_sorted.txt

cat Calann_genome_chromosomes.txt | grep -v "CM012146" | grep -v "CM012142" > Calann_genome_chromosomes_nosex.txt
cat UCE_Calann_sorted.txt | grep -v "CM012146" | grep -v "CM012142" > UCE_Calann_sorted_nosex.txt
cat ddRAD_kraken_liftover_Pufmau2Calann.bed | grep -v "CM012146" | grep -v "CM012142" > ddRAD_kraken_liftover_Pufmau2Calann_nosex.bed
cat CDS_Calann_sorted.txt | grep -v "CM012146" | grep -v "CM012142" > CDS_Calann_sorted_nosex.txt
```
### Calculate the mean distance to the closest marker of the same type ###
```
closest-features --closest --no-ref --no-overlaps --dist ddRAD_kraken_liftover_Pufmau2Calann_nosex.bed ddRAD_kraken_liftover_Pufmau2Calann_nosex.bed | cut -d'|' -f2 | grep -v NA | awk '{ if($1<=0){ $1*= -1;} print $1;}' > ddRAD_distances_nosex.txt
closest-features --closest --no-ref --no-overlaps --dist UCE_Calann_sorted_nosex.txt UCE_Calann_sorted_nosex.txt | cut -d'|' -f2 | grep -v NA | awk '{ if($1<=0){ $1*= -1;} print $1;}' > UCE_distances_nosex.txt
```
#############################################
## 8) regioneR permutation tests and plots ##
#############################################
```
#!/bin/R

library(BiocManager)
library(regioneR)
library(ggplot2)
library(scales)
library(ggpubr)

#Set directory and seed
setwd("~/Dropbox/Tesi/Data/Durham/Mapping_genome/regioneR/")
set.seed(12345)

#Load genome
chr <- read.delim("Calann_genome_chromosomes_nosex.txt", header = F)
calann <- getGenome(chr)

#Load region sets
UCE <- toGRanges("UCE_Calann_sorted_nosex.txt")
ddRAD <- toGRanges("ddRAD_kraken_liftover_Pufmau2Calann_nosex.bed")
CDS <- toGRanges("CDS_Calann_sorted_nosex.txt")

#Check overlaps and mean distances
over_UCE_ddRAD <- numOverlaps(UCE, ddRAD)
over_UCE_CDS <- numOverlaps(UCE, CDS)
over_ddRAD_CDS <- numOverlaps(ddRAD, CDS)
md_UCE_CDS <- meanDistance(UCE, CDS)
md_ddRAD_CDS <- meanDistance(ddRAD, CDS)
md_ddRAD_UCE <- meanDistance(ddRAD, UCE)

#Perform overlap tests
ddRAD_gene <- overlapPermTest(ddRAD, CDS, ntimes=5000, genome=calann)
UCE_gene <- overlapPermTest(UCE, CDS, ntimes=5000, genome=calann)
UCE_ddRAD <- overlapPermTest(UCE, ddRAD, ntimes=5000, genome=calann)

#Perform distance tests
ddRAD_gene_dist <- permTest(A=ddRAD, B=CDS, ntimes=5000,
                                randomize.function=randomizeRegions,
                                evaluate.function= meanDistance,
                                genome=calann, mc.set.seed=FALSE, mc.cores=4)

UCE_gene_dist <- permTest(A=UCE, B=CDS, ntimes=5000,
                            randomize.function=randomizeRegions,
                            evaluate.function= meanDistance,
                            genome=calann, mc.set.seed=FALSE, mc.cores=4)

UCE_ddRAD_dist <- permTest(A=UCE, B=ddRAD, ntimes=5000,
                            randomize.function=randomizeRegions,
                            evaluate.function= meanDistance,
                            genome=calann, mc.set.seed=FALSE, mc.cores=4)

#Perform distance test among the same data set to calculate the expected distances
ddRAD_ddRAD_dist <- permTest(A=ddRAD, B=ddRAD, ntimes=5000,
                           randomize.function=randomizeRegions,
                           evaluate.function= meanDistance,
                           genome=calann, mc.set.seed=FALSE, mc.cores=4)

UCE_UCE_dist <- permTest(A=UCE, B=UCE, ntimes=5000,
                             randomize.function=randomizeRegions,
                             evaluate.function= meanDistance,
                             genome=calann, mc.set.seed=FALSE, mc.cores=4)

#Histograms of distances among markers
dd <- read.delim("ddRAD_distances_nosex.txt", header=F)
uc <- read.delim("UCE_distances_nosex.txt", header=F)

ddp <- ggplot(dd, aes(x=V1)) +
  theme_minimal() +
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(), axis.text.x = element_text(angle=90)) +
  geom_histogram(binwidth = 10000, fill="#03A9F4") +
  scale_x_continuous(limits=c(0,500000), breaks=seq(0,500000,20000), labels = scales::comma) +
  scale_y_continuous(limits=c(0,900), breaks=seq(0,900,100)) +
  xlab("Distance to nearest ddRAD locus (bp)") + ylab("count") +
  geom_vline(aes(xintercept=mean(V1)), color="#ffab00", size=1) +
  geom_vline(aes(xintercept=median(V1)), color="#2e7d32", size=1) +
  geom_vline(aes(xintercept=quantile(V1, probs=0.25)), color="#2e7d32", alpha=0.5, size=1) +
  geom_vline(aes(xintercept=quantile(V1, probs=0.75)), color="#2e7d32", alpha=0.5, size=1)

ucp <- ggplot(uc, aes(x=V1)) +
  theme_minimal() +
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank(), axis.text.x = element_text(angle=90)) +
  geom_histogram(binwidth = 10000, fill="#FF5252") +
  scale_x_continuous(limits=c(0,500000), breaks=seq(0,500000,20000), labels = scales::comma) +
  scale_y_continuous(limits=c(0,900), breaks=seq(0,900,100)) +
  xlab("Distance to nearest UCE locus (bp)") + ylab("count") +
  geom_vline(aes(xintercept=mean(V1)), color="#ffab00", size=1) +
  geom_vline(aes(xintercept=median(V1)), color="#2e7d32", size=1) +
  geom_vline(aes(xintercept=quantile(V1, probs=0.25)), color="#2e7d32", alpha=0.5, size=1) +
  geom_vline(aes(xintercept=quantile(V1, probs=0.75)), color="#2e7d32", alpha=0.5, size=1)

ddp_ucp <- ggarrange(ddp, ucp, ncol=1, nrow=2, labels = "auto")
ggsave("distance_nearest_locus.pdf", ddp_ucp, device="pdf", units="cm", width=20, height=20, limitsize=FALSE)

## Plots of permutation tests ##
# ddRAD-ddRAD distance #
dd_dd_dist <- as.data.frame(ddRAD_ddRAD_dist$meanDistance$permuted)

#Calculate Z-score
dd_dd_dist_obs <- rbind(dd_dd_dist, mean(dd$V1))
zscore <- as.numeric(round(tail(scale(dd_dd_dist_obs),n=1),3))

dd_dd_distp <- ggplot(dd_dd_dist, aes(x=ddRAD_ddRAD_dist$meanDistance$permuted)) +
  theme_minimal() +
  ggtitle(paste0("p-value: 0.0002 \nZ-score: ",zscore," \nnperm: 5000")) +
  theme(panel.grid.minor.x=element_blank(), axis.text.x = element_text(angle=90), plot.title = element_text(size=10, hjust = 0.5)) +
  geom_histogram(binwidth = 300, fill="#03A9F4") +
  geom_vline(data=dd, aes(xintercept=mean(V1)), color="#ffab00", size=1) +
  geom_vline(aes(xintercept=mean(ddRAD_ddRAD_dist$meanDistance$permuted)), color="#424242", size=1) +
  geom_vline(aes(xintercept=quantile(ddRAD_ddRAD_dist$meanDistance$permuted, probs=0.05)), color="#d32f2f", size=1) +
  geom_segment(x = 82441.02, y = 375, xend = 97570.84, yend = 375,color="#424242", arrow=arrow(length = unit(0.03, "npc"))) +
  geom_segment(x = 97570.84, y = 375, xend = 82441.02, yend = 375,color="#424242", arrow=arrow(length = unit(0.03, "npc"))) +
  scale_x_continuous(name="Average distance between the closest ddRAD loci (bp)", limits=c(80000,105000), breaks=seq(80000,105000,2500)) +
  scale_y_continuous(limits=c(0,550), breaks=seq(0,550,100))

# UCE-UCE distance #

uc_uc_dist <- as.data.frame(UCE_UCE_dist$meanDistance$permuted)

#Calculate Z-score
uc_uc_dist_obs <- rbind(uc_uc_dist, mean(uc$V1))
zscore <- as.numeric(round(tail(scale(uc_uc_dist_obs),n=1),3))

uc_uc_distp <- ggplot(uc_uc_dist, aes(x=UCE_UCE_dist$meanDistance$permuted)) +
  theme_minimal() +
  ggtitle(paste0("p-value: 0.0002 \nZ-score: ",zscore," \nnperm: 5000")) +
  theme(panel.grid.minor.x=element_blank(), axis.text.x = element_text(angle=90), plot.title = element_text(size=10, hjust = 0.5)) +
  geom_histogram(binwidth = 3000, fill="#FF5252") +
  geom_vline(data=uc, aes(xintercept=mean(V1)), color="#ffab00", size=1) +
  geom_vline(aes(xintercept=mean(UCE_UCE_dist$meanDistance$permuted)), color="#424242", size=1) +
  geom_vline(aes(xintercept=quantile(UCE_UCE_dist$meanDistance$permuted, probs=0.05)), color="#d32f2f", size=1) +
  geom_segment(x = 522846.3, y = 500, xend = 98111.14, yend = 500,color="#424242", arrow=arrow(length = unit(0.03, "npc"))) +
  geom_segment(x = 98111.14, y = 500, xend = 522846.3, yend = 500,color="#424242", arrow=arrow(length = unit(0.03, "npc"))) +
  scale_x_continuous(name="Average distance between the closest UCE loci (bp)", limits=c(80000,600000), breaks=seq(80000,600000,50000)) +
  scale_y_continuous(limits=c(0,550), breaks=seq(0,550,100))

own_dist <- ggarrange(dd_dd_distp, uc_uc_distp, ncol=1, nrow=2, labels = "auto")
ggsave("perm_nearest_dist_own.pdf", own_dist, device="pdf", units="cm", width=20, height=20, limitsize=FALSE)

# Overlap ddRAD-gene #

ddRAD_CDS_over <- as.data.frame(ddRAD_gene$numOverlaps$permuted)

#Calculate Z-score
zscore <- round(ddRAD_gene$numOverlaps$zscore, 3)

ddRAD_CDS_overp <- ggplot(ddRAD_CDS_over, aes(x=ddRAD_gene$numOverlaps$permuted)) +
  theme_minimal() +
  ggtitle(paste0("p-value: ",format(round(ddRAD_gene$numOverlaps$pval,4)),"\nZ-score: ",zscore," \nnperm: 5000")) +
  theme(panel.grid.minor.x=element_blank(), axis.text.x = element_text(angle=90), plot.title = element_text(size=10, hjust = 0.5)) +
  geom_histogram(binwidth = 5, fill="#03A9F4") +
  geom_vline(aes(xintercept=ddRAD_gene$numOverlaps$observed), color="#ffab00", size=1) +
  geom_vline(aes(xintercept=mean(ddRAD_gene$numOverlaps$permuted)), color="#424242", size=1) +
  geom_vline(aes(xintercept=quantile(ddRAD_gene$numOverlaps$permuted, probs=0.95)), color="#d32f2f", size=1) +
  geom_segment(x = ddRAD_gene$numOverlaps$observed, y = 325, xend = mean(ddRAD_gene$numOverlaps$permuted), yend = 325,color="#424242", arrow=arrow(length = unit(0.03, "npc"))) +
  geom_segment(x = mean(ddRAD_gene$numOverlaps$permuted), y = 325, xend = ddRAD_gene$numOverlaps$observed, yend = 325,color="#424242", arrow=arrow(length = unit(0.03, "npc"))) +
  scale_x_continuous(name="Number of overlaps between ddRAD loci and CDS", limits=c(1290,1600), breaks=seq(1300,1600,50))

# Overlap UCE-gene #

UCE_CDS_over <- as.data.frame(UCE_gene$numOverlaps$permuted)

#Calculate Z-score
zscore <- round(UCE_gene$numOverlaps$zscore, 3)

UCE_CDS_overp <- ggplot(UCE_CDS_over, aes(x=UCE_gene$numOverlaps$permuted)) +
  theme_minimal() +
  ggtitle(paste0("p-value: ",format(round(UCE_gene$numOverlaps$pval,4), scientific=F),"\nZ-score: ",zscore," \nnperm: 5000")) +
  theme(panel.grid.minor.x=element_blank(), axis.text.x = element_text(angle=90), plot.title = element_text(size=10, hjust = 0.5)) +
  geom_histogram(binwidth = 5, fill="#FF5252") +
  geom_vline(aes(xintercept=UCE_gene$numOverlaps$observed), color="#ffab00", size=1) +
  geom_vline(aes(xintercept=mean(UCE_gene$numOverlaps$permuted)), color="#424242", size=1) +
  geom_vline(aes(xintercept=quantile(UCE_gene$numOverlaps$permuted, probs=0.95)), color="#d32f2f", size=1) +
  geom_segment(x = UCE_gene$numOverlaps$observed, y = 325, xend = mean(UCE_gene$numOverlaps$permuted), yend = 325,color="#424242", arrow=arrow(length = unit(0.03, "npc"))) +
  geom_segment(x = mean(UCE_gene$numOverlaps$permuted), y = 325, xend = UCE_gene$numOverlaps$observed, yend = 325,color="#424242", arrow=arrow(length = unit(0.03, "npc"))) +
  scale_x_continuous(name="Number of overlaps between UCE loci and CDS", limits=c(800,1050), breaks=seq(800,1050,50))

# Overlap UCE-ddRAD #

UCE_ddRAD_over <- as.data.frame(UCE_ddRAD$numOverlaps$permuted)

#Calculate Z-score
zscore <- round(UCE_ddRAD$numOverlaps$zscore, 3)

UCE_ddRAD_overp <- ggplot(UCE_ddRAD_over, aes(x=UCE_ddRAD$numOverlaps$permuted)) +
  theme_minimal() +
  ggtitle(paste0("p-value: ",format(round(UCE_ddRAD$numOverlaps$pval,4), scientific=F),"\nZ-score: ",zscore," \nnperm: 5000")) +
  theme(panel.grid.minor.x=element_blank(), axis.text.x = element_text(angle=90), plot.title = element_text(size=10, hjust = 0.5)) +
  geom_histogram(binwidth = 1, fill="#817EA3") +
  geom_vline(aes(xintercept=UCE_ddRAD$numOverlaps$observed), color="#ffab00", size=1) +
  geom_vline(aes(xintercept=mean(UCE_ddRAD$numOverlaps$permuted)), color="#424242", size=1) +
  geom_vline(aes(xintercept=quantile(UCE_ddRAD$numOverlaps$permuted, probs=0.95)), color="#d32f2f", size=1) +
  geom_segment(x = UCE_ddRAD$numOverlaps$observed, y = 325, xend = mean(UCE_ddRAD$numOverlaps$permuted), yend = 325,color="#424242", arrow=arrow(length = unit(0.03, "npc"))) +
  geom_segment(x = mean(UCE_ddRAD$numOverlaps$permuted), y = 325, xend = UCE_ddRAD$numOverlaps$observed, yend = 325,color="#424242", arrow=arrow(length = unit(0.03, "npc"))) +
  scale_x_continuous(name="Number of overlaps between UCE and ddRAD loci", limits=c(0,35), breaks=seq(0,35,5))

over <- ggarrange(ddRAD_CDS_overp, UCE_CDS_overp, UCE_ddRAD_overp, ncol=1, nrow=3, labels = "auto")
ggsave("overlap.pdf", over, device="pdf", units="cm", width=20, height=30, limitsize=FALSE)

# Distance ddRAD-gene #

ddRAD_CDS_dist <- as.data.frame(ddRAD_gene_dist$meanDistance$permuted)

#Calculate Z-score
zscore <- round(ddRAD_gene_dist$meanDistance$zscore, 3)

ddRAD_CDS_distp <- ggplot(ddRAD_CDS_dist, aes(x=ddRAD_gene_dist$meanDistance$permuted)) +
  theme_minimal() +
  ggtitle(paste0("p-value: ",format(round(ddRAD_gene_dist$meanDistance$pval,4), scientific=F),"\nZ-score: ",zscore," \nnperm: 5000")) +
  theme(panel.grid.minor.x=element_blank(), axis.text.x = element_text(angle=90), plot.title = element_text(size=10, hjust = 0.5)) +
  geom_histogram(binwidth = 250, fill="#03A9F4") +
  geom_vline(aes(xintercept=ddRAD_gene_dist$meanDistance$observed), color="#ffab00", size=1) +
  geom_vline(aes(xintercept=mean(ddRAD_gene_dist$meanDistance$permuted)), color="#424242", size=1) +
  geom_vline(aes(xintercept=quantile(ddRAD_gene_dist$meanDistance$permuted, probs=0.05)), color="#d32f2f", size=1) +
  geom_segment(x = ddRAD_gene_dist$meanDistance$observed, y = 325, xend = mean(ddRAD_gene_dist$meanDistance$permuted), yend = 325,color="#424242", arrow=arrow(length = unit(0.03, "npc"))) +
  geom_segment(x = mean(ddRAD_gene_dist$meanDistance$permuted), y = 325, xend = ddRAD_gene_dist$meanDistance$observed, yend = 325,color="#424242", arrow=arrow(length = unit(0.03, "npc"))) +
  scale_x_continuous(name="Average distance between the closest ddRAD loci and CDS (bp)", limits=c(50000,70000), breaks=seq(50000,70000,2500))

# Distance UCE-gene #

UCE_CDS_dist <- as.data.frame(UCE_gene_dist$meanDistance$permuted)

#Calculate Z-score
zscore <- round(UCE_gene_dist$meanDistance$zscore, 3)

UCE_CDS_distp <- ggplot(UCE_CDS_dist, aes(x=UCE_gene_dist$meanDistance$permuted)) +
  theme_minimal() +
  ggtitle(paste0("p-value: ",format(round(UCE_gene_dist$meanDistance$pval,4), scientific=F),"\nZ-score: ",zscore," \nnperm: 5000")) +
  theme(panel.grid.minor.x=element_blank(), axis.text.x = element_text(angle=90), plot.title = element_text(size=10, hjust = 0.5)) +
  geom_histogram(binwidth = 250, fill="#FF5252") +
  geom_vline(aes(xintercept=UCE_gene_dist$meanDistance$observed), color="#ffab00", size=1) +
  geom_vline(aes(xintercept=mean(UCE_gene_dist$meanDistance$permuted)), color="#424242", size=1) +
  geom_vline(aes(xintercept=quantile(UCE_gene_dist$meanDistance$permuted, probs=0.95)), color="#d32f2f", size=1) +
  geom_segment(x = UCE_gene_dist$meanDistance$observed, y = 325, xend = mean(UCE_gene_dist$meanDistance$permuted), yend = 325,color="#424242", arrow=arrow(length = unit(0.03, "npc"))) +
  geom_segment(x = mean(UCE_gene_dist$meanDistance$permuted), y = 325, xend = UCE_gene_dist$meanDistance$observed, yend = 325,color="#424242", arrow=arrow(length = unit(0.03, "npc"))) +
  scale_x_continuous(name="Average distance between the closest UCE loci and CDS (bp)", limits=c(55000,75000), breaks=seq(55000,75000,2500))

# Distance UCE-ddRAD #

UCE_dd_dist <- as.data.frame(UCE_ddRAD_dist$meanDistance$permuted)

#Calculate Z-score
zscore <- round(UCE_ddRAD_dist$meanDistance$zscore, 3)

UCE_ddRAD_distp <- ggplot(UCE_dd_dist, aes(x=UCE_ddRAD_dist$meanDistance$permuted)) +
  theme_minimal() +
  ggtitle(paste0("p-value: ",format(round(UCE_ddRAD_dist$meanDistance$pval,4), scientific=F),"\nZ-score: ",zscore," \nnperm: 5000")) +
  theme(panel.grid.minor.x=element_blank(), axis.text.x = element_text(angle=90), plot.title = element_text(size=10, hjust = 0.5)) +
  geom_histogram(binwidth = 250, fill="#817EA3") +
  geom_vline(aes(xintercept=UCE_ddRAD_dist$meanDistance$observed), color="#ffab00", size=1) +
  geom_vline(aes(xintercept=mean(UCE_ddRAD_dist$meanDistance$permuted)), color="#424242", size=1) +
  geom_vline(aes(xintercept=quantile(UCE_ddRAD_dist$meanDistance$permuted, probs=0.05)), color="#d32f2f", size=1) +
  geom_segment(x = UCE_ddRAD_dist$meanDistance$observed, y = 275, xend = mean(UCE_ddRAD_dist$meanDistance$permuted), yend = 275,color="#424242", arrow=arrow(length = unit(0.03, "npc"))) +
  geom_segment(x = mean(UCE_ddRAD_dist$meanDistance$permuted), y = 275, xend = UCE_ddRAD_dist$meanDistance$observed, yend = 275,color="#424242", arrow=arrow(length = unit(0.03, "npc"))) +
  scale_x_continuous(name="Average distance between the closest ddRAD and UCE loci (bp)", limits=c(70000,105000), breaks=seq(70000,105000,2500))

dist <- ggarrange(ddRAD_CDS_distp, UCE_CDS_distp, UCE_ddRAD_distp, ncol=1, nrow=3, labels = "auto")
ggsave("mean_distance.pdf", dist, device="pdf", units="cm", width=20, height=30, limitsize=FALSE)
```
