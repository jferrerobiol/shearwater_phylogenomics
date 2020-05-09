########################################
## 1) Calculate NeighbourNet networks ##
########################################
```
#!/bin/R

### Load libraries ###
library(ape)
library(phangorn)

### Read nexus ddRAD_75 sequence file and convert it to DNAbin object ###
phydat <- read.phyDat(file="ddRAD_75.nex", format = "nexus", type = "DNA")
al <- as.DNAbin(phydat)

### Calculate genetic distances with pairwise deletion and write distance for input in Splitstree ###
dist <- dist.dna(al, model = "K80", variance = FALSE,
                 gamma = FALSE, pairwise.deletion = TRUE,
                 base.freq = NULL, as.matrix = FALSE)
write.nexus.dist(dist, file = "ddRAD_75_dist_noN.nex", append = FALSE, upper = FALSE,
                 diag = TRUE, digits = getOption("digits"))
```
#######################################################
## 2) Calculate Patterson's D statistics with Dsuite ##
#######################################################

### Run pyRAD to obtain a VCF file including loci with data in at least 5 taxa ###
`/users/jferrer/programari/pyrad-3.0.66/pyrad/pyRAD.py -p params_clust89_min5.txt -s 7`

##### params_clust89_min5.txt file
```
==** parameter inputs for pyRAD version 3.0.66  **======================== affected step ==
./clust89                 ## 1. Working directory                                 (all)
                          ## 2. Loc. of non-demultiplexed files (if not line 18)  (s1)
                          ## 3. Loc. of barcode file (if not line 18)             (s1)
vsearch                   ## 4. command (or path) to call vsearch (or usearch)    (s3,s6)
muscle                    ## 5. command (or path) to call muscle                  (s3,s7)
AATT,CGG                  ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG)     (s1,s2)
10                        ## 7. N processors (parallel)                           (all)
7                         ## 8. Mindepth: min coverage for a cluster              (s4,s5)
8                         ## 9. NQual: max # sites with qual < 20 (or see line 20)(s2)
.89                       ## 10. Wclust: clustering threshold as a decimal        (s3,s6)
ddrad                 ## 11. Datatype: rad,gbs,pairgbs,pairddrad,(others:see docs)(all)
5                        ## 12. MinCov: min samples in a final locus             (s7)
5                         ## 13. MaxSH: max inds with shared hetero site          (s7)
Shearwaters_clust89_min5           ## 14. Prefix name for final output (no spaces)         (s7)
==== optional params below this line ===================================  affected step ==
                       ## 15.opt.: select subset (prefix* only selector)            (s2-s7)
                       ## 16.opt.: add-on (outgroup) taxa (list or prefix*)         (s6,s7)
                       ## 17.opt.: exclude taxa (list or prefix*)                   (s7)
./raw/merged/*.gz      ## 18.opt.: loc. of de-multiplexed data                      (s2)
                       ## 19.opt.: maxM: N mismatches in barcodes (def= 1)          (s1)
                       ## 20.opt.: phred Qscore offset (def= 33)                    (s2)
                       ## 21.opt.: filter: def=0=NQual 1=NQual+adapters. 2=strict   (s2)
                       ## 22.opt.: a priori E,H (def= 0.001,0.01, if not estimated) (s5)
8                      ## 23.opt.: maxN: max Ns in a cons seq (def=5)               (s5)
                       ## 24.opt.: maxH: max heterozyg. sites in cons seq (def=5)   (s5)
                       ## 25.opt.: ploidy: max alleles in cons seq (def=2;see docs) (s4,s5)
35                     ## 26.opt.: maxSNPs: (def=100). Paired (def=100,100)         (s7)
6,10                   ## 27.opt.: maxIndels: within-clust,across-clust (def. 3,99) (s3,s7)
                       ## 28.opt.: random number seed (def. 112233)              (s3,s6,s7)
1,1                    ## 29.opt.: trim overhang left,right on final loci, def(0,0) (s7)
a,p,v                  ## 30.opt.: output formats: p,n,a,s,v,u,t,m,k,g,* (see docs) (s7)
                       ## 31.opt.: maj. base call at depth>x<mindepth (def.x=mindepth) (s5)
                       ## 32.opt.: keep trimmed reads (def=0). Enter min length.    (s2)
                       ## 33.opt.: max stack size (int), def= max(500,mean+2*SD)    (s3)
                       ## 34.opt.: minDerep: exclude dereps with <= N copies, def=1 (s3)
1                      ## 35.opt.: use hierarchical clustering (def.=0, 1=yes)      (s6)
                       ## 36.opt.: repeat masking (def.=1='dust' method, 0=no)      (s3,s6)
                       ## 37.opt.: vsearch max threads per job (def.=6; see docs)   (s3,s6)
==== optional: list group/clade assignments below this line (see docs) ==================
Ardenna	2	A*
Calonectris	2	C*
Puffinus	2	P*,X*,aPEle
```
### Run Dsuite for Ardenna (the same for Calonectris and Puffinus) ### 

#### Create the input files for Dsuite ##
```
cat Shearwaters_clust89_min5.vcf | grep "^#C" | cut -f10- | tr '\t' '\n' | grep -vE '^a|^C|^P|^X' > ard_list # create a list with Ardenna individuals

vcftools --vcf Shearwaters_clust89_min5.vcf --non-ref-ac-any 1 --keep ard_list --recode --recode-INFO-all --stdout > Ardenna_clust89_min5_dsuite.vcf # generate a vcf file with only biallelic SNPs and Ardenna species
```
#### Modify the ard_list file and save it as popmap_dsuite_ardenna.tsv for input to Dsuite ##

#### Write tree (ardenna.nwk) from scratch based on relationships from the TENT 75 analyses ##

### Run Dsuite Dtrios ##

`/Users/apple/Documents/soft/Dsuite/Build/Dsuite Dtrios -j 300 -t ardenna.nwk Ardenna_clust89_min5_dsuite.vcf popmap_dsuite_ardenna.tsv`

################################################################
## 3) Evaluate ABBA SNPs: example A. tenuirostris - A. grisea ##
################################################################

### Extract ABBA SNPs fixed within species ###
`cat Shearwaters_clust89_min5.vcf | grep -v "^##" | cut -f1,2,5,10-22 | grep -v "1|0" | grep -v "0|1" | grep -v "0|2" | grep -v "2|0" | grep -v "1|2" | grep -v "2|1" | grep -vE "\,.\," | grep -v "./." | awk '{ if (($4 == $5) && ($4 == $6) && ($4 == $7) && ($4 == $8) && ($4 == $9) && ($4 == $10) && ($4 == $13) && ($4 == $14) && ($4 != $11) && ($11 == $12) && ($11 == $15) && ($11 == $16)) { print $0} }'`

##### There are 20 cases of ABBA fragments completely fixed within species

### Count the number of haplotypes per locus: example for locus 56829 ###
```
cat Shearwaters_clust89_min5.vcf | grep -v "^##" | grep "^56829	" | cut -f10- | tr '|' '\t' | awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' | sort | uniq | wc -l

### Calculate the mean number of haplotypes per loci with the same missing data restrictions ###

list=($(cat Shearwaters_clust89_min5.vcf | grep -v "^#" | cut -f1,2,5,10-22 | grep -vE "\,.\," | grep -v "./." | cut -f1 | sort | uniq))
for locus in ${list[@]}; do
cat Shearwaters_clust89_min5.vcf | grep -v "^##" | grep "^$locus	" | cut -f10- | tr '|' '\t' | awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' | sort | uniq | wc -l | sed 's/ //g' >> mean_num_haplotypes_ATen-AGri.tsv # list of the number of haplotypes for each locus

done
```
## # Calculate mean number of haplotypes
```
a=$(cat mean_num_haplotypes_ATen-AGri.tsv | awk '{sum+=$1}END{print sum/NR}')

printf "ATen-AGri haplotypes mean = "$a
```
### Calculate the mean number of SNPs per loci with the same missing data restrictions ###
```
cat Shearwaters_clust89_min5.vcf | grep -v "^#" | cut -f1,2,5,10-22 | grep -vE "\,.\," | grep -v "./." | cut -f1 | sort | uniq -c | sed -E 's/^ +//g' | tr ' ' '\t' | awk '{sum+=$1}END{print sum/NR}' # print the count
cat Shearwaters_clust89_min5.vcf | grep -v "^#" | cut -f1,2,5,10-22 | grep -vE "\,.\," | grep -v "./." | cut -f1 | sort | uniq -c | sed -E 's/^ +//g' | tr ' ' '\t' | cut -f1 > mean_num_snps_ATen-AGri.tsv # print a file with number of SNPs per locus
```
### Are the numbers obtained for ABBA SNPs significantly different than the average? ### 
```
#!/bin/R

## Calculate p-value of the mean number of SNPs seen in ABBA-BABA loci compared to all loci with the same missing data constraints ##

# Read vector of number of snps per locus
setwd("~/Dropbox/Tesi/Data/Durham/definitive_analyses/dtests/ABBA-BABA_vs_het/")
snps <- scan("mean_num_snps_ATen-AGri.tsv")

# Run permutations sampling the number of ABBA loci from the vector of number of snps per locus and calculate the mean 
t=""
for (i in 1:100000) {
  t[i] <- mean(sample(snps, 20))
}

# Plot a histogram of the means of the permutations
hist(as.numeric(t))

# Calculate Zscores and p-values for each mean
t <- as.numeric(t)
zscore <- as.numeric(round(scale(t),3))
pvalue <- 2*pnorm(-abs(zscore))

# Put mean number of SNPs, zscores and p-values in a data frame
atenagri <- data.frame(value=t, z=zscore, p_value=pvalue)

# Get the p-value for the observed value
head(atenagri$p_value[atenagri$value==13.5],1)

## Calculate p-value of the mean number of haplotypes seen in ABBA-BABA loci compared to all loci with the same missing data constraints ##

# Read vector of number of haplotypes per locus
setwd("~/Dropbox/Tesi/Data/Durham/definitive_analyses/dtests/ABBA-BABA_vs_het/")
haplos <- scan("mean_num_haplotypes_ATen-AGri.tsv")

# Run permutations sampling the number of ABBA loci from the vector of number of haplotypes per locus and calculate the mean
t=""
for (i in 1:100000) {
  t[i] <- mean(sample(haplos, 20))
}

# Plot a histogram of the means of the permutations
hist(as.numeric(t))

# Calculate Zscores and p-values for each mean
t <- as.numeric(t)
zscore <- as.numeric(round(scale(t),3))
pvalue <- 2*pnorm(-abs(zscore))

# Put mean number of haplotypes, zscores and p-values in a data frame
atenagri <- data.frame(value=t, z=zscore, p_value=pvalue)

# Get the p-value for the observed value
head(atenagri$p_value[atenagri$value==15.85],1)

## Calculate permutation test to see if proportion of w_s in ABBA-BABA loci is higher than expected ##

# Run permutation tests
t=""
for (i in 1:100000) {
  t[i] <- sum(sample(c("w_s","no"), 20, replace=T, prob=c(0.43425,0.56575)) == "w_s")
}

# Plot a histogram of the means of the permutations
hist(as.numeric(t),breaks=15)

# Calculate Zscores and p-values for each count
t <- as.numeric(t)
zscore <- as.numeric(round(scale(t),3))
pvalue <- pnorm(-abs(zscore))

# Put count of W-to-S mutations, zscores and p-values in a data frame
atenagri <- data.frame(value=t, z=zscore, p_value=pvalue)

# Get the p-value for the observed value
head(atenagri$p_value[atenagri$value==13],1)
```
##############################################################
## 4) Phylogenetic networks: example with Puffinus analyses ##
##############################################################

### Preparing input files to calculate CF tables ###
#### ddRAD SNPs dataset 
```
echo PLBo2 PMau2 PLBa3 PPuf2 PYel2 PBNi01 PHut8 PLLh1 PNNe2 PGav96 PNat1 PBDi1 PBBa65 POpi2 PEle2 | tr ' ' '\n' > Puffinus.txt 

#create a list of the most complete individual for every species 
cat Shearwaters_clust89_percent75.vcf | sed 's/.assembled//g' | vcftools --vcf - --keep Puffinus.txt --non-ref-ac-any 1 --recode --recode-INFO-all --stdout > Puffinus_clust89_percent75_vcftools.vcf #create a vcf file with only data for the selected individuals and only biallelic sites
loci=($(cat Puffinus_clust89_percent75_vcftools.vcf | grep -v "#" | cut -f1 | uniq)) #create a list of the loci included in the vcf file
cat Puffinus_clust89_percent75_vcftools.vcf | grep "^#" > Puffinus_clust89_percent75_vcftools_1snpxlocus.vcf #create a file where we will store a vcf file with only one SNP per locus
for item in ${loci[@]}; do grep -P "^$item\t" Puffinus_clust89_percent75_vcftools.vcf | shuf | head -1 >> Puffinus_clust89_percent75_vcftools_1snpxlocus.vcf; done #Select a random SNP per locus and write it to the new vcf file
python ../../../../../programari/vcf2phylip.py -i Puffinus_clust89_percent75_vcftools_1snpxlocus.vcf #convert the vcf file to a phylip file
cp  Puffinus_clust89_percent75_vcftools_1snpxlocus.min4.phy ../../../../phylonetworks/definitive_analyses/Puffinus_ddRAD.phy #rename it and this will be the input for the SNP2CF function
```
#### UCE SNPs dataset ##
```
echo POpi2 PBBa65 PNNeD1 PLLo1 PBNiD1 PBDiD1 PGav96 PLBa2 PHut8 PYel1 PPufD1 PMau2 PNatD1 PLBo2 PEle2 | tr ' ' '\n' > Puffinus.txt 

#create a list of the most complete individual for every species

for file in *vcf; do cat $file | vcftools --vcf - --keep Puffinus.txt --non-ref-ac-any 1 --recode --recode-INFO-all --stdout > ../dataset_out_phased_iupac_trimmed_cleaned_selected_75_vcf_Puffinus/$file; done #create one vcf file per locus containing only the selected species
for file in *vcf; do cat $file | grep "^#" > ../dataset_out_phased_iupac_trimmed_cleaned_selected_75_vcf_Puffinus_1SNP/$file; cat $file | grep -v "^#" | shuf | head -1 >> ../dataset_out_phased_iupac_trimmed_cleaned_selected_75_vcf_Puffinus_1SNP/$file; done #from the last vcf select only one SNP per locus and save it to a new vcf per locus
for file in *vcf; do #convert the single SNP vcf to fasta
	zero=$(cat $file | grep "uce" | cut -f4)
  	alt=$(cat $file | grep "uce" | cut -f5)
  	two=$(echo $alt | grep ",")
  	if [ -z "$two" ]; then one=$alt; else one=$(echo $alt | cut -d "," -f1); two=$(echo $alt | cut -d "," -f2); fi  
  	sp=($(cat $file | grep "CHROM" | cut -f10- | tr '\t' '\n'))
  	value=($(cat $file | grep "uce" | cut -f10- | tr '\t' '\n'))
  	for ((i=0;i<${#sp[@]};i++)); do 
  		if [ ${value[i]}  -eq 0 ]; then printf \>${sp[i]}'\n'$zero'\n'>> $path/dataset_out_phased_iupac_trimmed_cleaned_selected_75_fasta_Puffinus_1SNP/$(echo $file | sed 's/.vcf/.fa/')
    		elif [ ${value[i]}  -eq 1 ]; then printf \>${sp[i]}'\n'$one'\n'>> $path/dataset_out_phased_iupac_trimmed_cleaned_selected_75_fasta_Puffinus_1SNP/$(echo $file | sed 's/.vcf/.fa/')
    		else printf \>${sp[i]}'\n'$two'\n'>> $path/dataset_out_phased_iupac_trimmed_cleaned_selected_75_fasta_Puffinus_1SNP/$(echo $file | sed 's/.vcf/.fa/')
    		fi
  	done
done
```

```
TriSeq -in dataset_out_phased_iupac_trimmed_cleaned_selected_75_fasta_Puffinus_1SNP/*fa -o UCE_75_phased_iupac_Puffinus -of phylip #Concatenate all the previous fasta into a concatenated phylip file
mv UCE_75_phased_iupac_Puffinus.phy Puffinus_UCE.phy #rename the phylip file to input it to SNP2CF
```

### ddRAD gene trees ##
```
nw_prune -v ml_best_supports_ddRAD.trees PLBo2 PMau2 PLBa3 PPuf2 PYel2 PBNi01 PHut8 PLLh1 PNNe2 PGav96 PNat1 PBDi1 PBBa65 POpi2 PEle2 > ml_best_supports_ddRAD_pruned.trees #prune the trees to only include the selected samples
nw_ed ml_best_supports_ddRAD_pruned.trees 'i & b<=10' o > ml_best_supports_ddRAD_pruned_BS10.trees #collapse nodes with BS <10
nw_ed ml_best_supports_ddRAD_pruned.trees 'i & b<=50' o > ml_best_supports_ddRAD_pruned_BS50.trees #collapse nodes with BS < 50
```

### UCE gene trees ##
```
nw_prune -v ml_best_supports_UCE.trees POpi2 PBBa65 PNNeD1 PLLo1 PBNiD1 PBDiD1 PGav96 PLBa2 PHut8 PYel1 PPufD1 PMau2 PNatD1 PLBo2 PEle2 > ml_best_supports_UCE_pruned.trees #prune the trees to only include the selected samples
nw_ed ml_best_supports_UCE_pruned.trees 'i & b<=10' o > ml_best_supports_UCE_pruned_BS10.trees #collapse nodes with BS < 10
nw_ed ml_best_supports_UCE_pruned.trees 'i & b<=50' o > ml_best_supports_UCE_pruned_BS50.trees #collapse nodes wit BS < 50
```

### Calculate CF tables ###

##### SNPs
```
#!/bin/R

setwd("/Users/apple/Dropbox/Tesi/Data/Durham/PhyloNetworks/definitive_analyses")
Puffinus_ddRAD <- SNPs2CF(seqMatrix = "Puffinus_ddRAD.phy", outputName = "Puffinus_ddRAD_SNPs2CF.csv", save.progress = T, cores = 4) #ddRAD data
Puffinus_UCE <- SNPs2CF(seqMatrix = "Puffinus_UCE.phy", outputName = "Puffinus_UCE_SNPs2CF.csv", save.progress = T, cores = 4) #UCE data

## Gene trees: example ddRAD data ##

#!/bin/env/julia

using PhyloNetworks
cd("/Users/apple/Dropbox/Tesi/Data/Durham/PhyloNetworks/definitive_analyses")
ddRADtrees = ("ml_best_supports_ddRAD_pruned.trees")
ddRADgenetrees = readMultiTopology(ddRADtrees) #read trees
ddRADCFind = readTrees2CF(ddRADtrees, CFfile="Puffinus_ddRAD_genetreesCF.csv"); #calculate CF table
ddRADtrees = ("ml_best_supports_ddRAD_pruned_BS10.trees")
ddRADgenetrees = readMultiTopology(ddRADtrees)
ddRADCFind = readTrees2CF(ddRADtrees, CFfile="Puffinus_ddRAD_BS10_genetreesCF.csv");
ddRADtrees = ("ml_best_supports_ddRAD_pruned_BS50.trees")
ddRADgenetrees = readMultiTopology(ddRADtrees)
ddRADCFind = readTrees2CF(ddRADtrees, CFfile="Puffinus_ddRAD_BS50_genetreesCF.csv");

### Infer hmax=0 coalescent trees: example ddRAD BS10 data ###

using Distributed
addprocs(10)
@everywhere using PhyloNetworks
cd("/users/jferrer/Durham/phylonetworks/definitive_analyses")
dastralfile = ("starting_tree_astral.tre")
dastraltree = readTopology(dastralfile)
using CSV
df_sp_dgBS10 = CSV.read("Puffinus_ddRAD_BS10_genetreesCF.csv", categorical=false);
d_sp_dgBS10 = readTableCF!(df_sp_dgBS10);
net0_dgBS10 = snaq!(dastraltree, d_sp_dgBS10, hmax=0, filename="net0_dgBS10", seed=1234)

### I do the same for hmax=1,2,3 and 4 changing the starting tree to be the hmax-1 optimized topology and changing hmax in snaq! ###

### TICR test: example ddRAD BS10 data ###

#!/bin/R

library(phylolm)

#Read and explore CF data and tree
setwd("~/Dropbox/Tesi/Data/Durham/PhyloNetworks/definitive_analyses")
quartetCF = read.csv("Puffinus_ddRAD_BS10_genetreesCF.csv")
dim(quartetCF)
head(quartetCF)
dat = quartetCF[, c(1:7)]
for (i in 1:4){ dat[,i] = factor(dat[,i])}
head(dat)
tree = read.tree("net0_dgBS10.nwk")
tree
plot(tree)
edgelabels(round(tree$edge.length,3))

#Run test
prelim = test.tree.preparation(dat, tree)
Ntaxa = length(tree$tip.label)
internal.edges = which(tree$edge[,2] > Ntaxa+1)
internal.edges = c(3,  4,  5,  6,  9, 10, 15, 16, 17, 21, 23, 26)
res <- test.one.species.tree(dat,tree,prelim,edge.keep=internal.edges)
res[1:6]

#Get outlier taxa
outlier.4taxa.01 <- which(res$outlier.pvalues < 0.01)
length(outlier.4taxa.01)
q01 = as.matrix(quartetCF[outlier.4taxa.01,1:4])
q01
sort(table(as.vector(q01)),decreasing=TRUE)

dgBS10_01 <- cbind(dat[outlier.4taxa.01,],res$cf.exp[outlier.4taxa.01,])
setwd("~/Dropbox/Tesi/Data/Durham/PhyloNetworks/definitive_analyses/TICR_test/")
write.csv(dgBS10_01, file = "dgBS10_01_outliers.csv")

outlier.4taxa.05 <- which(res$outlier.pvalues < 0.05)
length(outlier.4taxa.05)
q05 = as.matrix(quartetCF[outlier.4taxa.05,1:4])
head(q05)
sort(table(as.vector(q05)),decreasing=TRUE)

### Optimise candidate networks: example BS10 data ###

## Read CF table and network ##
cd("/Users/apple/Dropbox/Tesi/Data/Durham/PhyloNetworks/definitive_analyses/")
df_sp = CSV.read("Puffinus_ddRAD_BS10_genetreesCF.csv", categorical=false);
d_sp = readTableCF!(df_sp);
net1file = ("net1.nwk")
net1 = readTopology(networkfile)

## Optimise every network with the CF values from the table: example network lherminieri -> boydi (net1) ##

net1alt = topologyMaxQPseudolik!(net1, d_sp);

## Calculate its pseudo-deviance ##

net1alt.loglik

## Reroot, save and plot to see the inheritance probabilities ##

rootatnode!(net1alt, "P_nativitatis")
writeTopology(net1alt, "net1alt.nwk")
plot(net1alt, :R, showGamma=true);
```
