# Gene conversion

#############################################
## 1) Calculate GC content per ddRAD locus ##
#############################################
```
python /users/jferrer/programari/AMAS/amas/AMAS.py summary -i out_fasta/*fa -f fasta -d dna -o summary_ddRAD_75.txt -c 20

cat summary_ddRAD_75.txt | cut -f1,3,11,12 | sed 's/Alignment_name/CHROM/' | sed -E 's/^locus_([0-9]+)\.fa/\1/g' > GC_content_locus.txt
```
###############################################################################                                                      
## 2) Calculate reference and derived allele frequencies for biallelic sites ##
###############################################################################
```
vcftools --vcf populations.snps.vcf --min-alleles 2 --max-alleles 2 --freq --out populations.snps

cat populations.snps.frq | cut -f1,2,5,6 | sed -E 's/POS	.*/POS	ref	ref_freq	alt	alt_freq/' | sed 's/:/	/g' > snps_freq.txt
```
##############################################################################
## 3) Summary plots to evaluate the signatures of GC-biased gene conversion ##
##############################################################################
```
#!/bin/R

### Load libraries and set working directory ###
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(grid)
library(gridExtra)
library(readxl)

setwd("~/Dropbox/Tesi/Data/Durham/gene_conversion/")

### Read input files and modify data frames to plot the data ###
snps <- read.delim("snps_freq.txt")
snps$mut_type <- as.factor(paste(snps$ref, snps$alt, sep="_"))
gc_content <- read.delim("GC_content_locus.txt")

gc <- full_join(snps, gc_content, by= "CHROM")
gc <- gc[!is.na(gc$ref),]

gc_final <- as.data.frame(gc %>%
  group_by(CHROM, ref, alt, mut_type) %>% 
  summarise_each(funs(mean)))

gc_final <- gc_final %>% mutate(
  mut_cat = as.factor(case_when(mut_type == "A_C" ~ "w_s",
            mut_type == "A_G" ~ "w_s",
            mut_type == "A_T" ~ "w_w",
            mut_type == "C_A" ~ "s_w",
            mut_type == "C_G" ~ "s_s",
            mut_type == "C_T" ~ "s_w",
            mut_type == "G_A" ~ "s_w",
            mut_type == "G_C" ~ "s_s",
            mut_type == "G_T" ~ "s_w",
            mut_type == "T_A" ~ "w_w",
            mut_type == "T_C" ~ "w_s",
            mut_type == "T_G" ~ "w_s",)))

### Plot mutation category vs GC content with size denoting alt. allele freq. ###
gc_plot_int <- gc_final %>% mutate(GC_int=cut(GC_content, breaks=seq(0,0.76,0.02), labels=seq(0.02,0.76,0.02)))

gc_plot_1 <- gc_plot_int %>% group_by(GC_int) %>% count(mut_cat) %>% 
          mutate(mut_cat_freq = n/sum(n)) %>% rename(mut_cat_count=n)
gc_plot_2 <- gc_plot_int %>% group_by(GC_int, mut_cat) %>%
  summarise(alt_freq_mean = mean(alt_freq), alt_freq_n = n())

gc_plot <- full_join(gc_plot_1, gc_plot_2)
gc_plot$GC_int <- as.numeric(levels(gc_plot$GC_int))[gc_plot$GC_int]
gc_plot <- gc_plot[gc_plot$GC_int > 0.29,]
gc_plot <- gc_plot[gc_plot$GC_int < 0.67,]

gc_plot$mut_cat <- factor(gc_plot$mut_cat, levels = c("s_w", "w_s", "s_s", "w_w"))

mut_cat_vs_gccontent <- ggplot(gc_plot, aes(x=GC_int, y=mut_cat_freq, 
                        color=mut_cat, group=mut_cat)) +
  theme_bw() +
  geom_point(aes(size=alt_freq_mean)) +
  geom_line() +
  scale_size(range = c(1,5)) +
  guides(size=F) +
  scale_color_manual(values=c("#481567","#3CBB75","#607D8B","#9E9E9E"),
                     labels=c("S to W","W to S","S to S","W to W")) +
  labs(x="GC content", y="Proportion", color="")

ggsave("mut_cat_vs_GCcont_vs_altfreqmean.pdf", mut_cat_vs_gccontent, device="pdf", units="cm", width=20, height=11, limitsize=FALSE)

### Plot mutation category vs alt. allele freq ###
gc_plot_int <- gc_final %>% mutate(alt_freq_int=cut(alt_freq, breaks=seq(0.04,0.5,0.01), labels=seq(0.05,0.5,0.01)))
gc_plot_int <- gc_plot_int[!is.na(gc_plot_int$alt_freq_int),]
  
gc_plot_1 <- gc_plot_int %>% group_by(alt_freq_int) %>% count(mut_cat) %>% 
  mutate(mut_cat_freq = n/sum(n)) %>% rename(mut_cat_count=n)
gc_plot_2 <- gc_plot_int %>% group_by(alt_freq_int, mut_cat) %>%
  summarise(GC_mean = mean(GC_content), GC_n = n())

gc_plot <- full_join(gc_plot_1, gc_plot_2)
gc_plot$alt_freq_int <- as.numeric(levels(gc_plot$alt_freq_int))[gc_plot$alt_freq_int]

gc_plot$mut_cat <- factor(gc_plot$mut_cat, levels = c("s_w", "w_s", "s_s", "w_w"))

mut_cat_vs_alt_freq <- ggplot(gc_plot, aes(x=alt_freq_int, y=mut_cat_freq, 
                                            color=mut_cat, group=mut_cat)) +
  theme_bw() +
  geom_point() +
  geom_line() +
  scale_color_manual(values=c("#481567","#3CBB75","#607D8B","#9E9E9E"),
                     labels=c("S to W","W to S","S to S","W to W")) +
  labs(x="Frequency alternative allele", y="Proportion", color="") +
  xlim(0.08,0.29)

ggsave("mut_cat_vs_altfreqmean.pdf", mut_cat_vs_alt_freq, device="pdf", units="cm", width=20, height=11, limitsize=FALSE)
```
#####################################################################################
## 4) Prepare input file with genotype counts per mutation type and per individual ##
#####################################################################################
```
printf sp'\t'ref'\t'alt'\t'ref-ref'\t'ref-alt'\t'alt-alt'\n' > freq_mut_types_ind.tsv
nuc=(A C T G)
for s in $(seq 10 60); do
	sp=$(cat *snps.vcf | grep "^#" | grep -v "^##" | cut -f$s)
	for n in ${nuc[@]}; do
		for n2 in ${nuc[@]}; do
			if [ $n2 != $n ]; then
				freq=$(cat *snps.vcf | grep -v "^#" | cut -f4,5,$s | cut -d ":" -f1 |\
  				grep "^$n	$n2" | grep -v "\./\." | cut -f3 |\
   				tr "/" "\t" | awk '!(NR%2){print $1+$2}' | sort | uniq -c |\
    				sed -E 's/^ +//' | cut -d " " -f1 | tr '\n' '_' | sed -E 's/	$//' | sed -E 's/_$//')
			fi
		printf $sp'\t'$n'\t'$n2'\t'$freq'\n' >> freq_mut_types_ind.tsv
    	done
    done
done
cat freq_mut_types_ind.tsv | tr '_' '\t' > freq_mut_types_ind_bo.tsv
rm freq_mut_types_ind.tsv
mv freq_mut_types_ind_bo.tsv freq_mut_types_ind.tsv
```
####################################################################################
## 4) Plot the W-to-S mutation proportion against life-history traits per species ##
####################################################################################
```
#!/bin/R

### Load libraries ###
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(grid)
library(gridExtra)
library(readxl)

### Read input file and modify data frame to calculate the proportion of each mutation category ###
mut_ind <- read.delim("freq_mut_types_ind.tsv")
mut_ind$mut_type <- paste(mut_ind$ref, mut_ind$alt, sep="_")

## Define mutation categories ##
mut_ind <- mut_ind %>% mutate(
  mut_cat = as.factor(case_when(mut_type == "A_C" ~ "w_s",
                                mut_type == "A_G" ~ "w_s",
                                mut_type == "A_T" ~ "w_w",
                                mut_type == "C_A" ~ "s_w",
                                mut_type == "C_G" ~ "s_s",
                                mut_type == "C_T" ~ "s_w",
                                mut_type == "G_A" ~ "s_w",
                                mut_type == "G_C" ~ "s_s",
                                mut_type == "G_T" ~ "s_w",
                                mut_type == "T_A" ~ "w_w",
                                mut_type == "T_C" ~ "w_s",
                                mut_type == "T_G" ~ "w_s",)))

## Calculate allele frequencies from genotype frequencies ##
mut_ind$ref_freq <- (2*mut_ind$ref.ref + mut_ind$ref.alt) / (2*mut_ind$ref.ref + 2*mut_ind$ref.alt + 2*mut_ind$alt.alt)
mut_ind$alt_freq <- (2*mut_ind$alt.alt + mut_ind$ref.alt) / (2*mut_ind$ref.ref + 2*mut_ind$ref.alt + 2*mut_ind$alt.alt)

## Calculate alternative allele counts ##
mut_ind$mut_freq <- 2*mut_ind$alt.alt + mut_ind$ref.alt

## Calculate proportions of W-to-S mutation and of S alleles in W-to-S mutations ##
a=""
b=""
c=""
d=""
e=""
f=""
g=""
h=""
k=""

for (i in 1:length(levels(mut_ind$sp))) {
  a[i] = levels(mut_ind$sp)[i]
  b[i]=sum(mut_ind[mut_ind$sp==a[i],11])
  c[i]=sum(mut_ind[mut_ind$sp==a[i] & mut_ind$mut_cat == "w_s",11])
  d[i]=as.numeric(c[i])/as.numeric(b[i])
  
  e[i]=2*sum(mut_ind[mut_ind$sp==a[i] & mut_ind$mut_cat == "w_s",6])
  f[i]=sum(mut_ind[mut_ind$sp==a[i] & mut_ind$mut_cat == "w_s",5])
  g[i]=as.numeric(e[i]) + as.numeric(f[i])
  h[i]=2*sum(mut_ind[mut_ind$sp==a[i] & mut_ind$mut_cat == "w_s",4]) + 2*sum(mut_ind[mut_ind$sp==a[i] & mut_ind$mut_cat == "w_s",5]) + 2*sum(mut_ind[mut_ind$sp==a[i] & mut_ind$mut_cat == "w_s",6])
  k[i]=as.numeric(g[i])/as.numeric(h[i])
  }

freqs = data.frame(sp=as.factor(a), ws_freq=as.numeric(d), s_in_ws=as.numeric(k))

### Read life-history data and merge with W-to-S mutation proportions data ###
meta <- read_excel("../Corr_Het_others/corr_het_others.xlsx")
meta <- meta %>% rename(sp = Id)
meta$sp <- as.factor(meta$sp)
corr <- full_join(freqs, meta, by="sp")
corr <- corr[corr$sp!="APacLMG" & corr$sp!="PAHa2" & corr$sp!="ABul1",]
genus <- c(rep("A",9),"P","A","A",rep("C",8),rep("P",28))
corr$genus <- as.factor(genus)

### Plot the W-to-S mutation proportion against life-history traits per species ###
## log(N pairs) ##
m = lm(log10(Npairs_pop) ~ ws_freq, corr)
eq1 = paste0("rÂ² = ", format(summary(m)$r.squared, digits = 3))
eq2 = paste0("p-value = ", format(summary(m)$coefficients[2,4], digits = 3))
eq = c(eq1, eq2)
eq = paste(eq, collapse="\n")
lognpairs_wsfreq <- ggplot(corr, aes(x=Npairs_pop, y=ws_freq)) +
  theme_bw() +
  geom_point(aes(color=genus), size = 1.5) +
  scale_colour_manual(values=c("#328fa0","#5b2679","#feea31"), guide=F) +
  labs(x = "N of breeding pairs", y = "W-to-S mutations proportion") +
  geom_smooth(method = "lm", se = F, color="#424242") +
  scale_x_continuous(trans = 'log10') +
  annotate("text", x = 7000, y = 0.43, label = eq, color="black", size = 4, parse = FALSE)

## Body mass ##
m = lm(Body_mass_hbw ~ ws_freq, corr)
eq1 = paste0("rÂ² = ", format(summary(m)$r.squared, digits = 3))
eq2 = paste0("p-value = ", format(summary(m)$coefficients[2,4], digits = 3))
eq = c(eq1, eq2)
eq = paste(eq, collapse="\n")
bm_wsfreq <- ggplot(corr, aes(x=Body_mass_hbw, y=ws_freq)) +
  theme_bw() +
  geom_point(aes(color=genus), size = 1.5) +
  scale_colour_manual(values=c("#328fa0","#5b2679","#feea31"), guide=F) +
  labs(x = "Body mass (g)", y = "") +
  geom_smooth(method = "lm", se = F, color="#424242") +
  annotate("text", x = 250, y = 0.43, label = eq, color="black", size = 4, parse = FALSE)

lifehist_ws <- ggarrange(lognpairs_wsfreq, bm_wsfreq, ncol=2, nrow=1, labels = c("a)", "b)"))
ggsave("life_hist_vs_WtoS_mutations.pdf", lifehist_ws, device="pdf", units="cm", width=30, height=13, limitsize=FALSE)
```
