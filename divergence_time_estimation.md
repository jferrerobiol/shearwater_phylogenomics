##################################################################
## 1) Prepare the 4 subsets for clock model selection in mcmc3r ##
##################################################################
```
sel=$(ls | shuf | head -60 | tr '\n' ' ') # select 60 UCE loci randomly from the 95% dataset
TriSeq -in $sel -grep ACar1 ACre1 AGra1 AGri1 ATen1 ABul1 APac1 CLeu1 CEdw1 CDio1 CBor1 PNat1 PHut1 PGav1 PEle1 POpi1 PNNe1 PBDi1 PBNi1 PBBa1 PPuf1 PMau1 PYel1 PLLh1 PLBo1 PLBa1 -of phylip -o ../../mcmctree/UCE_95_sel2.mcmctree
sel=$(ls | shuf | head -120 | tr '\n' ' ') # select 120 ddRAD loci randomly from the 95% dataset
TriSeq -in $sel -grep ACar1 ACre1 AGra1 AGri1 ATen1 ABul1 APac1 CLeu1 CEdw1 CDio1 CBor1 PNat1 PHut1 PGav1 PEle1 POpi1 PNNe1 PBDi1 PBNi1 PBBa1 PPuf1 PMau1 PYel1 PLLh1 PLBo1 PLBa1 -of phylip -o ../../mcmctree/ddRAD_95_sel2.mcmctree
```

#####################################################
## 2) Prepare the input tree for MCMCtree analyses ##
#####################################################

### Prune the Exabayes tree ###
```
#!/bin/R

## Loading libraries
library(ape)
library(phangorn)
library(phytools)

## Read and prune tree so as it contains a single tip for each taxa
setwd("~/Dropbox/Tesi/Data/Durham/definitive_analyses/exabayes/")
t <- read.tree("ExaBayes_ConsensusExtendedMajorityRuleNewick.concat_ddRAD-UCE_75_2part")
t <- root(t,"FGla")
t <- drop.tip(t,c("FGla","CLeu2","CEdw2","CDio2","CBor2","APac2","ATen2","AGri2","AGra2","ACar2","ACre2","PNat2","PGav2","PNNe2","PBNi2","PBBa2","PBBa3","PPuf2","PMau2","PYel2","PLLh3","PLLh2","PLBo2","PLBa2","PLBa3"))

## Remove branch lengths and supports
t$edge.length<-NULL # get rid of edge lengths
t$node.label<-NULL # get rid of branch supports
write.tree(t,file = "mcmctree.tre") # write tree without edge lengths
```
### Add a first line containing the number of species and the number of trees (26 1) ###

### Fix the root adding a root calibration for mcmc3r analyses: '>0.999,<1.001' ###
```
26 1
(((((((ACar1,ACre1),AGra1),AGri1),ATen1),(ABul1,APac1)),(CLeu1,(CEdw1,(CDio1,CBor1)))),(PNat1,((PHut1,PGav1),((PEle1,((POpi1,PNNe1),((PBDi1,PBNi1),PBBa1))),((PPuf1,(PMau1,PYel1)),(PLLh1,(PLBo1,PLBa1)))))))'>.999<1.001';
```

## The other trees that will be used will be modifications of this tree specifying different node calibrations ##

###############################################
## 3) mcmc3r analyses to chose a clock model ##
###############################################

### Estimate the parameters for the birth and death priors ###

## Calculate birth and death rates for the whole tree
```
calib <- makeChronosCalib(t, node = "root", age.min = 22.9, age.max = 23.1, interactive = FALSE, soft.bounds = FALSE)
chronos.control()
t2 <- chronos(t, lambda = 1, model = "correlated", quiet = FALSE,
        calibration = calib,
        control = chronos.control())
fit.bd<-birthdeath(t2)
fit.bd
bd(fit.bd)
```
### For each data subset I followed the tutorial Estimating the marginal likelihood of a relaxed-clock model with MCMCTree from dos Reis lab:
```
https://dosreislab.github.io/2017/10/24/marginal-likelihood-mcmc3r.html
```
## An example of an MCMCTree control file used in the analysis ##
```
seed = -1
seqfile = ./UCE_95_sel1.mcmctree.phy
treefile = ./mcmc3r.tre
outfile = out

ndata = 1      * number of partitions
usedata = 1    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
clock = 1      * 1: global clock; 2: independent rates; 3: correlated rates

model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
alpha = 0.02295   * alpha for gamma rates at sites; estimated from Exabayes UCE95_nopart analyses
ncatG = 5    * No. categories in discrete gamma

BDparas = 0.172 0.083 0.75      * birth, death, sampling;
kappa_gamma = 6 2   * gamma prior for kappa
alpha_gamma = 1 1    * gamma prior for alpha

rgene_gamma = 2 20   * gamma prior for mean rates for genes
sigma2_gamma = 1 10  * gamma prior for sigma^2 (for clock=2 or 3)

print = 1
burnin = 4000
sampfreq = 6
nsample = 20000
```
## The results for the best models are shown in Supplementary Table A3 ##

#################################################################
## 4) Prepare the alignments for dating analyses with MCMCTree ##
#################################################################

### Remove individuals keeping the most complete per taxon ###
```
TriSeq -c -in ddRAD_95.phy -grep ACar1 ACre1 AGra1 AGri1 ATen1 ABul1 APac1 CLeu1 CEdw1 CDio1 CBor1 PNat1 PHut1 PGav1 PEle1 POpi1 PNNe1 PBDi1 PBNi1 PBBa1 PPuf1 PMau1 PYel1 PLLh1 PLBo1 PLBa1 -of mcmctree -o ../mcmctree/ddRAD_95_mcmctree_tot.phy ## all ddRAD 95%
TriSeq -c -in UCE_95.phy -grep ACar1 ACre1 AGra1 AGri1 ATen1 ABul1 APac1 CLeu1 CEdw1 CDio1 CBor1 PNat1 PHut1 PGav1 PEle1 POpi1 PNNe1 PBDi1 PBNi1 PBBa1 PPuf1 PMau1 PYel1 PLLh1 PLBo1 PLBa1 -of mcmctree -o ../mcmctree/UCE_95_mcmctree_tot.phy ## all UCE 95%
TriSeq -in UCE_95_mcmctree_tot.phy_mcmctree.phy ddRAD_95_mcmctree_tot.phy_mcmctree.phy -of mcmctree -o concat_95_mcmctree_tot.phy ## concat 95% with 2 partitions (UCE, ddRAD)
```
#######################################################################################
## 5) Calculate a rough substitution rate in baseml to inform the rgene_gamma prior ##
#######################################################################################

### I run baseml for UCE and ddRAD datasets with the following tree file and control file ###

## Tree file with root fixed at 23.03 Ma ##
```
26 1
(((((((ACar1,ACre1),AGra1),AGri1),ATen1),(ABul1,APac1)),(CLeu1,(CEdw1,(CDio1,CBor1)))),(PNat1,((PHut1,PGav1),((PEle1,((POpi1,PNNe1),((PBDi1,PBNi1),PBBa1))),((PPuf1,(PMau1,PYel1)),(PLLh1,(PLBo1,PLBa1)))))))'@.2303';
```
## Example control file ##
```
      seqfile = ../../ddRAD_95_tot.mcmctree.phy 
     treefile = ./baseml.tre

      outfile = ./baseml.out       * main result file
        noisy = 3   * 0,1,2,3: how much rubbish on the screen
      verbose = 1   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
        model = 7   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
        Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff

        clock = 1   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
    fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below; 2: kappa for branches
        kappa = 2  * initial or fixed kappa

    fix_alpha = 0   * 0: estimate alpha; 1: fix alpha at value below
        alpha = 0.5   * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * 1: different alpha's for genes, 0: one alpha
        ncatG = 5   * # of categories in the dG, AdG, or nparK models of rates
        nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK 

      fix_rho = 1  
          rho = 0. * initial or given rho, 0:no correlation
        nhomo = 1   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0   * (0,1,2): rates (alpha>0) or ancestral states
    cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?
  fix_blength = 0
```
#################################################################################
## 6) MCMCTree analyses with different datasets, root calibrations and  bounds ##
#################################################################################

### For both, 4 calibrations and 3 calibrations we run the following bash script to run all the MCMCTree analyses at once ###
```
cd /ddn/data/sbvd77/mcmctree/dating_analyses
list=(UCE_95_tot_3calib ddRAD_95_tot_3calib concat_95_tot_3calib)

for case in ${list[@]}; do

  mkdir $case

  cd $case

  mkdir approx_$case\_diff_calib
  mkdir approx_$case\_diff_calib/approx_only_low_bound_calib_23 approx_$case\_diff_calib/approx_only_low_bound_calib_38 approx_$case\_diff_calib/approx_rootmax_23 approx_$case\_diff_calib/approx_rootmax_38
  mkdir approx_$case\_diff_calib/approx_only_low_bound_calib_23/approx approx_$case\_diff_calib/approx_only_low_bound_calib_23/hessian approx_$case\_diff_calib/approx_only_low_bound_calib_23/prior approx_$case\_diff_calib/approx_only_low_bound_calib_38/approx approx_$case\_diff_calib/approx_only_low_bound_calib_38/hessian approx_$case\_diff_calib/approx_only_low_bound_calib_38/prior approx_$case\_diff_calib/approx_rootmax_23/approx approx_$case\_diff_calib/approx_rootmax_23/hessian approx_$case\_diff_calib/approx_rootmax_23/prior approx_$case\_diff_calib/approx_rootmax_38/approx approx_$case\_diff_calib/approx_rootmax_38/hessian approx_$case\_diff_calib/approx_rootmax_38/prior

  cp ../template.ctl approx_$case\_diff_calib/approx_only_low_bound_calib_23/hessian/hessian.ctl; sed -i "s/UCE_95_sel1/$case/" approx_$case\_diff_calib/approx_only_low_bound_calib_23/hessian/hessian.ctl; sed -i 's/mcmctree.tre/mcmctree\_low\_23\_3calib.tre/' approx_$case\_diff_calib/approx_only_low_bound_calib_23/hessian/hessian.ctl; sed -i 's/.38/.23/' approx_$case\_diff_calib/approx_only_low_bound_calib_23/hessian/hessian.ctl
  cp ../template.ctl approx_$case\_diff_calib/approx_only_low_bound_calib_38/hessian/hessian.ctl; sed -i "s/UCE_95_sel1/$case/" approx_$case\_diff_calib/approx_only_low_bound_calib_38/hessian/hessian.ctl; sed -i 's/mcmctree.tre/mcmctree\_low\_38\_3calib.tre/' approx_$case\_diff_calib/approx_only_low_bound_calib_38/hessian/hessian.ctl
  cp ../template.ctl approx_$case\_diff_calib/approx_rootmax_23/hessian/hessian.ctl; sed -i "s/UCE_95_sel1/$case/" approx_$case\_diff_calib/approx_rootmax_23/hessian/hessian.ctl; sed -i 's/mcmctree.tre/mcmctree\_23\_3calib.tre/' approx_$case\_diff_calib/approx_rootmax_23/hessian/hessian.ctl; sed -i 's/.38/.23/' approx_$case\_diff_calib/approx_rootmax_23/hessian/hessian.ctl
  cp ../template.ctl approx_$case\_diff_calib/approx_rootmax_38/hessian/hessian.ctl; sed -i "s/UCE_95_sel1/$case/" approx_$case\_diff_calib/approx_rootmax_38/hessian/hessian.ctl; sed -i 's/mcmctree.tre/mcmctree\_38\_3calib.tre/' approx_$case\_diff_calib/approx_rootmax_38/hessian/hessian.ctl

  dire=($(find ./ -name "hessian" -type d))

  for i in ${dire[@]}; do

    printf "cd "$i"; mcmctree hessian.ctl; cp out.BV ../approx/in.BV; cp hessian.ctl ../approx/approx.ctl; cp hessian.ctl ../prior/prior.ctl; cd ../prior; sed -i 's/usedata\ \=\ 3/usedata\ \=\ 0/' prior.ctl; sed -i 's/hessian.txt/prior.txt/' prior.ctl; mcmctree prior.ctl; cd ../approx; sed -i 's/usedata\ \=\ 3/usedata\ \=\ 2/' approx.ctl; sed -i 's/hessian.txt/approx.txt/' approx.ctl; mcmctree approx.ctl"'\n' >> parallel_scripts

  done

  parallel < parallel_scripts

  cd ..

done
```
### Script to collect the results from all the runs ###
```
cd /ddn/data/sbvd77/mcmctree/dating_analyses/3calib

tests=($(find -name "*txt" | grep -vE "hessian|mcmc|tests"))

for i in ${tests[@]}; do

  var=$(echo $i | sed -E 's/^\.\///' | tr '/' '-' | sed -E 's/(\-[a-z]+\.txt)$//')
  cat $i | tail -40 | grep "t_n" | cut -d "_" -f2 | sed -E 's/\ +/	/g' | cut -f1,2,5,6 | sed -E 's/\((.*)\,(.*)\)/\1\2/' | sed -E "s/^/$var    /g" >> results_tests_definitive.tsv

done
```
### Example MCMCTree control file from the UCE analysis ###
```
seed = -1
seqfile = /ddn/data/sbvd77/mcmctree/UCE_95_tot_3calib.mcmctree.phy
treefile = /ddn/data/sbvd77/mcmctree/dating_analyses/mcmctree_23_3calib.tre
outfile = ./approx.txt
mcmcfile = ./mcmc.txt

ndata = 1      * number of partitions
seqtype = 0    * 0: nucleotides; 1:codons; 2:AAs
usedata = 2    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
clock = 2      * 1: global clock; 2: independent rates; 3: correlated rates
RootAge = '>0.152<0.23'  * safe constraint on root age, used if no fossil for root.

model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
alpha = 0.02295   * alpha for gamma rates at sites; estimated from Exabayes UCE95_nopart analyses
ncatG = 5    * No. categories in discrete gamma

cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

BDparas = 0.172 0.083 0.75      * birth, death, sampling; birth and death calculated with R ape; sampling =24/32 which is the coverage of shearwater species we have in the dataset
kappa_gamma = 6 2   * gamma prior for kappa
alpha_gamma = 1 1    * gamma prior for alpha

rgene_gamma = 2 376 1  * gammaDir prior for rate for genes  2/200 = 0.00532 overall substitution rate
sigma2_gamma = 1 10 1  * gamma prior for sigma^2 (for clock=2 or 3)

finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1): times, musigma2, rates, mixing, paras, FossilErr

print = 1
burnin = 5000
sampfreq = 500
nsample = 10000
```
######################################################
## 7) Plots for divergence time estimation analyses ##
######################################################
```
#!/bin/R

### Tree with rootmax_23 parameterization and posterior distributions plotted ###

t23 <- readMCMCtree("~/Dropbox/Tesi/Data/Durham/definitive_analyses/mcmctree/definitive/concat_95_tot/approx_concat_95_tot_diff_calib/approx_rootmax_23/approx/FigTree.tre")
t23_distrib <- read.table("~/Dropbox/Tesi/Data/Durham/definitive_analyses/mcmctree/definitive/concat_95_tot/approx_concat_95_tot_diff_calib/approx_rootmax_23/approx/calib_23_posterior.tsv", header=T, sep="\t") 
t23$apePhy$tip.label <- c("A.carneipes","A.creatopus","A.gravis","A.grisea","A.tenuirostris","A.bulleri","A.pacifica",
                          "C.leucomelas","C.edwardsii","C.diomedea","C.borealis","P.nativitatis","P.huttoni","P.gavia",
                          "P.elegans","P.opisthomelas","P.newelli","P.bailloni_dichrous","P.bailloni_nicolae",
                          "P.bailloni_bailloni","P.puffinus","P.mauretanicus","P.yelkouan","P.lherminieri",
                          "P.boydi","P.baroli")

pdf("tree_rootmax23.pdf", width=7, height=7)
MCMC.tree.plot(t23, MCMC.chain=t23_distrib, cex.tips=0.7, time.correction=100, scale.res=c("Period", "Epoch"),
               plot.type="distributions", cex.age=0.7, cex.labels=0.7, relative.height=0.08, col.tree="#607d8b", 
               label.offset=0.1, node.method="bar", no.margin=TRUE, lwd.bar = 2, density.col = "#3f51b540", 
               density.border.col = "#3f51b5")
dev.off()

### Plots comparing posterior means and CI for 2 different conditions ###

## ddRAD vs UCE ##
m <- lm(post_mean_uce ~ post_mean_ddrad, UCEvsddRAD_max23)
r2 = paste0("RÂ² = ", format(summary(m)$r.squared, digits = 4))

ddrad_uce_23 <- ggplot(data=UCEvsddRAD_max23, aes(x=post_mean_ddrad, y=post_mean_uce)) +
  geom_pointrange(aes(ymin=min_95_uce, ymax=max_95_uce), col="#455A64", shape=20, size=.2) +
  geom_pointrangeh(aes(xmin=min_95_ddrad, xmax=max_95_ddrad), col="#455A64", shape=20, size=.2) +
  geom_abline(intercept = 0, slope = 1, col="#455A64", size=0.5) +
  labs(x="Time estimates (Ma): ddRAD, root max. bound = 23 Ma", y="Time estimates (Ma): UCE, root max. bound = 23 Ma") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  annotate("text", x = 5, y = 25, label = r2, color="black", size = 4, parse = FALSE) +
  xlim(0, 30) + ylim(0, 30)

## 3 calibration points vs 4 ##
m <- lm(post_mean_4 ~ post_mean_3, difcalib23vs)
r2 = paste0("RÂ² = ", format(summary(m)$r.squared, digits = 4))

p3vs4_23 <- ggplot(data=difcalib23vs, aes(x=post_mean_3, y=post_mean_4)) +
  geom_pointrange(aes(ymin=min_95_4, ymax=max_95_4), col="#455A64", shape=20, size=.2) +
  geom_pointrangeh(aes(xmin=min_95_3, xmax=max_95_3), col="#455A64", shape=20, size=.2) +
  geom_abline(intercept = 0, slope = 1, col="#455A64", size=0.5) +
  labs(x="Time estimates (Ma): 3 calib. points, root max. bound = 23 Ma", y="Time estimates (Ma): 4 calib. points, root max. bound = 23 Ma") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  annotate("text", x = 5, y = 25, label = r2, color="black", size = 4, parse = FALSE) +
  xlim(0, 30) + ylim(0, 30)

## min. and max. bounds vs only min. ##
p23_low23 <- ggplot(data=max23vslowbound, aes(x=post_mean_23, y=post_mean_lowbound)) +
  geom_pointrange(aes(ymin=min_95_lowbound, ymax=max_95_lowbound), col="#455A64", shape=20, size=.2) +
  geom_pointrangeh(aes(xmin=min_95_23, xmax=max_95_23), col="#455A64", shape=20, size=.2) +
  geom_abline(intercept = 0, slope = 1, col="#455A64", size=0.5) +
  labs(x="Time estimates (Ma): max. and min. bounds, root max. bound = 23 Ma", y="Time estimates (Ma): only min. bounds, root max. bound = 23 Ma") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  xlim(0, 30) + ylim(0, 30)

## root max. bound = 38 vs 23 ##
p38_23 <- ggplot(data=max38vsmax23, aes(x=post_mean_38, y=post_mean_23)) +
  geom_pointrange(aes(ymin=min_95_23, ymax=max_95_23), col="#455A64", shape=20, size=.2) +
  geom_pointrangeh(aes(xmin=min_95_38, xmax=max_95_38), col="#455A64", shape=20, size=.2) +
  geom_abline(intercept = 0, slope = 1, col="#455A64", size=0.5) +
  labs(x="Posterior mean time estimates (Ma): root max. bound = 38 Ma", y="Posterior mean time estimates (Ma): root max. bound = 23 Ma") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  xlim(0, 40) + ylim(0, 40)

p_sel_cor <- ggarrange(ddrad_uce_23,p3vs4_23,p23_low23,p38_23,ncol=2, nrow=2,
                       labels=c("a)","b)","c)","d)"),
                       font.label=list(size=18,face="bold"), 
                       hjust=-.2, vjust=1)
ggsave("sel_dating_comparisons.pdf", p_sel_cor, device="pdf", units="cm",
       width=28, height=28, limitsize=FALSE)

### Posterior estimates and 95% CI for selected nodes with diff. param. and diff. datasets ###

p_key_nodes_dif_cal <- ggplot(data=dselnodes, aes(x=test_tot, y=post_mean, ymin=min_95, ymax=max_95, col=test_tot)) +
  geom_pointrange(fatten=1) +
  facet_grid(node_sp ~ dataset) +
  theme(panel.background = element_rect(fill = "#fafafa",colour = "#fafafa",size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2,linetype = 'solid',colour = "#212121")) +
  theme(legend.title=element_blank()) +
  labs(y = "Ma") +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank()) +
  guides(col = guide_legend(nrow =22)) +
  scale_colour_manual(values=c("#00000440","#000004","#33106840","#331068",
                               "#7D248240","#7D2482","#C83E7340","#C83E73",
                               "#F97C5D40","#F97C5D","#FED39540","#FED395"))

ggsave("key_nodes_dif_cal_prior_post.pdf", p_key_nodes_dif_cal, device="pdf", units="cm", width=18, height=17)

### Infinite-sites plots ###

d_infi_plot_23_ul <- subset(d_whole_nopriors, dataset == "concat_95_tot" & value == "rootmax_23")
d_infi_plot_23_ul$width_CI <- d_infi_plot_23_ul$max_95 - d_infi_plot_23_ul$min_95

d_infi_plot_23_ul_uce <- subset(d_whole_nopriors, dataset == "UCE_95_tot" & value == "rootmax_23")
d_infi_plot_23_ul_uce$width_CI <- d_infi_plot_23_ul_uce$max_95 - d_infi_plot_23_ul_uce$min_95

d_infi_plot_23_ul_ddRAD <- subset(d_whole_nopriors, dataset == "ddRAD_95_tot" & value == "rootmax_23")
d_infi_plot_23_ul_ddRAD$width_CI <- d_infi_plot_23_ul_ddRAD$max_95 - d_infi_plot_23_ul_ddRAD$min_95

m <- lm(width_CI ~ 0 + post_mean, d_infi_plot_23_ul)
w = paste0("w = ", format(coef(m)[1], digits = 2), "t")
r2 = paste0("RÂ² = ", format(summary(m)$r.squared, digits = 4))
eq <- c(w, r2)
eq = paste(eq, collapse="\n")

infi_plot_23 <- ggplot(data=d_infi_plot_23_ul, aes(x=post_mean, y=width_CI)) +
  geom_point(col="#455A64", size=2) +
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, col="#455A64", size=0.5, fullrange=T) +
  labs(x="Posterior mean times (Ma)", y="Posterior mean CI width (Myr)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  annotate("text", x = 5, y = 20, label = eq, color="black", size = 4, parse = FALSE) +
  xlim(0, 33) + ylim(0, 24)

m <- lm(width_CI ~ 0 + post_mean, d_infi_plot_23_ul_uce)
w = paste0("w = ", format(coef(m)[1], digits = 2), "t")
r2 = paste0("RÂ² = ", format(summary(m)$r.squared, digits = 4))
eq <- c(w, r2)
eq = paste(eq, collapse="\n")

infi_plot_23_uce <- ggplot(data=d_infi_plot_23_ul_uce, aes(x=post_mean, y=width_CI)) +
  geom_point(col="#455A64", size=2) +
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, col="#455A64", size=0.5, fullrange=T) +
  labs(x="Posterior mean times (Ma)", y="Posterior mean CI width (Myr)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  annotate("text", x = 5, y = 20, label = eq, color="black", size = 4, parse = FALSE) +
  xlim(0, 33) + ylim(0, 24)

m <- lm(width_CI ~ 0 + post_mean, d_infi_plot_23_ul_ddRAD)
w = paste0("w = ", format(coef(m)[1], digits = 2), "t")
r2 = paste0("RÂ² = ", format(summary(m)$r.squared, digits = 4))
eq <- c(w, r2)
eq = paste(eq, collapse="\n")

infi_plot_23_ddRAD <- ggplot(data=d_infi_plot_23_ul_ddRAD, aes(x=post_mean, y=width_CI)) +
  geom_point(col="#455A64", size=2) +
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, col="#455A64", size=0.5, fullrange=T) +
  labs(x="Posterior mean times (Ma)", y="Posterior mean CI width (Myr)") +
  theme_bw() +
  theme(panel.grid=element_blank()) +
  annotate("text", x = 5, y = 20, label = eq, color="black", size = 4, parse = FALSE) +
  xlim(0, 33) + ylim(0, 24)

p_infi_plot_dif_dat <- ggarrange(infi_plot_23_ddRAD, infi_plot_23_uce, infi_plot_23,
                                 labels=c("a)","b)","c)"),
                                 font.label=list(size=18,face="bold"), 
                                 hjust=0, vjust=1.3, ncol=3, nrow=1)
ggsave("infini_plots_dif_dat.pdf", p_infi_plot_dif_dat, device="pdf", units="cm", 
       width=42, height=14, limitsize=FALSE)

```
