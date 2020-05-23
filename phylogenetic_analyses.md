# Phylogenetic analyses

####################################
## 1) Partition scheme generation ##
####################################

### Script to generate a nexus file with each UCE in the 75% dataset defined as a CHARSET #

`cd $path/alignments/dataset_out_phased_iupac_trimmed_cleaned_selected_75`

### Convert nexus alignments for each locus to fasta
```
for file in ./*nexus; do cat $file | grep -vE "; | NEXUS | matrix" | sed -E 's/^/>/' | tr ' ' '\n' | awk 'NF' > ../dataset_out_ph$
```
### Merge all fasta files to generate a nexus concatenated file
`TriSeq -in *.fa -of nexus -o ../UCE_75_phased_iupac`
 
### Script to run Entropy-based automatic partitioning of UCE alignments (SWSC)
`python SWSCEN.py UCE_75_phased_iupac.nex $path/phylo_analyses/partition_finder/UCE_phased_iupac_SWSC_75`

### Run PartitionFinder2
`python2 /ddn/data/sbvd77/soft/partitionfinder-2.1.1/PartitionFinder.py $path/phylo_analyses/partition_finder/UCE_phased_iupac_SWS$`

########################################
## 2) Concatenated analyses: Exabayes ##
########################################

## For each dataset: ddRAD_75, ddRAD_95, UCE_75_contig, UCE_95_contig, UCE_75_iupac, UCE_95_iupac, UCE_75_iupac_partitioned, TENT_75_iupac ###
### Run Exabayes ##
```
cd $HOME/Durham/scripts
mkdir $HOME/Durham/phylo_analyses/ExaBayes/UCE_75_nopart/
output=$HOME/Durham/phylo_analyses/ExaBayes/UCE_75_nopart/
mpirun -np 64 exabayes -f /users/jferrer/Durham/scripts/UCE_75.phylip -m DNA -n UCE_75_nopart -q UCE_75.part -s $RANDOM -c config.nexus -w $output
```
### Create a consensus tree ##
`consense -f $output/ExaBayes_topologies.run-*.UCE_75_nopart -n UCE_75_nopart`

### Generate the 50% credible set of trees ##
`credibleSet -f $output/ExaBayes_topologies.run-*.UCE_75_nopart -n UCE_75_nopart`

### Check how well the parameters were sampled ##
`postProcParam -f $output/ExaBayes_parameters.run-*.UCE_75_nopart -n UCE_75_nopart`

####*.part
`DNA, p1=1-length alignment (except for the partitions analysis for which each partition is in a new line`

####config.nexus
####NEXUS
```
begin run;
 numruns 2
 numgen 1e6
 parsimonystart false
 numcoupledchains 4
 burninproportion 0.25
end;
```
########################################
## 3) Concatenated analyses: raxml-ng ##
########################################

## For each dataset: ddRAD_75, ddRAD_95, UCE_75_contig, UCE_95_contig, UCE_75_iupac, UCE_95_iupac, UCE_75_iupac_partitioned, TENT_75_iupac ###

### Change directory ##
`cd $HOME/Durham/phylo_analyses/raxml-ng/UCE_75_phased_iupac_75post/`

### Perform 50 tree searches using 25 random and 25 parsimony-based starting trees ##
```
raxml-ng --msa UCE_75_phased_iupac_partitions_75post.phy --model GTR+G --prefix UCE_75_phased_iupac_75post --seed 2 --threads 50 --brlen scaled --tree pars{25},rand{25}
```
### Perform 500 standard non-parametric bootstrap with 500 replicates ##
```
raxml-ng --bootstrap --msa UCE_75_phased_iupac_partitions_75post.phy --model GTR+G --prefix UCE_75_phased_iupac_75post  --seed 333 --threads 50 --bs-trees 500
```
### Check bootstrap convergence ##
```
raxml-ng --bsconverge --bs-trees UCE_75_phased_iupac_75post.raxml.bootstraps --prefix UCE_75_phased_iupac_75post --seed 2 --threads 32  --bs-cutoff 0.03
```
### Map bootstrap support values onto the best-scoring ML tree ##
```
raxml-ng --support --tree UCE_75_phased_iupac_75post.raxml.bestTree --bs-trees UCE_75_phased_iupac_75post.raxml.bootstraps --prefix UCE_75_phased_iupac_75post --threads 50
```
######################################
## 4) Species tree analyses: astral ##
######################################

## For ddRAD and UCE ###

### Generate a file with all the scripts to run gene trees with RAxML with the parallel function ##
```
for i in *phy; do gene_id=$(echo $i | sed 's/.phy//'); printf raxmlHPC\-PTHREADS\-SSE3\ \-s\ $i\ \-T\ 10\ \-n\ $gene_id\ \-m\ GTRGAMMA\ \-p\ $RANDOM\ \-f\ a\ \-x\ $RANDOM\ \-N\ 100\ \-w\ \/ddn\/data\/sbvd77\/UCE\/phylo\_analyses\/raxml\/gene\_trees'\n' >> ../../scripts/raxml_gene_trees_parallel; done
```
### Run the gene trees with RAxML in parallel ##
`parallel -j5 < raxml_gene_trees_parallel`

### Put all trees in a file for 75% dataset ##
`cat RAxML_bipartitions.* > ../../astral_new/ml_best_ddRAD_75.trees`

### Make a list of 95% loci ##
`cat PIS_in_95 | grep -v "locus" | cut -d "." -f1 | sed -E 's/$/\$/g' > ddRAD_95.tsv`

### Create the input files for the 95% dataset ##
`cat $(ls ../raxml/gene_trees/RAxML_bipartitions.uce* | grep -E -f ddRAD_95.tsv) > ml_best_ddRAD_95.trees`

### Collapse nodes with BS < 10 ##
`nw_ed ml_best_ddRAD_75.trees 'i & b<=10' o > ml_best_ddRAD_75_BS10.trees`
`nw_ed ml_best_ddRAD_95.trees 'i & b<=10' o > ml_best_ddRAD_95_BS10.trees`

### Run default for the non-contracted and contracted trees ##
```
java -jar ../../../programari/Astral/astral.5.6.3.jar -i ./ml_best_ddRAD_75.trees -o ./astral_default_ddRAD_75.tre 2> ./astral_default_ddRAD_75.log
java -jar ../../../programari/Astral/astral.5.6.3.jar -i ./ml_best_ddRAD_75_BS10.trees -o ./astral_default_ddRAD_75_BS10.tre 2> ./astral_default_ddRAD_75_BS10.log

java -jar ../../../programari/Astral/astral.5.6.3.jar -i ./ml_best_ddRAD_95.trees -o ./astral_default_ddRAD_95.tre 2> ./astral_default_ddRAD_95.log
java -jar ../../../programari/Astral/astral.5.6.3.jar -i ./ml_best_ddRAD_95_BS10.trees -o ./astral_default_ddRAD_95_BS10.tre 2> ./astral_default_ddRAD_95_BS10.log
```
### Run the scoring scripts with -t 2 (quartet scores) and -t 10 (polytomy test) ##
```
java -jar ../../../programari/Astral/astral.5.6.3.jar -q astral_default_ddRAD_75.tre -i ./ml_best_ddRAD_75.trees -t 2 -o astral_scored_t2_ddRAD_75.tre 2> astral_scored_t2_ddRAD_75.log
java -jar ../../../programari/Astral/astral.5.6.3.jar -q astral_default_ddRAD_75_BS10.tre -i ./ml_best_ddRAD_75_BS10.trees -t 2 -o astral_scored_t2_ddRAD_75_BS10.tre 2> astral_scored_t2_ddRAD_75_BS10.log

java -jar ../../../programari/Astral/astral.5.6.3.jar -q astral_default_ddRAD_95.tre -i ./ml_best_ddRAD_95.trees -t 2 -o astral_scored_t2_ddRAD_95.tre 2> astral_scored_t2_ddRAD_95.log
java -jar ../../../programari/Astral/astral.5.6.3.jar -q astral_default_ddRAD_95_BS10.tre -i ./ml_best_ddRAD_95_BS10.trees -t 2 -o astral_scored_t2_ddRAD_95_BS10.tre 2> astral_scored_t2_ddRAD_95_BS10.log

java -jar ../../../programari/Astral/astral.5.6.3.jar -q astral_default_ddRAD_75.tre -i ./ml_best_ddRAD_75.trees -t 10 -o astral_scored_t10_ddRAD_75.tre 2> astral_scored_t10_ddRAD_75.log
java -jar ../../../programari/Astral/astral.5.6.3.jar -q astral_default_ddRAD_75_BS10.tre -i ./ml_best_ddRAD_75_BS10.trees -t 10 -o astral_scored_t10_ddRAD_75_BS10.tre 2> astral_scored_t10_ddRAD_75_BS10.log

java -jar ../../../programari/Astral/astral.5.6.3.jar -q astral_default_ddRAD_95.tre -i ./ml_best_ddRAD_95.trees -t 10 -o astral_scored_t10_ddRAD_95.tre 2> astral_scored_t10_ddRAD_95.log
java -jar ../../../programari/Astral/astral.5.6.3.jar -q astral_default_ddRAD_95_BS10.tre -i ./ml_best_ddRAD_95_BS10.trees -t 10 -o astral_scored_t10_ddRAD_95_BS10.tre 2> astral_scored_t10_ddRAD_95_BS10.log
```
## For the UCE analyses ###
### Remove samples with high levels of missing data from the alignments and redo gene trees ##
```
TriSeq -c -in *.phy -rm PGav1 ACar1 AGri1 CDio2 -o ../dataset_out_phased_iupac_trimmed_cleaned_selected_75_pruned_phylip -of phylip
```
## For the TENT analyses ###
### Concatenate the tree files ##
```
cat ml_best_ddRAD_75.trees > ml_best_ddRAD-UCE_75.trees
cat ml_best_UCE_75.trees >> ml_best_ddRAD-UCE_75.trees

cat ml_best_ddRAD_75_BS10.trees > ml_best_ddRAD-UCE_75_BS10.trees
cat ml_best_UCE_75_BS10.trees >> ml_best_ddRAD-UCE_75_BS10.trees
```
#### Run the script to rename them according to the new names ##
```
ddRAD=($(cat $HOME//Durham/concat_UCE-ddRAD/correspondence_ddRAD-UCE-concat.tsv | grep -v "ddRAD" | cut -f1))
concat=($(cat $HOME/Durham/concat_UCE-ddRAD/correspondence_ddRAD-UCE-concat.tsv | grep -v "ddRAD" | cut -f3))
UCE=($(cat $HOME/Durham/concat_UCE-ddRAD/correspondence_ddRAD-UCE-concat.tsv | grep -v "ddRAD" | cut -f2))
for ((i = 0; i < ${#ddRAD[@]}; i++)); do sed -iE "s/${ddRAD[i]}/${concat[i]}/g" ml_best_ddRAD-UCE_75i.trees; done
for ((i = 0; i < ${#UCE[@]}; i++)); do sed -iE "s/${UCE[i]}/${concat[i]}/g" ml_best_ddRAD-UCE_75i.trees; done
sed -i 's/.assembled//g' ml_best_ddRAD-UCE_75i.trees
rm ml_best_ddRAD-UCE_75i.trees

for ((i = 0; i < ${#ddRAD[@]}; i++)); do sed -iE "s/${ddRAD[i]}/${concat[i]}/g" ml_best_ddRAD-UCE_75_BS10i.trees; done
for ((i = 0; i < ${#UCE[@]}; i++)); do sed -iE "s/${UCE[i]}/${concat[i]}/g" ml_best_ddRAD-UCE_75_BS10i.trees; done
sed -i 's/.assembled//g' ml_best_ddRAD-UCE_75_BS10i.trees
nw_prune ml_best_ddRAD-UCE_75_BS10i.trees PAAsD1 PAHa2 > ml_best_ddRAD-UCE_75_BS10.trees
rm ml_best_ddRAD-UCE_75_BS10i.trees
```
### Run default for the non-contracted and contracted trees ##
```
java -jar ../../../programari/Astral/astral.5.6.3.jar -i ./ml_best_ddRAD-UCE_75.trees -o ./astral_default_ddRAD-UCE_75.tre 2> ./astral_default_ddRAD-UCE_75.log
java -jar ../../../programari/Astral/astral.5.6.3.jar -i ./ml_best_ddRAD-UCE_75_BS10.trees -o ./astral_default_ddRAD-UCE_75_BS10.tre 2> ./astral_default_ddRAD-UCE_75_BS10.log
```
### Run the scoring scripts with -t 2 (quartet scores) and -t 10 (polytomy test) ##
```
java -jar ../../../programari/Astral/astral.5.6.3.jar -q astral_default_ddRAD-UCE_75.tre -i ./ml_best_ddRAD-UCE_75.trees -t 2 -o astral_scored_t2_ddRAD-UCE_75.tre 2> astral_scored_t2_ddRAD-UCE_75.log
java -jar ../../../programari/Astral/astral.5.6.3.jar -q astral_default_ddRAD-UCE_75_BS10.tre -i ./ml_best_ddRAD-UCE_75_BS10.trees -t 2 -o astral_scored_t2_ddRAD-UCE_75_BS10.tre 2> astral_scored_t2_ddRAD-UCE_75_BS10.log

java -jar ../../../programari/Astral/astral.5.6.3.jar -q astral_default_ddRAD-UCE_75.tre -i ./ml_best_ddRAD-UCE_75.trees -t 10 -o astral_scored_t10_ddRAD-UCE_75.tre 2> astral_scored_t10_ddRAD-UCE_75.log
java -jar ../../../programari/Astral/astral.5.6.3.jar -q astral_default_ddRAD-UCE_75_BS10.tre -i ./ml_best_ddRAD-UCE_75_BS10.trees -t 10 -o astral_scored_t10_ddRAD-UCE_75_BS10.tre 2> astral_scored_t10_ddRAD-UCE_75_BS10.log
```
#####################################
## 5) Species tree analyses: SNAPP ##
#####################################

### Build xml files with beauti ##
`cl89_pc75_allspecies_1.xml => ddRAD`
`snapp_UCE_phased_default.xml => UCE`

### Run the analyses with BEAST2 ##
`beast -threads 25 -seed $RANDOM cl89_pc75_allspecies_1.xml`

`beast -threads 25 -seed $RANDOM snapp_UCE_phased_default.xml`

### Summarize the sample of trees with Treeannotator ##

###############
## End of script ##
###############
