# Tent data assembly

#########################################################################################################################
## 1) Select a representative sequence of the genus Puffinus for each ddRAD and UCE 75% datasets to map ddRADs to UCEs ##
#########################################################################################################################
```
for file in *.phy; do cat $file | grep -E "^P" | shuf | head -1 | sed 's/^/>/' | sed -E "s/(assembled)\ +([a-z]+.*)/\1\_$file\n\2/" | sed 's/.phy//' >> ../ddRAD_mapping_UCEs.fa; done # ddRAD fragments

sed -i 's/-//g' ddRAD_mapping_UCEs.fa # Take out the gap symbols

for file in *phy; do cat $file | grep -E "^P" | shuf | head -1 | sed 's/^/>/' | sed -E "s/(^>[a-zA-Z0-9]+)\ +([a-z]+.*)/\1\_$file\n\2/" | sed 's/.phy//' | sed -E 's/^n+//' | sed -E 's/n+$//' | sed 's/-//g' >> ../UCE_to_map_ddRAD.fa; done # UCE fragments
```
########################################
## 2) Blast the ddRAD to the UCE loci ##
########################################

### Build the database for blasat ###
```
makeblastdb -in UCE_to_map_ddRAD.fa -dbtype nucl -out blast/index/UCE_db
```
### Run blast ###
```
blastn -query ddRAD_mapping_UCEs.fa -db blast/index/UCE_db -evalue 1e-3 -max_target_seqs 1 -out blast/ddRAD2UCE.tsv -outfmt 6
```
######################################################
## 3) Remove the ddRAD loci that mapped to UCE loci ##
######################################################
```
cat ddRAD2UCE.tsv | cut -f1 | cut -d "_" -f2 | sed -E 's/$/.phy/g' | tr '\n' ' '

rm <output from the previous command>
```
###################################
## 4) Prepare the TENT alignment ##
###################################

### ddRAD phylip alignment ###
```
TriSeq -in p*phy -of phylip -rm PAHa2.assembled --missing-filter 25 25 -o ddRAD_75
```
### UCE phylip alignment ###
```
TriSeq -in uce*phy -of phylip -rm PAAsD1 --missing-filter 25 25 -o UCE_75
```
### TENT phylip alignment ###
```
TriSeq -in ddRAD_75.phy UCE_75.phy -of phylip -o concat_UCE-ddRAD
```
