# Extracting Raw Sequences for Specific OTUs from Selected Samples in mothur  

### Author: Conor O'Kane  
### Tools Used: mothur, R  
---

## Overview  
This document provides a step-by-step guide for extracting raw sequences from specific **OTUs** in selected samples using **mothur** and **R**. 

For example, in my blocking primer project, the same dietary sample was generating different OTUs depending on which blocking primer was included (lake trout reads in some, salmonid reads in others). I wanted to take a look at the raw amplicon sequences that were being grouped into those OTUs for each sample to assess their similarity. This was the way I did so. 

---

## 1. Setting Up mothur
### Load mothur
Before running mothur, ensure that the correct version is loaded. Use `module spider` to see the current version on the HPCC, then copy/paste to load.
```bash
module spider mothur
module load Mothur/1.48.0-foss-2023a-Python-3.11.3  # Current version at the time of writing
```
We are then read to open up mothur:
```bash
mothur
```
### Select the Sample Groups
Run `get.groups` in mothur to select the relevant samples (groups) that you want to include from a .shared file. Use the sample names separated by "-". This uses the .shared file generated from `make.shared` in previous mothur protocols.

In my case, my three dietary samples were:  
40_PMBlock_B_C3_CA14.12S  
40_PMBlock_A_dT_H_CA14.12S  
40_PMBlock_B_C3_H_CA14.12S
```bash
get.groups(shared=stability.trim.contigs.good.unique.good.precluster.pick.opti_mcc.shared, groups=40_PMBlock_B_C3_CA14.12S-40_PMBlock_A_dT_H_CA14.12S-40_PMBlock_B_C3_H_CA14.12S)
```
**Output Files (in my case):**
```
stability.trim.contigs.good.unique.good.precluster.pick.opti_mcc.0.01.pick.shared
```
Rename the output file with for easier reference. Copy the previous output for the `shared` argument:
```bash
rename.file(shared=stability.trim.contigs.good.unique.good.precluster.pick.opti_mcc.0.01.pick.shared, new=CA14.Blockers2.3.5.shared)
```

---

## 2. Identify Top OTUs in R

To better visualize the top OTUs in your samples, we can organize them with a bit of R. This can all still be done in the terminal window, and R should already be loaded on the HPCC. So `quit` out of **mothur** and open up **R**. 

```bash
R
```

### Read in the Shared File
This is the renamed .shared file from the previous mothur step.
```r
shared <- read.table("CA14.Blockers2.3.5.shared", header=TRUE, sep="\t", stringsAsFactors=FALSE)
```

### Extract the Top 5 OTUs per Sample
In my case, I was interested in only the top few OTUs per sample. Feel free to adjust the number of "top" OTUs that you want to inspect. This is done through the `head()` command in the final line of the loop. 
```r
otu_counts <- shared[ , -(1:3)]
top5_OTUs_list <- list()

for(i in 1:nrow(otu_counts)) {
  row_counts <- as.numeric(otu_counts[i, ])
  names(row_counts) <- names(otu_counts)
  sorted <- sort(row_counts, decreasing=TRUE)
  top5_OTUs_list[[shared$Group[i]]] <- head(sorted, 5)
}

top5_OTUs_list
```

**Output:**
You should get something that looks like this, depending on how many samples and OTUs you choose to look at. 
```
$`40_PMBlock_A_dT_H_CA14.12S`
Otu0005 Otu0004 Otu0229 Otu0365 Otu0243
  12754    1141      61      45      44

$`40_PMBlock_B_C3_CA14.12S`
Otu0006 Otu0008 Otu0001 Otu0004 Otu0147
   5074    4353    2746     541     361

$`40_PMBlock_B_C3_H_CA14.12S`
Otu0001 Otu0004 Otu0028 Otu0010 Otu0020 
  13991    2507      63      62      54 
```
In my case, I'll take a look at OTUs 0005, 0004, 0006, 0008, and 0001.

### Save the Selected OTUs to a File
Modify the list below to match the OTUs you want to extract:
```r
selected_otus <- c("Otu0005","Otu0004","Otu0006","Otu0008","Otu0001")
write.table(
  x = selected_otus,
  file = "CA14.blockers235.accnos", # modify to your own file name
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = "\t"
)
```
The output should look like this:
```bash
Otu0005
Otu0004
Otu0006
Otu0008
Otu0001
```
You can `quit()` out of R

---

## 3. Extract Sequences with mothur
Hop back into **mothur** to finish the sequence extractions
```bash
mothur
```


Use `get.otus` to filter the dataset to only the selected OTUs from a list file. This list file would have been generated from `cluster` with the lab mothur protocol. Use the .accnos file you just saved for OTU selection:

```bash
get.otus(list=stability.trim.contigs.good.unique.good.precluster.pick.opti_mcc.list, accnos=CA14.blockers235.accnos)
```

Use `list.seqs` to list all the sequence names (accession numbers) from selected OTUs. This uses the output from the previous step:

```bash
list.seqs(list=stability.trim.contigs.good.unique.good.precluster.pick.opti_mcc.0.01.pick.list)
```
Generate filtered FASTA and count files with `get.seqs`. This extracts only the sequences corresponding to the selected OTUs from the original FASTA file. Use the FASTA and count table outputs from previous mothur protocol, and the accnos file from the previous step:

```bash
get.seqs(fasta=stability.trim.contigs.good.unique.fasta, count=stability.trim.contigs.good.unique.good.precluster.count_table, accnos=stability.trim.contigs.good.unique.good.precluster.pick.opti_mcc.0.01.pick.accnos)
```

Split the extracted sequences into separate FASTA and count table files for each sample with `split.groups`. Use the same groups as used with `get.groups` from before:
```bash
split.groups(fasta=stability.trim.contigs.good.unique.pick.fasta, count=stability.trim.contigs.good.unique.good.precluster.pick.count_table, groups=40_PMBlock_B_C3_CA14.12S-40_PMBlock_A_dT_H_CA14.12S-40_PMBlock_B_C3_H_CA14.12S)
```
**Output Files:**
These are the output files, split by each of your selected samples that contain the raw sequences (.fasta) and count numbers (.count_table) for your selected OTUs. 
```
stability.trim.contigs.good.unique.good.precluster.pick.40_PMBlock_B_C3_CA14.12S.count_table
stability.trim.contigs.good.unique.pick.40_PMBlock_B_C3_CA14.12S.fasta
stability.trim.contigs.good.unique.good.precluster.pick.40_PMBlock_A_dT_H_CA14.12S.count_table
stability.trim.contigs.good.unique.pick.40_PMBlock_A_dT_H_CA14.12S.fasta
stability.trim.contigs.good.unique.good.precluster.pick.40_PMBlock_B_C3_H_CA14.12S.count_table
stability.trim.contigs.good.unique.pick.40_PMBlock_B_C3_H_CA14.12S.fasta
```

---

## 4. Classify Sequences (Optional)
For myself, I wanted to verify the classification chain to see which sequences were being grouped into the different OTUs (e.g. what was being classified as "lake trout" vs "salmonid").

To do so, I used `classify.seqs` with the generated count and fasta files along with my reference database and taxonomy files. This generates wang.taxonomy and wang.tax.summary files to further examine how sequences were classified.
```bash
classify.seqs(count=stability.trim.contigs.good.unique.good.precluster.pick.40_PMBlock_B_C3_CA14.12S.count_table, fasta=stability.trim.contigs.good.unique.pick.40_PMBlock_B_C3_CA14.12S.fasta, reference=FishOnly_12S_align_JK_noprimers_091922.fas, taxonomy=FishOnly_12S_rDNA_taxonomy_091922.txt, cutoff=80)
```
Repeat for other sample groups:
```bash
classify.seqs(count=stability.trim.contigs.good.unique.good.precluster.pick.40_PMBlock_A_dT_H_CA14.12S.count_table, fasta=stability.trim.contigs.good.unique.pick.40_PMBlock_A_dT_H_CA14.12S.fasta, reference=FishOnly_12S_align_JK_noprimers_091922.fas, taxonomy=FishOnly_12S_rDNA_taxonomy_091922.txt, cutoff=80)
classify.seqs(count=stability.trim.contigs.good.unique.good.precluster.pick.40_PMBlock_B_C3_H_CA14.12S.count_table, fasta=stability.trim.contigs.good.unique.pick.40_PMBlock_B_C3_H_CA14.12S.fasta, reference=FishOnly_12S_align_JK_noprimers_091922.fas, taxonomy=FishOnly_12S_rDNA_taxonomy_091922.txt,cutoff=80)
classify.seqs(count=stability.trim.contigs.good.unique.good.precluster.pick.40_PMBlock_B_C3_H_CA14.12S.count_table, fasta=stability.trim.contigs.good.unique.pick.40_PMBlock_B_C3_H_CA14.12S.fasta, reference=FishOnly_12S_align_JK_noprimers_091922.fas, taxonomy=FishOnly_12S_rDNA_taxonomy_091922.txt,cutoff=80)
```

