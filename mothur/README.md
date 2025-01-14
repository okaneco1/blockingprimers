## mothur Resources

This folder contains the resources for demultiplexing and bioinformatic filtering via mothur.


### Files

- **FASTA File**:  
  This is the alignment database used for this project—a list of Great Lakes fish. Primers are not included in the alignment.

- **MakeContigs**:  
  A job submission script for a SLURM workload manager. It runs the `make.contigs` command to assemble paired-end reads into contigs, handle demultiplexing, and perform quality control. This step prepares the sequencing data for downstream analysis.

- **mothur_job_submission**:  
  Another job submission script that calls `SL_12S_BP_mothur_script.sh`.

- **mothur_script.sh**:  
  This script performs the bioinformatic workflow using **mothur** to process the sequencing data. Tasks include quality control, sequence alignment, clustering, and taxonomic classification.



### mothur Script Breakdown

**Input and Filtering:**  
`summary.seqs`: Generates sequence statistics (`stability.trim.contigs.fasta`).  
`screen.seqs`: Filters sequences by ambiguity (`maxambig=0`), length (`maxlength=107`), and quality metrics.  
`unique.seqs`: Removes duplicate sequences to retain unique ones.  
`count.seqs`: Updates sequence counts based on group information.  

**Alignment and Preprocessing:**  
`align.seqs`: Aligns sequences to a reference (`FishOnly_12S_align_JK_noprimers_091922.fas`).  
`screen.seqs`: Further filters alignments (`minlength=85`, `maxlength=110`, `maxhomop=7`).  
`pre.cluster`: Groups sequences with up to 0 mismatches (`diffs=0`).  

**Chimera Detection:**  
`chimera.vsearch`: Identifies chimeric sequences (`dereplicate=t`).  
`remove.seqs`: Removes detected chimeric sequences.  

**Taxonomic Classification:**  
`classify.seqs`: Assigns taxonomy using a reference and taxonomy file (`cutoff=80`).  

**Clustering and OTU Generation:**  
`dist.seqs`: Calculates pairwise distances (`cutoff=0.03`).  
`cluster`: Clusters sequences into OTUs (`cutoff=0.01`).  

**Output Generation:**  
`make.shared`: Creates a shared file for community analysis.  
`classify.otu`: Assigns taxonomy to OTUs.  
`count.groups`: Counts the number of sequences per group.  
