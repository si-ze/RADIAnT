# RADIAnT - Unified identification of statistically-robust RNA-DNA interactions from diverse data types

RADIAnT is a reads-to-interactions pipeline for analyzing RNA-DNA ligation data. Currently, RADIAnT supports the analysis of reads from following sources: 
- RADICL-seq
- GRID-seq
- RedC

RADIAnT calls interactions against a dataset-specific, unified background which considers RNA binding site-TSS distance and genomic region bias. By scaling the background with RNA abundance, RADIAnT is sensitive enough to detect specific interactions of lowly expressed transcripts, while remaining specific enough to discount false positive interactions of highly abundant RNAs.


## Getting started
### Requirements

RADIAnT is provided as a snakemake pipeline, which makes use of several bioinformatic tools. Please make sure that the following software is installed:
- [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [Python](https://www.python.org/downloads/)
- [R](https://cran.r-project.org/doc/manuals/r-release/R-admin.html)
- [STAR](https://github.com/alexdobin/STAR)
- [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)
- [SAMtools](https://www.htslib.org/download/)
- [BAMcoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html)
- [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/)

### Download and Install

To obtain RADIAnT from GitHub, you can clone the directory in the following way:

```git clone https://github.com/si-ze/RADIAnT.gitt```

or using pip:

```pip install git+git://github.com/si-ze/RADIAnT.git```


### Run

Navigate to the directory to which you have downloaded the RADIAnT repository and issue the snakemake command. Example using the provided test files 

```snakemake -s /mnt/d/RADIAnT/example/Snakefile.smk --configfile=example/test_run/config_RADICL_mESCs.yaml --cores 32```

##  Git repository structure


This repository follows the [recommended snakemake repository structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html)

```
RADIAnT/
├── LICENSE
├── README.md
├── Rplot001.svg
├── Rplots.pdf
├── config
│   ├── config_GRID_mESCs.yaml
│   ├── config_RADICL_mESCs.yaml
│   └── config_RedC_K562.yaml
├── example
│   ├── GRID
│   │   ├── GSM2396700_DNA.fastq.gz
│   │   ├── GSM2396700_RNA.fastq.gz
│   │   └── README.md
│   ├── RADICL
│   │   ├── README.md
│   │   ├── mESC_2FA_n1_DNA.fastq.gz
│   │   ├── mESC_2FA_n1_RNA.fastq.gz
│   │   ├── mESC_2FA_n1_RNA_depleted.fastq
│   │   ├── mESC_2FA_n1_rRNA.fastq
│   │   └── mESC_2FA_n1_rRNA_removal_stats.txt
│   ├── README.md
│   └── RedC
│       ├── README.md
│       ├── Red-C_K562_N2_DNA.fastq.gz
│       ├── Red-C_K562_N2_RNA3.fastq.gz
│       └── Red-C_K562_N2_RNA5.fastq.gz
├── resources
│   ├── human
│   │   ├── README.md
│   │   ├── Silva_rRNA_Homo_sapiens.fasta.gz
│   │   ├── gencode.v45.annotation.genes.bed.gz
│   │   ├── gencode.v45.annotation.gtf.gz
│   │   ├── hg38-blacklist.v2.bed.gz
│   │   └── hg38_genome_5kb_bins_named.bed.gz
│   └── mouse
│       ├── README.md
│       ├── arb-silva.de_2024-03-25_id1313737_tax_silva_trunc.fasta.gz
│       ├── gencode.vM29.annotation.gtf.gz
│       ├── mm39.05kb.bed.gz
│       ├── mm39.excluderanges.bed.gz
│       ├── mm39.ncbiRefSeq.genes.bed.gz
│       └── mm39.ncbiRefSeq.gtf.gz
└── workflow
    ├── RADIAnT.smk
    └── scripts
        ├── RADIAnT_command_line.R
        ├── intersect_processing.R
        ├── plot_read_stats.R
        └── redc_intersect_max.R
```

## Example runs 

The subdirectory ```config``` holds an exemplary configuration file for a RADIAnT analysis of GRID-sequencing data (```config_GRID_mESCs.yaml```),  RADICL-sequencing data (```config_RADICL_mESCs.yaml```),  and RedC data (```config_RedC_K562.yaml```), respectively. Input reads, sampled from publicly available data, can be found in the subdirectory  ```example``` . \\\

To e.g. run the RADICL test case, issue 

```snakemake -s /mnt/d/github_test/RADIAnT/workflow/RADIAnT.smk --configfile=/mnt/d/github_test/RADIAnT/config/config_RADICL_mESCs.yaml --cores 32```

Any of the provided config files will create a new subdirectory in the ```RADIAnT``` directory, called ```results```. For each of the different sequencing methods, another subdirectory will be created in ```results```, holding the outputs of the RADIAnT pipeline. The location of the output files can be defined in the configuration file (parameter ```output_directory```).

```
RADIAnT/results/RADICL/
├── bam
├── counts
├── fastq
├── interactions
├── intersects
├── logs
└── merge
```

The main results are located in the ```interactions``` direcotry and described in the following subsection. 



### Output figures

RADIAnT outputs several logs and figures.

The Sankey plot illustrates the information flow through the pipeline, providing insight into the efficiency of the sequencing method at hand. It provides insight into the proportion of initial input reads considered valid post preprocessing and how many of those are linked to significant interactions. The corresponding counts can be found in the ```*_read_stats.txt```

<p align="center">
  <img src="https://github.com/user-attachments/assets/2a68de4a-b972-4efb-aada-b7efdac537f9" width=75% height=75%>
</p>


```genes_barplot.all_interactions.png``` and ```genes_barplot.significant_interactions.png``` show the top 15 genes by total number of interactions and significant number of interactions, respectively. 

<p float="left">
  <img src="https://github.com/user-attachments/assets/8cbae406-9b54-404d-9e34-441d893020e9" width=45% height=45%>
  <img src="https://github.com/user-attachments/assets/2ff51a62-40c6-41ea-9c9c-11527aa706de" width=45% height=45%>
</p>




```genes_barplot.all_reads.png``` and ```genes_barplot.significant_reads.png``` show the top 15 genes by total number of reads and reads originating from significant interactions, respectively. 

<p float="left">
  <img src="https://github.com/user-attachments/assets/3408380a-698f-4146-bc72-02f8cd7c2a8e" width=45% height=45%>
  <img src="https://github.com/user-attachments/assets/e89670ef-0ec4-44a6-a6bf-b3b9dc84f943" width=45% height=45%>
</p>

