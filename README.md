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

```git clone https://github.com/SchulzLab/STARE.git```

or using pip:

```pip install git+git://github.com/si-ze/RADIAnT.git```


### Run

Navigate to the directory to which you have downloaded the RADIAnT repository and issue the snakemake command. Example using the provided test files 

```snakemake -s /mnt/d/RADIAnT/example/Snakefile.smk --configfile=example/test_run/config_RADICL_mESCs.yaml --cores 32```

