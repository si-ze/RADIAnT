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

```git clone https://github.com/si-ze/RADIAnT.git```

or using pip:

```pip install git+git://github.com/si-ze/RADIAnT.git```


### Configure the run 

In ```/path/to/RADIAnT/config/``` example configuration files have been added for the analysis of RADICL-seq, GRID-seq and Red-C, respectively. Edit the adequate file to configure your snakemake run of RADIAnT. 
* Replace ```/path/to/RADIAnT/``` with the absolute path of the local directory you have cloned the repository to. One way to do so is to replace  ```/path/to/RADIAnT/``` in an text editor. It could also be achieved using sed:

```
cd config
mv config_RADICL_mESCs.yaml backup_config_RADICL_mESCs.yaml
sed -e 's|/path/to/RADIAnT/|/my/actual/path/RADIAnT/|g' backup_config_RADICL_mESCs.yaml > config_RADICL_mESCs.yaml
```

where  ```/my/actual/path/RADIAnT/``` is the path output by running the  ```pwd``` command in the ```RADIAnT``` directory.

* Specify sequencing method used to produce reads (```method``` variable) 
* Provide RADIAnT with the aboslute path to the RADIANT ```workflow``` directory (```workflow_directory``` variable)
* Provide RADIAnT with the aboslute path to the fastq files of your experiment (```fastq_directory``` variable)
* All fastq files located in ```fastq_directory```which belong to the experiment to be analysed must start with the same prefix. Specify this in the ```sample_base``` variable.
* Also specify the suffixes of the DNA and RNA reads. If for example your reads are named ```my_run.DNA.fastq``` and ```my_run.RNA.fastq```, set
    * ```sample_base: my_run.```
    * ```dna_fastq_suffix: DNA.fastq``` # This variable *must* end in ```fastq```. Even if the files are compressed: please *omit* the .gz extension.
    * ```rna_fastq_suffix: RNA.fastq``` # This variable *must* end in ```fastq```. Even if the files are compressed: please *omit* the .gz extension.
* Please note that the pipeline will automatically detect if the files are gzipped. Please *do not* add the .gz extension, the pipeline implicitly checks for that.
* Specify which directory the results should be written to (```output_directory``` variable)
* Provide RADIAnT with information on the organism sequenced (```species``` variable, currently supoported: ```mouse``` and ```human```) and appropriate annotation files
* Provide RADIAnT with information on where required bioinformatic tools have been installed to. 


### Run RADIAnT

Navigate to the directory to which you have downloaded the RADIAnT repository and issue the snakemake command. Example using the provided test files 

```snakemake -s /path/to/RADIAnT/workflow/RADIAnT.smk --configfile=/path/to/RADIAnT/config/config_RADICL_mESCs.yaml --cores 8```

Please make sure to **replace** ```/path/to/RADIAnT/``` with the **absolute path of the local directory** you have cloned the repository to. **Both** in 
* the **command** issued to run snakemake and
* in **all paths referenced in the config file**.






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

```snakemake -s /path/to/RADIAnT/workflow/RADIAnT.smk --configfile=/path/to/RADIAnT/config/config_RADICL_mESCs.yaml --cores 32```

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

-------------------------------------------------------------------------------

```
RADIAnT/results/RADICL/
├── ...
├── logs 📂
    └── *_read_stats.txt
    └── *_sankey.*
├── ...
```


The Sankey plot illustrates the information flow through the pipeline, providing insight into the efficiency of the sequencing method at hand. It provides insight into the proportion of initial input reads considered valid post preprocessing and how many of those are linked to significant interactions. The corresponding counts can be found in the ```*_read_stats.txt```


<p align="center">
  <img src="https://github.com/user-attachments/assets/2a68de4a-b972-4efb-aada-b7efdac537f9" width=55% height=55%>
</p>



-------------------------------------------------------------------------------
```
RADIAnT/results/RADICL/
├── ...
├── interactions 📂
    └── ...
    └── *_RADIAnT_results.txt
    └── ...
├── ...
```

```*_RADIAnT_results.txt``` is the **main output file** hodling ALL (also not significant) interactions, information on the RNA and DNA parts involved on the information, as well as significance values. The results are provided in a tab-delimited table.  Examplary exerpt: 

```
InteractionID      Symbol  Bin     BinChr  BinCentre       GeneChr GeneLeft        GeneRight       ReadCount       ExpectedCount  Method  OriginalMethod  P       Padj
...
Neat1::chr19_1839       Neat1   chr19_1839      chr19   9192500 chr19   5867500 5902500 2       0.000423059069622596  Distance    Cis      8.94642526389563e-08    1.65994756850035e-07
Neat1::chr19_3788       Neat1   chr19_3788      chr19   18937500        chr19   5867500 5902500 1       7.05098449370993e-05    Distance Cis     1       1
Necab2::chr8_24172      Necab2  chr8_24172      chr8    120857500       chr8    120167500       120202500       1     0.000111640587817074    Distance Cis     1       1
Necab3::chr2_10018      Necab3  chr2_10018      chr2    50087500        chr2    154382500       154407500       1     2.35032816456998e-05    Bin      Cis     1       1
Necab3::chr2_31558      Necab3  chr2_31558      chr2    157787500       chr2    154382500       154407500       1     0.000117516408228499    Distance Cis     1       1
...
```

The columns hold the following information: 


| Column        | Description                                                                                                                                 |
|---------------|---------------------------------------------------------------------------------------------------------------------------------------------|
| InteractionID | Unique interaction identifier                                                                                             |
| Symbol        | Gene name (HGNC symbol) of interacting RNA part                                                                                             |
| Bin           | Bin identifier of interacting DNA part<br>(format: BinChr_BinNum where BinNum is the enumerator of the 5kb DNA bin on the given chromosome) |
| BinChr        | Chromosome which the interacting DNA bin is located on                                                                                      |
| BinCentre     | Centre of interacting DNA bin (for distance calculation)                                                                                    |
| GeneChr       | Chromosome which the interacting RNA part is located on                                                                                     |
| GeneLeft      | Leftmost gene coordinate ("start") of interacting RNA part                                                                                  |
| GeneRight     | Rightmost gene coordinate ("end") of interacting RNA part                                                                                   |
| ReadCount     | Number of RNA-DNA reads supporting this interaction                                                                                         |
| ExpectedCount | Expected count modelled by background frequency x RNA abundance                                                                             |
| Method        | Method selected for background construction<br>(maximum value arising from distance-based or bin-based approach?)                           |
| P             | Poisson test significance value (ReadCount > ExpectedCount?)                                                                                |
| Padj          | Poisson test significance value after correction for multiple testing                                                                       |
|               |                                                                                                                                             |


-------------------------------------------------------------------------------
```
RADIAnT/results/RADICL/
├── ...
├── interactions 📂
    └── ...
    └── genes_venn.total_intra_trans.*
    └── ...
├── ...
```

The Venn diagram shows the proportions of genes with significant intrachromosomal interactions, with significant transchromsomal interactions, genes with both kinds of interactions, and genes simply annotated not interacting with chromatin in a significant manner. Counts can be found in the file ```*_genes.number_of_interactions.txt```, with the columns ```Symbol``` (gene name), ```Type``` (of count: transchromosomal or intrachromosomal), ```Total_Interactions``` (count) and ```Significant_Interactions``` (count). 

<p align="center">
  <img src="https://github.com/user-attachments/assets/7aadf29b-8717-45a0-96c0-cbf79893a8a9" width=50% height=50%>
</p>

-------------------------------------------------------------------------------
```
RADIAnT/results/RADICL/
├── ...
├── interactions 📂
    └── ...
    └── genes_barplot.all_interactions.*
    └── genes_barplot.all_reads.*
    └── genes_barplot.significant_interactions.*
    └── genes_barplot.significant_reads.*
    └── genes_barplot.significant_reads.*
    └── ...
├── ...
```


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

-------------------------------------------------------------------------------





## Individual GOI plots

While the above data gets generated automaticall with each RADIAnT run, two types of plots are supported which have to be customised by the user and produced running standalone R scripts. 

### Viewpoint plot: GOI-chromatin interaction at a given region of a genome

The viewpoint plot provides insight into significant interactions of a gene of interest (GOI) with a specific region of the genome. This is an example plot for the interactions of the lncRNA _Malat1_ with the chromatin around its own gene locus in mouse embryonic stem cells: 

![Viewpoint Malat1 chr19_3500000_1 1e+07](https://github.com/user-attachments/assets/5fc7d60e-0dae-4446-8765-2f03ced27af4)

Explanation: 
* The experiment-specific **background** is plotted in **grey** . It reflects the expected interaction frequency of the GOI with each position in the genome (either defined by bin or by distance to the gene locus) as calculated from the experiment, and scaled to fit the abundance of the GOI. The gap in the grey background is located at the gene body as intragenic interactions are excluded to filter out transcriptional clouding
* The **red** frequency plot shows the actual **observed interaction counts** of the GOI with each position in the genome.
* Significant interactions are highlighted by grey dots. The size of the dot indicates the "significance" of the interaction.

To create a viewpoint plot, run the R script provided with RADIAnT ```plot_GOI_viewpoint.r``` in the following manner:

```
Rscript /path/to/RADIAnT/workflow/scripts/plot_GOI_viewpoint.r --results /path/to/my_RADIAnT_results.txt --goi myGOI --binAnnotation /path/to/RADIAnT/resources/myorganism/my_bins_named.bed.gz --chr chromosomeName --start coordinate --end coordinate 
```

The script takes the following parameters: 


| Parameter       | Description                                                                                            |
|-----------------|--------------------------------------------------------------------------------------------------------|
| --results       | Absolute path to RADIAnT results file (*_RADIAnT_results.txt                                           |
| --goi           | Gene symbol of gene of interest                                                                        |
| --genome        | Genome version (currently supported: hg38 and mm39)                                                    |
| --binAnnotation | Absolute path to the bin annotation file. Also provided in ```resources/```.                           |
| --chr           | Chromosome name of the region of the genome of which interactions with the GOI should be plotted for.  |
| --start         | Start coordinate of the region to be plotted.                                                          |
| --end           | End coordinate of the region to be plotted.                                                            |
| --outdir        | Optional parameter. Desired output directory. Default is the directory plot_GOI_viewpoint.r is run in. |
| --outformat     | Optional parameter. Desired output format. Supported: svg / png / both. Default: both.                 |


### Chromosome plot: visualise how a GOI interacts with different chromosomes

This plot provides insight into interactions sites of a gene of interest (GOI) across the whole genome. The chromosomes are visualised in a circular manner. This is an example plot for significant interactions of the lncRNA _NEAT1_ across chromosomes of the human genome (HUVEC data set):

![NEAT1 hg38 chromosome_plot](https://github.com/user-attachments/assets/f6c60c5d-79d9-4d4f-82a2-3b0211edc8fe)


## Common issues

### Error: FATAL INPUT ERROR: unrecognized parameter name "genomeType" in input "genomeParameters.txt"

Please make sure that the mapping is run with a STAR version which is **equal to or newer than** the one used to generate the index. You can check your STAR version by running ```/path/to/STAR --version```. You can check the STAR version used to generate the index by checking the file at ```/path/to/STAR_index/Log.out``` (first line). 




 

