# Snakemake for processing of RNA-DNA interaction data from split fastq files to Gene-Bin interaction counts

# TO ADD?:
# FASTQ deduplication? (CZID) ideally performed by the user on unsplit reads though
# FASTQ quality trimming (trimmomatic), actually not required because STAR employs soft clipping of read ends
# Add interaction calling Rscript from SZ (probably need a list of required R packages for user to pre-install)
# Logs/metrics for Sankey plot and general quality overview
# Multi-mapping handling for Red-C data
# Blacklist regions excluded from counting/alignment for Red-C
# Limit to only canonical chromosomes sooner in the workflow
# Multi-mapping as a parameter in config file
# Multi-mapping for DNA reads?

# Config file

#configfile: "config.yaml"


# script directory
workflow_dir = config["workflow_directory"] if config["workflow_directory"].endswith("/") else config["workflow_directory"] + "/"



# method

method = config["method"][0]

# Directory holding the fastq files of the experiment to be analysed
fq_dir = config["fastq_directory"] if config["fastq_directory"].endswith("/") else config["fastq_directory"] + "/"


# Samples basename

samples = config["sample_base"]

print(samples)

# Mapping 

star_index = config["star_index"] if config["star_index"].endswith("/") else config["star_index"] + "/"
gtf = config["gtf"]
genome = config["genome_fasta"]

# output directory

outdir_base = config["output_directory"] if config["output_directory"].endswith("/") else config["output_directory"] + "/"

outdir_bam = outdir_base + "bam/"

outdir_bw = outdir_base + "bw/"

outdir_fastq = outdir_base + "fastq/"

outdir_intersects = outdir_base + "intersects/"

outdir_merge = outdir_base + "merge/"

outdir_counts = outdir_base + "counts/"

outdir_interactions = outdir_base + "interactions/"

outdir_logs = outdir_base + "logs/"

# all

if method == 'Red-C':

    rule all:
        input:
            #expand("{sample}DNA_dedup.bam", sample=config["sample_base"]),
            #expand("{sample}RNA_dedup.bam", sample=config["sample_base"]),
            #expand(outdir_bam + "{sample}DNA_sorted.bam", sample=config["sample_base"]),
            #expand(outdir_bam + "{sample}RNA_sorted.bam", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}DNA_bin_intersect.txt", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}RNA5_sorted.bed", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}RNA3_sorted.bed", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}RNA5_gene_intersect.txt", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}RNA3_gene_intersect.txt", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}RNA5_gene_intersect_maximums.txt", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}RNA3_gene_intersect_maximums.txt", sample=config["sample_base"]),
            #expand("{sample}DNA_dedup_cpm.bw", sample=config["sample_base"]),
            #expand("{sample}RNA_dedup_cpm.bw", sample=config["sample_base"]),
            #expand(outdir_bw + "{sample}DNA_sorted_cpm.bw", sample=config["sample_base"]),
            #expand(outdir_bw + "{sample}RNA_sorted_cpm.bw", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}DNA_bin_intersect_maximums.txt", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}RNA_gene_intersect_maximums.txt", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}RNA5_gene_intersect_maximums_cut.txt", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}RNA3_gene_intersect_maximums_cut.txt", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}DNA_bin_intersect_maximums_cut.txt", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}RNA5_gene_unstranded_intersect.txt", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}RNA3_gene_unstranded_intersect.txt", sample=config["sample_base"]),
            #expand(outdir_counts + "{sample}RNA-bin_counts.txt", sample=config["sample_base"]),
            expand(outdir_interactions + "{sample}RADIAnT_results.txt", sample=config["sample_base"]),
            expand(outdir_logs + "{sample}sankey.svg", sample=config["sample_base"]), 
            expand(outdir_logs + "{sample}sankey.png", sample=config["sample_base"]), 
            expand(outdir_logs + "{sample}read_stats.txt", sample=config["sample_base"]),
            expand(outdir_interactions + "{sample}genes.number_of_interactions.txt", sample=config["sample_base"]),
            star_index + "Log.out"

else:
    rule all:
        input:
            #expand("{sample}DNA_dedup.bam", sample=config["sample_base"]),
            #expand("{sample}RNA_dedup.bam", sample=config["sample_base"]),
            #expand(outdir_bam + "{sample}DNA_sorted.bam", sample=config["sample_base"]),
            #expand(outdir_bam + "{sample}RNA_sorted.bam", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}DNA_bin_intersect.txt", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}RNA_gene_intersect.txt", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}RNA_multi_gene_intersect.txt", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}RNA_gene_intersect_unstranded.txt", sample=config["sample_base"]),
            #expand("{sample}DNA_dedup_cpm.bw", sample=config["sample_base"]),
            #expand("{sample}RNA_dedup_cpm.bw", sample=config["sample_base"]),
            #expand(outdir_bw + "{sample}DNA_sorted_cpm.bw", sample=config["sample_base"]),
            #expand(outdir_bw + "{sample}RNA_sorted_cpm.bw", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}DNA_bin_intersect_maximums.txt", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}RNA_gene_intersect_maximums.txt", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}RNA_gene_intersect_maximums_cut.txt", sample=config["sample_base"]),
            #expand(outdir_intersects + "{sample}DNA_bin_intersect_maximums_cut.txt", sample=config["sample_base"]),
            #expand(outdir_counts + "{sample}RNA-bin_counts.txt", sample=config["sample_base"]),
            expand(outdir_interactions + "{sample}RADIAnT_results.txt", sample=config["sample_base"]),
            expand(outdir_logs + "{sample}sankey.svg", sample=config["sample_base"]), 
            expand(outdir_logs + "{sample}sankey.png", sample=config["sample_base"]), 
            expand(outdir_logs + "{sample}read_stats.txt", sample=config["sample_base"]),
            expand(outdir_interactions + "{sample}genes.number_of_interactions.txt", sample=config["sample_base"]),
            star_index + "Log.out"


# If no STAR index provided, decompress provided genome FASTA to build index
rule gunzip_genome_fasta:
    input: 
        genome_fasta = config["genome_fasta"] if config["genome_fasta"].endswith(".gz") else config["genome_fasta"] + ".gz"
    output: 
        decompressed_fasta = temporary(config["genome_fasta"][:-3] if config["genome_fasta"].endswith(".gz") else config["genome_fasta"])
    run: 
        shell("pigz -k -d -p {threads} {input.genome_fasta}")

# If no STAR index provided, decompress provided genome annotation to build index
rule gunzip_gtf:
    input: 
        gtf = config["gtf"] if config["gtf"].endswith(".gz") else config["gtf"] + ".gz"
    output: 
        decompressed_gtf = temporary(config["gtf"][:-3] if config["gtf"].endswith(".gz") else config["gtf"])
    run: 
        shell("pigz -k -d -p {threads} {input.gtf}")

# If no STAR index provided, build index
rule build_star_index: 
    input:
        gtf = re.sub(r"\.gz$", "", config["gtf"]),
        genome_fasta = re.sub(r"\.gz$", "", config["genome_fasta"])
    params:
        star_binary = config["star_binary"],
        star_index = config["star_index"]
    output: 
        star_log = star_index + "Log.out"
        # star_done = star_index + "index.done"
    run:
        shell("{params.star_binary} \
        --runMode genomeGenerate \
        --genomeDir {params.star_index} \
        --genomeFastaFiles {input.genome_fasta} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 50")
        # shell("touch {output.star_done}")


# Decompress FASTQs
rule gunzip_dna:
    input:
        dna_gz = fq_dir + "{sample}" + (config["dna_fastq_suffix"] if config["dna_fastq_suffix"].endswith(".gz") else config["dna_fastq_suffix"] + ".gz")
    threads:
        config["threads"]
    output:
        dna_fastq = temporary(fq_dir + "{sample}" + (config["dna_fastq_suffix"][:-3] if config["dna_fastq_suffix"].endswith(".gz") else config["dna_fastq_suffix"]))
    run:
        shell("pigz -k -d -p {threads} {input.dna_gz}")



# DNA alignment


rule align_dna:
    input:
        dna_fastq = fq_dir + "{sample}"+config["dna_fastq_suffix"],
        genome_parameters = star_index + "Log.out"
    threads:
        config["threads"]
    params:
        # star_index = config["star_index"] if os.path.exists(config["star_index"]) else outdir_base + "resources/" + config["species"] + "/star_index",
        star_index = config["star_index"], 
        star_binary = config["star_binary"],
        base_name = outdir_bam + "{sample}" + 'DNA_'
    output:
        aligned_dna = outdir_bam + "{sample}DNA_Aligned.out.bam"
    run:
        shell("{params.star_binary} \
               --runThreadN {threads} \
               --genomeDir {params.star_index} \
               --genomeLoad NoSharedMemory \
               --limitBAMsortRAM 30000000000 \
               --readFilesIn {input.dna_fastq} \
               --outFileNamePrefix {params.base_name} \
               --outSAMtype BAM Unsorted \
               --alignIntronMax 1 \
               --alignMatesGapMax 1 \
               --outFilterScoreMinOverLread 0 \
               --outFilterMatchNminOverLread 0 \
               --outFilterMatchNmin 0")


# rule align_dna:
#     input:
#         dna_fastq = fq_dir + "{sample}"+config["dna_fastq_suffix"]
#     threads:
#         config["threads"]
#     params:
#         star_binary = config["star_binary"],
#         star_index = config["star_index"],
#         base_name = outdir_bam + "{sample}" + 'DNA_'
#     output:
#         aligned_dna = outdir_bam + "{sample}DNA_Aligned.out.bam"
#     run:
#         shell("{params.star_binary} \
#                --runThreadN {threads} \
#                --genomeDir {params.star_index} \
#                --genomeLoad NoSharedMemory \
#                --limitBAMsortRAM 30000000000 \
#                --readFilesIn {input.dna_fastq} \
#                --outFileNamePrefix {params.base_name} \
#                --outSAMtype BAM Unsorted \
#                --alignIntronMax 1 \
#                --alignMatesGapMax 1 \
#                --outFilterScoreMinOverLread 0 \
#                --outFilterMatchNminOverLread 0 \
#                --outFilterMatchNmin 0")

# Remove blacklisted regions from DNA

rule no_blacklist_dna:
    input:
        aligned_dna = outdir_bam + "{sample}DNA_Aligned.out.bam"
    params:
        blacklist = config["blacklist"]
    output:
        no_blacklist_dna = temporary(outdir_bam + "{sample}DNA_Aligned.out.bl.bam")
    run:
        shell("bedtools intersect -v -a {input.aligned_dna} -b {params.blacklist} > {output.no_blacklist_dna}")


# extract uniquely mapping reads

rule unique_dna:
    input:
        no_blacklist_dna = outdir_bam + "{sample}DNA_Aligned.out.bl.bam"
    threads:
        config["threads"]
    params:
        samtools_binary = config["samtools_binary"]
    output:
        unique_dna = temporary(outdir_bam + "{sample}DNA_unique.bam")
    run:
        shell("{params.samtools_binary} view -@ {threads} -q 255 -o {output.unique_dna} {input.no_blacklist_dna}")

# Collate bam (samtools collate)

rule collate_dna:
    input:
        unique_dna = outdir_bam + "{sample}DNA_unique.bam"
    threads:
        config["threads"]
    params:
        samtools_binary = config["samtools_binary"]
    output:
        collated_dna = temporary(outdir_bam + "{sample}DNA_collated.bam")
    run:
        shell("{params.samtools_binary} collate -@ {threads} -o {output.collated_dna} {input.unique_dna}")

# Fixmate (samtool fixmate)

rule fixmate_dna:
    input:
        collated_dna = outdir_bam + "{sample}DNA_collated.bam"
    threads:
        config["threads"]
    params:
        samtools_binary = config["samtools_binary"]
    output:
        fixmate_dna = temporary(outdir_bam + "{sample}DNA_fixmate.bam")
    run:
        shell("{params.samtools_binary} fixmate -@ {threads} -m {input.collated_dna} {output.fixmate_dna}")


# Sort by coordinate (samtools sort)

rule sort_dna:
    input:
        fixmate_dna = outdir_bam + "{sample}DNA_fixmate.bam"
    threads:
        config["threads"]
    params:
        samtools_binary = config["samtools_binary"]
    output:
        coord_sorted_dna = outdir_bam + "{sample}DNA_sorted.bam"
    run:
        shell("{params.samtools_binary} sort -@ {threads} -o {output.coord_sorted_dna} {input.fixmate_dna}")

# remove duplicates (samtools markdup -r)

# rule dedup_dna:
#     input:
#         coord_sorted_dna = "{sample}DNA_sorted.bam"
#     threads:
#         config["threads"]
#     output:
#         dedup_dna = "{sample}DNA_dedup.bam"
#     run:
#         shell("samtools markdup -@ {threads} -r {input.coord_sorted_dna} {output.dedup_dna}")

# index deduplicated DNA bam

# rule index_dna_bam:
#     input:
#         dedup_dna = "{sample}DNA_dedup.bam"
#     threads:
#         config["threads"]
#     output:
#         dedup_dna_index = "{sample}DNA_dedup.bam.bai"
#     run:
#         shell("samtools index -@ {threads} {input.dedup_dna}")

rule index_dna_bam:
    input:
        dedup_dna = outdir_bam + "{sample}DNA_sorted.bam"
    threads:
        config["threads"]
    params:
        samtools_binary = config["samtools_binary"]
    output:
        dedup_dna_index = outdir_bam + "{sample}DNA_sorted.bam.bai"
    run:
        shell("{params.samtools_binary} index -@ {threads} {input.dedup_dna}")

# get genome coverage for DNA

# rule dna_coverage:
#     input:
#         dedup_dna = "{sample}DNA_dedup.bam",
#         dedup_dna_index = "{sample}DNA_dedup.bam.bai"
#     threads:
#         config["threads"]
#     params:
#         blacklist = config["blacklist"]
#     output:
#         dedup_dna_bw = "{sample}DNA_dedup_cpm.bw"
#     run:
#         shell("bamCoverage --bam {input.dedup_dna} \
#                -o {output.dedup_dna_bw} \
#                -of bigwig \
#                -bs 1 \
#                --blackListFileName {params.blacklist} \
#                --normalizeUsing CPM \
#                -p {threads}")

rule dna_coverage:
    input:
        dedup_dna = outdir_bam + "{sample}DNA_sorted.bam",
        dedup_dna_index = outdir_bam + "{sample}DNA_sorted.bam.bai"
    threads:
        config["threads"]
    params:
        bamCoverage_binary = config["bamCoverage_binary"],
        blacklist = config["blacklist"]
    output:
        dedup_dna_bw = outdir_bw + "{sample}DNA_sorted_cpm.bw"
    run:
        shell("{params.bamCoverage_binary} --bam {input.dedup_dna} \
               -o {output.dedup_dna_bw} \
               -of bigwig \
               -bs 1 \
               --blackListFileName {params.blacklist} \
               --normalizeUsing CPM \
               -p {threads}")

# intersect DNA and bins

# rule dna_bin_intersect:
#     input:
#         dedup_dna = "{sample}DNA_dedup.bam"
#     params:
#         genome_bins = config["genome_bins"]
#     output:
#         dna_bin_intersect = "{sample}DNA_bin_intersect.txt"
#     run:
#         shell("bedtools intersect -bed -wo -a {input.dedup_dna} -b {params.genome_bins} > {output.dna_bin_intersect}")

rule dna_bin_intersect:
    input:
        dedup_dna = outdir_bam + "{sample}DNA_sorted.bam"
    params:
        genome_bins = config["genome_bins"],
        bedtools_binary = config["bedtools_binary"]
    output:
        dna_bin_intersect = outdir_intersects + "{sample}DNA_bin_intersect.txt"
    run:
        shell("{params.bedtools_binary} intersect -bed -wo -a {input.dedup_dna} -b {params.genome_bins} > {output.dna_bin_intersect}")

# Add intersect as proportion of gene width column to RNA-gene intersect

rule intersect_DNA_proportion:
    input:
        dna_bin_intersect = outdir_intersects + "{sample}DNA_bin_intersect.txt"
    output:
        bin_width = outdir_intersects + "{sample}DNA_bin_intersect_proportions.txt"
    run:
        shell("awk 'BEGIN{{OFS=\"\t\"}} {{$18 = $17/($15-$14); print}}' {input.dna_bin_intersect} > {output.bin_width}")

# select maximum proportion alignments per read

rule maximum_DNA_proportion:
    input:
        dna_bin_proportions = outdir_intersects + "{sample}DNA_bin_intersect_proportions.txt"
    params:
        bedtools_binary = config["bedtools_binary"]
    output:
        dna_bin_maximums = outdir_intersects + "{sample}DNA_bin_intersect_maximums_cut.txt"
    run:
        shell("sort -r -k 4,4 -k 18,18 {input.dna_bin_proportions} | {params.bedtools_binary} groupby -g 4 -c 18 -o first -full | cut -f 4,16 | sort -k 1,1 > {output.dna_bin_maximums}")

# cut DNA intersection to only read and bin for joining with RNA

# rule cut_and_sort_DNA:
#     input:
#         dna_bin_maximums = outdir_intersects + "{sample}DNA_bin_intersect_maximums.txt"
#     output:
#         dna_bin_cut = outdir_intersects + "{sample}DNA_bin_intersect_maximums_cut.txt"
#     run:
#         shell("cut -f 4,16 {input.dna_bin_maximums} | sort -k 4,4 > {output.dna_bin_cut}")

# RNA processing for Red-C

if method == "Red-C":

    # 5' RNA (actually 3' I think)
    
    rule gunzip_rna5:
        input:
            rna_gz = fq_dir + "{sample}" + config["rna_fastq_suffix"] + ".gz"
        threads:
            config["threads"]
        output:
            rna_fastq = temporary(fq_dir + "{sample}" + config["rna_fastq_suffix"])
        run:
            shell("pigz -k -d -p {threads} {input.rna_gz}")
    
    rule remove_rrna5:
        input:
            rna_fastq = fq_dir + "{sample}"+config["rna_fastq_suffix"]
        threads:
            config["threads"]
        params:
            rrna_index = config["rrna_fasta"],
            bbduk = config["bbduk_script"]
        output:
            rna_clean = outdir_fastq + "{sample}RNA5_depleted.fastq",
            rna_ribo = outdir_fastq + "{sample}rRNA5.fastq",
            stats= outdir_fastq + "{sample}rRNA5_removal_stats.txt"
        run:
            shell("{params.bbduk} in={input.rna_fastq} out={output.rna_clean} outm={output.rna_ribo} ref={params.rrna_index} k=13 hdist=1 stats={output.stats}")
            #shell("ribodetector_cpu -t {threads} -i {input.rna_fastq} -l 21 -e rrna --chunk_size 256 -o {output.rna_clean}")
    
    # rule remove_rrna_rna5:
    #     input:
    #         rna_fastq = fq_dir + "{sample}"+config["rna_fastq_suffix"]
    #     threads:
    #         config["threads"]
    #     output:
    #         rna_clean = outdir_fastq + "{sample}RNA5_depleted.fastq"
    #     run:
    #         shell("ribodetector_cpu -t {threads} -i {input.rna_fastq} -l 30 -e rrna --chunk_size 256 -o {output.rna_clean}")
    
    rule align_rna5:
        input:
            rna_clean = outdir_fastq + "{sample}RNA5_depleted.fastq"
        threads:
            config["threads"]
        params:
            star_binary = config["star_binary"],
            star_index = config["star_index"],
            base_name = outdir_bam + "{sample}"+'RNA5_'
        output:
            aligned_rna = outdir_bam + "{sample}RNA5_Aligned.out.bam",
            rna_map_log = outdir_bam + "{sample}RNA5_Log.final.out"
        run:
            shell("{params.star_binary} \
                --runThreadN {threads} \
                --genomeDir {params.star_index} \
                --genomeLoad NoSharedMemory \
                --limitBAMsortRAM 30000000000 \
                --readFilesIn {input.rna_clean} \
                --outFileNamePrefix {params.base_name} \
                --outSAMtype BAM Unsorted \
                --alignIntronMax 1 \
                --alignMatesGapMax 1 \
                --outFilterScoreMinOverLread 0 \
                --outFilterMatchNminOverLread 0 \
                --outFilterMatchNmin 0")

    rule unique_rna5:
        input:
            aligned_rna = outdir_bam + "{sample}RNA5_Aligned.out.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            unique_rna = temporary(outdir_bam + "{sample}RNA5_unique.bam"),
            multi_rna = temporary(outdir_bam + "{sample}RNA5_multi.bam")
        run:
            shell("{params.samtools_binary} view -@ {threads} -q 255 -U {output.multi_rna} -o {output.unique_rna} {input.aligned_rna}")

    # Collate bam (samtools collate)

    rule collate_rna5:
        input:
            unique_rna = outdir_bam + "{sample}RNA5_unique.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            collated_rna = temporary(outdir_bam + "{sample}RNA5_collated.bam")
        run:
            shell("{params.samtools_binary} collate -@ {threads} -o {output.collated_rna} {input.unique_rna}")

    rule collate_multi_rna5:
        input:
            unique_rna = outdir_bam + "{sample}RNA5_multi.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            collated_rna = temporary(outdir_bam + "{sample}RNA5_multi_collated.bam")
        run:
            shell("{params.samtools_binary} collate -@ {threads} -o {output.collated_rna} {input.unique_rna}")

    # Fixmate (samtool fixmate)

    rule fixmate_rna5:
        input:
            collated_rna = outdir_bam + "{sample}RNA5_collated.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            fixmate_rna = temporary(outdir_bam + "{sample}RNA5_fixmate.bam")
        run:
            shell("{params.samtools_binary} fixmate -@ {threads} -m {input.collated_rna} {output.fixmate_rna}")

    rule fixmate_multi_rna5:
        input:
            collated_rna = outdir_bam + "{sample}RNA5_multi_collated.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            fixmate_rna = temporary(outdir_bam + "{sample}RNA5_multi_fixmate.bam")
        run:
            shell("{params.samtools_binary} fixmate -@ {threads} -m {input.collated_rna} {output.fixmate_rna}")

    # Sort by coordinate (samtools sort)

    rule sort_rna5:
        input:
            fixmate_rna = outdir_bam + "{sample}RNA5_fixmate.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            coord_sorted_rna = outdir_bam + "{sample}RNA5_sorted.bam"
        run:
            shell("{params.samtools_binary} sort -@ {threads} -o {output.coord_sorted_rna} {input.fixmate_rna}")

    rule sort_multi_rna5:
        input:
            fixmate_rna = outdir_bam + "{sample}RNA5_multi_fixmate.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            coord_sorted_rna = outdir_bam + "{sample}RNA5_multi_sorted.bam"
        run:
            shell("{params.samtools_binary} sort -@ {threads} -o {output.coord_sorted_rna} {input.fixmate_rna}")

    # Index bam file

    rule index_rna5_bam:
        input:
            dedup_rna = outdir_bam + "{sample}RNA5_sorted.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            dedup_rna_index = outdir_bam + "{sample}RNA5_sorted.bam.bai"
        run:
            shell("{params.samtools_binary} index -@ {threads} {input.dedup_rna}")

    rule index_multi_rna5_bam:
        input:
            dedup_rna = outdir_bam + "{sample}RNA5_multi_sorted.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            dedup_rna_index = outdir_bam + "{sample}RNA5_multi_sorted.bam.bai"
        run:
            shell("{params.samtools_binary} index -@ {threads} {input.dedup_rna}")

    # Intersect RNA alignments with gene loci

    rule rna5_gene_intersect:
        input:
            dedup_rna = outdir_bam + "{sample}RNA5_sorted.bam"
        params:
            genes = config["genes"],
            bedtools_binary = config["bedtools_binary"]
        output:
            rna_gene_intersect = outdir_intersects + "{sample}RNA5_gene_intersect.txt"
        run:
            shell("{params.bedtools_binary} intersect -S -bed -wo -a {input.dedup_rna} -b {params.genes} > {output.rna_gene_intersect}")

    rule rna5_multi_gene_intersect:
        input:
            dedup_rna = outdir_bam + "{sample}RNA5_multi_sorted.bam"
        params:
            genes = config["genes"],
            bedtools_binary = config["bedtools_binary"]
        output:
            rna_gene_intersect = outdir_intersects + "{sample}RNA5_multi_gene_intersect.txt"
        run:
            shell("{params.bedtools_binary} intersect -S -bed -wo -a {input.dedup_rna} -b {params.genes} > {output.rna_gene_intersect}")

    rule rna5_unstranded_gene_intersect:
        input:
            dedup_rna = outdir_bam + "{sample}RNA5_sorted.bam"
        params:
            genes = config["genes"],
            bedtools_binary = config["bedtools_binary"]
        output:
            rna_gene_unstranded_intersect = outdir_intersects + "{sample}RNA5_gene_unstranded_intersect.txt"
        run:
            shell("{params.bedtools_binary} intersect -bed -wo -a {input.dedup_rna} -b {params.genes} > {output.rna_gene_unstranded_intersect}")

    rule rna5_multi_unstranded_gene_intersect:
        input:
            dedup_rna = outdir_bam + "{sample}RNA5_multi_sorted.bam"
        params:
            genes = config["genes"],
            bedtools_binary = config["bedtools_binary"]
        output:
            rna_gene_unstranded_intersect = outdir_intersects + "{sample}RNA5_gene_multi_unstranded_intersect.txt"
        run:
            shell("{params.bedtools_binary} intersect -bed -wo -a {input.dedup_rna} -b {params.genes} > {output.rna_gene_unstranded_intersect}")

    # rule process_intersections_rna5_unstranded:
    #     input:
    #         unique_intersect = outdir_intersects + "{sample}RNA5_gene_unstranded_intersect.txt",
    #         multi_intersect = outdir_intersects + "{sample}RNA5_gene_multi_unstranded_intersect.txt"
    #     output:
    #         maximum_intersects = outdir_intersects + "{sample}RNA5_gene_unstranded_intersect_maximums_cut.txt"
    #     run:
    #         shell("Rscript /mnt/y/Fileserver_NGS/Current/RNA_DNA_interactions/Method/snakemake/utility_scripts/intersect_processing.R --unique {input.unique_intersect} --multi {input.multi_intersect} --output {output.maximum_intersects}")

    # rule process_intersections_rna5:
    #     input:
    #         unique_intersect = outdir_intersects + "{sample}RNA5_gene_intersect.txt",
    #         multi_intersect = outdir_intersects + "{sample}RNA5_multi_gene_intersect.txt"
    #     output:
    #         maximum_intersects = outdir_intersects + "{sample}RNA5_gene_intersect_maximums_cut.txt"
    #     run:
    #         shell("Rscript /mnt/y/Fileserver_NGS/Current/RNA_DNA_interactions/Method/snakemake/utility_scripts/intersect_processing.R --unique {input.unique_intersect} --multi {input.multi_intersect} --output {output.maximum_intersects}")


    # 3' RNA (actually 5' I think)

    rule gunzip_rna3:
        input:
            rna_gz = fq_dir + "{sample}" + config["rna_2_fastq_suffix"] + ".gz"
        threads:
            config["threads"]
        output:
            rna_fastq = temporary(fq_dir + "{sample}" + config["rna_2_fastq_suffix"])
        run:
            shell("pigz -k -d -p {threads} {input.rna_gz}")
    
    rule remove_rrna3:
        input:
            rna_fastq = fq_dir + "{sample}"+config["rna_2_fastq_suffix"]
        threads:
            config["threads"]
        params:
            rrna_index = config["rrna_fasta"],
            bbduk = config["bbduk_script"]
        output:
            rna_clean = outdir_fastq + "{sample}RNA3_depleted.fastq",
            rna_ribo = outdir_fastq + "{sample}rRNA3.fastq",
            stats= outdir_fastq + "{sample}rRNA3_removal_stats.txt"
        run:
            shell("{params.bbduk} in={input.rna_fastq} out={output.rna_clean} outm={output.rna_ribo} ref={params.rrna_index} k=13 hdist=1 stats={output.stats}")
            #shell("ribodetector_cpu -t {threads} -i {input.rna_fastq} -l 21 -e rrna --chunk_size 256 -o {output.rna_clean}")
    
    # rule remove_rrna_rna3:
    #     input:
    #         rna_fastq = fq_dir + "{sample}"+config["rna_2_fastq_suffix"]
    #     threads:
    #         config["threads"]
    #     output:
    #         rna_clean = outdir_fastq + "{sample}RNA3_depleted.fastq"
    #     run:
    #         shell("ribodetector_cpu -t {threads} -i {input.rna_fastq} -l 30 -e rrna --chunk_size 256 -o {output.rna_clean}")
    
    rule align_rna3:
        input:
            rna_clean = outdir_fastq + "{sample}RNA3_depleted.fastq"
        threads:
            config["threads"]
        params:
            star_binary = config["star_binary"],
            star_index = config["star_index"],
            base_name = outdir_bam + "{sample}"+'RNA3_'
        output:
            aligned_rna = outdir_bam + "{sample}RNA3_Aligned.out.bam", 
            rna_map_log = outdir_bam + "{sample}RNA3_Log.final.out"
        run:
            shell("{params.star_binary} \
                --runThreadN {threads} \
                --genomeDir {params.star_index} \
                --genomeLoad NoSharedMemory \
                --limitBAMsortRAM 30000000000 \
                --readFilesIn {input.rna_clean} \
                --outFileNamePrefix {params.base_name} \
                --outSAMtype BAM Unsorted \
                --alignIntronMax 1 \
                --alignMatesGapMax 1 \
                --outFilterScoreMinOverLread 0 \
                --outFilterMatchNminOverLread 0 \
                --outFilterMatchNmin 0")

    rule unique_rna3:
        input:
            aligned_rna = outdir_bam + "{sample}RNA3_Aligned.out.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            unique_rna = temporary(outdir_bam + "{sample}RNA3_unique.bam"),
            multi_rna = temporary(outdir_bam + "{sample}RNA3_multi.bam")
        run:
            shell("{params.samtools_binary} view -@ {threads} -q 255 -U {output.multi_rna} -o {output.unique_rna} {input.aligned_rna}")

    # Collate bam (samtools collate)

    rule collate_rna3:
        input:
            unique_rna = outdir_bam + "{sample}RNA3_unique.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            collated_rna = temporary(outdir_bam + "{sample}RNA3_collated.bam")
        run:
            shell("{params.samtools_binary} collate -@ {threads} -o {output.collated_rna} {input.unique_rna}")

    rule collate_multi_rna3:
        input:
            unique_rna = outdir_bam + "{sample}RNA3_multi.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            collated_rna = temporary(outdir_bam + "{sample}RNA3_multi_collated.bam")
        run:
            shell("{params.samtools_binary} collate -@ {threads} -o {output.collated_rna} {input.unique_rna}")

    # Fixmate (samtool fixmate)

    rule fixmate_rna3:
        input:
            collated_rna = outdir_bam + "{sample}RNA3_collated.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            fixmate_rna = temporary(outdir_bam + "{sample}RNA3_fixmate.bam")
        run:
            shell("{params.samtools_binary} fixmate -@ {threads} -m {input.collated_rna} {output.fixmate_rna}")

    rule fixmate_multi_rna3:
        input:
            collated_rna = outdir_bam + "{sample}RNA3_multi_collated.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            fixmate_rna = temporary(outdir_bam + "{sample}RNA3_multi_fixmate.bam")
        run:
            shell("{params.samtools_binary} fixmate -@ {threads} -m {input.collated_rna} {output.fixmate_rna}")

    # Sort by coordinate (samtools sort)

    rule sort_rna3:
        input:
            fixmate_rna = outdir_bam + "{sample}RNA3_fixmate.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            coord_sorted_rna = outdir_bam + "{sample}RNA3_sorted.bam"
        run:
            shell("{params.samtools_binary} sort -@ {threads} -o {output.coord_sorted_rna} {input.fixmate_rna}")

    rule sort_multi_rna3:
        input:
            fixmate_rna = outdir_bam + "{sample}RNA3_multi_fixmate.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            coord_sorted_rna = outdir_bam + "{sample}RNA3_multi_sorted.bam"
        run:
            shell("{params.samtools_binary} sort -@ {threads} -o {output.coord_sorted_rna} {input.fixmate_rna}")

    # Index bam file

    rule index_rna3_bam:
        input:
            dedup_rna = outdir_bam + "{sample}RNA3_sorted.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            dedup_rna_index = outdir_bam + "{sample}RNA3_sorted.bam.bai"
        run:
            shell("{params.samtools_binary} index -@ {threads} {input.dedup_rna}")

    rule index_multi_rna3_bam:
        input:
            dedup_rna = outdir_bam + "{sample}RNA3_multi_sorted.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            dedup_rna_index = outdir_bam + "{sample}RNA3_multi_sorted.bam.bai"
        run:
            shell("{params.samtools_binary} index -@ {threads} {input.dedup_rna}")

    # Intersect RNA alignments with gene loci

    rule rna3_gene_intersect:
        input:
            dedup_rna = outdir_bam + "{sample}RNA3_sorted.bam"
        params:
            genes = config["genes"],
            bedtools_binary = config["bedtools_binary"]
        output:
            rna_gene_intersect = outdir_intersects + "{sample}RNA3_gene_intersect.txt"
        run:
            shell("{params.bedtools_binary} intersect -s -bed -wo -a {input.dedup_rna} -b {params.genes} > {output.rna_gene_intersect}")

    rule rna3_multi_gene_intersect:
        input:
            dedup_rna = outdir_bam + "{sample}RNA3_multi_sorted.bam"
        params:
            genes = config["genes"],
            bedtools_binary = config["bedtools_binary"]
        output:
            rna_gene_intersect = outdir_intersects + "{sample}RNA3_multi_gene_intersect.txt"
        run:
            shell("{params.bedtools_binary} intersect -s -bed -wo -a {input.dedup_rna} -b {params.genes} > {output.rna_gene_intersect}")

    rule rna3_unstranded_gene_intersect:
        input:
            dedup_rna = outdir_bam + "{sample}RNA3_sorted.bam"
        params:
            genes = config["genes"],
            bedtools_binary = config["bedtools_binary"]
        output:
            rna_gene_unstranded_intersect = outdir_intersects + "{sample}RNA3_gene_unstranded_intersect.txt"
        run:
            shell("{params.bedtools_binary} intersect -bed -wo -a {input.dedup_rna} -b {params.genes} > {output.rna_gene_unstranded_intersect}")

    rule process_intersections_rna3:
        input:
            unique_intersect = outdir_intersects + "{sample}RNA3_gene_intersect.txt",
            multi_intersect = outdir_intersects + "{sample}RNA3_multi_gene_intersect.txt"
        output:
            maximum_intersects = outdir_intersects + "{sample}RNA3_gene_intersect_maximums_cut.txt"
        params: 
            intersect_processing = workflow_dir + "scripts/intersect_processing.R"
        run:
            shell("Rscript {params.intersect_processing} --unique {input.unique_intersect} --multi {input.multi_intersect} --output {output.maximum_intersects}")

    # # Add intersect as proportion of gene width column to RNA-gene intersect

    # rule intersect_rna5_proportion:
    #     input:
    #         rna_gene_intersect = outdir_intersects + "{sample}RNA5_gene_intersect.txt"
    #     output:
    #         gene_width = outdir_intersects + "{sample}RNA5_gene_intersect_proportions.txt"
    #     run:
    #         shell("awk 'BEGIN{{OFS=\"\t\"}} {{$20 = $19/($15-$14); print}}' {input.rna_gene_intersect} > {output.gene_width}")

    # # select maximum proportion alignments per read

    # rule maximum_rna5_proportion:
    #     input:
    #         rna_gene_proportions = outdir_intersects + "{sample}RNA5_gene_intersect_proportions.txt"
    #     output:
    #         rna_gene_maximums = outdir_intersects + "{sample}RNA5_gene_intersect_maximums.txt"
    #     run:
    #         shell("sort -r -k 4,4 -k 20,20 {input.rna_gene_proportions} | groupBy -g 4 -c 20 -o first -full | cut -f 4,17 | sort -k 1,1 > {output.rna_gene_maximums}")

    # # 3' RNA processing

    # rule gunzip_rna3:
    #     input:
    #         rna_gz = fq_dir + "{sample}" + config["rna_2_fastq_suffix"] + ".gz"
    #     threads:
    #         config["threads"]
    #     output:
    #         rna_fastq = temporary(fq_dir + "{sample}" + config["rna_2_fastq_suffix"])
    #     run:
    #         shell("pigz -k -d -p {threads} {input.rna_gz}")
    
    # rule remove_rrna3:
    #     input:
    #         rna_fastq = fq_dir + "{sample}"+config["rna_2_fastq_suffix"]
    #     threads:
    #         config["threads"]
    #     params:
    #         rrna_index = config["rrna_fasta"],
    #         bbduk = config["bbduk_script"]
    #     output:
    #         rna_clean = outdir_fastq + "{sample}RNA3_depleted.fastq",
    #         rna_ribo = outdir_fastq + "{sample}rRNA3.fastq",
    #         stats= outdir_fastq + "{sample}rRNA3_removal_stats.txt"
    #     run:
    #         shell("{params.bbduk} in={input.rna_fastq} out={output.rna_clean} outm={output.rna_ribo} ref={params.rrna_index} k=13 hdist=1 stats={output.stats}")
    #         #shell("ribodetector_cpu -t {threads} -i {input.rna_fastq} -l 21 -e rrna --chunk_size 256 -o {output.rna_clean}")

    # # rule remove_rrna_rna3:
    # #     input:
    # #         rna_fastq = fq_dir + "{sample}"+config["rna_2_fastq_suffix"]
    # #     threads:
    # #         config["threads"]
    # #     output:
    # #         rna_clean = outdir_fastq + "{sample}RNA3_depleted.fastq"
    # #     run:
    # #         shell("ribodetector_cpu -t {threads} -i {input.rna_fastq} -l 20 -e rrna --chunk_size 4096 -o {output.rna_clean}")

    # rule align_rna3:
    #     input:
    #         rna_clean = outdir_fastq + "{sample}RNA3_depleted.fastq"
    #     threads:
    #         config["threads"]
    #     params:
    #         star_binary = config["star_binary"],
    #         star_index = config["star_index"],
    #         base_name = outdir_bam + "{sample}"+'RNA3_'
    #     output:
    #         aligned_rna = outdir_bam + "{sample}RNA3_Aligned.out.bam"
    #     run:
    #         shell("{config.star_binary} \
    #             --runThreadN {threads} \
    #             --genomeDir {params.star_index} \
    #             --genomeLoad NoSharedMemory \
    #             --limitBAMsortRAM 30000000000 \
    #             --readFilesIn {input.rna_clean} \
    #             --outFileNamePrefix {params.base_name} \
    #             --outSAMtype BAM Unsorted \
    #             --alignIntronMax 1 \
    #             --alignMatesGapMax 1 \
    #             --outFilterScoreMinOverLread 0 \
    #             --outFilterMatchNminOverLread 0 \
    #             --outFilterMatchNmin 0")
            
    # rule unique_rna3:
    #     input:
    #         aligned_rna = outdir_bam + "{sample}RNA3_Aligned.out.bam"
    #     threads:
    #         config["threads"]
    #     params:
    #         samtools_binary = config["samtools_binary"]
    #     output:
    #         unique_rna = temporary(outdir_bam + "{sample}RNA3_unique.bam")
    #     run:
    #         shell("{params.samtools_binary} view -@ {threads} -q 255 -o {output.unique_rna} {input.aligned_rna}")

    # # Collate bam (samtools collate)

    # rule collate_rna3:
    #     input:
    #         unique_rna = outdir_bam + "{sample}RNA3_unique.bam"
    #     threads:
    #         config["threads"]
    #     params:
    #         samtools_binary = config["samtools_binary"]
    #     output:
    #         collated_rna = temporary(outdir_bam + "{sample}RNA3_collated.bam")
    #     run:
    #         shell("{params.samtools_binary} collate -@ {threads} -o {output.collated_rna} {input.unique_rna}")

    # # Fixmate (samtool fixmate)

    # rule fixmate_rna3:
    #     input:
    #         collated_rna = outdir_bam + "{sample}RNA3_collated.bam"
    #     threads:
    #         config["threads"]
    #     params:
    #         samtools_binary = config["samtools_binary"]
    #     output:
    #         fixmate_rna = temporary(outdir_bam + "{sample}RNA3_fixmate.bam")
    #     run:
    #         shell("{params.samtools_binary} fixmate -@ {threads} -m {input.collated_rna} {output.fixmate_rna}")

    # # Sort by coordinate (samtools sort)

    # rule sort_rna3:
    #     input:
    #         fixmate_rna = outdir_bam + "{sample}RNA3_fixmate.bam"
    #     threads:
    #         config["threads"]
    #     params:
    #         samtools_binary = config["samtools_binary"]
    #     output:
    #         coord_sorted_rna = outdir_bam + "{sample}RNA3_sorted.bam"
    #     run:
    #         shell("{params.samtools_binary} sort -@ {threads} -o {output.coord_sorted_rna} {input.fixmate_rna}")

    # rule index_rna3_bam:
    #     input:
    #         dedup_rna = outdir_bam + "{sample}RNA3_sorted.bam"
    #     threads:
    #         config["threads"]
    #     params:
    #         samtools_binary = config["samtools_binary"]
    #     output:
    #         dedup_rna_index = outdir_bam + "{sample}RNA3_sorted.bam.bai"
    #     run:
    #         shell("{params.samtools_binary} index -@ {threads} {input.dedup_rna}")

    # # rule rna3_bam2bed:
    # #     input:
    # #         dedup_rna = outdir_bam + "{sample}RNA3_sorted.bam"
    # #     output:
    # #         dedup_rna_bed = outdir_intersects + "{sample}RNA3_sorted.bed"
    # #     run:
    # #         shell("bam2bed -R < {input.dedup_rna} | sort -k 4,4 > {output.dedup_rna_bed}")

    # rule rna3_gene_intersect:
    #     input:
    #         dedup_rna = outdir_bam + "{sample}RNA3_sorted.bam"
    #     params:
    #         genes = config["genes"],
    #         bedtools_binary = config["bedtools_binary"]
    #     output:
    #         rna_gene_intersect = outdir_intersects + "{sample}RNA3_gene_intersect.txt"
    #     run:
    #         shell("{params.bedtools_binary} intersect -s -bed -wo -a {input.dedup_rna} -b {params.genes} > {output.rna_gene_intersect}")

    # rule cat_rna_intersects:
    #     input:
    #         rna5_intersect = outdir_intersects + "{sample}RNA5_gene_intersect.txt",
    #         rna3_intersect = outdir_intersects + "{sample}RNA3_gene_intersect.txt"
    #     output:
    #         cat_intersect = outdir_intersects + "{sample}RNA_gene_intersect_cat.txt"
    #     run:
    #         shell("cat {input.rna5_intersect} {input.rna3_intersect} > {output.cat_intersect}")

    # rule permissive_rna_gene_intersect:
    #     input:
    #         cat_intersect = outdir_intersects + "{sample}RNA_gene_intersect_cat.txt"
    #     params:
    #         genes = config["genes"],
    #         bedtools_binary = config["bedtools_binary"]
    #     output:
    #         perm_rna_intersect = outdir_intersects + "{sample}RNA_gene_intersect_maximums_cut.txt"
    #     run:
    #         shell("sort -k 4,4 {input.cat_intersect} | \
    #         {params.bedtools_binary} groupby -g 4 -c 1,2,3,6 -o mode,min,max,mode | \
    #         awk '{{OFS=\"\t\"}} {{($2-$1) < 5000}} {{print $2,$3,$4,$1,\".\",$5}}' | \
    #         {params.bedtools_binary} intersect -s -wo -a stdin -b {params.genes} | \
    #         sort -k 4,4 | \
    #         {params.bedtools_binary} groupby -g 4 -c 13 -o max -full | \
    #         cut -f 4,11 > {output.perm_rna_intersect}")

    rule red_max_intersects:
        input:
            rna3_max = outdir_intersects + "{sample}RNA3_gene_intersect_maximums_cut.txt",
            rna5_unique = outdir_intersects + "{sample}RNA5_gene_unstranded_intersect.txt",
            rna5_multi = outdir_intersects + "{sample}RNA5_gene_multi_unstranded_intersect.txt",
        output:
            max_intersects = outdir_intersects + "{sample}RNA_gene_intersect_maximums_cut.txt"
        params: 
            redc_intersect_max = workflow_dir + "scripts/redc_intersect_max.R"
        run:
            shell("Rscript {params.redc_intersect_max} --three {input.rna3_max} --uniquefive {input.rna5_unique} --multifive {input.rna5_multi} --output {output.max_intersects}")

    # if config["redc_mode"] == "strict":

    #     rule strict_rna_gene_combine:
    #         input:
    #             rna5 = outdir_intersects + "{sample}RNA5_gene_intersect_maximums_cut.txt",
    #             rna3 = outdir_intersects + "{sample}RNA3_gene_intersect_maximums_cut.txt"
    #         output:
    #             combined = outdir_intersects + "{sample}RNA_gene_intersect_maximums_cut.txt"
    #         run:
    #             shell("join -j 1 {input.rna5} {input.rna3} | sed -e 's/ /\t/g' | awk '$2==$3' | cut -f 1,2 > {output.combined}")

    # elif config["redc_mode"] == "permissive":

    #     rule permissive_rna_gene_combine:
    #         input:
    #             #rna5 = outdir_intersects + "{sample}RNA5_gene_intersect_maximums_cut.txt",
    #             rna3 = outdir_intersects + "{sample}RNA3_gene_intersect_maximums_cut.txt",
    #             rna5 = outdir_intersects + "{sample}RNA5_gene_unstranded_intersect_maximums_cut.txt"
    #         output:
    #             combined = outdir_intersects + "{sample}RNA_gene_intersect_maximums_cut.txt"
    #         run:
    #             shell("cat {input.rna5} {input.rna3} | bedtools groupby -g 1 -c 2 -o first | sort -k1 -k2 > {output.combined}")

    # # Add intersect as proportion of gene width column to RNA-gene intersect

    # rule rna_gene_intersect_redc:
    #     input:
    #         rna5_read_bed = outdir_intersects + "{sample}RNA5_sorted.bed",
    #         rna3_read_bed = outdir_intersects + "{sample}RNA3_sorted.bed"
    #     params:
    #         gene_bed = config["genes"]
    #     output:
    #         rna_gene_intersect = outdir_intersects + "{sample}RNA_gene_intersect.txt"
    #     run:
    #         shell("join -t$'\t' -j 4 {input.rna5_read_bed} {input.rna3_read_bed} | \
    #         awk '$2==$7' | \
    #         bedtools groupby -i stdin -g 1 -c 2,8,9,11 -o distinct,min,max,mode | \
    #         awk '{{OFS=\"\t\"}} {{print $2,$3,$4,$1,$5}}' | \
    #         sortBed -i stdin | \
    #         bedtools intersect -wo -a stdin -b {params.gene_bed} > {output.rna_gene_intersect}")

    # # Proportion of intersect relative to gene size

    # rule intersect_rna_proportion_redc:
    #     input:
    #         rna_gene_intersect = outdir_intersects + "{sample}RNA_gene_intersect.txt"
    #     output:
    #         gene_width = outdir_intersects + "{sample}RNA_gene_intersect_proportions.txt"
    #     run:
    #         shell("awk 'BEGIN{{OFS=\"\t\"}} {{$13 = $12/($8-$7); print}}' {input.rna_gene_intersect} > {output.gene_width}")

    # # select maximum proportion alignments per read

    # rule maximum_rna_proportion_redc:
    #     input:
    #         rna_gene_proportions = outdir_intersects + "{sample}RNA_gene_intersect_proportions.txt"
    #     output:
    #         rna_gene_maximums = outdir_intersects + "{sample}RNA_gene_intersect_maximums_cut.txt"
    #     run:
    #         shell("sort -r -k 4,4 -k 13,13 {input.rna_gene_proportions} | groupBy -g 4 -c 13 -o first -full | cut -f 4,10 | sort -k 1,1 > {output.rna_gene_maximums}")

    # rule intersect_rna3_proportion:
    #     input:
    #         rna_gene_intersect = outdir_intersects + "{sample}RNA3_gene_intersect.txt"
    #     output:
    #         gene_width = outdir_intersects + "{sample}RNA3_gene_intersect_proportions.txt"
    #     run:
    #         shell("awk 'BEGIN{{OFS=\"\t\"}} {{$20 = $19/($15-$14); print}}' {input.rna_gene_intersect} > {output.gene_width}")

    # # select maximum proportion alignments per read

    # rule maximum_rna3_proportion:
    #     input:
    #         rna_gene_proportions = outdir_intersects + "{sample}RNA3_gene_intersect_proportions.txt"
    #     output:
    #         rna_gene_maximums = outdir_intersects + "{sample}RNA3_gene_intersect_maximums.txt"
    #     run:
    #         shell("sort -r -k 4,4 -k 20,20 {input.rna_gene_proportions} | groupBy -g 4 -c 20 -o first -full | cut -f 4,17 | sort -k 1,1 > {output.rna_gene_maximums}")

    # rule combine_rna_intersects:
    #     input:
    #         rna5_intersect = outdir_intersects + "{sample}RNA5_gene_intersect_maximums.txt",
    #         rna3_intersect = outdir_intersects + "{sample}RNA3_gene_intersect_maximums.txt"
    #     output:
    #         rna_intersect = outdir_intersects + "{sample}RNA_gene_intersect_maximums_cut.txt"
    #     run:
    #         shell("join -t $'\t' -j 1 -o 1.1,1.2,2.1,2.2 {input.rna5_intersect} {input.rna3_intersect} | awk '$2==$4' | cut -f 1,2 > {output.rna_intersect}")

else:
        
    # remove ribosomal RNA reads

    rule gunzip_rna:
        input:
            rna_gz = fq_dir + "{sample}" + (config["rna_fastq_suffix"] if config["rna_fastq_suffix"].endswith(".gz") else config["rna_fastq_suffix"] + ".gz")
        threads:
            config["threads"]
        output:
            rna_fastq = temporary(fq_dir + "{sample}" + (config["rna_fastq_suffix"][:-3] if config["rna_fastq_suffix"].endswith(".gz") else config["rna_fastq_suffix"]))
        run:
            shell("pigz -k -d -p {threads} {input.rna_gz}")
    
    rule remove_rrna:
        input:
            rna_fastq = fq_dir + "{sample}"+config["rna_fastq_suffix"]
        threads:
            config["threads"]
        params:
            rrna_index = config["rrna_fasta"],
            bbduk = config["bbduk_script"]
        output:
            rna_clean = outdir_fastq + "{sample}RNA_depleted.fastq",
            rna_ribo = outdir_fastq + "{sample}rRNA.fastq",
            stats= outdir_fastq + "{sample}rRNA_removal_stats.txt"
        run:
            shell("{params.bbduk} in={input.rna_fastq} out={output.rna_clean} outm={output.rna_ribo} ref={params.rrna_index} k=13 hdist=1 stats={output.stats}")
            #shell("ribodetector_cpu -t {threads} -i {input.rna_fastq} -l 21 -e rrna --chunk_size 256 -o {output.rna_clean}")

    # align RNA

    rule align_rna:
        input:
            rna_clean = outdir_fastq + "{sample}RNA_depleted.fastq",
            genome_parameters = star_index + "Log.out"
        threads:
            config["threads"]
        params:
            # star_index = config["star_index"] if os.path.exists(config["star_index"]) else outdir_base + "resources/" + config["species"] + "/star_index",
            star_index = config["star_index"], 
            star_binary = config["star_binary"],
            base_name = outdir_bam + "{sample}"+'RNA_'
        output:
            aligned_rna = outdir_bam + "{sample}RNA_Aligned.out.bam",
            rna_map_log = outdir_bam + "{sample}RNA_Log.final.out"
        run:
            shell("{params.star_binary} \
                --runThreadN {threads} \
                --genomeDir {params.star_index} \
                --genomeLoad NoSharedMemory \
                --limitBAMsortRAM 30000000000 \
                --readFilesIn {input.rna_clean} \
                --outFileNamePrefix {params.base_name} \
                --outSAMtype BAM Unsorted \
                --alignIntronMax 1 \
                --alignMatesGapMax 1 \
                --outFilterScoreMinOverLread 0 \
                --outFilterMatchNminOverLread 0 \
                --outFilterMatchNmin 0")

    # rule align_rna:
    #     input:
    #         rna_clean = outdir_fastq + "{sample}RNA_depleted.fastq"
    #     threads:
    #         config["threads"]
    #     params:
    #         star_binary = config["star_binary"],
    #         star_index = config["star_index"],
    #         base_name = outdir_bam + "{sample}"+'RNA_'
    #     output:
    #         aligned_rna = outdir_bam + "{sample}RNA_Aligned.out.bam",
    #         rna_map_log = outdir_bam + "{sample}RNA_Log.final.out"
    #     run:
    #         shell("{params.star_binary} \
    #             --runThreadN {threads} \
    #             --genomeDir {params.star_index} \
    #             --genomeLoad NoSharedMemory \
    #             --limitBAMsortRAM 30000000000 \
    #             --readFilesIn {input.rna_clean} \
    #             --outFileNamePrefix {params.base_name} \
    #             --outSAMtype BAM Unsorted \
    #             --alignIntronMax 1 \
    #             --alignMatesGapMax 1 \
    #             --outFilterScoreMinOverLread 0 \
    #             --outFilterMatchNminOverLread 0 \
    #             --outFilterMatchNmin 0")










    # Remove blacklisted regions from RNA

    rule no_blacklist_rna:
        input:
            aligned_rna = outdir_bam + "{sample}RNA_Aligned.out.bam"
        params:
            blacklist = config["blacklist"]
        output:
            no_blacklist_rna = temporary(outdir_bam + "{sample}RNA_Aligned.out.bl.bam")
        run:
            shell("bedtools intersect -v -a {input.aligned_rna} -b {params.blacklist} > {output.no_blacklist_rna}")

    # Collate bam (samtools collate)

    rule unique_rna:
        input:
            aligned_rna = outdir_bam + "{sample}RNA_Aligned.out.bl.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            unique_rna = temporary(outdir_bam + "{sample}RNA_unique.bam"),
            multi_rna = temporary(outdir_bam + "{sample}RNA_multi.bam")
        run:
            shell("{params.samtools_binary} view -@ {threads} -q 255 -U {output.multi_rna} -o {output.unique_rna} {input.aligned_rna}")

    # Collate bam (samtools collate)

    rule collate_unique_rna:
        input:
            unique_rna = outdir_bam + "{sample}RNA_unique.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            collated_rna = temporary(outdir_bam + "{sample}RNA_unique_collated.bam")
        run:
            shell("{params.samtools_binary} collate -@ {threads} -o {output.collated_rna} {input.unique_rna}")

    rule collate_multi_rna:
        input:
            multi_rna = outdir_bam + "{sample}RNA_multi.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            collated_rna = temporary(outdir_bam + "{sample}RNA_multi_collated.bam")
        run:
            shell("{params.samtools_binary} collate -@ {threads} -o {output.collated_rna} {input.multi_rna}")            

    # Fixmate (samtool fixmate)

    rule fixmate_rna:
        input:
            collated_rna = outdir_bam + "{sample}RNA_unique_collated.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            fixmate_rna = temporary(outdir_bam + "{sample}RNA_fixmate.bam")
        run:
            shell("{params.samtools_binary} fixmate -@ {threads} -m {input.collated_rna} {output.fixmate_rna}")

    rule fixmate_multi_rna:
        input:
            collated_rna = outdir_bam + "{sample}RNA_multi_collated.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            fixmate_rna = temporary(outdir_bam + "{sample}RNA_multi_fixmate.bam")
        run:
            shell("{params.samtools_binary} fixmate -@ {threads} -m {input.collated_rna} {output.fixmate_rna}")

    # Sort by coordinate (samtools sort)

    rule sort_rna:
        input:
            fixmate_rna = outdir_bam + "{sample}RNA_fixmate.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            coord_sorted_rna = outdir_bam + "{sample}RNA_sorted.bam"
        run:
            shell("{params.samtools_binary} sort -@ {threads} -o {output.coord_sorted_rna} {input.fixmate_rna}")

    rule sort_multi_rna:
        input:
            fixmate_rna = outdir_bam + "{sample}RNA_multi_fixmate.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            coord_sorted_rna = outdir_bam + "{sample}RNA_multi_sorted.bam"
        run:
            shell("{params.samtools_binary} sort -@ {threads} -o {output.coord_sorted_rna} {input.fixmate_rna}")

    # # remove duplicates (samtools markdup -r)

    # rule dedup_rna:
    #     input:
    #         coord_sorted_rna = "{sample}RNA_sorted.bam"
    #     threads:
    #         config["threads"]
    #     output:
    #         dedup_rna = "{sample}RNA_dedup.bam"
    #     run:
    #         shell("samtools markdup -@ {threads} -r {input.coord_sorted_rna} {output.dedup_rna}")

    # # index deduplicated RNA bam

    # rule index_rna_bam:
    #     input:
    #         dedup_rna = "{sample}RNA_dedup.bam"
    #     threads:
    #         config["threads"]
    #     output:
    #         dedup_rna_index = "{sample}RNA_dedup.bam.bai"
    #     run:
    #         shell("samtools index -@ {threads} {input.dedup_rna}")

    # index deduplicated RNA bam

    rule index_rna_bam:
        input:
            dedup_rna = outdir_bam + "{sample}RNA_sorted.bam"
        threads:
            config["threads"]
        params:
            samtools_binary = config["samtools_binary"]
        output:
            dedup_rna_index = outdir_bam + "{sample}RNA_sorted.bam.bai"
        run:
            shell("{params.samtools_binary} index -@ {threads} {input.dedup_rna}")

    # # get genome coverage for RNA

    # rule rna_coverage:
    #     input:
    #         dedup_rna = "{sample}RNA_dedup.bam",
    #         dedup_rna_index = "{sample}RNA_dedup.bam.bai"
    #     threads:
    #         config["threads"]
    #     params:
    #         blacklist = config["blacklist"]
    #     output:
    #         dedup_rna_bw = "{sample}RNA_dedup_cpm.bw"
    #     run:
    #         shell("bamCoverage --bam {input.dedup_rna} \
    #                -o {output.dedup_rna_bw} \
    #                -of bigwig \
    #                -bs 1 \
    #                --blackListFileName {params.blacklist} \
    #                --normalizeUsing CPM \
    #                -p {threads}")

    # get genome coverage for RNA

    rule rna_coverage:
        input:
            dedup_rna = outdir_bam + "{sample}RNA_sorted.bam",
            dedup_rna_index = outdir_bam + "{sample}RNA_sorted.bam.bai"
        threads:
            config["threads"]
        params:
            bamCoverage_binary = config["bamCoverage_binary"],
            blacklist = config["blacklist"]
        output:
            dedup_rna_bw = outdir_bw + "{sample}RNA_sorted_cpm.bw"
        run:
            shell("{params.bamCoverage_binary} --bam {input.dedup_rna} \
                -o {output.dedup_rna_bw} \
                -of bigwig \
                -bs 1 \
                --blackListFileName {params.blacklist} \
                --normalizeUsing CPM \
                -p {threads}")

    # # RNA-gene intersect

    # rule rna_gene_intersect:
    #     input:
    #         dedup_rna = "{sample}RNA_dedup.bam"
    #     params:
    #         genes = config["genes"]
    #     output:
    #         rna_gene_intersect = "{sample}RNA_gene_intersect.txt"
    #     run:
    #         shell("bedtools intersect -bed -wo -a {input.dedup_rna} -b {params.genes} > {output.rna_gene_intersect}")

    if method == "RADICL":
        
        rule rna_gene_intersect_stranded_radicl:
            input:
                dedup_rna = outdir_bam + "{sample}RNA_sorted.bam"
            params:
                genes = config["genes"],
                bedtools_binary = config["bedtools_binary"]
            output:
                rna_gene_intersect = outdir_intersects + "{sample}RNA_gene_intersect.txt"
            run:
                shell("{params.bedtools_binary} intersect -s -bed -wo -abam {input.dedup_rna} -b {params.genes} > {output.rna_gene_intersect}")

        rule rna_gene_intersect_stranded_radicl_multi:
            input:
                dedup_rna = outdir_bam + "{sample}RNA_multi_sorted.bam"
            params:
                genes = config["genes"],
                bedtools_binary = config["bedtools_binary"]
            output:
                rna_gene_intersect = outdir_intersects + "{sample}RNA_multi_gene_intersect.txt"
            run:
                shell("{params.bedtools_binary} intersect -s -bed -wo -abam {input.dedup_rna} -b {params.genes} > {output.rna_gene_intersect}")

    elif method == "GRID" or method == "iMARGI":

        rule rna_gene_intersect_stranded_grid:
            input:
                dedup_rna = outdir_bam + "{sample}RNA_sorted.bam"
            params:
                genes = config["genes"],
                bedtools_binary = config["bedtools_binary"]
            output:
                rna_gene_intersect = outdir_intersects + "{sample}RNA_gene_intersect.txt"
            run:
                shell("{params.bedtools_binary} intersect -S -bed -wo -abam {input.dedup_rna} -b {params.genes} > {output.rna_gene_intersect}")

        rule rna_gene_intersect_stranded_grid_multi:
            input:
                dedup_rna = outdir_bam + "{sample}RNA_multi_sorted.bam"
            params:
                genes = config["genes"],
                bedtools_binary = config["bedtools_binary"]
            output:
                rna_gene_intersect = outdir_intersects + "{sample}RNA_multi_gene_intersect.txt"
            run:
                shell("{params.bedtools_binary} intersect -S -bed -wo -abam {input.dedup_rna} -b {params.genes} > {output.rna_gene_intersect}")


    rule rna_gene_intersect_unstranded:
        input:
            dedup_rna = outdir_bam + "{sample}RNA_sorted.bam"
        params:
            genes = config["genes"],
            bedtools_binary = config["bedtools_binary"]
        output:
            rna_gene_intersect = outdir_intersects + "{sample}RNA_gene_intersect_unstranded.txt"
        run:
            shell("{params.bedtools_binary} intersect -bed -wo -abam {input.dedup_rna} -b {params.genes} > {output.rna_gene_intersect}")

    # rule intersect_rna_proportion:
    #     input:
    #         rna_gene_intersect = outdir_intersects + "{sample}RNA_gene_intersect.txt"
    #     output:
    #         gene_width = outdir_intersects + "{sample}RNA_gene_intersect_proportions.txt"
    #     run:
    #         shell("awk 'BEGIN{{OFS=\"\t\"}} {{$20 = $19/($15-$14); print}}' {input.rna_gene_intersect} > {output.gene_width}")

    # # select maximum proportion alignments per read

    # rule maximum_rna_proportion:
    #     input:
    #         rna_gene_proportions = outdir_intersects + "{sample}RNA_gene_intersect_proportions.txt"
    #     params:
    #         bedtools_binary = config["bedtools_binary"]
    #     output:
    #         rna_gene_maximums = outdir_intersects + "{sample}RNA_gene_intersect_maximums_cut.txt"
    #     run:
    #         shell("sort -r -k 4,4 -k 20,20 {input.rna_gene_proportions} | {params.bedtools_binary} groupby -g 4 -c 20 -o first -full | cut -f 4,17 | sort -k 1,1 > {output.rna_gene_maximums}")

    rule process_intersections:
        input:
            unique_intersect = outdir_intersects + "{sample}RNA_gene_intersect.txt",
            multi_intersect = outdir_intersects + "{sample}RNA_multi_gene_intersect.txt"
        output:
            maximum_intersects = outdir_intersects + "{sample}RNA_gene_intersect_maximums_cut.txt"
        params: 
            intersect_processing = workflow_dir + "scripts/intersect_processing.R"
        run:
            shell("Rscript {params.intersect_processing} --unique {input.unique_intersect} --multi {input.multi_intersect} --output {output.maximum_intersects}")

# join RNA and DNA reads

rule join_RNA_and_DNA:
    input:
        cut_rna = outdir_intersects + "{sample}RNA_gene_intersect_maximums_cut.txt",
        cut_dna = outdir_intersects + "{sample}DNA_bin_intersect_maximums_cut.txt"
    output:
        joined_rna_dna = outdir_merge + "{sample}RNA-bin_pairs.txt"
    run:
        shell("join -t$'\t' -j 1 -o 1.1,1.2,2.1,2.2 {input.cut_rna} {input.cut_dna} > {output.joined_rna_dna}")

# count RNA-bin pairs

rule count_joins:
    input:
        joined_rna_dna = outdir_merge + "{sample}RNA-bin_pairs.txt"
    output:
        rna_bin_counts = outdir_counts + "{sample}RNA-bin_counts.txt"
    run:
        shell("sort -k 2,2 -k 4,4 {input.joined_rna_dna} | cut -f 2,4 | uniq -c | awk '{{print $2\"\t\"$3\"\t\"$1}}' > {output.rna_bin_counts}")

# Call interactions with RADIAnT

rule radiant:
    input:
        counts = outdir_counts + "{sample}RNA-bin_counts.txt"    
    params:
        gtf = config["gtf"],
        counts = outdir_counts + "{sample}RNA-bin_counts.txt",
        bins = config["genome_bins"],
        species = config["species"],
        outdir = outdir_interactions,
        name = "{sample}",
        RADIAnT_command_line = workflow_dir + "scripts/RADIAnT_command_line.R"
    output:
        rna_bin_interactions = outdir_interactions + "{sample}RADIAnT_results.txt"
    run:
        shell("Rscript {params.RADIAnT_command_line} \
               --gtf {params.gtf} \
               --counts {input.counts} \
               --bins {params.bins} \
               --species {params.species} \
               --outdir {params.outdir} \
               --name {params.name}")


if method == 'Red-C':
    rule sankey: 
        input: 
            ribo_stats= outdir_fastq + "{sample}rRNA5_removal_stats.txt",
            rna_map_log = outdir_bam + "{sample}RNA5_Log.final.out",
            rna_bin_interactions = outdir_interactions + "{sample}RADIAnT_results.txt"
        output:
            svg = outdir_logs + "{sample}sankey.svg",
            png = outdir_logs + "{sample}sankey.png",
            txt = outdir_logs + "{sample}read_stats.txt"
        params: 
            plot_read_stats = workflow_dir + "scripts/plot_read_stats.R"
        run:
            shell("Rscript {params.plot_read_stats} {input.ribo_stats} {input.rna_map_log} {input.rna_bin_interactions} {output.svg} {output.png} {output.txt}")
else:
    rule sankey: 
        input: 
            ribo_stats= outdir_fastq + "{sample}rRNA_removal_stats.txt",
            rna_map_log = outdir_bam + "{sample}RNA_Log.final.out",
            rna_bin_interactions = outdir_interactions + "{sample}RADIAnT_results.txt"
        output:
            svg = outdir_logs + "{sample}sankey.svg",
            png = outdir_logs + "{sample}sankey.png",
            txt = outdir_logs + "{sample}read_stats.txt"
        params: 
            plot_read_stats = workflow_dir + "scripts/plot_read_stats.R"
        run:
            shell("Rscript {params.plot_read_stats} {input.ribo_stats} {input.rna_map_log} {input.rna_bin_interactions} {output.svg} {output.png} {output.txt}")




rule gene_int_stats: 
    input: 
        rna_bin_interactions = outdir_interactions + "{sample}RADIAnT_results.txt"
    params: 
        plot_gene_stats = workflow_dir + "scripts/plot_gene_stats.R",
        gtf = config["gtf"],
        outdir_interactions = outdir_interactions
    output:
        txt = outdir_interactions + "{sample}genes.number_of_interactions.txt"
    run:
        shell("Rscript {params.plot_gene_stats} {input.rna_bin_interactions} {params.gtf} {params.outdir_interactions} {output.txt}")



