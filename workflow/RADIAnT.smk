# Snakemake for processing of RNA-DNA interaction data from split fastq files to Gene-Bin interaction counts


#################################################################################################################################
# Config file ===================================================================================================================
#################################################################################################################################

# script directory ==============================================================================================================

workflow_dir = config["workflow_directory"] if config["workflow_directory"].endswith("/") else config["workflow_directory"] + "/"

# resource directory ============================================================================================================

resource_dir = config["resource_directory"] if config["resource_directory"].endswith("/") else config["resource_directory"] + "/"

# species =======================================================================================================================

species = config["species"]

# bin_base =======================================================================================================================

bin_base = config["bin_base"] 

# method ========================================================================================================================

method = config["method"][0]

# Directory holding the fastq files of the experiment to be analysed ============================================================

fq_dir = config["fastq_directory"] if config["fastq_directory"].endswith("/") else config["fastq_directory"] + "/"

# Samples basename ==============================================================================================================

samples = config["sample_base"]

print(samples)

# Bin sizes for analysis ========================================================================================================

bin_sizes = [str(x) for x in config["bin_sizes"]] 

print(bin_sizes)

# Mapping  ======================================================================================================================

star_index = config["star_index"] if config["star_index"].endswith("/") else config["star_index"] + "/"

# output directories ============================================================================================================

outdir_base = config["output_directory"] if config["output_directory"].endswith("/") else config["output_directory"] + "/"

outdir_bam = outdir_base + "bam/"

outdir_bw = outdir_base + "bw/"

outdir_fastq = outdir_base + "fastq/"

outdir_intersects = outdir_base + "intersects/"

outdir_merge = outdir_base + "merge/"

outdir_counts = outdir_base + "counts/"

outdir_interactions = outdir_base + "interactions/"

outdir_logs = outdir_base + "logs/"


#################################################################################################################################
# Rules =========================================================================================================================
#################################################################################################################################

# constrain wildcards ===========================================================================================================

wildcard_constraints:
	sample = "(" + "|".join(map(str,samples)) + ")",
	bin_size = "(" + "|".join(map(str,bin_sizes)) + ")"





# all ===========================================================================================================================


rule all:
    input:
        expand(outdir_interactions + "{sample}{bin_size}_RADIAnT_results.txt", sample=samples, bin_size = bin_sizes),
        star_index + "Log.out"


# build blacklist  ==============================================================================================================

rule build_effective_blacklist:
  input:
    blacklist = config["blacklist"],
    gtf = config["gtf"],
    script = workflow_dir + "scripts/build_effective_blacklist.R"
  output:
    effective_blacklist = temporary(outdir_base + "effective_blacklist.bed")
  params:
    biotypes = ",".join(config.get("blacklist_biotypes", [])),
    bedtools_binary = config["bedtools_binary"]
  run:
    shell("""Rscript {input.script} \
    --blacklist {input.blacklist} \
    --gtf {input.gtf} \
    --biotypes "{params.biotypes}" \
    --outbed {output.effective_blacklist}
    """)



# If no STAR index provided, decompress provided genome FASTA to build index ====================================================

rule gunzip_genome_fasta:
    input: 
        genome_fasta = config["genome_fasta"] if config["genome_fasta"].endswith(".gz") else config["genome_fasta"] + ".gz"
    output: 
        decompressed_fasta = temporary(config["genome_fasta"][:-3] if config["genome_fasta"].endswith(".gz") else config["genome_fasta"])
    run: 
        shell("pigz -k -d -p {threads} {input.genome_fasta}")

# If no STAR index provided, decompress provided genome annotation to build index ===============================================

rule gunzip_gtf:
    input: 
        gtf = config["gtf"] if config["gtf"].endswith(".gz") else config["gtf"] + ".gz"
    output: 
        decompressed_gtf = temporary(config["gtf"][:-3] if config["gtf"].endswith(".gz") else config["gtf"])
    run: 
        shell("pigz -k -d -p {threads} {input.gtf}")

# If no STAR index provided, build index ========================================================================================

rule build_star_index: 
    input:
        gtf = re.sub(r"\.gz$", "", config["gtf"]),
        genome_fasta = re.sub(r"\.gz$", "", config["genome_fasta"])
    params:
        star_binary = config["star_binary"],
        star_index = config["star_index"]
    threads:
        config["threads"]
    output: 
        star_log = star_index + "Log.out"
    run:
        shell("{params.star_binary} \
        --runMode genomeGenerate \
	--runThreadN {threads} \
        --genomeDir {params.star_index} \
        --genomeFastaFiles {input.genome_fasta} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 50")


# Decompress FASTQs =============================================================================================================

rule gunzip_dna:
    input:
        dna_gz = fq_dir + "{sample}" + (config["dna_fastq_suffix"] if config["dna_fastq_suffix"].endswith(".gz") else config["dna_fastq_suffix"] + ".gz")
    threads:
        config["threads"]
    output:
        dna_fastq = temporary(fq_dir + "{sample}" + (config["dna_fastq_suffix"][:-3] if config["dna_fastq_suffix"].endswith(".gz") else config["dna_fastq_suffix"]))
    run:
        shell("pigz -k -d -p {threads} {input.dna_gz}")

# DNA alignment =================================================================================================================

rule align_dna:
    input:
        dna_fastq = fq_dir + "{sample}" + (config["dna_fastq_suffix"][:-3] if config["dna_fastq_suffix"].endswith(".gz") else config["dna_fastq_suffix"]), # fq_dir + "{sample}"+config["dna_fastq_suffix"],
        genome_parameters = star_index + "Log.out"
    threads:
        config["threads"]
    params:
        # star_index = config["star_index"] if os.path.exists(config["star_index"]) else outdir_base + "resources/" + config["species"] + "/star_index",
        star_index = config["star_index"], 
        star_binary = config["star_binary"],
        base_name = outdir_bam + "{sample}" + 'DNA_'
    output:
        aligned_dna = outdir_bam + "{sample}DNA_Aligned.out.bam",
        dna_log = outdir_bam + "{sample}DNA_Log.final.out"
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

# Remove blacklisted regions from DNA ===========================================================================================

rule blacklist_filter_dna:
    input:
      aligned_dna = outdir_bam + "{sample}DNA_Aligned.out.bam",
      blacklist = outdir_base + "effective_blacklist.bed"
    output:
        blacklist_filtered_dna = temporary(outdir_bam + "{sample}DNA_Aligned.out.bl_filt.bam")
    run:
        shell("bedtools intersect -v -a {input.aligned_dna} -b {input.blacklist} > {output.blacklist_filtered_dna}")


# extract uniquely mapping reads ================================================================================================

rule unique_dna:
    input:
        blacklist_filtered_dna = outdir_bam + "{sample}DNA_Aligned.out.bl_filt.bam"
    threads:
        config["threads"]
    params:
        samtools_binary = config["samtools_binary"]
    output:
        unique_dna = temporary(outdir_bam + "{sample}DNA_unique.bam")
    run:
        shell("{params.samtools_binary} view -@ {threads} -q 255 -o {output.unique_dna} {input.blacklist_filtered_dna}")

# Collate bam (samtools collate) ================================================================================================

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

# Fixmate (samtool fixmate) =====================================================================================================

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


# Sort by coordinate (samtools sort) ============================================================================================

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

# Index DNA BAM file ============================================================================================================

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

# Normalized DNA coverage =======================================================================================================

rule dna_coverage:
    input:
        dedup_dna = outdir_bam + "{sample}DNA_sorted.bam",
        dedup_dna_index = outdir_bam + "{sample}DNA_sorted.bam.bai"
    threads:
        config["threads"]
    params:
        bamCoverage_binary = config["bamCoverage_binary"],
        blacklist = outdir_base + "effective_blacklist.bed"
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

# Intersect DNA reads with genomic bins =========================================================================================

rule dna_bin_intersect:
    input:
        dedup_dna = outdir_bam + "{sample}DNA_sorted.bam", 
        genome_bins = resource_dir + species + "/" + bin_base + ".{bin_size}_bins.bed.gz"
    params:
        bedtools_binary = config["bedtools_binary"]
    output: 
        dna_bin_intersect = outdir_intersects + "{sample}DNA_bin_intersect_{bin_size}.txt"
    run:
        shell("{params.bedtools_binary} intersect -bed -f 0.51 -wo -a {input.dedup_dna} -b {input.genome_bins} > {output.dna_bin_intersect}")

# # Add intersect as proportion of gene width column to RNA-gene intersect

# rule intersect_DNA_proportion:
#     input:
#         dna_bin_intersect = outdir_intersects + "{sample}DNA_bin_intersect_{bin_size}.txt"
#     output:
#         bin_width = outdir_intersects + "{sample}DNA_bin_intersect_proportions_{bin_size}.txt"
#     run:
#         shell("awk 'BEGIN{{OFS=\"\t\"}} {{$18 = $17/($15-$14); print}}' {input.dna_bin_intersect} > {output.bin_width}")

# # Cut DNA intersect proportions to only required columns

# rule cut_DNA_proportion:
#     input:
#         dna_bin_proportions = outdir_intersects + "{sample}DNA_bin_intersect_proportions_{bin_size}.txt"
#     output:
#         dna_props_cut = temporary(outdir_intersects + "{sample}DNA_bin_intersect_proportions_cut_{bin_size}.txt")
#     run:
#         shell("cut -f 4,16,18 {input.dna_bin_proportions} > {output.dna_props_cut}")

# # Sort DNA proportions by proportion of overlap

# rule sort_dna_props:
#     input:
#         dna_props_cut = outdir_intersects + "{sample}DNA_bin_intersect_proportions_cut_{bin_size}.txt"
#     output:
#         dna_props_cut_sorted = temporary(outdir_intersects + "{sample}DNA_bin_intersect_proportions_cut_sorted_{bin_size}.txt")
#     threads:
#         config["threads"]
#     run:
#         shell("sort -S 50% --parallel={threads} -k 1,1 -k 3,3r {input.dna_props_cut} > {output.dna_props_cut_sorted}")

# rule unique_sorted_dna_props:
#     input:
#         dna_props_cut_sorted = outdir_intersects + "{sample}DNA_bin_intersect_proportions_cut_sorted_{bin_size}.txt"
#     output:
#         dna_props_cut_unique = temporary(outdir_intersects + "{sample}DNA_bin_intersect_proportions_cut_sorted_unique_{bin_size}.txt")
#     threads:
#         config["threads"]
#     run:
#         shell("sort -S 50% --parallel={threads} -u -k 1b,1 {input.dna_props_cut_sorted} > {output.dna_props_cut_unique}")

rule max_cut_dna:
    input:
        dna_props_cut_unique = outdir_intersects + "{sample}DNA_bin_intersect_{bin_size}.txt"
    output:
        dna_bin_maximums = outdir_intersects + "{sample}DNA_bin_intersect_maximums_cut_{bin_size}.txt"
    threads:
        config["threads"]
    run:
        shell("cut -f 4,16 {input.dna_props_cut_unique} > {output.dna_bin_maximums}") 

rule sort_max_cut_dna:
    input:
        dna_bin_maximums = outdir_intersects + "{sample}DNA_bin_intersect_maximums_cut_{bin_size}.txt"
    output:
        dna_bin_maximums_sorted = outdir_intersects + "{sample}DNA_bin_intersect_maximums_cut_sorted_{bin_size}.txt"
    threads:
        config["threads"]
    run:
        shell("sort --parallel {threads} -S 50% -k 1,1 {input.dna_bin_maximums} > {output.dna_bin_maximums_sorted}")

# rule max_cut_dna:
#     input:
#         dna_props_cut_unique = outdir_intersects + "{sample}DNA_bin_intersect_proportions_cut_sorted_unique_{bin_size}.txt"
#     output:
#         dna_bin_maximums = outdir_intersects + "{sample}DNA_bin_intersect_maximums_cut_{bin_size}.txt"
#     run:
#         shell("cut -f 1,2 {input.dna_props_cut_unique} > {output.dna_bin_maximums}") 


# select maximum proportion alignments per read

# rule maximum_DNA_proportion:
#     input:
#         dna_bin_proportions = outdir_intersects + "{sample}DNA_bin_intersect_proportions.txt"
#     params:
#         bedtools_binary = config["bedtools_binary"]
#     output:
#         dna_bin_maximums = outdir_intersects + "{sample}DNA_bin_intersect_maximums_cut.txt"
#     run:
#         #shell("sort -r -k 4,4 -k 18,18 {input.dna_bin_proportions} | {params.bedtools_binary} groupby -g 4 -c 18 -o first -full | cut -f 4,16 | sort -k 1,1 > {output.dna_bin_maximums}")
#         shell("cut -f 4,16,18 {input.dna_bin_proportions} | sort -S 50% --parallel={threads} -k 1,1 -k 3,3r | sort -S 50% --parallel={threads} -u -k 1b,1 | cut -f 1,2 > {output.dna_bin_maximums}")

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
                

    rule blacklist_filter_rna5:
        input:
            aligned_rna = outdir_bam + "{sample}RNA5_Aligned.out.bam",
            blacklist = outdir_base + "effective_blacklist.bed"
        output:
            blacklist_filtered_rna = temporary(outdir_bam + "{sample}RNA5_Aligned.out.bl_filt.bam") 
        run:
            shell("bedtools intersect -v -a {input.aligned_rna} -b {input.blacklist} > {output.blacklist_filtered_rna}")

    rule unique_rna5:
        input:
            aligned_rna = outdir_bam + "{sample}RNA5_Aligned.out.bl_filt.bam"
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
                
    rule blacklist_filter_rna3:
        input:
            aligned_rna = outdir_bam + "{sample}RNA3_Aligned.out.bam",
            blacklist = outdir_base + "effective_blacklist.bed"
        output:
            blacklist_filtered_rna = temporary(outdir_bam + "{sample}RNA3_Aligned.out.bl_filt.bam") 
        run:
            shell("bedtools intersect -v -a {input.aligned_rna} -b {input.blacklist} > {output.blacklist_filtered_rna}")
    
    rule unique_rna3:
        input:
            aligned_rna = outdir_bam + "{sample}RNA3_Aligned.out.bl_filt.bam"
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

    rule sort_red_rna_intersects:
        input:
            intersects = outdir_intersects + "{sample}RNA_gene_intersect_maximums_cut.txt"
        output:
            sorted_intersects = outdir_intersects + "{sample}RNA_gene_intersect_maximums_cut_sorted.txt"
        threads:
            config["threads"]
        run:
            shell("sort -S 50% --parallel={threads} -k 1b,1 {input.intersects} > {output.sorted_intersects}")

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
            genome_parameters = config["star_index"] + "Log.out"
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

    # Remove blacklisted regions from RNA

    rule blacklist_filter_rna:
        input:
            aligned_rna = outdir_bam + "{sample}RNA_Aligned.out.bam",
            blacklist = outdir_base + "effective_blacklist.bed"
        output:
            blacklist_filtered_rna = temporary(outdir_bam + "{sample}RNA_Aligned.out.bl_filt.bam") 
        run:
            shell("bedtools intersect -v -a {input.aligned_rna} -b {input.blacklist} > {output.blacklist_filtered_rna}")

    # Collate bam (samtools collate)

    rule unique_rna:
        input:
            aligned_rna = outdir_bam + "{sample}RNA_Aligned.out.bl_filt.bam" 
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
    #         blacklist = outdir_base + "effective_blacklist.bed"
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
            blacklist = outdir_base + "effective_blacklist.bed"
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

  
    rule process_intersections:
        input:
            unique_intersect = outdir_intersects + "{sample}RNA_gene_intersect.txt",
            multi_intersect = outdir_intersects + "{sample}RNA_multi_gene_intersect.txt"
        output:
            maximum_intersects = outdir_intersects + "{sample}RNA_gene_intersect_maximums_cut.txt"
        params: 
            intersect_processing = workflow_dir + "scripts/intersect_processing.R",
            intersect_rule_unique = config.get("intersect_rule_unique", "all"),
            intersect_rule_multi = config.get("intersect_rule_multi", "all")
        run:
            shell("Rscript {params.intersect_processing} "
              "--unique {input.unique_intersect} "
              "--multi {input.multi_intersect} "
              "--output {output.maximum_intersects} "
              "--intersect-rule-unique {params.intersect_rule_unique} "
              "--intersect-rule-multi {params.intersect_rule_multi}")

    rule sort_rna_intersects:
        input:
            intersects = outdir_intersects + "{sample}RNA_gene_intersect_maximums_cut.txt"
        output:
            sorted_intersects = outdir_intersects + "{sample}RNA_gene_intersect_maximums_cut_sorted.txt"
        threads:
            config["threads"]
        run:
            shell("sort -S 50% --parallel={threads} -k 1b,1 {input.intersects} > {output.sorted_intersects}")

# join RNA and DNA reads

rule join_RNA_and_DNA:
    input:
        cut_rna = outdir_intersects + "{sample}RNA_gene_intersect_maximums_cut_sorted.txt",
        cut_dna = outdir_intersects + "{sample}DNA_bin_intersect_maximums_cut_sorted_{bin_size}.txt"
    output:
        joined_rna_dna = outdir_merge + "{sample}RNA-bin_pairs_{bin_size}.txt"
    run:
        shell("join -t$'\t' -j 1 -o 1.1,1.2,2.1,2.2 {input.cut_rna} {input.cut_dna} > {output.joined_rna_dna}")

# sort joined reads by RNA/bin

rule sort_join:
    input:
        joined_rna_dna = outdir_merge + "{sample}RNA-bin_pairs_{bin_size}.txt"
    output:
        sorted_rna_dna = temporary(outdir_merge + "{sample}RNA-bin_pairs_sorted_{bin_size}.txt")
    threads:
        config["threads"]
    run:
        shell("sort -S 50% --parallel {threads} -k 2,2 -k 4,4 {input.joined_rna_dna} > {output.sorted_rna_dna}")

rule cut_join:
    input:
        sorted_rna_dna = outdir_merge + "{sample}RNA-bin_pairs_sorted_{bin_size}.txt"
    output:
        cut_rna_dna = temporary(outdir_merge + "{sample}RNA-bin_pairs_cut_{bin_size}.txt")
    run:
        shell("cut -f 2,4 {input.sorted_rna_dna} > {output.cut_rna_dna}")

rule count_joins:
    input:
        cut_rna_dna = outdir_merge + "{sample}RNA-bin_pairs_cut_{bin_size}.txt"
    output:
        rna_bin_counts = temporary(outdir_counts + "{sample}RNA-bin_counts_raw_{bin_size}.txt")
    run:
        shell("uniq -c {input.cut_rna_dna} > {output.rna_bin_counts}")

rule format_counts:
    input:
        rna_bin_counts = outdir_counts + "{sample}RNA-bin_counts_raw_{bin_size}.txt"
    output:
        formatted_counts = outdir_counts + "{sample}RNA-bin_counts_{bin_size}.txt"
    run:
        shell("awk 'OFS=\"\t\" {{print $2,$3,$1}}' {input.rna_bin_counts} > {output.formatted_counts}")

# count RNA-bin pairs

# rule count_joins:
#     input:
#         joined_rna_dna = outdir_merge + "{sample}RNA-bin_pairs_{bin_size}.txt"
#     output:
#         rna_bin_counts = outdir_counts + "{sample}RNA-bin_counts_{bin_size}.txt"
#     threads:
#         config["threads"]
#     run:
#         shell("sort -S 50% --parallel {threads} -k 2,2 -k 4,4 {input.joined_rna_dna} | cut -f 2,4 | uniq -c | awk '{{print $2\"\t\"$3\"\t\"$1}}' > {output.rna_bin_counts}")

# Call interactions with RADIAnT

rule radiant:
    input:
        counts = outdir_counts + "{sample}RNA-bin_counts_{bin_size}.txt"
    params:
        gtf = config["gtf"],
        counts = outdir_counts + "{sample}RNA-bin_counts_{bin_size}.txt",
        bins =  resource_dir + species + "/" + bin_base + ".{bin_size}_bins.bed.gz", 
        species = species,
        outdir = outdir_interactions,
        name = "{sample}{bin_size}_",
        RADIAnT_command_line = workflow_dir + "scripts/RADIAnT_command_line.R"
    output:
        rna_bin_interactions = outdir_interactions + "{sample}{bin_size}_RADIAnT_results.txt"
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

    rule updated_sankey: 
        input: 
            unsplit = fq_dir + "",
            rna_fastq = fq_dir + "{sample}"+config["rna_fastq_suffix"],
            rna_fastq_depl = outdir_fastq + "{sample}RNA_depleted.fastq",
            dna_fastq = fq_dir + "{sample}"+config["dna_fastq_suffix"],
            rna_map_log = outdir_bam + "{sample}RNA_Log.final.out",
            dna_map_log = outdir_bam + "{sample}DNA_Log.final.out",
            rna_unique_intersect = outdir_intersects + "{sample}RNA_gene_intersect.txt",
            rna_multi_intersect = outdir_intersects + "{sample}RNA_multi_gene_intersect.txt",
            dna_bin_intersect = outdir_intersects + "{sample}DNA_bin_intersect_maximums_cut.txt",
            pairs = outdir_merge + "{sample}RNA-bin_pairs.txt",
            radiant = outdir_interactions + "{sample}RADIAnT_results.txt"
        output:
            svg = outdir_logs + "{sample}Sankey_updated.svg",
            #png = outdir_logs + "{sample}Sankey_updated.png"
            #txt = outdir_logs + "{sample}read_stats.txt"
        params: 
            plot_read_stats = workflow_dir + "scripts/combined_Sankey.R",
            outdir = outdir_logs
        run:
            shell("Rscript {params.plot_read_stats} --unsplit {input.unsplit} --rnafastq {input.rna_fastq} --deplrnafastq {input.rna_fastq_depl} --dnafastq {input.dna_fastq} --rnalog {input.rna_map_log} --dnalog {input.dna_map_log} --rnauniqueintersect {input.rna_unique_intersect} --rnamultiintersect {input.rna_multi_intersect} --dnaintersect {input.dna_bin_intersect} --pairs {input.pairs} --radiant {input.radiant} --output {params.outdir}{wildcards.sample}")

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



