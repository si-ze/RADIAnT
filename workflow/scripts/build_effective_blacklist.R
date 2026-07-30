#!/usr/bin/env Rscript

message("Creating effective blacklist")

library(data.table)
library(argparser)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)


parser = arg_parser(description='')

parser = add_argument(parser, '--blacklist', type = 'character', help = 'BED-file holding already blacklisted regions.')
parser = add_argument(parser, '--gtf', type = 'character', help = 'GTF holding gene annotations')
parser = add_argument(parser, '--biotypes', type = 'character', default = "", help = 'Comma-separated list of biotypes to blacklist, i.e. --biotypes snRNA,snoRNA,rRNA')
parser = add_argument(parser, '--outbed', type = 'character', help = 'Provide the name for the BED file to output the effective blacklist to.')
arg_vector = parse_args(parser)

blacklist_file = arg_vector$blacklist
gtf_file = arg_vector$gtf
biotypes_string = arg_vector$biotypes



# TESTING
# gtf_file = "/projects/tim/RADICL/cat_RADICL_endo/config/config_cat_RADICL_EMT.yaml"
# blacklist_file = "/projects/bioinformatics/genome_annotations/t2t/T2T_mock_blacklist.bed"


# gtf_file = "/projects/bioinformatics/tools/RADIAnT/resources/mouse/gencode.vM36.annotation.gtf.gz"
# blacklist_file = "/projects/bioinformatics/tools/RADIAnT/resources/mouse/mm39.excluderanges.bed.gz"
# biotypes_string="snRNA,        snoRNA"
# TESTING



# Get the existing blacklist
message("Used blacklist: ", blacklist_file)
blacklist <- data.table::fread(blacklist_file)
blacklist_gr <- GRanges(
  seqnames = blacklist$V1,
  ranges= IRanges(
    start=blacklist$V2 + 1, # since BED is 0-based
    end=blacklist$V3
  )
)



# Get all biotypes which should be added as blacklist-regions
message("Biotypes selected to be blacklisted: ", biotypes_string)
if (biotypes_string == "") {
  biotypes = character(0)
} else {
  biotypes = strsplit(biotypes_string, ",")[[1]]
  biotypes = gsub(" ", "", biotypes)
  biotypes = biotypes[biotypes != ""]
}


# Get all gene annotations
genes <- rtracklayer::readGFF(gtf_file) %>% 
  dplyr::filter(type=="gene")

# Search for the biotype column in the GTF (naming differs, depending on the origin of the file, e.g. )
# and then filter for only those rows which contain the biotype to be blacklisted
if ("gene_biotype" %in% colnames(genes)) {
  biotype_rows <- genes %>% dplyr::filter(gene_biotype %in% biotypes)
} else if ("gene_type" %in% colnames(genes)) {
  biotype_rows <- genes %>% dplyr::filter(gene_type %in% biotypes)
} else if ("transcript_biotype" %in% colnames(genes)) {
  biotype_rows <- genes %>% dplyr::filter(transcript_biotype %in% biotypes)
} else if ("transcript_type" %in% colnames(genes)) {
  biotype_rows <- genes %>% dplyr::filter(transcript_type %in% biotypes)
} else {
  stop("No biotype column found in GTF. Checked gene_biotype, gene_type, transcript_biotype, transcript_type.")
}

# get the GRanges of the biotypes of interest
biotypes_gr <- GRanges(
  seqnames = biotype_rows$seqid,
  ranges = IRanges(start=biotype_rows$start, end=biotype_rows$end)
)

message("Adding biotypes to blacklist ...")
# combine the pre-existing blacklist with the biotypes to exclude --> this is the effective blacklist
effective_blacklist_gr <- GenomicRanges::reduce(c(blacklist_gr, biotypes_gr))
# write to BED (should be found in /resources/effective_blacklist.bed)
rtracklayer::export.bed(effective_blacklist_gr, con=arg_vector$outbed)
message("Effective blacklist generated.")

