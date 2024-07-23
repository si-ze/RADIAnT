library(data.table)
library(dplyr)
library(argparser)

# Set up argument parser

message('Parsing command line arguments')

parser = arg_parser(description='Retrieve maximal RNA-gene intersects from RNA-DNA ligation data')

parser = add_argument(parser, '--unique', type = 'character', help = 'RNA-gene intersects of uniquely mapping reads')

parser = add_argument(parser, '--multi', type = 'character', help = 'RNA-gene intersects of multi-mapping reads')

parser = add_argument(parser, '--output', type = 'character', help = 'Maximal intersects outfile')

# Parse command-line arguments

arg_vector = parse_args(parser)

# Read in and process unique intersections

message('Reading unique intersects')

unique_intersect = data.table::fread(arg_vector$unique)

column_names = c('read_chr', 'read_start', 'read_stop', 'read_id', 'mapping_quality', 'read_strand',
                 'read_start_2', 'read_stop_2', 'rgb', 'misc', 'align_length', 'misc_2', 'gene_chr',
                 'gene_start', 'gene_stop', 'score', 'gene_id', 'gene_strand', 'intersect_length')

colnames(unique_intersect) = column_names

unique_intersect$gene_size = unique_intersect$gene_stop - unique_intersect$gene_start

unique_intersect$intersect_proportion = unique_intersect$intersect_length / unique_intersect$gene_size

setkey(unique_intersect, 'read_id')

data.table::setorder(unique_intersect, 'read_id')

message('Processing unique intersects')

max_unique_intersect = unique_intersect[,.(gene_id = gene_id[which.max(intersect_proportion)]), by = 'read_id']

# dplyr_start = Sys.time()
# max_unique_intersect_dplyr = unique_intersect %>%
#     group_by(read_id) %>%
#     arrange(desc(intersect_proportion)) %>%
#     slice(1L) %>%
#     setDT()
# dplyr_end = Sys.time()
# dplyr_duration = dplyr_end - dplyr_start

cts_per_gene = max_unique_intersect[,.(sum_cts = length(read_id)), by = 'gene_id']

# Read in and process multi-mapping intersects

message('Reading multi-mapping intersects')

multi_intersect = data.table::fread(arg_vector$multi)

colnames(multi_intersect) = column_names

multi_intersect$unique_cts_per_gene = cts_per_gene$sum_cts[match(multi_intersect$gene_id, cts_per_gene$gene_id)]

multi_intersect$unique_cts_per_gene[is.na(multi_intersect$unique_cts_per_gene)] = 0

multi_intersect$gene_length = multi_intersect$gene_stop - multi_intersect$gene_start

multi_intersect$intersect_proportion = multi_intersect$intersect_length/multi_intersect$gene_length

message('Processing multi-mapping intersects')

max_multi_intersect = multi_intersect[,.(gene_id = gene_id[which.max(unique_cts_per_gene*intersect_proportion)]), by = 'read_id']

# max_multi_intersect = multi_intersect %>%
#     group_by(read_id) %>%
#     arrange(desc(unique_cts_per_gene), desc(intersect_proportion)) %>%
#     slice(1L) %>%
#     setDT()

max_intersect_bind = rbind(max_unique_intersect[,c('read_id', 'gene_id')],
                           max_multi_intersect[,c('read_id', 'gene_id')])

max_intersect_bind = max_intersect_bind[order(max_intersect_bind$read_id),]                           

message('Writing max intersects')

data.table::fwrite(max_intersect_bind, file = arg_vector$output, col.names = FALSE, row.names = FALSE, sep = '\t', eol = '\n')
