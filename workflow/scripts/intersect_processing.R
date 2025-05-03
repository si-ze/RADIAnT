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

intersect_type = "all"

multi_intersect_type = "all"

# Troubleshooting

#arg_vector$unique = "M:/RADICL/Data/GRID/MDA-MB231/radiant_out/intersects/GRID_MDA-MB-231_cat_RNA_gene_intersect.txt"

#arg_vector$multi = "M:/RADICL/Data/GRID/MDA-MB231/radiant_out/intersects/GRID_MDA-MB-231_cat_RNA_multi_gene_intersect.txt"

#arg_vector$output = "M:/RADICL/Data/GRID/MDA-MB231/radiant_out/intersects/GRID_MDA-MB-231_cat_RNA_gene_intersect_maximums_cut_test.txt"

# Read in and process unique intersections

message('Reading unique intersects')

unique_intersect = data.table::fread(arg_vector$unique)

column_names = c('read_chr', 'read_start', 'read_stop', 'read_id', 'mapping_quality', 'read_strand',
                 'read_start_2', 'read_stop_2', 'rgb', 'misc', 'align_length', 'misc_2', 'gene_chr',
                 'gene_start', 'gene_stop', 'score', 'gene_id', 'gene_strand', 'intersect_length')

colnames(unique_intersect) = column_names

if(intersect_type == "proportion"){
  
  unique_intersect$gene_size = unique_intersect$gene_stop - unique_intersect$gene_start
  
  unique_intersect$intersect_proportion = unique_intersect$intersect_length / unique_intersect$gene_size
  
  setkey(unique_intersect, 'read_id')
  
  data.table::setorder(unique_intersect, 'read_id')
  
  max_intersect = unique_intersect[,.(gene_id = gene_id[which.max(intersect_proportion)]), by = 'read_id']

} else if (intersect_type == "expression"){
  
  cts_per_gene = unique_intersect[,.(sum_cts = length(read_id)), by = 'gene_id']
  
  unique_intersect$gene_count = cts_per_gene$sum_cts[match(unique_intersect$gene_id, cts_per_gene$gene_id)]
  
  max_intersect = unique_intersect[,.(gene_id = gene_id[which.max(gene_count)]), by = "read_id"]
  
} else if (intersect_type == "all"){
  
  max_intersect = unique_intersect
  
}

message('Processing unique intersects')


# dplyr_start = Sys.time()
# max_unique_intersect_dplyr = unique_intersect %>%
#     group_by(read_id) %>%
#     arrange(desc(intersect_proportion)) %>%
#     slice(1L) %>%
#     setDT()
# dplyr_end = Sys.time()
# dplyr_duration = dplyr_end - dplyr_start

cts_per_gene = max_intersect[,.(sum_cts = length(read_id)), by = 'gene_id']

# Read in and process multi-mapping intersects

message('Reading multi-mapping intersects')

multi_intersect = data.table::fread(arg_vector$multi)

colnames(multi_intersect) = column_names

if(multi_intersect_type == "propcount"){
  
  multi_intersect$unique_cts_per_gene = cts_per_gene$sum_cts[match(multi_intersect$gene_id, cts_per_gene$gene_id)]

  multi_intersect$unique_cts_per_gene[is.na(multi_intersect$unique_cts_per_gene)] = 0
  
  multi_intersect$gene_length = multi_intersect$gene_stop - multi_intersect$gene_start
  
  multi_intersect$intersect_proportion = multi_intersect$intersect_length/multi_intersect$gene_length
  
  message('Processing multi-mapping intersects')
  
  max_multi_intersect = multi_intersect[,.(gene_id = gene_id[which.max(unique_cts_per_gene*intersect_proportion)]), by = 'read_id']

} else if(multi_intersect_type == "all"){
  
  max_multi_intersect = multi_intersect
  
}

# max_multi_intersect = multi_intersect %>%
#     group_by(read_id) %>%
#     arrange(desc(unique_cts_per_gene), desc(intersect_proportion)) %>%
#     slice(1L) %>%
#     setDT()

max_intersect_bind = rbind(max_intersect[,c('read_id', 'gene_id')],
                           max_multi_intersect[,c('read_id', 'gene_id')])

max_intersect_bind = max_intersect_bind[order(max_intersect_bind$read_id, max_intersect_bind$gene_id),]                           

message('Writing max intersects')

data.table::fwrite(max_intersect_bind, file = arg_vector$output, col.names = FALSE, row.names = FALSE, sep = '\t', eol = '\n')
