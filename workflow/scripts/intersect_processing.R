library(data.table)
library(dplyr)
library(argparser)

# Set up argument parser

message('Parsing command line arguments')

parser = arg_parser(description='Retrieve maximal RNA-gene intersects from RNA-DNA ligation data')

# Per read, one entry describing its mapping to a gene
parser = add_argument(parser, '--unique', type = 'character', help = 'RNA-gene intersects of uniquely mapping reads')

# Per read, one or more entries describing its mapping to a gene
parser = add_argument(parser, '--multi', type = 'character', help = 'RNA-gene intersects of multi-mapping reads')

# Output file 
parser = add_argument(parser, '--output', type = 'character', help = 'Maximal intersects outfile')
# format: two-column read_id, gene_id ("max" / best intersect)


# Parse command-line arguments

arg_vector = parse_args(parser)


# Read in and process unique intersections
message('Reading unique intersects')

unique_intersect = data.table::fread(arg_vector$unique)

column_names = c('read_chr', 'read_start', 'read_stop', 'read_id', 'mapping_quality', 'read_strand',
                 'read_start_2', 'read_stop_2', 'rgb', 'misc', 'align_length', 'misc_2', 'gene_chr',
                 'gene_start', 'gene_stop', 'score', 'gene_id', 'gene_strand', 'intersect_length')
colnames(unique_intersect) = column_names
# -> table with 1 entry per read, giving information about the read, as well as the quality of mapping and the characteristics of the overlap with each respective gene 


if(intersect_type == "proportion"){
   # Compute gene length (as 100%) 
  unique_intersect$gene_size = unique_intersect$gene_stop - unique_intersect$gene_start

  # Compute intersect length as proportion of gene length
  unique_intersect$intersect_proportion = unique_intersect$intersect_length / unique_intersect$gene_size

  # Order by read ID
  setkey(unique_intersect, 'read_id')
  data.table::setorder(unique_intersect, 'read_id')

  # Per read, identify the gene with the maximum overlap (% gene length): return table ofgene - read pairs (two-columned table with col.names =  read_id, gene_id - of the maximum overlap)  
  max_intersect = unique_intersect[,.(gene_id = gene_id[which.max(intersect_proportion)]), by = 'read_id']

} else if (intersect_type == "expression"){

  # Get total number of reads mapping to gene (two-columned table with col.names =  gene_id, sum_cts, where sum_cts stands for the summed read counts per gene)  
  cts_per_gene = unique_intersect[,.(sum_cts = length(read_id)), by = 'gene_id']
  
  # Add read count per gene info to each unique intersect
  unique_intersect$gene_count = cts_per_gene$sum_cts[match(unique_intersect$gene_id, cts_per_gene$gene_id)]
  
  # Identify the gene with the highest read count support (correlate of expression) and assign the read to this gene
  max_intersect = unique_intersect[,.(gene_id = gene_id[which.max(gene_count)]), by = "read_id"]
  
} else if (intersect_type == "all"){
  
  max_intersect = unique_intersect
  
}

message('Processing unique intersects')


# Get total number of reads mapping to gene (two-columned table with col.names =  gene_id, sum_cts, where sum_cts stands for the summed read counts per gene)
cts_per_gene = max_intersect[,.(sum_cts = length(read_id)), by = 'gene_id']

# Read in and process multi-mapping intersects

message('Reading multi-mapping intersects')

multi_intersect = data.table::fread(arg_vector$multi)

colnames(multi_intersect) = column_names

if(multi_intersect_type == "propcount"){
  
  # Refer to the read counts from the unique mapping to assign total number of reads per gene 
  multi_intersect$unique_cts_per_gene = cts_per_gene$sum_cts[match(multi_intersect$gene_id, cts_per_gene$gene_id)]

  # For genes which might not have been featured in the unique counts, set NA to 0 counts
  multi_intersect$unique_cts_per_gene[is.na(multi_intersect$unique_cts_per_gene)] = 0
  
  # Compute gene length (as 100%)
  multi_intersect$gene_length = multi_intersect$gene_stop - multi_intersect$gene_start
  
  # Compute intersect length as proportion of gene length   
  multi_intersect$intersect_proportion = multi_intersect$intersect_length/multi_intersect$gene_length
  
  message('Processing multi-mapping intersects')
  
  max_multi_intersect = multi_intersect[,.(gene_id = gene_id[which.max(unique_cts_per_gene*intersect_proportion)]), by = 'read_id']

} else if(multi_intersect_type == "all"){
  
  max_multi_intersect = multi_intersect
  
}


max_intersect_bind = rbind(max_intersect[,c('read_id', 'gene_id')],
                           max_multi_intersect[,c('read_id', 'gene_id')])

max_intersect_bind = max_intersect_bind[order(max_intersect_bind$read_id, max_intersect_bind$gene_id),]                           

message('Writing max intersects')

data.table::fwrite(max_intersect_bind, file = arg_vector$output, col.names = FALSE, row.names = FALSE, sep = '\t', eol = '\n')
