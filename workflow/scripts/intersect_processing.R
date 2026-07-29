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

# Tie-breaking rules
parser = add_argument(parser, '--intersect-rule-unique', type = 'character', help = 'Which criterion to define the main gene by. "proportion" = longest overlap with read; "count" = highest read count (correlate of expression); "all" = return all overlaps)')
parser = add_argument(parser, '--intersect-rule-multi', type = 'character', help = 'Which criterion to define the main gene by. "proportion" = longest overlap with read; "count" = highest read count (correlate of expression); "all" = return all overlaps)')

# Parse command-line arguments
arg_vector = parse_args(parser)

###############################################################################################################


# Read in and process unique intersections
message('Reading unique intersects')
unique_intersect = data.table::fread(arg_vector$unique)
column_names = c('read_chr', 'read_start', 'read_stop', 'read_id', 'mapping_quality', 'read_strand',
                 'read_start_2', 'read_stop_2', 'rgb', 'misc', 'align_length', 'misc_2', 'gene_chr',
                 'gene_start', 'gene_stop', 'score', 'gene_id', 'gene_strand', 'intersect_length')
colnames(unique_intersect) = column_names
# -> table with 1 entry per read, giving information about the read, as well as the quality of mapping and the characteristics of the overlap with each respective gene 


message('Processing unique intersects')
if(arg_vector$intersect_rule_unique == "proportion"){
  # ----------- Per uniquely mapping read, assign it to the gene for which this read -----------
  # covers the largest fraction of the gene body. ----------------------------------------------
  
  
  # Compute gene length (as 100%)
  unique_intersect$gene_length = unique_intersect$gene_stop - unique_intersect$gene_start
  
  # Compute intersect length as proportion of gene length
  unique_intersect$intersect_proportion = unique_intersect$intersect_length / unique_intersect$gene_length
  
  # Order by read ID
  setkey(unique_intersect, 'read_id')  
  data.table::setorder(unique_intersect, 'read_id')
  # Per read, identify the gene with the maximum overlap (% gene length): return table ofgene - read pairs (two-columned table with col.names =  read_id, gene_id - of the maximum overlap)
  max_unique_intersect = unique_intersect[,.(gene_id = gene_id[which.max(intersect_proportion)]), by = 'read_id']

} else if (arg_vector$intersect_rule_unique == "count"){
  # -----------  Use only unique-read support per gene as an evidence of expression/interaction strength. ----------
  # Per read, assign it to the gene that has the highest total number of uniquely mapped reads, --------------------
  # regardless of the exact overlap proportion for this specific read. ---------------------------------------------
  
  # Get total number of reads mapping to gene (two-columned table with col.names =  gene_id, sum_cts, where sum_cts stands for the summed read counts per gene)
  cts_per_gene = unique_intersect[,.(sum_cts = length(read_id)), by = 'gene_id']
 
  # Add read count per gene info to each unique intersect
  unique_intersect$gene_count = cts_per_gene$sum_cts[match(unique_intersect$gene_id, cts_per_gene$gene_id)]
  
  # Order by read ID
  setkey(unique_intersect, 'read_id')  
  data.table::setorder(unique_intersect, 'read_id')
  # Identify the gene with the highest read count support (correlate of expression) and assign the read to this gene
  max_unique_intersect = unique_intersect[,.(gene_id = gene_id[which.max(gene_count)]), by = "read_id"]
  
} else if (arg_vector$intersect_rule_unique == "propcount") {
  
  # Hybrid rule for unique-mapping reads:
  # combine global unique-read support per gene with local overlap proportion.
  # Per read, assign to the gene maximizing (gene_count * intersect_proportion).
  
  # Compute gene length and overlap proportion as in "proportion"
  unique_intersect$gene_length = unique_intersect$gene_stop - unique_intersect$gene_start
  unique_intersect$intersect_proportion = unique_intersect$intersect_length / unique_intersect$gene_length
  
  # Total unique-read counts per gene as in "count"
  cts_per_gene = unique_intersect[
    , .(sum_cts = length(read_id)),
    by = 'gene_id'
  ]
  unique_intersect$gene_count = cts_per_gene$sum_cts[
    match(unique_intersect$gene_id, cts_per_gene$gene_id)
  ]
  
  # Per read, pick gene with maximal combined score
  setkey(unique_intersect, 'read_id')
  data.table::setorder(unique_intersect, 'read_id')
  max_unique_intersect = unique_intersect[
    , .(gene_id = gene_id[which.max(gene_count * intersect_proportion)]),
    by = 'read_id'
  ] } else if (arg_vector$intersect_rule_unique == "all"){
    
  # -------------------- Keep all gene overlaps for each multi-mapping read ------------
  # (no tie-breaking). -----------------------------------------------------------------
   

  # Order by read ID
  setkey(unique_intersect, 'read_id')
  data.table::setorder(unique_intersect, 'read_id')
  # Keep all intersections of each read with any genes it maps to
  max_unique_intersect = unique_intersect[,c("read_id", "gene_id")]
}

###############################################################################################################



cts_per_gene = max_unique_intersect[,.(sum_cts = length(read_id)), by = 'gene_id']

# Read in and process multi-mapping intersects

message('Reading multi-mapping intersects')

multi_intersect = data.table::fread(arg_vector$multi)
colnames(multi_intersect) = column_names

message('Processing multi-mapping intersects')

if (arg_vector$intersect_rule_multi == "proportion") {
  # ----------- For multi-mapping reads, pick per read the gene where this read ---------
  # covers the largest fraction of the gene body. --------------------------------------- 
  
  # compute gene length for each intersect
  multi_intersect$gene_length = multi_intersect$gene_stop - multi_intersect$gene_start
  
  # intersect length as proportion of gene length
  multi_intersect$intersect_proportion = multi_intersect$intersect_length / multi_intersect$gene_length
  
  # per read, select the gene with the max overlap proportion
  setkey(multi_intersect, 'read_id')
  data.table::setorder(multi_intersect, 'read_id')
  max_multi_intersect = multi_intersect[
    intersect_proportion, .(gene_id = gene_id[which.max(intersect_proportion)]),
  by = 'read_id'
  ]
} else if (arg_vector$intersect_rule_multi == "count"){
  # -----------  For multi-mapping reads, use unique-read support per gene as a prior. ---------
  # Per read, assign to the gene that has the highest number of uniquely mapped reads. ---------
  
  # Total unique-read counts per gene (precomputed from unique_intersect)
  cts_per_gene = unique_intersect[
    , .(sum_cts = length(read_id)),
    by = 'gene_id'
  ]
  
  # Attach unique counts to each multi intersect
  multi_intersect$unique_cts_per_gene = cts_per_gene$sum_cts[
    match(multi_intersect$gene_id, cts_per_gene$gene_id)
  ]
  
  # Genes not seen in the unique table get 0 counts
  multi_intersect$unique_cts_per_gene[is.na(multi_intersect$unique_cts_per_gene)] = 0
  
  # Per read, pick the gene with maximal unique support
  setkey(multi_intersect, 'read_id')
  data.table::setorder(multi_intersect, 'read_id')
  max_multi_intersect = multi_intersect[
    , .(gene_id = gene_id[which.max(unique_cts_per_gene)]),
    by = 'read_id'
  ]
  
} else if(arg_vector$intersect_rule_multi == "propcount"){
  
  # -----------  Hybrid rule for multi-mapping reads: ------------------------------------------------------
  # combine global unique-read support per gene with local overlap proportion. -----------------------------
  # Per read, assign to the gene maximizing (unique_cts_per_gene * intersect_proportion). ------------------
  
  multi_intersect$unique_cts_per_gene = cts_per_gene$sum_cts[match(multi_intersect$gene_id, cts_per_gene$gene_id)]
  
  multi_intersect$unique_cts_per_gene[is.na(multi_intersect$unique_cts_per_gene)] = 0
  
  multi_intersect$gene_length = multi_intersect$gene_stop - multi_intersect$gene_start
  
  multi_intersect$intersect_proportion = multi_intersect$intersect_length/multi_intersect$gene_length
  
  message('Processing multi-mapping intersects')
  
  max_multi_intersect = multi_intersect[,.(gene_id = gene_id[which.max(unique_cts_per_gene*intersect_proportion)]), by = 'read_id']
  
} else if(arg_vector$intersect_rule_multi == "all"){
  # Keep all gene overlaps for each multi-mapping read (no tie-breaking).
  max_multi_intersect = multi_intersect[, c("read_id", "gene_id")]
  
}

# max_multi_intersect = multi_intersect %>%
#     group_by(read_id) %>%
#     arrange(desc(unique_cts_per_gene), desc(intersect_proportion)) %>%
#     slice(1L) %>%
#     setDT()

max_intersect_bind = rbind(max_unique_intersect[,c('read_id', 'gene_id')],
                           max_multi_intersect[,c('read_id', 'gene_id')])

max_intersect_bind = max_intersect_bind[order(max_intersect_bind$read_id, max_intersect_bind$gene_id),]                           

message('Writing max intersects')

data.table::fwrite(max_intersect_bind, file = arg_vector$output, col.names = FALSE, row.names = FALSE, sep = '\t', eol = '\n')



