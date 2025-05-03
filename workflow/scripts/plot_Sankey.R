# Combined Sankey for RNA & DNA reads -------------------------------------

# Unsplit fastq
# DNA fastq
# RNA fastq
# RNA alignment (unique mapping, multi-mapping, no mapping)
# DNA alignment (unique mapping, multi-mapping, no mapping)
# DNA-bin intersection (multi-rescued & unique & no intersection)
# RNA-bin intersection (multi-rescued & unique & no intersection)
# Valid RNA-DNA pair or no pair
# Sig interaction or non-sig interaction

# Library -----------------------------------------------------------------

library(ShortRead)
library(R.utils)
library(argparser)
library(data.table)
library(ggplot2)
library(ggsankey)

# Argument parsing setup -------------------------------------------------------

message('Parsing command line arguments')

parser = arg_parser(description='Plot RADIAnT pipeline overview as Sankey plot')

parser = add_argument(parser, '--unsplit', type = 'character', help = 'Unsplit FASTQ file (if available)')

parser = add_argument(parser, '--rnafastq', type = 'character', help = 'RNA FASTQ file')

parser = add_argument(parser, '--deplrnafastq', type = 'character', help = 'rRNA-depleted RNA FASTQ file')

parser = add_argument(parser, '--dnafastq', type = 'character', help = 'DNA FASTQ file')

parser = add_argument(parser, '--rnalog', type = 'character', help = 'RNA STAR log')

parser = add_argument(parser, '--dnalog', type = 'character', help = 'DNA STAR log')

parser = add_argument(parser, '--rnauniqueintersect', type = 'character', help = 'Unique RNA-gene intersect')

parser = add_argument(parser, '--rnamultiintersect', type = 'character', help = 'Multi-mapping RNA-gene intersect')

parser = add_argument(parser, '--dnaintersect', type = 'character', help = 'DNA-gene intersect')

parser = add_argument(parser, '--pairs', type = 'character', help = 'RNA-DNA intersect pairs')

parser = add_argument(parser, '--radiant', type = 'character', help = 'RADIAnT results file')

parser = add_argument(parser, '--output', type = 'character', help = 'Outfile name prefix')

# Parse arguments ---------------------------------------------------------

arg_list = parse_args(parser)

# Test files --------------------------------------------------------------

unsplit_fq = arg_list$unsplit

rna_fastq = arg_list$rnafastq

depleted_rna_fastq = arg_list$deplrnafastq

dna_fastq = arg_list$dnafastq

rna_star_log = arg_list$rnalog

dna_star_log = arg_list$dnalog

rna_gene_unique_intersect = arg_list$rnauniqueintersect

rna_gene_multi_intersect = arg_list$rnamultiintersect

dna_bin_unique_intersect = arg_list$dnaintersect

valid_pairs = arg_list$pairs

radiant_interactions = arg_list$radiant

output = arg_list$output

# Number of unsplit reads ----------------------------------------

if(!file.exists(unsplit_fq)){
  unsplit_reads = countFastq(unsplit_fq)$records
} else {
  unsplit_reads = NA
}

# Number of RNA reads -----------------------------------------------------

rna_reads = countFastq(rna_fastq)$records

# Number of depleted RNA reads --------------------------------------------

depleted_rna_reads = countFastq(depleted_rna_fastq)$records

# Number of DNA reads -----------------------------------------------------

# dna_reads = countFastq(dna_fastq)$records

# Alignment numbers DNA -------------------------------------------------------

dna_alignment = data.table::fread(input = dna_star_log, sep = '|', fill = T, header = F, skip = 4, col.names = c('Metric', 'Total'))

dna_reads = as.numeric(dna_alignment$Total[dna_alignment$Metric == 'Number of input reads'])

dna_unique = as.numeric(dna_alignment$Total[dna_alignment$Metric == 'Uniquely mapped reads number'])

dna_multi = as.numeric(dna_alignment$Total[dna_alignment$Metric == 'Number of reads mapped to multiple loci'])

dna_unmapped = sum(as.numeric(dna_alignment$Total[grepl('Number of reads unmapped', dna_alignment$Metric)]))

# Alignment numbers RNA ---------------------------------------------------

rna_alignment = data.table::fread(input = rna_star_log, sep = '|', fill = T, header = F, skip = 4, col.names = c('Metric', 'Total'))

rna_reads = as.numeric(rna_alignment$Total[rna_alignment$Metric == 'Number of input reads'])

rna_unique = as.numeric(rna_alignment$Total[rna_alignment$Metric == 'Uniquely mapped reads number'])

rna_multi = as.numeric(rna_alignment$Total[rna_alignment$Metric == 'Number of reads mapped to multiple loci'])

rna_unmapped = sum(as.numeric(rna_alignment$Total[grepl('Number of reads unmapped', rna_alignment$Metric)]))

# DNA intersects ----------------------------------------------------------

#dna_bin_intersect_reads = nrow(data.table::fread(dna_bin_unique_intersect))

dna_bin_intersect_reads = countLines(dna_bin_unique_intersect)[1]

dna_bin_non_intersects = dna_unique - dna_bin_intersect_reads

# RNA intersects ----------------------------------------------------------

rna_gene_unique_reads = length(unique(data.table::fread(rna_gene_unique_intersect, header = F)$V4))

rna_gene_multi_reads = length(unique(data.table::fread(rna_gene_multi_intersect, header = F)$V4))

rna_gene_multi_non_intersects = rna_multi - rna_gene_multi_reads

rna_gene_unique_non_intersects = rna_unique - rna_gene_unique_reads

# Valid pairs -------------------------------------------------------------

valid_read_pairs = countLines(valid_pairs)

# Interactions ------------------------------------------------------------

interactions = data.table::fread(radiant_interactions)

sig_interactions = sum(interactions$ReadCount[interactions$Padj < 0.05])

nonsig_interactions = sum(interactions$ReadCount[interactions$Padj >= 0.05])

# Plot Sankey -------------------------------------------------------------

# Splitting 

unsplit_to_dna = data.frame(
  x = 'Unsplit reads',
  Node = 'All',
  next_x = 'Split reads',
  next_node = 'DNA',
  nreads = rep(dna_reads, round(dna_reads/100000))
)

unsplit_to_rna = data.frame(
  x = 'Unsplit reads',
  Node = 'All',
  next_x = 'Split reads',
  next_node = 'RNA',
  nreads = rep(rna_reads+depleted_rna_reads, round((rna_reads+depleted_rna_reads)/100000))
)

# Alignment

dna_to_unique = data.frame(
  x = 'Split reads',
  Node = 'DNA',
  next_x = 'Aligned reads',
  next_node = 'DNA unique',
  nreads = rep(dna_unique, round(dna_unique/100000))
)

dna_to_multi = data.frame(
  x = 'Split reads',
  Node = 'DNA',
  next_x = 'Aligned reads',
  next_node = 'DNA multi',
  nreads = rep(dna_multi, round(dna_multi/100000))
)

dna_multi_to_na = data.frame(
  x = 'Aligned reads',
  Node = 'DNA multi',
  next_x = NA,
  next_node = NA,
  nreads = rep(dna_multi, round(dna_multi/100000))
)

dna_to_unmapped = data.frame(
  x = 'Split reads',
  Node = 'DNA',
  next_x = 'Aligned reads',
  next_node = 'DNA unmapped',
  nreads = rep(dna_unmapped, round(dna_unmapped/100000))
)

dna_unmapped_to_na = data.frame(
  x = 'Aligned reads',
  Node = 'DNA unmapped',
  next_x = NA,
  next_node = NA,
  nreads = rep(dna_unmapped, round(dna_unmapped/100000))
)

rna_to_unique = data.frame(
  x = 'Split reads',
  Node = 'RNA',
  next_x = 'Aligned reads',
  next_node = 'RNA unique',
  nreads = rep(rna_unique, round(rna_unique/100000))
)

rna_to_multi = data.frame(
  x = 'Split reads',
  Node = 'RNA',
  next_x = 'Aligned reads',
  next_node = 'RNA multi',
  nreads = rep(rna_multi, round(rna_multi/100000))
)

rna_to_none = data.frame(
  x = 'Split reads',
  Node = 'RNA',
  next_x = 'Aligned reads',
  next_node = 'RNA unmapped',
  nreads = rep(rna_unmapped, round(rna_unmapped/100000))
)

rna_unmapped_to_na = data.frame(
  x = 'Aligned reads',
  Node = 'RNA unmapped',
  next_x = NA,
  next_node = NA,
  nreads = rep(rna_unmapped, round(rna_unmapped/100000))
)

rna_to_rrna = data.frame(
  x = 'Split reads',
  Node = 'RNA',
  next_x = 'Aligned reads',
  next_node = 'rRNA',
  nreads = rep(depleted_rna_reads, round(depleted_rna_reads/100000))
)

rrna_to_na = data.frame(
  x = 'Aligned reads',
  Node = 'rRNA',
  next_x = NA,
  next_node = NA,
  nreads = rep(depleted_rna_reads, round(depleted_rna_reads/100000))
)

# Intersects

dna_unique_to_intersect = data.frame(
  x = 'Aligned reads',
  Node = 'DNA unique',
  next_x = 'Intersect reads',
  next_node = 'DNA-bin',
  nreads = rep(dna_bin_intersect_reads, round(dna_bin_intersect_reads/100000))
)

dna_unique_to_na = data.frame(
  x = 'Aligned reads',
  Node = 'DNA unique',
  next_x = 'Intersect reads',
  next_node = 'No intersect',
  nreads = rep(dna_unique - dna_bin_intersect_reads[1], round((dna_unique - dna_bin_intersect_reads[1])/100000))
)

rna_unique_to_intersect = data.frame(
  x = 'Aligned reads',
  Node = 'RNA unique',
  next_x = 'Intersect reads',
  next_node = 'RNA-gene',
  nreads = rep(rna_gene_unique_reads, round(rna_gene_unique_reads/100000))
)

rna_unique_to_non_intersect = data.frame(
  x = 'Aligned reads',
  Node = 'RNA unique',
  next_x = 'Intersect reads',
  next_node = 'No RNA-gene intersect',
  nreads = rep(rna_gene_unique_non_intersects, round(rna_gene_unique_non_intersects/100000))
)

rna_multi_to_intersect = data.frame(
  x = 'Aligned reads',
  Node = 'RNA multi',
  next_x = 'Intersect reads',
  next_node = 'RNA-gene',
  nreads = rep(rna_gene_multi_reads, round(rna_gene_multi_reads/100000))
)

rna_multi_to_non_intersect = data.frame(
  x = 'Aligned reads',
  Node = 'RNA multi',
  next_x = 'Intersect reads',
  next_node = 'No RNA-gene intersect',
  nreads = rep(rna_gene_multi_non_intersects, round(rna_gene_multi_non_intersects/100000))
)

rna_non_intersect_to_na = data.frame(
  x = 'Intersect reads',
  Node = 'No RNA-gene intersect',
  next_x = NA,
  next_node = NA,
  nreads = rep(rna_gene_unique_non_intersects + rna_gene_multi_non_intersects, round((rna_gene_unique_non_intersects + rna_gene_multi_non_intersects)/100000))
)

# Valid pairs

dna_intersect_to_pairs = data.frame(
  x = 'Intersect reads',
  Node = 'DNA-bin',
  next_x = 'Valid pairs',
  next_node = 'Valid pairs',
  nreads = rep(valid_read_pairs[1], round(valid_read_pairs[1]/100000))
)

dna_intersect_to_no_pair = data.frame(
  x = 'Intersect reads',
  Node = 'DNA-bin',
  next_x = 'Valid pairs',
  next_node = 'No RNA partner',
  nreads = rep(dna_bin_intersect_reads[1] - valid_read_pairs[1], round((dna_bin_intersect_reads[1] - valid_read_pairs[1])/100000))
)

dna_no_pair_to_na = data.frame(
  x = 'Valid pairs',
  Node = 'No RNA partner',
  next_x = NA,
  next_node = NA,
  nreads = rep(dna_bin_intersect_reads[1] - valid_read_pairs[1], round((dna_bin_intersect_reads[1] - valid_read_pairs[1])/100000))
)

rna_intersect_to_pairs = data.frame(
  x = 'Intersect reads',
  Node = 'RNA-gene',
  next_x = 'Valid pairs',
  next_node = 'Valid pairs',
  nreads = rep(valid_read_pairs[1], round(valid_read_pairs[1]/100000))
)

rna_intersect_to_no_pair = data.frame(
  x = 'Intersect reads',
  Node = 'RNA-gene',
  next_x = 'Valid pairs',
  next_node = 'No DNA partner',
  nreads = rep((rna_gene_multi_reads + rna_gene_unique_reads) - valid_read_pairs[1], round(((rna_gene_multi_reads + rna_gene_unique_reads) - valid_read_pairs[1])/100000))
)

rna_no_pair_to_na = data.frame(
  x = 'Valid pairs',
  Node = 'No DNA partner',
  next_x = NA,
  next_node = NA,
  nreads = rep((rna_gene_multi_reads + rna_gene_unique_reads) - valid_read_pairs[1], round(((rna_gene_multi_reads + rna_gene_unique_reads) - valid_read_pairs[1])/100000))
)

# Interactions

pairs_to_sig_interactions = data.frame(
  x = 'Valid pairs',
  Node = 'Valid pairs',
  next_x = 'Interactions',
  next_node = 'Significant interactions',
  nreads = rep(sig_interactions, round(sig_interactions/100000))
)

pairs_to_nonsig_interactions = data.frame(
  x = 'Valid pairs',
  Node = 'Valid pairs',
  next_x = 'Interactions',
  next_node = 'Non-significant interactions',
  nreads = rep(nonsig_interactions, round(nonsig_interactions/100000))
)

# Interactions to NA

sig_interactions_to_na = data.frame(
  x = 'Interactions',
  Node = 'Significant interactions',
  next_x = NA,
  next_node = NA,
  nreads = rep(sig_interactions, round(sig_interactions/100000))
)

nonsig_interactions_to_na = data.frame(
  x = 'Interactions',
  Node = 'Non-significant interactions',
  next_x = NA,
  next_node = NA,
  nreads = rep(nonsig_interactions, round(nonsig_interactions/100000))
)
# Combine to data frame

sankey_list = list(
  unsplit_to_dna,
  unsplit_to_rna,
  dna_to_unique,
  dna_to_multi,
  dna_multi_to_na,
  dna_to_unmapped,
  dna_unmapped_to_na,
  rna_to_unique, 
  rna_to_multi,
  rna_to_none,
  rna_to_rrna,
  rrna_to_na,
  rna_unmapped_to_na,
  dna_unique_to_intersect,
  dna_unique_to_na,
  rna_unique_to_intersect,
  rna_unique_to_non_intersect,
  rna_multi_to_intersect,
  rna_multi_to_non_intersect,
  rna_non_intersect_to_na,
  dna_intersect_to_pairs,
  dna_intersect_to_no_pair,
  dna_no_pair_to_na,
  rna_intersect_to_pairs,
  rna_intersect_to_no_pair,
  rna_no_pair_to_na,
  pairs_to_sig_interactions,
  pairs_to_nonsig_interactions,
  sig_interactions_to_na,
  nonsig_interactions_to_na
)

sankey_df = data.table::rbindlist(sankey_list)

# Set factor levels

sankey_df$x = factor(sankey_df$x, levels = c('Unsplit reads', 'Split reads', 'Aligned reads', 'Intersect reads', 'Valid pairs', 'Interactions'))

sankey_df$next_x = factor(sankey_df$next_x, levels = c('Unsplit reads', 'Split reads', 'Aligned reads', 'Intersect reads', 'Valid pairs', 'Interactions'))

node_levels = c('All', 'RNA', 'DNA', 'rRNA', 'RNA unmapped', 'RNA multi', 'RNA unique', 'DNA unique', 'DNA multi', 'DNA unmapped', 'No RNA-gene intersect', 'RNA-gene', 'DNA-bin', 'No DNA partner', 'Valid pairs', 'No RNA partner', 'Significant interactions', 'Non-significant interactions')

sankey_df$Node = factor(sankey_df$Node, levels = rev(node_levels))

sankey_df$next_node = factor(sankey_df$next_node, levels = rev(node_levels))

# Plot Sankey itself

sankey_plot = ggplot(data = sankey_df, mapping = aes(x = x, next_x = next_x, node = Node, next_node = next_node, fill = Node)) +
  geom_sankey(flow.alpha = 0.2, flow.fill = 'black', node.alpha = 0.5) +
  #geom_sankey_text() +
  theme_sankey() +
  geom_sankey_text(mapping = aes(group = Node, label = Node), size = 7/.pt) +
  scale_fill_viridis_d() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = 'black', size = 9))

# Output

# Plots

ggsave(plot = sankey_plot, filename = paste0(output, 'Sankey_updated.svg'), width = 6, height = 4)
ggsave(plot = sankey_plot, filename = paste0(output, 'Sankey_updated.png'), width = 6, height = 4, units = 'in', dpi = 600)

# Underlying numbers

sankey_stats = na.omit(unique(sankey_df))

sankey_totals = sankey_stats[,.(Total_reads = sum(nreads)), by = 'next_node']

0.5*(sankey_totals$Total_reads[sankey_totals$next_node=='Valid pairs']) / sankey_totals$Total_reads[sankey_totals$next_node=='DNA']

data.table::fwrite(x = na.omit(unique(sankey_df)), file = paste0(output, 'Sankey_updated_stats.txt'), sep = '\t', col.names = T, row.names = F)
