# Library -----------------------------------------------------------------

library(dplyr)
library(data.table)
library(argparser)

# Set up argument parser

message('Parsing command line arguments')

parser = arg_parser(description='Retrieve maximal RNA-gene intersects from Red-C data')

parser = add_argument(parser, '--three', type = 'character', help = 'Maximum cut RNA3-gene intersects')

parser = add_argument(parser, '--uniquefive', type = 'character', help = 'RNA-gene intersects of unique-mapping RNA5 reads')

parser = add_argument(parser, '--multifive', type = 'character', help = 'RNA-gene intersects of multi-mapping RNA5 reads')

parser = add_argument(parser, '--output', type = 'character', help = 'Maximal intersects outfile')

# Parse command-line arguments

arg_vector = parse_args(parser)

# test assignment of Red-C 5' RNA reads to genes --------------------------

#rna3_file = 'M:/RADICL/Data/Red-C_JO/snakemake/results/intersects/JO500_Red-C_RNA3_gene_intersect_maximums_cut.txt'

message(paste("3' file is", arg_vector$three))

rna3_file = arg_vector$three

rna3_intersect = data.table::fread(rna3_file, header = F) %>%
  magrittr::set_colnames(c('Read', 'Gene'))

#rna5_unique_intersect_file = 'M:/RADICL/Data/Red-C_JO/snakemake/results/intersects/JO500_Red-C_RNA5_gene_unstranded_intersect.txt'

message("Reading 5' unique RNA intersects")

rna5_unique_intersect_file = arg_vector$uniquefive

rna5_unique_intersect = data.table::fread(rna5_unique_intersect_file, header = F)

#rna5_multi_intersect_file = 'M:/RADICL/Data/Red-C_JO/snakemake/results/intersects/JO500_Red-C_RNA5_gene_multi_unstranded_intersect.txt'

message("Reading 5' multi RNA intersects")

rna5_multi_intersect_file = arg_vector$multifive

rna5_multi_intersect_file = data.table::fread(rna5_multi_intersect_file, header = F)

rna5_intersect_bind = rbind(rna5_unique_intersect, rna5_multi_intersect_file)

rna3_counts = as.data.frame(table(rna3_intersect$Gene)) %>%
  magrittr::set_colnames(c('RNA', 'Count'))

rna5_intersect_bind$rna3_count = rna3_counts$Count[match(rna5_intersect_bind$V17, rna3_counts$RNA)]

rna5_maximums = rna5_intersect_bind[,.(Gene = V17[which.max(rna3_count)]), by = 'V4'] %>%
  dplyr::rename(Read = V4)

final_bind = rbind(rna3_intersect, rna5_maximums)

final_bind = final_bind[order(final_bind$Read, final_bind$Gene),]

message("Writing output")

data.table::fwrite(x = final_bind, file = arg_vector$output, sep = '\t', row.names = F, col.names = F, eol = '\n')

