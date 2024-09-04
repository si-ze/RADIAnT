# Usage: 
# Rscript workflow/scripts/plot_circos.r --results /path/to/my_RADIAnT_results.txt --goi myGOI --genome mygenome --cytoband resources/myorganism/my.cytoBand.txt.gz

# Library -----------------------------------------------------------------

library(ggplot2)
library(GenomicRanges)
library(ggforce)
library(dplyr)


# Parameters --------------------------------------------------------------


library(argparser)

parser = arg_parser(description='Visualise GOI interactions across chromosomes in circular plot')

parser = add_argument(parser, '--results', type = 'character', help = 'Absolute path to RADIAnT results file (*_RADIAnT_results.txt)')

parser = add_argument(parser, '--goi', type = 'character', help = 'Gene symbol of gene of interest')

parser = add_argument(parser, '--genome', type = 'character', help = 'Genome version (currently supported: hg38 and mm39)')

parser = add_argument(parser, '--cytoband', type = 'character', help = 'Absolute path to cytoband file. Also provided in resources/ directory together with RADIAnT')

parser = add_argument(parser, '--outdir', type = 'character', help = 'Desired output directory', default='.')

parser = add_argument(parser, '--outformat', type = 'character', help = 'Desired output file format. Supported: svg / png / both. Default (also if omitted): both', default='both')

parser = add_argument(parser, '--spacer', type = 'integer', help = 'Spacer size. Default: 20000000', default=20000000)

arg_vector = parse_args(parser)



res = data.table::fread(arg_vector$results) %>% as.data.frame()
goi = arg_vector$goi
genome = arg_vector$genome
cytoband_file = arg_vector$cytoband

outdir = arg_vector$outdir
output_format = arg_vector$outformat
spacer_size = arg_vector$spacer


# args = commandArgs(trailingOnly = TRUE)
# res=data.table::fread(args[1]) # %>% as.data.frame()
# goi=args[2]
# genome=args[3]
# cytoband_file=args[4]

# output_format="both"
# if (length(args)>4) {
#   output_format=args[5]
# }

# spacer_size=20000000

# genome = "mm39"
# cytoband_file = '/mnt/d/RADIAnT/resources/mouse/mm10.cytoBand.txt.gz'
# # cytoband_file = '/mnt/d/RADIAnT/resources/mouse/mm39.cytoBandIdeo.txt.gz'

# experiment_name = 'test_mm10'

# res_file = '/mnt/y/Fileserver_NGS/Current/RNA_DNA_interactions/RADICL_mESC/RADICL_mESC_2FA/interactions/mESC_2FA_n1_RADIAnT_results.txt'

# rna = 'Malat1'






# Get chromosome sizes ----------------------------------------------------

chroms = getChromInfoFromUCSC(genome = genome, assembled.molecules.only = T)[,1:2]

# Create spacer DF --------------------------------------------------------

spacer_df = data.frame(
  chrom = 'spacer',
  size = spacer_size
)

# Assemble into geom_rect compatible --------------------------------------

rect_df = do.call('rbind', lapply(1:nrow(chroms), function(chrom){
  return(rbind(chroms[chrom,], spacer_df))
}))

rect_df$end = cumsum(as.numeric(rect_df$size))

rect_df$start = c(0, rect_df$end[1:(nrow(rect_df)-1)])

rect_df$class = ifelse(rect_df$chrom == 'spacer', yes = 'spacer', no = 'chr')

# convert to degrees

rect_df$degrees_start = 360*(rect_df$start/max(rect_df$end))

rect_df$degrees_end = 360*(rect_df$end/max(rect_df$end))

# convert to radians

rect_df$radians_start = 2*pi*rect_df$degrees_start/360

rect_df$radians_end = 2*pi*rect_df$degrees_end/360

# Create label df ---------------------------------------------------------

rect_df$label_radians = rowMeans(rect_df[,c('radians_start', 'radians_end')])

rect_df$label_x = 0.685 * cos(rect_df$label_radians)

rect_df$label_y = 0.685 * sin(rect_df$label_radians)

# Add cytoband information ------------------------------------------------

cytoband = read.table(cytoband_file, fill = T)

colnames(cytoband) = c('chr', 'start', 'end', 'name', 'type')

cytoband$type = sub('pos.*', 'pos', cytoband$type)

cytoband$type[cytoband$type == ""] = 'gneg'

cytoband = cytoband[cytoband$chr %in% rect_df$chrom,]

cytoband$Cum_chrom_start = rect_df$start[match(cytoband$chr, rect_df$chrom)]

cytoband$Cum_band_start = cytoband$Cum_chrom_start + cytoband$start

cytoband$Cum_band_end = cytoband$Cum_chrom_start + cytoband$end

cytoband$Degrees_band_start = 360 * (cytoband$Cum_band_start / max(rect_df$end))

cytoband$Degrees_band_end = 360 * (cytoband$Cum_band_end / max(rect_df$end))

cytoband$Radians_band_start = 2 * pi * cytoband$Degrees_band_start / 360

cytoband$Radians_band_end = 2 * pi * cytoband$Degrees_band_end / 360

# Arcs from genomic coordinates -------------------------------------------



interactions = res[res$Symbol == goi & res$BinChr != res$GeneChr & res$Padj < 0.05,]

arc_df = data.frame(
  Start_chr = interactions$GeneChr,
  Start_coord = rowMeans(interactions[,c('GeneLeft', 'GeneRight')]),
  End_chr = interactions$BinChr,
  End_coord = interactions$BinCentre,
  LogPadj = -log10(interactions$Padj)
)

arc_df$LogPadj[is.infinite(arc_df$LogPadj)] = max(arc_df$LogPadj[is.finite(arc_df$LogPadj)])

arc_df$Cum_chr_start = rect_df$start[match(arc_df$Start_chr, rect_df$chrom)]

arc_df$Cum_chr_end = rect_df$start[match(arc_df$End_chr, rect_df$chrom)]

arc_df$Cum_arc_start = arc_df$Cum_chr_start + arc_df$Start_coord

arc_df$Cum_arc_end = arc_df$Cum_chr_end + arc_df$End_coord

arc_df$Degrees_start = 360 * (arc_df$Cum_arc_start/max(rect_df$end))

arc_df$Degrees_end = 360 * (arc_df$Cum_arc_end/max(rect_df$end))

arc_df$Curve_threshold = arc_df$Degrees_start - 180

arc_df$Degrees_difference = (arc_df$Degrees_start - arc_df$Degrees_end) %% 360

arc_df$Curvature_weight = abs(arc_df$Degrees_difference - 180)

arc_df$Curvature = scales::rescale(arc_df$Curvature_weight, to = c(0,0.5))

arc_df$Curvature = ifelse(arc_df$Degrees_difference < 180, yes = -arc_df$Curvature, no = arc_df$Curvature)

arc_df$Radians_start = 2 * pi * arc_df$Degrees_start / 360

arc_df$Radians_end = 2 * pi * arc_df$Degrees_end / 360

arc_df$X = 0.55 * cos(arc_df$Radians_start)

arc_df$Xend = 0.55 * cos(arc_df$Radians_end)

arc_df$Y = 0.55 * sin(arc_df$Radians_start)

arc_df$Yend = 0.55 * sin(arc_df$Radians_end)

arc_df = arc_df[arc_df$LogPadj > 10,]


# Plot --------------------------------------------------------------------

plot = ggplot() +
  lapply(1:nrow(arc_df), function(arc){
    #geom_curve(arc_df[arc,], mapping = aes(x = Y, xend = Yend, y = X, yend = Xend, colour = Condition, alpha = LogPadj, size = LogPadj), curvature = arc_df$Curvature[arc])
    geom_curve(arc_df[arc,], mapping = aes(x = Y, xend = Yend, y = X, yend = Xend, alpha = LogPadj, size = LogPadj), curvature = arc_df$Curvature[arc])
  }) +
  geom_arc_bar(data = cytoband, mapping = aes(x0 = 0, y0 = 0, r0 = 0.55, r = 0.6, start = Radians_band_start, end = Radians_band_end, fill = type), colour = NA) +
  geom_arc_bar(data = rect_df[rect_df$class=='chr',], mapping = aes(x0 = 0, y0 = 0, r0 = 0.55, r = 0.6, start = radians_start, end = radians_end), fill = NA) +
  coord_fixed() +
  geom_text(data = rect_df[rect_df$class=='chr',], mapping = aes(x = label_y, y = label_x, label = chrom), size = 7/.pt) +
  #facet_wrap(~Condition) +
  scale_fill_manual(breaks = c('gneg', 'gpos', 'gvar', 'acen', 'stalk'), values = c('white', 'grey30', 'grey', 'darkred', 'steelblue4')) +
  scale_colour_manual(values = c('black', 'darkred')) +
  scale_size_continuous(range = c(.1,1)) +
  theme_void() +
  theme(legend.position = 'none',
        panel.spacing = unit(2,'lines'))

# Export ------------------------------------------------------------------


if(!dir.exists(outdir)) {dir.create(outdir, recursive = TRUE)}

my_filename=paste0(outdir, "/", goi,  ".", genome, ".chromosome_plot")



if (stringr::str_to_lower(output_format) =="svg") {
  ggsave(plot = plot, filename=paste0(my_filename, ".svg"), width = 2.5, height = 2.5, device=svg)
} else if (stringr::str_to_lower(output_format)=="png") {
  ggsave(plot = plot, filename=paste0(my_filename, ".png"), width = 2.5, height = 2.5)
} else {
  ggsave(plot = plot, filename=paste0(my_filename, ".svg"), width = 2.5, height = 2.5, device=svg)
  ggsave(plot = plot, filename=paste0(my_filename, ".png"), width = 2.5, height = 2.5)
}
