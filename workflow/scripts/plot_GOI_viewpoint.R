# Usage: 
# Rscript workflow/scripts/plot_gene_viewpoint.r --results /path/to/my_RADIAnT_results.txt --goi myGOI --binAnnotation resources/myorganism/my_bins_named.bed.gz --chr chromosomeName --start coordinate --end coordinate 


library(ggplot2)
library(slider)
library(stringr)
library(dplyr)
library(data.table)
library(GenomicRanges)


# res=data.table::fread("/mnt/d/Emma/JUN22_12/interactions/JUN22_12_RADIAnT_results.txt")
# goi = "MALAT1"
# binAnnotation = data.table:fread("/mnt/d/RADIAnT/resources/human/hg38_5kb_bins_named.bed.gz")
# # MALAT1 on human Chromosome 11: 65,497,640-65,508,073 forward strand. (8Mb to the left and to the right)
# chr = 'chr11'
# start=57500000
# stop=74000000


# res=data.table::fread("/mnt/y/Fileserver_NGS/Current/RNA_DNA_interactions/RADICL_mESC/RADICL_mESC_2FA/interactions/mESC_2FA_n1_RADIAnT_results.txt")
# res=data.table::fread("/mnt/y/Fileserver_NGS/Current/RNA_DNA_interactions/RADICL_mESC/RADICL_mESC_2FA/interactions/RADICL_mESC_2FA_pooled_RADIAnT_results.txt")
# binAnnotation = data.table:fread("/mnt/d/RADIAnT/resources/mouse/mm39_5kb_bins_named.bed.gz")
# goi = "Malat1"
# chr = 'chr19'
# start = 3500000
# stop = 11000000






library(argparser)

parser = arg_parser(description='Visualise GOI interactions across chromosomes in circular plot')

parser = add_argument(parser, '--results', type = 'character', help = 'Absolute path to RADIAnT results file (*_RADIAnT_results.txt)')

parser = add_argument(parser, '--goi', type = 'character', help = 'Gene symbol of gene of interest')

parser = add_argument(parser, '--genome', type = 'character', help = 'Genome version (currently supported: hg38 and mm39)')

parser = add_argument(parser, '--binAnnotation', type = 'character', help = 'Absolute path to bin annotation file. Also provided in resources/ directory together with RADIAnT')

parser = add_argument(parser, '--chr', type = 'character', help = 'Which chromosome to plot')

parser = add_argument(parser, '--start', type = 'integer', help = 'Start coordinate of region to be plotted')

parser = add_argument(parser, '--end', type = 'integer', help = 'End coordinate of region to be plotted')

parser = add_argument(parser, '--outdir', type = 'character', help = 'Desired output directory', default='.')

parser = add_argument(parser, '--outformat', type = 'character', help = 'Desired output file format. Supported: svg / png / both. Default (also if omitted): both', default='both')

arg_vector = parse_args(parser)



results = data.table::fread(arg_vector$results) %>% as.data.frame()
goi = arg_vector$goi
genome = arg_vector$genome
binAnnotation = data.table::fread(arg_vector$binAnnotation) %>% as.data.frame() %>% setNames( c('BinChr', 'BinStart', 'BinEnd', 'Bin')) %>% dplyr::mutate(Centre = rowMeans(.[,c('BinStart', 'BinEnd')]))
chr = arg_vector$chr
outdir = arg_vector$outdir
output_format = arg_vector$outformat
start=as.numeric(arg_vector$start)
stop=as.numeric(arg_vector$end)




theme_TW = function(){
  theme_bw() +
    theme(axis.text = element_text(colour = 'black', size=12),
          axis.line = element_line(),
          panel.border = element_blank(),
          axis.ticks = element_line(colour = 'black'),
          strip.background = element_rect(colour = NA))
}


theme_top = function(){
  theme_void() +
    theme(text=element_text(size=7), 
          axis.line.y.left = element_line(colour = 'black', lineend = 'square'),
          axis.ticks.y.left = element_line(colour = 'black'),
          axis.text.y.left = element_text(colour = 'black'),
          axis.title.y.left = element_text(colour = 'black'),
          legend.text = element_text(colour = 'black'),
          legend.title = element_text(colour = 'black'),
          plot.margin = unit(c(0,0,0,0), 'lines'), 
          plot.background = element_rect(fill="white", colour = NA   ))
}



plotResults = function(
  results, 
  chr, 
  start, 
  stop, 
  goi, 
  threshold=0.05, 
  binAnnotation, 
  additionalTrackFiles = NULL,
  additionalTrackNames = NULL,
  additionalTrackTypes = NULL,
  outFile = 'viewpoint.svg', 
  smoothing_windows = 0, 
  span = 0.01, 
  ylim = NULL, 
  name = ""){

  # subset results to provided region
  
  goiData = setDT(results[results$Symbol==goi,])
  
  message(paste('Effective', goi, 'sequencing depth is:', sum(goiData$ReadCount)))
  
  effectiveSeqDepth = sum(goiData$ReadCount)
  
  # goiData$Method
  
  goiCisCount = sum(goiData$ReadCount[goiData$OriginalMethod=='Cis'])
  
  goiTransCount = sum(goiData$ReadCount[goiData$OriginalMethod=='Trans'])
  
  goiData$TotalCount = ifelse(goiData$Method == 'Bin', yes = goiTransCount, no = goiCisCount)
  

  regionData = goiData %>%
    dplyr::filter(BinChr == chr & dplyr::between(BinCentre, start, stop)) %>%
    dplyr::select(Symbol, Bin, BinChr, BinCentre, GeneChr, GeneLeft, GeneRight, ReadCount, ExpectedCount, TotalCount, Padj)
  
  # identify uncovered bins and add as zero values
  
  coveredBins = unique(regionData$Bin)
  

  uncoveredBins = binAnnotation %>%
    dplyr::filter(BinChr == chr & dplyr::between(Centre, start, stop) & !(Bin %in% coveredBins)) %>%
    dplyr::rename(BinCentre = Centre) %>%
    dplyr::mutate(Symbol = goi, 
                  GeneChr = unique(regionData$GeneChr),
                  GeneLeft = unique(regionData$GeneLeft),
                  GeneRight = unique(regionData$GeneRight),
                  ReadCount = 0,
                  ExpectedCount = NA,
                  TotalCount = NA,
                  Padj = 1) %>%
    dplyr::select(Symbol, Bin, BinChr, BinCentre, GeneChr, GeneLeft, GeneRight, ReadCount, ExpectedCount, TotalCount, Padj)
  
  regionData = rbind(regionData, uncoveredBins)
  

  ggplot(regionData, mapping = aes(x = BinCentre, y = ExpectedCount)) +
    geom_line()
  
  regionData$InteractionFreq = regionData$ReadCount / regionData$TotalCount
  
  regionData$InteractionFreq_bg = regionData$ExpectedCount/ regionData$TotalCount
  
  regionData$InteractionFreq_bg = unlist(slider::slide_dbl(.x = regionData$InteractionFreq_bg, .f = mean, .before = smoothing_windows/2, .after = smoothing_windows/2))
  
  plotFreqs = reshape2::melt(regionData[,c('BinCentre', 'InteractionFreq', 'InteractionFreq_bg')], id = 'BinCentre', value.name = 'Frequency', variable.name = 'Type')
  
  vp_data = setDT(ggplot_build(
    ggplot(data = plotFreqs, mapping = aes(x = as.numeric(BinCentre), y = Frequency)) + # first plot to retrieve coordinates of smoothed frequency lines
      stat_smooth(mapping = aes(colour = Type), geom = 'line', n = 2000, span = span, method = 'loess') +
      scale_colour_manual(values = c('black', 'darkred'), breaks = c('InteractionFreq', 'InteractionFreq_bg'))
  )$data[[1]])

  
  vpRanges = GRanges(seqnames = unique(regionData$BinChr), ranges = IRanges(start = vp_data$x, end = vp_data$x))
  
  binRanges = GRanges(seqnames = unique(regionData$BinChr), ranges = IRanges(start = regionData$BinCentre-2500, end = regionData$BinCentre+2500, names = regionData$Bin))
  
  ovl = as.data.frame(findOverlaps(vpRanges, binRanges))
  
  vp_data$Bin = NA
  
  vp_data$Bin[ovl$queryHits] = names(binRanges)[ovl$subjectHits]
  
  vp_data = vp_data[,.(Max_y = max(y)), by = 'Bin']

  
  sigInts = merge(vp_data, regionData, by = 'Bin')
  
  sigInts = sigInts[sigInts$Padj < threshold, c('BinCentre', 'Max_y', 'Padj')]
  
  sigInts$Position = factor(ifelse(sigInts$BinCentre<unique(regionData$GeneLeft), yes = 'Upstream', no = 'Downstream'), levels = c('Upstream', 'Downstream'))
  
  plotFreqs$Type = factor(plotFreqs$Type, levels = c('InteractionFreq', 'InteractionFreq_bg'))
  
  plotFreqs$Distance = regionData$Distance[match(as.numeric(plotFreqs$BinCentre),regionData$BinCentre)]
  
  plotFreqs$Position = factor(ifelse(plotFreqs$BinCentre<unique(regionData$GeneLeft), yes = 'Upstream', no = 'Downstream'), levels = c('Upstream', 'Downstream'))
  
  labels_fun = function(x){
    paste(x/1000000, 'Mb')
  }
  
  plotFreqs$DataType = 'RADICL'
  
  if(is.null(ylim)){
    ylim = max(plotFreqs$Frequency, na.rm=TRUE)*1000000
  }
  
  depth_df = data.frame(
    x = max(plotFreqs$BinCentre),
    y = ylim,
    label = paste('Effective', goi, 'seq. depth:', effectiveSeqDepth),
    Position = factor('Downstream', levels = c('Upstream', 'Downstream'))
  )
  
  if(!is.null(name)){
  
    name_df = data.frame(
      x = min(plotFreqs$BinCentre),
      y = ylim,
      label = name,
      Position = factor('Upstream', levels = c('Upstream', 'Downstream'))
    )
  
  }
  
  vp = ggplot(data = plotFreqs, mapping = aes(x = BinCentre, y = Frequency*1000000)) + # final viewpoint plot
    #stat_smooth(mapping = aes(colour = Type), geom = 'line', n = 2000, span = 0.01, method = 'loess', lwd = 0.75, alpha = 0.5) +
    stat_smooth(mapping = aes(fill = Type), geom = 'area', n = 2000, span = span, method = 'loess', lwd = 0.75, alpha = 0.5) +
    scale_colour_manual(values = c('black', 'darkred'), breaks = c('InteractionFreq_bg', 'InteractionFreq'), labels = c('Background', bquote(italic(.(goi)))), aesthetics = c('fill', 'colour')) +
    geom_point(data = sigInts[sigInts$Padj<threshold,], mapping = aes(x = BinCentre, y = (Max_y*1.5)*1000000, stroke = NA, size = -log10(Padj)), pch = 21, fill = 'black', alpha = 0.5) +
    scale_size_continuous(range=c(0.5,3)) +
    coord_cartesian(ylim = c(0,ylim), expand = F) +
    scale_x_continuous(breaks = pretty, labels = labels_fun, expand = c(0,0)) +
    #scale_x_continuous(breaks = pretty, labels = labels_fun) +
    facet_grid(cols = vars(Position), scales = 'free_x', space = 'free_x') +
    #facet_wrap(~Position, scales = 'free_x') +
    xlab(unique(regionData$BinChr)) +
    ylab('Interactions per\nmillion transcripts') +
    labs(size = expression(paste('-log'[10]*'(Adj.'~italic('P')~'value)')), colour = 'RNA', fill = 'RNA') +
    geom_text(data = depth_df, mapping = aes(x = x, y = y, label = label), colour = 'black', hjust = 1, vjust = 1.05, size = 7/.pt) +
    geom_text(data = name_df, mapping = aes(x = x, y = y, label = label), colour = 'black', hjust = -0.05, vjust = 1.05, size = 7/.pt) +
    theme_top() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          legend.box = 'horizontal',
          legend.position = 'top',
          legend.direction = 'horizontal', 
          axis.line.x=element_line(), 
          axis.text.x=element_text(), 
          axis.title.x=element_text(), 
          axis.title.y.left=element_text(angle=90, margin=margin(r=1, unit="lines")))


  
  if(!is.null(additionalTrackFiles)){
    
    additionalPlots = plotAdditionalTracks(
      regionData = regionData, 
      additionalTrackFiles = additionalTrackFiles, 
      additionalTrackNames = additionalTrackNames, 
      additionalTrackTypes = additionalTrackTypes
    )
    
    vp = vp + theme(axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank())
    
    combinedPlotList = c(list(vp), additionalPlots)
    
    # combinedPlotList = lapply(combinedPlotList, function(plot){
    #   newPlot = plot + theme(plot.margin = unit(c(0,0,0,0), 'lines'),
    #                          axis.text = element_text(size = 7),
    #                          axis.title = element_text(size = 9))
    #   return(newPlot)
    # })
    
    combinedPlot = cowplot::plot_grid(plotlist = combinedPlotList, align = 'v', axis = 'lr', rel_heights = c(1,rep(0.2, length(additionalPlots))), ncol = 1)

    return(combinedPlot)
    
  } else {
    
    return(vp)
    
  }
  
  return(message('Complete'))
  
}








vp <- plotResults(
  results, 
  chr, 
  start, 
  stop, 
  goi, 
  threshold=0.05, 
  binAnnotation
)

vp <- vp + 
    theme_top() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          legend.box = 'horizontal',
          legend.position = 'top',
          legend.direction = 'horizontal', 
          axis.line.x=element_line(), 
          axis.text.x=element_text(), 
          axis.title.x=element_text(), 
          axis.title.y.left=element_text(angle=90, margin=margin(r=1, unit="lines")))





if(!dir.exists(outdir)) {dir.create(outdir, recursive = TRUE)}

my_filename=paste0(outdir, "/", goi, ".viewpoint.",  genome, ".", chr, "_", start, "_", stop)


if (stringr::str_to_lower(output_format) =="svg") {
  ggsave(plot=vp, filename=paste0(my_filename, ".svg"), height=2, width=5, device=svg)
} else if (stringr::str_to_lower(output_format)=="png") {
  ggsave(plot=vp, filename=paste0(my_filename, ".png"), height=2, width=5)
} else {
  ggsave(plot=vp, filename=paste0(my_filename, ".svg"), height=2, width=5, device=svg)
  ggsave(plot=vp, filename=paste0(my_filename, ".png"), height=2, width=5)
}





