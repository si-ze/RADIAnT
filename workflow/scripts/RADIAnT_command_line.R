library(data.table)
library(ggplot2)
library(ggarchery)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#library(EnsDb.Hsapiens.v86)
#library(EnsDb.Mmusculus.v79)
#library(org.Mm.eg.db)
#library(org.Hs.eg.db)
library(ChIPpeakAnno)
library(GenomicRanges)
library(GenomicFeatures)
library(IRanges)
#library(rcartocolor)
#library(ggrepel)
library(dplyr)
library(purrr)
library(pbapply)

#library(RColorBrewer)
#library(ggVennDiagram)
#library(ComplexUpset)
#library(viridis)

# --------------------------

# plotting functions

theme_TW = function(){
  theme_bw() +
    theme(axis.text = element_text(colour = 'black', size=12),
          axis.line = element_line(),
          panel.border = element_blank(),
          axis.ticks = element_line(colour = 'black'),
          strip.background = element_rect(colour = NA))
}

theme_heatmap = function(){
  theme_bw() +
    theme(axis.text = element_text(colour = 'black'),
          axis.line = element_blank(),
          panel.border = element_rect(colour = 'black', fill = NA),
          axis.ticks = element_blank(),
          axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))
}

colourbar_TW = function(){
  guide_colorbar(frame.colour = 'black', ticks.colour = 'black', draw.llim = F, draw.ulim = F, title.hjust = 0.5, title.vjust = 0.5, title.theme = element_text(angle = 270), title.position = 'right')
}

ggGenes = function(gtf, chr, start, stop){
  #gtf = data.table::fread(gtfFile, sep = '\t', header = F)
  conds = gtf$V3 == 'gene' & gtf$V1 == chr & gtf$V4 > start & gtf$V5 < stop
  genes = gtf[conds,]
  geneDF = data.table::rbindlist(lapply(1:nrow(genes), function(geneIndex){
    start = ifelse(genes$V7[geneIndex] == '+', yes = genes$V4[geneIndex], no = genes$V5[geneIndex])
    stop = ifelse(genes$V7[geneIndex] == '+', yes = genes$V5[geneIndex], no = genes$V4[geneIndex])
    symbol = sub('.*gene_name "', '', sub('"; level.*', '', genes$V9[geneIndex]))
    return(data.frame(start = start, stop = stop, symbol = symbol))
  }))
  geneDF = removeGeneOverlaps(geneDF)
  return(geneDF)
}

removeGeneOverlaps = function(geneDF){
  geneDF$Size = abs(geneDF$stop - geneDF$start)
  geneDF = geneDF[order(geneDF$Size, decreasing = T),]
  geneDF$y = 0
  for(geneIndex in 2:nrow(geneDF)){
    geneCoord = sort(c(geneDF$start[geneIndex], geneDF$stop[geneIndex]))
    overlaps = sapply(1:nrow(geneDF), function(overlapIndex){
      overlapCoord = sort(c(geneDF$start[overlapIndex], geneDF$stop[overlapIndex]))
      overlapCondStart = geneCoord[1] > overlapCoord[1] & geneCoord[1] < overlapCoord[2]
      overlapCondEnd = geneCoord[2] > overlapCoord[1] & geneCoord[2] < overlapCoord[2]
      overlap = overlapCondStart | overlapCondEnd
      return(overlap)
    })
    geneDF$y[geneIndex] =0-sum(overlaps)
  }
  return(geneDF)
}

# RADIAnT functions -------------------------------------------------------

# Function to find the closest value in Bin Centre for snapping gene to right/left

snapToBinCentre <- function(binAnnotationList, chr_name, gene_coord, left=TRUE) {
  
  # Extract the corresponding dataframe from binAnnotationList
  my_chr_df <- binAnnotationList[[chr_name]]
  
  if (left) {
    # for left gene ends, consider only bins left from gene    
    my_chr_df <- my_chr_df[as.numeric(unlist(my_chr_df$Centre))<=gene_coord,]
  } else { # otherwise only right bin centres
    my_chr_df <- my_chr_df[as.numeric(unlist(my_chr_df$Centre))>=gene_coord,]
  }
  
  # Find the index of the closest value in the "Bin Centre" column
  closest_index <- which.min(abs(my_chr_df$Centre - gene_coord))
  
  
  # Get the closest value in "Bin Centre"
  closest_value <- my_chr_df$Centre[closest_index]
  
  return(as.numeric(unlist(closest_value)))
}

# Snapping function altered to use dplyr::round_any, faster

snapToBinCentreTW = function(genes){
  
  geneSnapped = genes %>%
    dplyr::filter(type == 'gene') %>%
    dplyr::select('seqid', 'start', 'end', 'gene_name') %>%
    dplyr::mutate(LeftBinEnd = plyr::round_any(start, 5000, floor)) %>%
    dplyr::mutate(RightBinStart = plyr::round_any(end, 5000, ceiling)) %>%
    dplyr::mutate(GeneLeftSnapped = LeftBinEnd - 2500, GeneRightSnapped = RightBinStart + 2500) %>%
    select(seqid, start, end, gene_name, GeneLeftSnapped, GeneRightSnapped)
  
  return(geneSnapped)
  
}

# Function to read in required input data

getInputDataSZ <- function(countsFile, binAnnotationFile, gtfFile) {
  
  # message("Reading GTF file ...")
  # gtf = data.table::fread(gtfFile, sep = '\t', header = F)
  # message("Constructing taxonomic database ...")
  # myTxDB  =  makeTxDbFromGFF(file=gtfFile, format = "gtf")
  # 
  # 
  # if(species == 'human'){
  #     # ensdb = ensdb.Hsapiens.v86
  #     ensdb = org.Hs.eg.db
  # } else if(species == 'mouse'){
  #     ensdb = org.Mm.eg.db
  # }
  
  message("Reading in count data ...")
  # read in interaction counts
  counts <- data.table::fread(file = countFile) 
  counts <- counts[,c("geneSymbol", "binID", "intctCounts")] 
  colnames(counts) <- c("Symbol", "Bin", "ReadCount")
  counts <- counts %>% group_by(Symbol, Bin) %>% summarise(ReadCount = sum(ReadCount))
  # > head(counts)
  #     Symbol        Bin ReadCount
  # 1: Gm14508 chr5_22952    1
  # 2:   C2cd5 chr6_28529    3
  # 3:    Dbnl chr8_13600    1
  
  message("Reading in genome data ...")
  binAnnotation = data.table::fread(file = binAnnotationFile)
  colnames(binAnnotation) = c("BinChr", "BinStart", "BinEnd", "Bin")
  # > head(binAnnotation)
  #    BinChr BinStart BinEnd    Bin
  # 1:   chr1        0   5000 chr1_1
  # 2:   chr1     5000  10000 chr1_2
  # 3:   chr1    10000  15000 chr1_3
  
  # remove bins with rare chrom annotations like "chr1_MU069434v1_random"
  binAnnotation = binAnnotation[-grep('_', binAnnotation$BinChr)]
  
  # compute bin centre ...
  binAnnotation$Centre = rowMeans(binAnnotation[,2:3])
  # and add to interaction count info
  counts[,c('BinChr', 'BinCentre')] = binAnnotation[match(counts$Bin, binAnnotation$Bin), c('BinChr', 'Centre')]
  # > head(counts)
  #     Symbol        Bin ReadCount BinChr BinCentre
  # 1: Gm14508 chr5_22952    1   chr5 114757500
  # 2:   C2cd5 chr6_28529    3   chr6 142642500
  # 3:    Dbnl chr8_13600    1   chr8  67997500
  
  message("Getting info on genes ...")
  
  # read in gene info
  message("Reading GTF...")
  
  genes = rtracklayer::readGFF(gtfFile)
  
  message('Snapping gene coordinates...')
  
  geneData = snapToBinCentreTW(genes)
  
  colnames(geneData) = c('seqnames', 'start', 'end', 'Symbol', 'GeneLeftSnapped', 'GeneRightSnapped')
  
  geneData$Symbol = make.unique(geneData$Symbol)
  
  # geneRanges = genes(myTxDB) # note: this function automatically assigns gene start as lesser than gene end
  # names(geneRanges) = sub('\\..*', '', names(geneRanges))
  # entrez2symbol = AnnotationDbi::select(x = ensdb, keys = names(geneRanges), columns = c('ENTREZID', 'SYMBOL'), keytype = 'ENSEMBL') # delete ENTREZID
  # geneRanges$Symbol = entrez2symbol$SYMBOL[match(names(geneRanges), entrez2symbol$ENSEMBL)]
  # geneRanges$Symbol[is.na(geneRanges$Symbol)] = names(geneRanges)[is.na(geneRanges$Symbol)]
  # geneData = as.data.frame(geneRanges, row.names=NULL)
  # #geneData$GeneStart = as.numeric(ifelse(geneData$strand == '-', yes = geneData$end, no = geneData$start))
  # #geneData$GeneEnd = as.numeric(ifelse(geneData$strand == '-', yes = geneData$start, no = geneData$end))
  # 
  # 
  # message("Correcting gene positions ...")
  # # create seperate bin annotation tables per chromosome
  # binAnnotationList = split(binAnnotation, binAnnotation$BinChr)
  # 
  # geneData$GeneCentre = as.numeric(rowMeans(geneData[,c('start', 'end')]))
  # # snap left gene coordinates to next bin centre to the left
  # geneData$GeneLeftSnapped <- apply(geneData, 1, function(gene) {
  # snapToBinCentre(binAnnotationList, as.character(gene["seqnames"]), as.numeric(gene["start"]), left=TRUE)
  # })
  # geneData$GeneLeftSnapped <- as.numeric(geneData$GeneLeftSnapped)
  # 
  # # snap right gene coordinates to next bin centre to the right
  # geneData$GeneRightSnapped <- apply(geneData, 1, function(gene) {
  # snapToBinCentre(binAnnotationList, as.character(gene["seqnames"]), as.numeric(gene["end"]), left=FALSE)
  # })
  # geneData$GeneRightSnapped <- as.numeric(geneData$GeneRightSnapped)
  # > head(geneData)
  #                    seqnames     start       end  width strand  gene_id                Symbol GeneCentre GeneLeftSnapped GeneRightSnapped
  # ENSMUSG00000000001     chr3 108014596 108053462  38867      -  ENSMUSG00000000001.5   Gnai3  108034029       108012500  108057500
  # ENSMUSG00000000003     chrX  76881507  76897229  15723      -  ENSMUSG00000000003.16  Pbsn   76889368        76877500   76897500
  # ENSMUSG00000000028    chr16  18599197  18630737  31541      -  ENSMUSG00000000028.16  Cdc45   18614967        18597500   18632500
  
  
  inputData <- list(
    "gtf"=genes, 
    #"myTxDB"=myTxDB,
    "counts"=counts,
    "binAnnotation"=binAnnotation,
    #"geneRanges"=geneRanges,
    "geneData"=geneData
    #"ensdb"=ensdb
  )
  
  gc()
  
  return(inputData)
}

# Function to read in required input data updated

getInputDataTW <- function(countsFile, binAnnotationFile, gtfFile) {
  
  # Read in interaction counts
  
  message("Reading in count data ...")
  
  counts <- data.table::fread(file = countFile) 
  
  #counts <- counts[,c("geneSymbol", "binID", "intctCounts")] 
  
  colnames(counts) <- c("Symbol", "Bin", "ReadCount")
  
  counts <- counts %>% 
    dplyr::group_by(Symbol, Bin) %>% 
    dplyr::summarise(ReadCount = sum(ReadCount)) %>%
    dplyr::ungroup() # handles multiple loci/samples if files are concatenated
  
  # Read in bin annotations
  
  message("Reading in genome data ...")
  
  binAnnotation = data.table::fread(file = binAnnotationFile)
  
  colnames(binAnnotation) = c("BinChr", "BinStart", "BinEnd", "Bin")
  
  #binAnnotation = binAnnotation[-grepl('_', binAnnotation$BinChr)]
  
  binAnnotation$Centre = rowMeans(binAnnotation[,2:3]) # centres for snapping
  
  counts[,c('BinChr', 'BinCentre')] = binAnnotation[match(counts$Bin, binAnnotation$Bin), c('BinChr', 'Centre')]
  
  # Read in gene data
  
  message("Getting info on genes ...")
  
  message("Reading GTF...")
  
  genes = rtracklayer::readGFF(gtfFile) # reads in a better format than data.table::fread, but slower
  
  genes = genes %>%
    filter(type == 'gene')
  
  genes$gene_name = make.unique(genes$gene_name)
  
  message('Snapping gene coordinates...')
  
  geneData = snapToBinCentreTW(genes) # snaps using dplyr::round_any, faster than previous method
  
  colnames(geneData) = c('seqnames', 'start', 'end', 'Symbol', 'GeneLeftSnapped', 'GeneRightSnapped')
  
  # Convert to gene symbol if needed ----------------------------------------
  
  if(startsWith(counts$Symbol[1], 'ENS')){
    
    geneid2symbol = genes %>%
      filter(type == 'gene') %>%
      dplyr::select(gene_id, gene_name)
    
    counts$Symbol = geneid2symbol$gene_name[match(counts$Symbol, geneid2symbol$gene_id)]
    
  }
  
  inputData <- list(
    "gtf"=genes, 
    "counts"=counts,
    "binAnnotation"=binAnnotation,
    "geneData"=geneData
  )
  
  gc()
  
  return(inputData)
}

geneArrows = function(geneDF, distance){
  geneDF$width = abs(geneDF$stop - geneDF$start)
  arrowDF = data.table::rbindlist(lapply(1:nrow(geneDF), function(geneIndex){
    sense = geneDF$start[geneIndex] < geneDF$stop[geneIndex]
    if(sense){
      arrowStart = seq(from = geneDF$stop[geneIndex], to = geneDF$start[geneIndex], by = -distance)
      arrowDF = data.frame(start = arrowStart, end = arrowStart+1)
    } else {
      arrowEnd = seq(from = geneDF$stop[geneIndex], to = geneDF$start[geneIndex], by = distance)
      arrowDF = data.frame(start = arrowEnd+1, end = arrowEnd)
    }
    arrowDF$y = geneDF$y[geneIndex]
    return(arrowDF)
  }))
  return(arrowDF)
}

getCisResults <- function(myCisCounts, cisDistanceCounts, binAnnotation){
  message("Computing results for cis interactions ...")
  
  myCisCounts$InteractionFreq_bg <- cisDistanceCounts$InteractionFreq[match(myCisCounts$Distance, cisDistanceCounts$Distance)]
  
  myCisCounts$ExpectedCounts <- myCisCounts$InteractionFreq_bg*myCisCounts$TotalCisCount
  
  threshCisCounts <- myCisCounts[myCisCounts$ReadCount >= count_thresh & myCisCounts$ReadCount > myCisCounts$ExpectedCounts,]
  
  threshCisCounts$P <- pbsapply(1:nrow(threshCisCounts), function(row){
    poisson.test(x = threshCisCounts$ReadCount[row],
                 r = threshCisCounts$ExpectedCounts[row], 
                 alternative = 'greater')$p.value
  })
  
  threshCisCounts$Padj_global <- p.adjust(p = threshCisCounts$P, 
                                          method = 'fdr', 
                                          n = nrow(threshCisCounts))
  
  tmpCisRes <- merge(myCisCounts, threshCisCounts[,c('Symbol', 'Distance', 'P', 'Padj_global')], by = c('Symbol', 'Distance'), all = T)
  
  tmpCisRes[,c("BinStart", "BinEnd")] <- binAnnotation[match(tmpCisRes$Bin, binAnnotation$Bin), c("BinStart", "BinEnd")]
  
  output_table_cis <- tmpCisRes %>% select("Symbol", "Bin", "P", "Padj_global", "ExpectedCounts", "ReadCount", "GeneChr", "GeneLeft", "GeneRight", "BinChr","BinStart", "BinCentre", "BinEnd", "Distance")
  
  rm(tmpCisRes)
  
  gc()
  
  output_table_cis %>% group_by(Symbol) %>% mutate(Padj_gene = p.adjust(P, method = 'fdr')) %>% as.data.frame() %>% fwrite(file=paste0(outDir, "/cisResults.csv"), sep=",")
  
  return(output_table_cis)
  
}

getTransResults <- function(transCounts, transBinMeans, binAnnotation){
  message("Computing results for trans interactions ...")
  myTransCounts <- merge(transBinMeans, transCounts[,c('Bin', 'Symbol', 'GeneLeft', 'GeneRight', 'GeneChr', 'ReadCount')], by = 'Bin')
  
  myTransCounts$ExpectedCounts <- myTransCounts$IntFreqByBin * transCounts$TotalTransCount 
  
  threshTransCounts <- myTransCounts[myTransCounts$ReadCount >= count_thresh & myTransCounts$ReadCount > myTransCounts$ExpectedCounts,]
  
  threshTransCounts$P <- pbsapply(1:nrow(threshTransCounts), function(row){
    poisson.test(x = threshTransCounts$ReadCount[row],
                 r = threshTransCounts$ExpectedCounts[row], 
                 alternative = 'greater')$p.value
  })
  
  threshTransCounts$Padj_global <- p.adjust(p = threshTransCounts$P, 
                                            method = 'fdr', 
                                            n = nrow(threshTransCounts))
  
  tmpTransRes <- merge(myTransCounts, threshTransCounts[,c('Symbol', 'P', 'Padj_global', 'Bin')], by = c('Symbol', 'Bin'), all = T)
  
  tmpTransRes[,c("BinStart", "BinCentre", "BinEnd", "BinChr")] <- binAnnotation[match(tmpTransRes$Bin, binAnnotation$Bin), c("BinStart", "Centre", "BinEnd", "BinChr")]
  
  output_table_trans <- tmpTransRes %>% select("Symbol", "Bin", "P", "Padj_global", "ExpectedCounts", "ReadCount", "GeneChr", "GeneLeft", "GeneRight", "BinChr","BinStart", "BinCentre", "BinEnd")
  
  output_table_trans$Distance <- rep(NA, nrow(output_table_trans))
  
  rm(tmpTransRes)
  
  gc()
  
  output_table_trans %>% group_by(Symbol) %>% mutate(Padj_gene = p.adjust(P, method = 'fdr')) %>% as.data.frame() %>% fwrite(file=paste0(outDir, "/transResults.csv"), sep=",")
  
  return(output_table_trans)
  # NAs occur in P & Padj when threshold isn't met 
}

getCisCounts <- function(bgCounts, binAnnotation, geneData) {
  # subset to only cis frequencies
  cisCounts = bgCounts[bgCounts$BinChr==bgCounts$GeneChr,] 
  # compute distance
  cisCounts$Distance = ifelse(cisCounts$BinCentre < cisCounts$GeneLeft, # if bin centre outside gene body (to the left)
                              yes = cisCounts$BinCentre - cisCounts$GeneLeft, # compute distance between bins and gene left/right
                              no =  ifelse(cisCounts$BinCentre > cisCounts$GeneRight, # if bin centre outside gene body (to the right)
                                           yes = cisCounts$BinCentre - cisCounts$GeneRight, # compute distance between bins and gene left/right
                                           no = NA)) # not interested bin centres falling within gene body
  
  # # ensure frequencies lie within window L
  # cisCounts = cisCounts[abs(cisCounts$Distance) <= myL,] 
  
  
  
  
  cisCounts = na.omit(cisCounts) # remove any NA values (likely gene body)
  # > head(cisCounts)
  #     Symbol        Bin ReadCount BinChr BinCentre GeneChr  GeneLeft GeneRight Distance
  # 1:   C2cd5 chr6_28529    3   chr6 142642500    chr6 142952500 143047500 -310000
  # 2: Macrod1 chr19_2102    1  chr19  10507500   chr19   7032500   7177500 3330000
  # 3:  Lrrc8d chr5_21159   11   chr5 105792500    chr5 105847500 105982500 -55000
  
  cisCounts <- as.data.table(cisCounts)
  
  
  # list all distances which are covered for each RNA
  coveredDistances = paste0(cisCounts$Symbol, '::', cisCounts$Distance) 
  # > head(coveredDistances)
  # [1] "C2cd5::-310000"   "Macrod1::3330000" "Lrrc8d::-55000"   "Agap1::255000"   
  # [5] "Lrrc8d::735000"   "Pex14::-685000"  
  
  # list all possible distances for every RNA
  possibleDistances = unlist(pblapply(unique(cisCounts$Symbol), function(rna){ 
    paste0(rna, '::', seq(-5000000, 5000000, by = 5000))
  }))
  
  # subset to only those which have no frequency detected (i.e. zero bins)
  uncoveredDistances = possibleDistances[!(possibleDistances %in% coveredDistances)] 
  
  
  # create data.frame which contains all zero counts
  zeroCounts = data.frame( 
    Symbol = sub('::.*', '', uncoveredDistances),
    Distance = as.numeric(sub('.*::', '', uncoveredDistances))
  )
  # > head(zeroCounts)
  #   Symbol Distance ReadCount
  # 1  C2cd5 -5000000    0
  # 2  C2cd5 -4995000    0
  # 3  C2cd5 -4990000    0
  
  # zeroCounts$GeneLeftSnapped <- inputData$geneData$GeneLeftSnapped[match(zeroCounts$Symbol, inputData$geneData$Symbol)]
  # zeroCounts$GeneRightSnapped <- inputData$geneData$GeneRightSnapped[match(zeroCounts$Symbol, inputData$geneData$Symbol)]
  # zeroCounts$GeneChr <- as.character(inputData$geneData$seqnames[match(zeroCounts$Symbol, inputData$geneData$Symbol)])
  
  zeroCounts[,c('GeneLeftSnapped', 'GeneRightSnapped', 'GeneChr')] = geneData[match(zeroCounts$Symbol, inputData$geneData$Symbol), c('GeneLeftSnapped', 'GeneRightSnapped', 'seqnames')]
  
  zeroCounts$BinCentre <- ifelse(zeroCounts$Distance<0, 
                                 yes=zeroCounts$GeneLeftSnapped+zeroCounts$Distance,
                                 no=zeroCounts$GeneRightSnapped+zeroCounts$Distance)
  
  chrMaxBins = binAnnotation[,.(ChrMax = max(BinEnd)), by = 'BinChr']
  
  zeroCounts$ChrMax = chrMaxBins$ChrMax[match(zeroCounts$GeneChr, chrMaxBins$BinChr)]
  
  zeroCounts = zeroCounts[zeroCounts$BinCentre >= 0 & zeroCounts$BinCentre <= zeroCounts$ChrMax,]
  
  zeroCounts$Bin = binAnnotation$Bin[match(paste0(zeroCounts$GeneChr,'::',zeroCounts$BinCentre), paste0(binAnnotation$BinChr, '::', binAnnotation$Centre))]
  
  #zeroCounts$Bin <- binAnnotation$Bin[binAnnotation$BinChr==zeroCounts$GeneChr & binAnnotation$Centre==zeroCounts$BinCentre]
  
  zeroCounts$ReadCount <- 0
  
  zeroCounts[,c("GeneLeftSnapped", "GeneRightSnapped", "GeneChr", "ChrMax", "BinCentre")] <- NULL
  
  # combine non-zero and zero interaction counts
  # allCisCounts = rbind(zeroCounts, cisCounts[,c("Symbol","Distance","Bin","ReadCount")]) 
  allCisCounts = rbind(zeroCounts, cisCounts[,.(Symbol,Distance,Bin,ReadCount)]) 
  
  
  # compute total counts within L per RNA
  rnaTotals = cisCounts[,.(TotalCisCount = sum(ReadCount)), by = Symbol] 
  # > head(rnaTotals)
  #     Symbol TotalCisCount
  # 1:   C2cd5       622
  # 2: Macrod1       479
  # 3:  Lrrc8d      1591
  
  # add total counts to interaction counts
  allCisCounts$TotalCisCount = rnaTotals$TotalCisCount[match(allCisCounts$Symbol, rnaTotals$Symbol)] 
  # compute normalised interaction frequencies
  allCisCounts$FreqNormByCisCount = allCisCounts$ReadCount/allCisCounts$TotalCisCount 
  allCisCounts$FreqNormByCisCount[is.na(allCisCounts$FreqNormByCisCount)] = 0 # replace any NA values with zero
  # > tail(allCisCounts)
  #            Symbol Distance ReadCount TotalCisCount FreqNormByCisCount
  # 39123457   Nfkbid -1600000    1        33     0.030303030
  # 39123458    Kif3a -3180000    1        41     0.024390244
  # 39123459    Barx2   150000    1         6     0.166666667
  
  
  # calculate mean frequency per distance
  cisDistanceCounts = setDT(allCisCounts)[,.(FreqNormByCisCount = mean(FreqNormByCisCount)), by = Distance] 
  cisDistanceCounts$Type = 'Background' # add type column to background frequencies
  # > head(cisDistanceCounts)
  #    Distance FreqNormByCisCount       Type
  # 1: -5000000    6.416667e-05 Background
  # 2: -4995000    6.975770e-05 Background
  # 3: -4990000    7.306042e-05 Background
  
  myCisCounts <- allCisCounts 
  
  # merge later 
  # myCisCounts[,c("GeneChr", "GeneLeft", "GeneRight")] <- cisCounts[match(myCisCounts$Symbol, cisCounts$Symbol), c("GeneChr", "GeneLeft", "GeneRight")]
  # myCisCounts$BinCentre <- ifelse(myCisCounts$Distance<0, (myCisCounts$GeneLeft+myCisCounts$Distance), (myCisCounts$GeneRight+myCisCounts$Distance))
  # myCisCounts$BinChr <- myCisCounts$GeneChr
  # tmpBinAnnotation <- binAnnotation
  # names(tmpBinAnnotation)[5] <- "BinCentre"
  # myCisCounts <- merge(myCisCounts, tmpBinAnnotation, by=c("BinChr", "BinCentre"))
  
  gc()
  
  return(list(
    "cisDistanceCounts"=cisDistanceCounts,
    "cisCounts"=myCisCounts
  ))
}

getCisCountsTW <- function(counts, binAnnotation, geneData, window) {
  
  # Subset counts to only those in window minus the gene body
  
  cisCounts = counts %>%
    dplyr::filter(GeneChr == BinChr) %>% # filter for interactions on the same chromosome
    dplyr::filter(dplyr::between(BinCentre, GeneLeft-window, GeneRight+window)) %>% # filter for interactions within the specified window
    dplyr::filter(!dplyr::between(BinCentre, GeneLeft, GeneRight)) # filter to remove counts within the gene body
  
  # Compute bin - gene distances (left & right)
  
  cisCounts$Distance = ifelse(cisCounts$BinCentre < cisCounts$GeneLeft, yes = cisCounts$BinCentre - cisCounts$GeneLeft, no = cisCounts$BinCentre - cisCounts$GeneRight)
  
  cisCounts = data.table::setDT(cisCounts)
  
  # Total read counts per distance in window
  
  rnasByDistance = cisCounts[,.(SumCount = sum(ReadCount), CoveredRNAs = length(unique(Symbol))), by = 'Distance']
  
  # Total read counts within window
  
  rnasByDistance$OverallRNACount = sum(counts$ReadCount)
  
  # Interaction frequency per distance
  
  rnasByDistance$InteractionFrequency = rnasByDistance$SumCount/rnasByDistance$OverallRNACount
  
  # Add to window counts
  
  cisCounts$ExpInteractionFrequencyWindow = rnasByDistance$InteractionFrequency[match(cisCounts$Distance, rnasByDistance$Distance)]
  
  # Compute total window counts per RNA
  
  rnaWindowCounts = data.table::setDT(cisCounts)[,.(WindowReadCount = sum(ReadCount)), by = 'Symbol']
  
  cisCounts$WindowReadCount = rnaWindowCounts$WindowReadCount[match(cisCounts$Symbol, rnaWindowCounts$Symbol)]
  
  cisCounts$ExpInteractionCountWindow = cisCounts$WindowReadCount * cisCounts$ExpInteractionFrequencyWindow
  
  # Return counts and generalised background
  
  return(list(
    "cisDistanceCounts"=rnasByDistance,
    "cisCounts"=cisCounts
  ))
}


getTransCounts <- function(bgCounts) {
  
  # define trans conditions
  trans_conds = bgCounts$BinChr!=bgCounts$GeneChr | (bgCounts$BinCentre<(bgCounts$GeneLeft-myL)) | (bgCounts$BinCentre>(bgCounts$GeneRight+myL))
  # different chromosomes AND interactions on same chromosome outside L
  transCounts = bgCounts[trans_conds,]
  transCounts <- as.data.table(transCounts)
  #     GeneID        Bin Count BinChr BinCentre GeneChr  GeneLeft GeneRight
  # 2:   C2cd5 chr6_28529    3   chr6 142642500    chr6 142952500 143047500
  # 3:    Dbnl chr8_13600    1   chr8  67997500   chr11   5737500   5752500
  # 4: Map3k20 chr2_14437    7   chr2  72182500    chr2  72112500  72277500
  
  # get total number of RNAs
  totalTransRNAs = length(unique(transCounts$Symbol))
  # > totalTransRNAs
  # [1] 30338
  
  # get number of genes interacting per bin
  rnasPerBin = transCounts[,.(TotalTransRNAsPerBin = length(unique(Symbol)), TotalZerosPerBin = totalTransRNAs-length(unique(Symbol))), by = 'Bin']
  # > head(rnasPerBin)
  #                 Bin TotalTransRNAsPerBin TotalZerosPerBin
  #      1:  chr8_13600         23       30315
  #      2:  chr18_7965         20       30318
  #      3: chr14_11104         34       30304
  
  # add total zeros to trans interaction counts
  transCounts$TotalZerosPerBin = rnasPerBin$TotalZerosPerBin[match(transCounts$Bin, rnasPerBin$Bin)]
  
  transRNAtotals = transCounts[,.(TotalTransReadsPerRNA = sum(ReadCount)), by = 'Symbol']
  # > head(transRNAtotals)
  #      Symbol TotalTransReadsPerRNA
  # 1:     Dbnl       405
  # 2:     Rtn4      3648
  # 3:    Ikbkb      1281
  
  
  transCounts$TotalTransCount = transRNAtotals$TotalTransReadsPerRNA[match(transCounts$Symbol, transRNAtotals$Symbol)]
  
  transCounts$FreqNormInTrans = transCounts$ReadCount/transCounts$TotalTransCount
  
  transCounts  <- as.data.table(transCounts)
  
  transBinMeans = transCounts[,.(IntFreqByBin = mean(c(FreqNormInTrans, rep(0,unique(TotalZerosPerBin))))), by = 'Bin']
  # > head(transBinMeans)
  #            Bin     IntFreqByBin
  # 1:  chr8_13600 4.197029e-06
  # 2:  chr18_7965 7.059736e-07
  # 3: chr14_11104 7.844734e-06
  
  
  
  return(list(
    "transCounts"=transCounts,
    "transBinMeans"=transBinMeans
  ))
}

getTransCountsTW <- function(counts, cisWindow, geneData) {
  
  counts[,c('GeneChr', 'GeneLeft', 'GeneRight')] = geneData[match(counts$Symbol, geneData$Symbol),c('seqnames', 'GeneLeftSnapped', 'GeneRightSnapped')] # add gene information to frequencies
  
  transCounts = counts %>%
    dplyr::filter(BinChr != GeneChr | !dplyr::between(BinCentre, GeneLeft-cisWindow, GeneRight+cisWindow)) %>%
    data.table::setDT()
  
  countsByBin = transCounts[,.(BinCounts = sum(ReadCount)), by = Bin]
  
  countsByBin$TotalCount = sum(counts$ReadCount)
  
  countsByBin$InteractionFrequency = countsByBin$BinCounts/countsByBin$TotalCount
  
  transCounts$ExpectedInteractionFrequencyBin = countsByBin$InteractionFrequency[match(transCounts$Bin, countsByBin$Bin)]
  
  nonWindowCounts = transCounts[,.(NonWindowCount = sum(ReadCount)), by = 'Symbol']
  
  transCounts$NonWindowCounts = nonWindowCounts$NonWindowCount[match(transCounts$Symbol, nonWindowCounts$Symbol)]
  
  transCounts$ExpectedInteractionCountBin = transCounts$NonWindowCounts*transCounts$ExpectedInteractionFrequencyBin
  
  return(list(
    "transCounts"=transCounts,
    "transBinMeans"=countsByBin
  ))
}

combinedDataTW = function(cisCounts, transCounts, cisBackground, transBackground){
  
  cisMetadata = cisCounts %>% 
    dplyr::select(Symbol, Bin, BinChr, BinCentre, GeneChr, GeneLeft, GeneRight, ReadCount, WindowReadCount, ExpInteractionFrequencyWindow, ExpInteractionCountWindow) %>%
    dplyr::mutate(JoinColumn = paste(Symbol, Bin, sep = '::')) %>%
    dplyr::mutate() %>%
    dplyr::relocate(JoinColumn, .before = Symbol)
  
  transMetadata = transCounts %>%
    dplyr::mutate(JoinColumn = paste(Symbol, Bin, sep = '::')) %>%
    dplyr::relocate(JoinColumn, .before = Symbol)
  
  cisMetadata$ExpectedInteractionFrequencyBin = transBackground$InteractionFrequency[match(cisMetadata$Bin, transBackground$Bin)]
  
  cisMetadata$ExpectedInteractionFrequencyBin[is.na(cisMetadata$ExpectedInteractionFrequencyBin)] = 0
  
  cisMetadata$ExpectedInteractionCountBin = cisMetadata$WindowReadCount*cisMetadata$ExpectedInteractionFrequencyBin
  
  #cisMetadata$ExpectedInteractionCountBin = transMetadata$ExpectedInteractionCountBin[match(cisMetadata$JoinColumn, transMetadata$JoinColumn)]
  
  cisMetadata$MaxExpectedCount = do.call(pmax, cisMetadata[,c('ExpInteractionCountWindow', 'ExpectedInteractionCountBin')])
  
  cisMetadata = cisMetadata %>%
    dplyr::mutate(Method = ifelse(ExpInteractionCountWindow > ExpectedInteractionCountBin, yes = 'Distance', no = 'Bin')) %>%
    dplyr::select(JoinColumn, Symbol, Bin, BinChr, BinCentre, GeneChr, GeneLeft, GeneRight, ReadCount, MaxExpectedCount, Method) %>%
    dplyr::rename(ExpectedCount = MaxExpectedCount) %>%
    dplyr::mutate(OriginalMethod = 'Cis')
  
  transMetadata = transCounts %>% 
    dplyr::select(Symbol, Bin, BinChr, BinCentre, GeneChr, GeneLeft, GeneRight, ReadCount, ExpectedInteractionCountBin) %>%
    dplyr::mutate(JoinColumn = paste(Symbol, Bin, sep = '::'), Method = 'Bin') %>%
    dplyr::relocate(JoinColumn, .before = Symbol) %>%
    dplyr::rename(ExpectedCount = ExpectedInteractionCountBin) %>%
    dplyr::mutate(OriginalMethod = 'Trans')
  
  combinedData = rbind(cisMetadata, transMetadata)
  
  return(combinedData)
  
}


getCountsAndFreqs  <- function(counts, binAnnotation, geneData) {
  # get each gene ID and respective counts for background building 
  message("Getting background counts ...")
  bgCounts = counts[counts$Symbol %in% unique(counts$Symbol),] 
  bgCounts[,c('GeneChr', 'GeneLeft', 'GeneRight')] = geneData[match(bgCounts$Symbol, geneData$Symbol),c('seqnames', 'GeneLeftSnapped', 'GeneRightSnapped')] # add gene information to frequencies
  # > head(bgCounts)
  #     Symbol        Bin ReadCount BinChr BinCentre GeneChr  GeneLeft GeneRight
  # 2:   C2cd5 chr6_28529    3   chr6 142642500    chr6 142952500 143047500
  # 3:    Dbnl chr8_13600    1   chr8  67997500   chr11   5737500   5752500
  # 4: Map3k20 chr2_14437    7   chr2  72182500    chr2  72112500  72277500
  message("Getting count/frequency data in cis ...")
  cisData <- getCisCounts(bgCounts = bgCounts, binAnnotation = binAnnotation, geneData = geneData)
  
  
  message("Getting count/frequency data in trans ...")
  transData <- getTransCounts(bgCounts)
  
  return(list(
    "cisDistanceCounts"=cisData$cisDistanceCounts,
    "cisCounts"=cisData$cisCounts,
    "transCounts"=transData$transCounts,
    "transBinMeans"=transData$transBinMeans
  ))
}

getCountsAndFreqsTW  <- function(counts, binAnnotation, geneData, cisWindow) {
  
  # get each gene ID and respective counts for background building 
  
  message("Getting background counts ...")
  
  counts[,c('GeneChr', 'GeneLeft', 'GeneRight')] = geneData[match(counts$Symbol, geneData$Symbol),c('seqnames', 'GeneLeftSnapped', 'GeneRightSnapped')] # add gene information to frequencies
  
  message("Getting count/frequency data in cis ...")
  
  cisData <- getCisCountsTW(counts = counts, binAnnotation = binAnnotation, geneData = geneData, window = cisWindow)
  
  message("Getting count/frequency data in trans ...")
  
  transData <- getTransCountsTW(counts = counts, geneData = geneData, cisWindow = cisWindow)
  
  return(list(
    "cisDistanceCounts"=cisData$cisDistanceCounts,
    "cisCounts"=cisData$cisCounts,
    "transCounts"=transData$transCounts,
    "transBinMeans"=transData$transBinMeans
  ))
}


plot_bed_track = function(bed_file, region_range, name){
  
  bed_data = data.table::fread(bed_file) %>%
    select(1:4) %>%
    rename(Chr = V1, Start = V2, Stop = V3, Name = V4) %>%
    makeGRangesFromDataFrame() %>%
    filter_by_overlaps(region_range) %>%
    mutate(Position = factor(ifelse(start < unique(goiData$GeneLeft), yes = 'Upstream', no = 'Downstream'), levels = c('Upstream', 'Downstream'))) %>%
    as.data.table()
  
  zero_data = data.table(
    seqnames = unique(bed_data$seqnames),
    start = c(start(region_range), unique(goiData$GeneLeft)-1, unique(goiData$GeneRight), end(region_range)-1),
    end = c(start(region_range)+1, unique(goiData$GeneLeft), unique(goiData$GeneRight)+1, end(region_range)),
    width = 1,
    strand = '*',
    Position = c('Upstream', 'Downstream')
  )
  
  plot_data = rbind(bed_data, zero_data)
  
  bed_plot = ggplot(data = plot_data) +
    geom_rect(mapping = aes(xmin = start, xmax = end, ymin = 0, ymax = 1), fill = 'black', colour = 'black', size = 0.1) +
    scale_x_continuous(expand = c(0,0)) +
    theme_TW() +
    facet_wrap(~Position, scales = 'free_x') +
    ylab(name) +
    theme(strip.background = element_blank(), 
          strip.text = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y.left = element_text(angle = 0))
  
  return(bed_plot)
  
}

# Plotting function -------------------------------------------------------

plotResults = function(results, cisCounts, transCounts, chr, start, stop, goi, threshold=0.05, binAnnotation, additionalTrackFile = NULL, outFile = 'viewpoint.svg', smoothing_windows = 0, span = 0.01, ylim = NULL){
  
  # subset results to provided region
  
  goiData = setDT(results[results$Symbol==goi,])
  
  message(paste('Effective', goi, 'sequencing depth is:', sum(goiData$ReadCount)))
  
  effectiveSeqDepth = sum(goiData$ReadCount)
  
  goiData$Method
  
  goiCisCount = data.table::setDT(cisCounts)[cisCounts$Symbol==goi,.(CisTotal = sum(ReadCount)), by = 'Symbol']
  
  goiTransCount = data.table::setDT(transCounts)[transCounts$Symbol==goi,.(TransTotal = sum(ReadCount)), by = 'Symbol']
  
  goiData$TotalCount = ifelse(goiData$Method == 'Bin', yes = max(goiTransCount$TransTotal), no = max(goiCisCount$CisTotal))
  
  regionData = goiData %>%
    dplyr::filter(BinChr == chr & dplyr::between(BinCentre, start, stop)) %>%
    select(Symbol, Bin, BinChr, BinCentre, GeneChr, GeneLeft, GeneRight, ReadCount, ExpectedCount, TotalCount, Padj)
  
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
    ylim = max(plotFreqs$Frequency)*1000000
  }
  
  depth_df = data.frame(
    x = max(plotFreqs$BinCentre),
    y = ylim,
    label = paste('Effective', goi, 'seq. depth:', effectiveSeqDepth),
    Position = factor('Downstream', levels = c('Upstream', 'Downstream'))
  )
  
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
    theme_TW() +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          legend.box = 'horizontal',
          legend.position = 'top',
          legend.direction = 'horizontal')
  
  if(!is.null(additionalTrackFile)){
    
    trackData = data.table::fread(additionalTrackFile, sep = '\t', header = F)
    
    subTrackData = trackData[trackData$V1 == unique(regionData$GeneChr) & (between(as.numeric(trackData$V2), min(regionData$BinStart),  unique(regionData$GeneLeft)) | between(as.numeric(trackData$V3), unique(regionData$GeneRight), max(regionData$BinEnd))),]
    
    subTrackData$MeanX = rowMeans(subTrackData[,c('V2', 'V3')])
    
    missingDataMin = data.frame(V1 = unique(subTrackData$V1),
                                V2 = min(regionData$BinCentre),
                                V3 = min(regionData$BinCentre),
                                V4 = 0,
                                MeanX = min(regionData$BinCentre))
    
    missingDataMax = data.frame(V1 = unique(subTrackData$V1),
                                V2 = max(regionData$BinCentre),
                                V3 = max(regionData$BinCentre),
                                V4 = 0,
                                MeanX = max(regionData$BinCentre))
    
    subTrackData = rbindlist(list(subTrackData, missingDataMin, missingDataMax))
    
    subTrackData$Position = factor(ifelse(subTrackData$MeanX < unique(regionData$GeneLeft), yes = 'Upstream', no = 'Downstream'), levels = c('Upstream', 'Downstream'))
    
    subTrackData$DataType = 'CAGE'
    
    cage_plot = ggplot(subTrackData, aes(x = MeanX, y = V4)) +
      #geom_bar(stat='identity', colour = 'black') +
      #scale_y_log10() +
      geom_line(stat = 'identity', linewidth = 0.1) +
      facet_wrap(~Position, scales = 'free_x') +
      scale_x_continuous(breaks = pretty, labels = labels_fun, expand = c(0,0)) +
      scale_y_continuous(limits = c(0,2), oob = scales::squish) +
      theme_TW() +
      ylab('log CAGE\nexpression') +
      xlab(unique(regionData$GeneChr)) +
      theme(strip.background = element_blank(),
            strip.text = element_blank(),
            legend.box = 'horizontal')
    
    vp = vp + theme(axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank())
    
    combinedPlot = cowplot::plot_grid(vp, cage_plot, align = 'v', axis = 'lr', rel_heights = c(1,0.5), ncol = 1)
    
    ggsave(plot = combinedPlot, filename = outFile, width = 8, height = 2.5)
    ggsave(plot = combinedPlot, filename = outFile, width = 8, height = 2.5, units = 'in', dpi = 600)
    
    return(vp)
    
  } else {
    
    #ggsave(plot = vp, filename = outFile, width = 8, height = 2)
    
    return(vp)
    
  }
  
  return(message('Complete'))
  
}

plotGenomicFeatures <- function(sigTable, extension) {
  goiRanges = GRanges(seqnames = sigTable$BinChr, 
                      ranges = IRanges(start = sigTable$BinStart, end = sigTable$BinEnd))
  
  # # obsolete precedence: 
  # precedence= c('Promoters', 'Introns', 'Exons', 'immediateDownstream', 'fiveUTRs', 'threeUTRs')
  # goiLocations1= as.data.frame(assignChromosomeRegion(peaks.RD = goiRanges, nucleotideLevel = F, TxDb = inputData$myTxDB, precedence = precedence))
  # goiLocations1$feature = rownames(goiLocations)
  
  # Marcel's suggestion:
  precedence =c('Promoters', 'fiveUTRs', 'threeUTRs', 'Introns', 'Exons', 'immediateDownstream' )
  goiLocations= as.data.frame(assignChromosomeRegion(peaks.RD = goiRanges, nucleotideLevel = F, TxDb = inputData$myTxDB, precedence = precedence))
  goiLocations <- goiLocations %>% rename(percentage = percentage.Freq, feature = percentage.subjectHits)
  
  
  ggplot(data = goiLocations, mapping = aes(y = feature, x = percentage, fill = feature))+ 
    geom_bar(stat = 'identity', show.legend = F, colour = 'black', width = 0.7) +
    rcartocolor::scale_fill_carto_d() +
    scale_x_continuous(expand = expansion(mult = c(0,0.05))) +
    theme_TW() +
    theme(axis.ticks.y = element_blank()) +
    ylab('Genomic feature') +
    xlab(bquote('Percentage of interactions in' ~ italic(.(extension))))
  
  ggsave(filename = paste0(outDir, '/genomic_feature_barplot.', extension, '.png'), width = 7, height = 4, units = 'in', dpi = 600)
}

# Run RADIAnT -------------------------------------------------------------

# Parse arguments

library(argparser)

parser = arg_parser(description='Run RADIAnT on RNA-Bin interaction counts')

parser = add_argument(parser, '--gtf', type = 'character', help = 'GTF file for genome of interest (ideally from GENCODE)')

parser = add_argument(parser, '--counts', type = 'character', help = 'RNA-Bin interaction counts from RADIAnT pre-processing pipeline')

parser = add_argument(parser, '--bins', type = 'character', help = 'BED file for bins used in the pre-processing pipeline')

parser = add_argument(parser, '--species', type = 'character', help = 'Species of interest')

parser = add_argument(parser, '--outdir', type = 'character', help = 'Desired output directory')

parser = add_argument(parser, '--name', type = 'character', help = 'Experiment/Sample name')

arg_vector = parse_args(parser)

# Establish required input data from command line arguments ---------------------------------------

gtfFile = arg_vector$gtf
countFile = arg_vector$counts
binAnnotationFile = arg_vector$bins
count_thresh = 2
myL = 200000000
species = arg_vector$species
outDir = arg_vector$outdir
experimentName = arg_vector$name

# Run RADIAnT

inputDataTW = getInputDataTW(countsFile = countFile, 
                             binAnnotationFile = binAnnotationFile, 
                             gtfFile = gtfFile)

countsAndFreqsTW = getCountsAndFreqsTW(counts = inputDataTW$counts, 
                                       binAnnotation = inputDataTW$binAnnotation, 
                                       geneData = inputDataTW$geneData, 
                                       cisWindow = myL)

combinedData = combinedDataTW(cisCounts = countsAndFreqsTW$cisCounts, 
                              transCounts = countsAndFreqsTW$transCounts, 
                              cisBackground = countsAndFreqsTW$cisDistanceCounts, 
                              transBackground = countsAndFreqsTW$transBinMeans)

testableData = combinedData[combinedData$ReadCount >= count_thresh,]

testableData$P = pbsapply(1:nrow(testableData), function(row){
  poisson.test(x = testableData$ReadCount[row],
               r = testableData$ExpectedCount[row],
               alternative = 'greater')$p.value
})

testableData$Padj = p.adjust(p = testableData$P, method = 'fdr', n = nrow(testableData))

combinedData[,c('P', 'Padj')] = testableData[match(combinedData$JoinColumn, testableData$JoinColumn), c('P', 'Padj')]

combinedData$P[is.na(combinedData$P)] = 1

combinedData$Padj[is.na(combinedData$Padj)] = 1

# Write output

combinedData <- combinedData %>% dplyr::rename(InteractionID = JoinColumn)

data.table::fwrite(x = combinedData, file = paste0(outDir, experimentName, 'RADIAnT_results.txt'), quote = F, sep = '\t', eol = '\n', row.names = F, col.names = T)
