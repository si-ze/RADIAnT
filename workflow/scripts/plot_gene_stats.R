# library(data.table)
# library(ggplot2)
# library(ggarchery)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# #library(EnsDb.Hsapiens.v86)
# #library(EnsDb.Mmusculus.v79)
# library(org.Mm.eg.db)
# library(org.Hs.eg.db)
# library(ChIPpeakAnno)
# library(GenomicRanges)
# library(GenomicFeatures)
# library(IRanges)
# library(rcartocolor)
# library(ggrepel)
# library(dplyr)
# library(purrr)
# library(pbapply)

# library(RColorBrewer)
# 
# library(ComplexUpset)
# library(viridis)

# library(GenomicRanges)
# library(GenomicFeatures)
# library(IRanges)


library(data.table)
library(dplyr)
library(ggplot2)
# library(ggrepel)
library(forcats)
library(GenomicRanges)
library(GenomicFeatures)
library(ggVennDiagram)


args = commandArgs(trailingOnly = TRUE)
res_file <- args[1] 
gtf_file <- args[2] 
interactions_dir <- args[3] 
txt_output_file <- args[4]

# # DELETE
# setwd("/mnt/d/sankey_tst/test_run/RADICL/")
# res_file <- "interactions/mESC_2FA_n1_RADIAnT_results.txt"
# gtf_file <- "/mnt/d/upload_test/RADIAnT/resources/mouse/gencode.vM29.annotation.gtf.gz"
# interactions_dir <- "/mnt/d/sankey_tst/test_run/RADICL/logs/"
# # DELETE


res <- data.table::fread(res_file) %>% as.data.frame()

#########################################################
# Plot significant interactions 
#########################################################



message("==================================================================")
message("========= Plotting top 15 genes (by # sig. interactions) =========")
message("==================================================================")

res$Type <- ifelse(res$BinChr==res$GeneChr, "Intrachromosomal", "Transchromosomal")


sig_counts_df <- res %>%
  filter(Padj<0.05) %>%
  group_by(Symbol, Type) %>%
  summarise(
    SigIntsDiff = n(),
    SigReadsDiff = sum(ReadCount)
)

sig_counts_df <- sig_counts_df %>% 
  group_by(Symbol) %>%
  mutate(SigInts = sum(SigIntsDiff), 
         SigReads=sum(SigReadsDiff))





tryCatch({
  sig_counts_df <- sig_counts_df %>%
    arrange(desc(SigInts))
  stackedBarData <- sig_counts_df %>% head(15)
  stackedBarData$Type <- factor(stackedBarData$Type, levels = c("Intrachromosomal", "Transchromosomal"))
  stackedBarData$label_y <- ifelse(stackedBarData$Type=="Intrachromosomal", stackedBarData$SigIntsDiff*0.05, stackedBarData$SigInts)
  my_cols <- c("Intrachromosomal" = "grey70", "Transchromosomal" = "grey90")
  my_order <- stackedBarData %>% arrange(-SigInts) %>% pull(Symbol) %>% unique() %>% as.vector()
  stackedBarData$Symbol <- factor(stackedBarData$Symbol, levels = my_order)
  ggplot(stackedBarData, aes(fill=forcats::fct_rev(Type), y=SigIntsDiff, x=Symbol)) + 
      geom_bar(color=NA, position="stack", stat="identity") +
      geom_text(aes(label=scales::comma(SigIntsDiff), y=label_y), size=12/.pt, position="stack", vjust=-0.25) +
      theme_bw() +
      theme(
        axis.title= element_text(colour = 'black', size=12),
        axis.text = element_text(colour = 'black', size=12),
        axis.ticks = element_line(colour = 'black'),
        axis.line = element_line(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(colour = NA), 
        legend.position="top", 
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) +
      scale_fill_manual("Interaction Type", values=my_cols) + 
      scale_y_continuous("Significant interactions", labels = function(x) format(x, scientific = TRUE)) +
      scale_x_discrete("Gene") +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
      ggsave(paste0(interactions_dir, "/genes_barplot.significant_interactions.png"), width = 12, height = 5, units = 'in', dpi = 600)
      ggsave(paste0(interactions_dir, "/genes_barplot.significant_interactions.svg"), width = 12, height = 5, dpi = 600, device=svg)

}, error = function(e) {
  # Code to handle errors goes here
  print(paste("An error occurred:", e))
})


message("==================================================================")
message("=========== Plotting top 15 genes (by # sig. reads) ==============")
message("==================================================================")
 
tryCatch({
  sig_counts_df <- sig_counts_df %>%
    arrange(desc(SigReads))
  stackedBarData <- sig_counts_df %>% head(15)
  stackedBarData$Type <- factor(stackedBarData$Type, levels = c("Intrachromosomal", "Transchromosomal"))
  stackedBarData$label_y <- ifelse(stackedBarData$Type=="Intrachromosomal", stackedBarData$SigReadsDiff*0.05, stackedBarData$SigReads)
  my_cols <- c("Intrachromosomal" = "grey70", "Transchromosomal" = "grey90")
  my_order <- stackedBarData %>% arrange(-SigReads) %>% pull(Symbol) %>% unique() %>% as.vector()
  stackedBarData$Symbol <- factor(stackedBarData$Symbol, levels = my_order)
  ggplot(stackedBarData, aes(fill=forcats::fct_rev(Type), y=SigReadsDiff, x=Symbol)) + 
      geom_bar(color=NA, position="stack", stat="identity") +
      geom_text(aes(label=scales::comma(SigReadsDiff), y=label_y), size=12/.pt, position="stack", vjust=-0.25) +
      theme_bw() +
      theme(
        axis.title= element_text(colour = 'black', size=12),
        axis.text = element_text(colour = 'black', size=12),
        axis.ticks = element_line(colour = 'black'),
        axis.line = element_line(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(colour = NA), 
        legend.position="top", 
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) +
      scale_fill_manual("Interaction Type", values=my_cols) + 
      scale_y_continuous("Significant reads", labels = function(x) format(x, scientific = TRUE)) +
      scale_x_discrete("Gene") +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
      ggsave(paste0(interactions_dir, "/genes_barplot.significant_reads.png"), width = 12, height = 5, units = 'in', dpi = 600)
      ggsave(paste0(interactions_dir, "/genes_barplot.significant_reads.svg"), width = 12, height = 5, units = 'in', device=svg)
}, error = function(e) {
  # Code to handle errors goes here
  print(paste("An error occurred:", e))
})



message("==================================================================")
message("============ Plotting top 15 genes (by # interactions) ===========")
message("==================================================================")

all_counts_df <- res %>%
  group_by(Symbol, Type) %>%
  summarise(
    IntsDiff = n(),
    ReadsDiff = sum(ReadCount)
)

all_counts_df <- all_counts_df %>% 
  group_by(Symbol) %>%
  mutate(Ints = sum(IntsDiff), 
  Reads=sum(ReadsDiff))




tryCatch({
  all_counts_df <- all_counts_df %>% arrange(desc(Ints))
  stackedBarData <- all_counts_df %>% head(15)
  stackedBarData$Type <- factor(stackedBarData$Type, levels = c("Intrachromosomal", "Transchromosomal"))
  stackedBarData$label_y <- ifelse(stackedBarData$Type=="Intrachromosomal", stackedBarData$IntsDiff*0.05, stackedBarData$Ints)
  my_cols <- c("Intrachromosomal" = "grey70", "Transchromosomal" = "grey90")
  my_order <- stackedBarData %>% arrange(-Ints) %>% pull(Symbol) %>% unique() %>% as.vector()
  stackedBarData$Symbol <- factor(stackedBarData$Symbol, levels = my_order)
  ggplot(stackedBarData, aes(fill=forcats::fct_rev(Type), y=IntsDiff, x=Symbol)) + 
      geom_bar(color=NA, position="stack", stat="identity") +
      geom_text(aes(label=scales::comma(IntsDiff), y=label_y), size=12/.pt, position="stack", vjust=-0.25) +
      theme_bw() +
      theme(
        axis.title= element_text(colour = 'black', size=12),
        axis.text = element_text(colour = 'black', size=12),
        axis.ticks = element_line(colour = 'black'),
        axis.line = element_line(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(colour = NA), 
        legend.position="top", 
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) +
      scale_fill_manual("Interaction Type", values=my_cols) + 
      scale_y_continuous("Total number of interactions", labels = function(x) format(x, scientific = TRUE)) +
      scale_x_discrete("Gene") +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
      ggsave(paste0(interactions_dir, "/genes_barplot.all_interactions.png"), width = 12, height = 5, units = 'in', dpi = 600)
      ggsave(paste0(interactions_dir, "/genes_barplot.all_interactions.svg"), width = 12, height = 5, units = 'in', device=svg)
}, error = function(e) {
  # Code to handle errors goes here
  print(paste("An error occurred:", e))
})


message("==================================================================")
message("============== Plotting top 15 genes (by # reads) ================")
message("==================================================================")
 
 
tryCatch({
  all_counts_df <- all_counts_df %>% arrange(desc(Reads))
  stackedBarData <- all_counts_df %>% head(15)
  stackedBarData$Type <- factor(stackedBarData$Type, levels = c("Intrachromosomal", "Transchromosomal"))
  stackedBarData$label_y <- ifelse(stackedBarData$Type=="Intrachromosomal", stackedBarData$ReadsDiff*0.05, stackedBarData$Reads)
  my_cols <- c("Intrachromosomal" = "grey70", "Transchromosomal" = "grey90")
  my_order <- stackedBarData %>% arrange(-Reads) %>% pull(Symbol) %>% unique() %>% as.vector()
  stackedBarData$Symbol <- factor(stackedBarData$Symbol, levels = my_order)
  ggplot(stackedBarData, aes(fill=forcats::fct_rev(Type), y=ReadsDiff, x=Symbol)) + 
      geom_bar(color=NA, position="stack", stat="identity") +
      geom_text(aes(label=scales::comma(ReadsDiff), y=label_y), size=12/.pt, position="stack", vjust=-0.25) +
      theme_bw() +
      theme(
        axis.title= element_text(colour = 'black', size=12),
        axis.text = element_text(colour = 'black', size=12),
        axis.ticks = element_line(colour = 'black'),
        axis.line = element_line(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(colour = NA), 
        legend.position="top", 
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) +
      scale_fill_manual("Interaction Type", values=my_cols) + 
      scale_y_continuous("Total number of reads", labels = function(x) format(x, scientific = TRUE)) +
      scale_x_discrete("Gene") +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
      ggsave(paste0(interactions_dir, "/genes_barplot.all_reads.png"), width = 12, height = 5, units = 'in', dpi = 600)
      ggsave(paste0(interactions_dir, "/genes_barplot.all_reads.svg"), width = 12, height = 5, units = 'in',device=svg)
}, error = function(e) {
  # Code to handle errors goes here
  print(paste("An error occurred:", e))
})

message("==================================================================")
message("============== Venn diagram  ================")
message("==================================================================")





gtf = data.table::fread(gtf_file, sep = '\t', header = F)
genes = gtf %>% filter(V3=="gene") %>% mutate(symbol=sub('.*gene_name "', '', sub('"; level.*', '', V9))) 
genes_annotated <- genes %>% dplyr::select(symbol) %>% unique() %>% unlist() %>% as.vector()

genes_with_counts <- res %>% filter(ReadCount>0) %>% pull(Symbol) %>% unique() %>% as.vector()

genes_with_sig_cis_int <- res %>% filter(Type=="Intrachromosomal" & Padj<0.05) %>% pull(Symbol) %>% unique() %>% as.vector()
genes_with_sig_trans_int <- res %>% filter(Type=="Transchromosomal" & Padj<0.05) %>% pull(Symbol) %>% unique() %>% as.vector()


format_labels <- function(l) {
  rounded <-  ifelse(l<=1e3, l, ifelse(l <= 1e6, round(l / 1e3), round(l / 1e6)))
  suffix <- ifelse(l<=1e3, "", ifelse(l <= 1e6, "K", "M"))
  label <- paste0(rounded, suffix)
  return(label)
}


set_lists = list(
   "Annotated"=unlist(genes_annotated),
#   "with interactions"=unlist(genes_with_int),
#  "int. in cis"=unlist(genes_with_cis_int),
  "Significant\nIntrachromosomal"=unlist(genes_with_sig_cis_int), 
#  "int. in trans"=unlist(genes_with_trans_int),
  "Significant\nTranschromosomal"=unlist(genes_with_sig_trans_int)
)


tryCatch({
  ggVennDiagram(set_lists, 
                  label_alpha=0, 
                  label_geom="label",
                  edge_size=2) + 
  scale_fill_gradient(low=NA,high =NA) +
  theme_void() + 
  theme(legend.position = "none")
  ggsave(paste0(interactions_dir, "/genes_venn.total_intra_trans.png"), width = 10, height = 10, units = 'in', dpi = 600)
  ggsave(paste0(interactions_dir, "/genes_venn.total_intra_trans.svg"), width = 10, height = 10, units = 'in', device=svg)
}, error = function(e) {
  # Code to handle errors goes here
  print(paste("An error occurred:", e))
})






message("==================================================================")
message("============== Output text  ================")
message("==================================================================")
txt_output <- merge(all_counts_df[,c("Symbol", "Type", "IntsDiff")], sig_counts_df[,c("Symbol", "Type", "SigIntsDiff")], by=c("Symbol", "Type"))
names(txt_output) <- c("Symbol", "Type", "Total_Interactions", "Significant_Interactions")
txt_output %>% write.table(file=txt_output_file, quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")










