library(ggsankey)
library(dplyr)
library(ggplot2)
library(scales)
library(data.table)



args = commandArgs(trailingOnly = TRUE)
rRNA_log_file = args[1] 
RNA_map_log_file =args[2]
RADIAnT_res_file =args[3]
svg_out = args[4]
png_out = args[5]
txt_out = args[6]



# rRNA_log_file <- "D:/test_files/RADICL_EndoMT/subs.SS_RADICL_N1_CTL_FA_rRNA_removal_stats.txt"
# RADIAnT_res_file <- "D:f
# RNA_map_log_file <- "D:/snakemake_test/bam/subs.SS_RADICL_N1_CTL_FA_RNA_Log.final.out"



stats_template <- data.frame(
  step1 = c("total", "total", "total", "total", "total", "total"),
  step2 = c("retained", "retained", "retained", "retained", "retained", "ribosomal"),
  step3 = c("unique", "unique", "unique", "multi", "unmapped", NA),
  step4 = c("in_interactions", "in_interactions", "not_in_ints", NA, NA, NA),
  step5 = c("in_sig_ints", "not_in_sig_ints", NA, NA, NA, NA),
  stringsAsFactors = FALSE
)


rRNA_log <- data.table::fread(cmd=(paste("sed 's/#//g'", rRNA_log_file)), fill=TRUE)
total_counts <- rRNA_log %>% filter(V1=="Total") %>% select(V2) %>% as.numeric()
ribosomal_counts <- rRNA_log %>% filter(V1=="Matched") %>% select(V2) %>% as.numeric()
retained_counts <- total_counts - ribosomal_counts

read_counts <- data.frame(type=c("total", "ribosomal", "retained"), 
                          count=c(total_counts, ribosomal_counts, retained_counts))

RNA_map_log <- data.table::fread(RNA_map_log_file, fill=TRUE, sep="|", strip.white=TRUE, header=FALSE)

uniquely_mapped <- RNA_map_log %>% filter(V1=="Uniquely mapped reads number") %>% select(V2) %>% as.numeric()

multi_mapped <- RNA_map_log %>% 
  filter(V1=="Number of reads mapped to multiple loci" | V1=="Number of reads mapped to too many loci") %>% 
  mutate(V2 = as.numeric(V2)) %>%
  summarise( sum(V2, na.rm = TRUE)) %>% as.numeric()

unmapped <- RNA_map_log %>% 
  filter(V1=="Number of reads unmapped: too many mismatches" | V1=="Number of reads unmapped: too short" | V1=="Number of reads unmapped: other") %>% 
  mutate(V2 = as.numeric(V2)) %>%
  summarise( sum(V2, na.rm = TRUE)) %>% as.numeric()

read_counts <- read_counts %>% 
                add_row(type="unique", count=uniquely_mapped) %>%
                add_row(type="multi", count=multi_mapped) %>%
                add_row(type="unmapped", count=unmapped)
  
  
RADIAnT_res <- data.table::fread(RADIAnT_res_file)
in_interactions <- RADIAnT_res %>% select(ReadCount) %>% sum()
not_in_ints <- uniquely_mapped - in_interactions
in_sig_ints <- RADIAnT_res %>% filter(Padj<0.05) %>% select(ReadCount) %>% sum()
not_in_sig_ints <- in_interactions - in_sig_ints
read_counts <- read_counts %>% 
                add_row(type="in_interactions", count=in_interactions) %>% 
                add_row(type="not_in_ints", count=not_in_ints) %>%
                add_row(type="in_sig_ints", count=in_sig_ints) %>% 
                add_row(type="not_in_sig_ints", count=not_in_sig_ints)             

to_multiply <- c("in_sig_ints", "not_in_sig_ints", "not_in_ints",  "multi", "unmapped", "ribosomal")





read_counts$rounded_count <- ifelse(round(read_counts$count) < 1e3, read_counts$count, ifelse(round(read_counts$count/ 1e3) < 1e3, round(read_counts$count / 1e3), round(read_counts$count / 1e6)))
read_counts$suffix <- ifelse(round(read_counts$count) < 1e3, "", ifelse(round(read_counts$count/ 1e3) < 1e3, "K", "M"))
read_counts$label <- paste0(read_counts$type, " (", read_counts$rounded_count, read_counts$suffix,  ")", sep = "")
read_counts$label <- gsub("_", " ", read_counts$label)
read_counts$label <- gsub("sig", "sig.", read_counts$label)
read_counts$label <- gsub("ints", "int.", read_counts$label)



# get all counts of the read sets for which rows in ... should be repeated
indices <- match(to_multiply, read_counts$type)
# get the times each corresponding row should be repeated
rep_times <- read_counts$count[indices] %>% rescale(to=c(1,1000)) %>% round()


my_stats <- stats_template[rep(row.names(stats_template), rep_times), ]


stats_long <- my_stats %>% make_long(step1, step2, step3, step4, step5)

stats_long$x <- recode(stats_long$x, step1="", step2="ribosomal cleaning", 
                      step3="mapping", step4="ident. of interactions", step5="sig. test")


read_levels <- c("total", "retained", "ribosomal",
                 "multi", "unmapped",  "unique", 
                 "in_interactions", "not_in_ints", 
                 "non_in_sig_ints", "in_sig_ints", NA)

stats_long$node <- factor(stats_long$node, 
                          levels=read_levels)


stats_long$label <- sapply(stats_long$node, function(node) {
  matching_label <- read_counts$label[grep(paste0("^", node), read_counts$type)]
  if (length(matching_label) > 0) {
    matching_label[1]
  } else {
    NA
  }
})


tryCatch({
stats_long <- stats_long  %>% 
  filter(!is.na(node))

sankey_plot <- ggplot(stats_long, aes(x = x,
               next_x = next_x,
               node = node,
               next_node = next_node,
               label=label, 
               fill=node)) +
  geom_sankey(fill="lightsalmon", flow.alpha = .3, color="black") +
  geom_sankey_label(size=5, color = "black", fill="white") +
  #scale_fill_viridis_d() +
  theme_void() +
  labs(x = NULL) +
  theme(legend.position = "none",
       plot.title = element_text(hjust = .5), 
       axis.text.x = element_text(size=14))

ggsave(filename = svg_out, plot=sankey_plot, device=svg, width = 12, height = 5)
ggsave(filename = png_out, plot=sankey_plot, width = 12, height = 5, units = "in", dpi = 300)
}, error = function(e) {
  # Code to handle errors goes here
  print(paste("An error occurred:", e))
})


read_counts %>% select(type, count) %>% write.table(file=txt_out, quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")


