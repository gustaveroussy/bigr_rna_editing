### script to summarise results from Editing analysis by SPRINT.

library(optparse)
library(ggplot2)
library(stringr)
library(patchwork)
set.seed(1234)

## parameters
option_list <- list(
  ### Project
  make_option("--SPRINT_path", help="Path to the SPRINT analysis"),
  make_option("--samples_order_for_ggplot", help="Samples order for ggplot"),
)
parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
ggplot_samples_x_order <- unlist(stringr::str_split(args$options$samples_order_for_ggplot, ","))
data_path <- args$options$SPRINT_path
files <- paste0(data_path,list.files(path = data_path, pattern = "SPRINT_identified_hyper.res", recursive = TRUE))
data_output <- data_path

## load and format data
df <- data.frame()
for(file in files){
  sample_name=strsplit(file,"/",fixed=T)[[1]][11] ### to do: to check if it works
  tmp_df <- read.table(file, sep="\t", header=FALSE)
  tmp_df$samples=sample_name
  df <- rbind(df,tmp_df)
}
colnames(df) <- c("chromonsome", "start", "end", "type", "supporting", "strand", "AD:DP", "samples")
df$samples <- factor(df$samples, level = ggplot_samples_x_order)
# change AC to A-to-C etc
first_base <- str_split_fixed(df$type, "", 2)[,1]
second_base <- str_split_fixed(df$type, "", 2)[,2]
df$type <- paste0(first_base,"-to-",second_base)

##compute percentage of read supporting the edition
df$total_reads <- str_split_fixed( df$`AD:DP`, ":", 2)[,2]
df$percent_reads_supporting <- apply(df, 1, function(x){ round((as.numeric(x[5])*100)/as.numeric(x[9]),2) })
sub_df <- table(df[,c("type","samples")])
write.table(sub_df, file=paste0(data_output, "Summary_table_counts_SPRINT.tsv"), sep="\t", quote=FALSE)
sub_df <- table(df[,c("type","samples","strand")])
write.table(sub_df, file=paste0(data_output, "Summary_table_counts_SPRINT_by_strand.tsv"), sep="\t", quote=FALSE)

## graphs
#number of all edition by sample
ggplot(data=df, aes(x=samples)) +
  geom_histogram(stat="count", position=position_dodge()) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Number of edition") +
  xlab("")
ggsave(paste0(data_output, "Number_of_all_edition_by_sample.png"), width = 8, height = 5)

#number of each Edition by sample
ggplot(data=df, aes(x=samples, fill=type)) +
  geom_histogram(stat="count", position=position_dodge()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Number of edition") +
  xlab("")
ggsave(paste0(data_output, "Number_of_each_edition_by_sample.png"), width = 16, height = 5)

#number of each Edition by sample (log y scale)
ggplot(data=df, aes(x=samples, fill=type)) +
  geom_histogram(stat="count", position=position_dodge()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(trans="log10") + 
  ylab("Number of edition (log10 scale)") +
  xlab("")
ggsave(paste0(data_output, "Number_of_each_edition_by_sample_log_y_scale.png"), width = 16, height = 5)

#number and percentage of read supporting the Edition by sample
plot_percent <- list()
plot_count <- list()
for (sample in unique(sort(df$samples))){
  sub_sample_df <- df[ which(df$samples==sample), ]
  plot_percent[[sample]] <- ggplot(data=sub_sample_df, aes(x=percent_reads_supporting, fill=type)) +
          geom_histogram() +
          theme_classic() +
          ggtitle(sample) +
          xlab("Percentages of reads supporting the edition")
  plot_count[[sample]] <- ggplot(data=sub_sample_df, aes(x=supporting, fill=type)) +
          geom_histogram() +
          scale_x_log10(breaks=c(1,2,3,4,5,10,max(sub_sample_df$supporting))) +
          theme_classic() +
          ggtitle(sample) +
          xlab("Number of reads supporting the edition")
}
png(paste0(data_output, 'Percentages_of_reads_supporting_the_edition_for_each_sample.png'), width = 500, height = 4500)
wrap_plots(plot_percent, nrow=length(plot_percent))
dev.off()
png(paste0(data_output, 'Number_of_reads_supporting_the_edition_for_each_sample.png'), width = 500, height = 4500)
wrap_plots(plot_count, nrow=length(plot_count))
dev.off()
