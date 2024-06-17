### script to summarise results from Editing analysis by RNAEditingIndexer.

pdf(NULL)
library(optparse)
library(ggplot2)
library(stringr)
library(tidyr)
library(dplyr)
library(patchwork)
set.seed(1234)

## parameters
option_list <- list(
  ### Project
  make_option("--EditingIndex", help="Path to the result table of RNAEditingIndexer analysis"),
  make_option("--samples_order_for_ggplot", help="Samples order for ggplot"),
)
parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
ggplot_samples_x_order <- unlist(stringr::str_split(args$options$samples_order_for_ggplot, ","))
files <- args$options$EditingIndex
data_output <- gsub("summary/EditingIndex.csv", "", file)


## load and format data
df <- read.table(file, sep=",", header=TRUE)
# Lines:
# RefSeqThenMMSites: signal of interest
# MMSitesThenRefSeq: opposé de RefSeqThenMMSites. Identification sur les sites de mismatch.
# Randomly lines: Il s'agit principalement d'un contrôle négatif. Idéalement, le signal serait deux fois moins élevé, le bruit toujours le même.
# Columns:
# <mismatch type>EditingIndex: score of interest
# NumOf<base type>PositionsCovered: Number of reference genome positions of it indexed
# TotalCoverageAt<base type>Positions: Coverage at these positions (including mismatches, without SNPs)
# So number of reads with edition = TotalCoverageAt<base type>Positions - NumOf<base type>PositionsCovered

## filtering table to keep only interesting columns
df <- subset(df, StrandDecidingMethod %in% c("Randomly","RefSeqThenMMSites"))
column_to_keep <- grep("^StrandDecidingMethod$|^Sample$|EditingIndex|PositionsCovered|TotalCoverageAt", colnames(df), value = TRUE)
df <- df[,column_to_keep]

## format
column_to_keep <- grep("^StrandDecidingMethod$|^Sample$|EditingIndex", colnames(df), value = TRUE)

df_pivot <- df %>% select(all_of(column_to_keep)) %>% pivot_longer(cols = A2CEditingIndex:C2TEditingIndex,
                                names_to = "Edition_types",
                                values_to = "Edition_score") %>% as.data.frame()

head(df_pivot)
df_pivot$Edition_types <- str_replace(str_replace(df_pivot$Edition_types, "2", "-to-"), "EditingIndex", "")
sub_df_pivot <- subset(df_pivot, StrandDecidingMethod  == "RefSeqThenMMSites")

## graphs
#number of each Edition by sample
ggplot(data=sub_df_pivot, aes(x=factor(Sample, levels = ggplot_samples_x_order), y = Edition_score, fill = Edition_types )) +
  geom_bar(stat="identity", position=position_dodge()) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Scores of edition") +
  xlab("")
ggsave(paste0(data_output, "Scores_of_all_edition_by_sample.png"), width = 8, height = 5)

#number of each Edition by sample (log y scale) -> impossible to see something log10(0)=1 and log10(0.5)=negatif, so graph is reverted in value and impossible to interpret easily.


#graphs of editing scores and background noise for each edition type
edition_types <- grep("EditingIndex", colnames(df), value = TRUE)
# plot_score_noise <- list()
# for (edition in edition_types){
#  print(edition)
#  edition_name <- str_replace(str_replace(edition, "2", "-to-"), "EditingIndex", "")
#  plot_score_noise[[edition]] <- ggplot(data=df, aes(x=factor(Sample, levels = ggplot_samples_x_order), y=get(edition),fill=StrandDecidingMethod)) +
#                              geom_bar(stat="identity", position=position_dodge()) +
#                              theme_classic() +
#                              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#                              ylab(paste0("Scores of ", edition_name," edition")) +
#                              xlab("") +
#                              scale_fill_brewer(palette="Paired")
#  print(plot_score_noise[[edition]])
# }
# 
# png(paste0(data_output, 'Scores_of_edition_and_background_noise_for_each_edition_type.png'), width = 500, height = 1500)
# wrap_plots(plot_score_noise, nrow=length(plot_score_noise))
# dev.off()
edition <- "A2CEditingIndex"
edition_name <- str_replace(str_replace(edition, "2", "-to-"), "EditingIndex", "")
p1 <- ggplot(data=df, aes(x=factor(Sample, levels = ggplot_samples_x_order), y=A2CEditingIndex,fill=StrandDecidingMethod)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab(paste0("Scores of ", edition_name," edition")) +
  xlab("") +
  scale_fill_brewer(palette="Paired")
edition <- "A2GEditingIndex"
edition_name <- str_replace(str_replace(edition, "2", "-to-"), "EditingIndex", "")
p2 <- ggplot(data=df, aes(x=factor(Sample, levels = ggplot_samples_x_order), y=A2GEditingIndex,fill=StrandDecidingMethod)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab(paste0("Scores of ", edition_name," edition")) +
  xlab("") +
  scale_fill_brewer(palette="Paired")
edition <- "A2TEditingIndex"
edition_name <- str_replace(str_replace(edition, "2", "-to-"), "EditingIndex", "")
p3 <- ggplot(data=df, aes(x=factor(Sample, levels = ggplot_samples_x_order), y=A2TEditingIndex,fill=StrandDecidingMethod)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab(paste0("Scores of ", edition_name," edition")) +
  xlab("") +
  scale_fill_brewer(palette="Paired")
edition <- "C2AEditingIndex"
edition_name <- str_replace(str_replace(edition, "2", "-to-"), "EditingIndex", "")
p4 <- ggplot(data=df, aes(x=factor(Sample, levels = ggplot_samples_x_order), y=C2AEditingIndex,fill=StrandDecidingMethod)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab(paste0("Scores of ", edition_name," edition")) +
  xlab("") +
  scale_fill_brewer(palette="Paired")
edition <- "C2GEditingIndex"
edition_name <- str_replace(str_replace(edition, "2", "-to-"), "EditingIndex", "")
p5 <- ggplot(data=df, aes(x=factor(Sample, levels = ggplot_samples_x_order), y=C2GEditingIndex,fill=StrandDecidingMethod)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab(paste0("Scores of ", edition_name," edition")) +
  xlab("") +
  scale_fill_brewer(palette="Paired")
edition <- "C2TEditingIndex"
edition_name <- str_replace(str_replace(edition, "2", "-to-"), "EditingIndex", "")
p6 <- ggplot(data=df, aes(x=factor(Sample, levels = ggplot_samples_x_order), y=C2TEditingIndex,fill=StrandDecidingMethod)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab(paste0("Scores of ", edition_name," edition")) +
  xlab("") +
  scale_fill_brewer(palette="Paired")
png(paste0(data_output, 'Scores_of_edition_and_background_noise_for_each_edition_type.png'), width = 500, height = 1500)
wrap_plots(list(p1,p2,p3,p4,p5,p6), nrow=6)
dev.off()

#coverage by sample
sub_df <- subset(df, StrandDecidingMethod  == "RefSeqThenMMSites") %>% 
              mutate(NumOfAPositionsCovered_Edited = TotalCoverageAtAPositions - NumOfAPositionsCovered) %>%
              mutate(NumOfCPositionsCovered_Edited = TotalCoverageAtCPositions - NumOfCPositionsCovered) %>%
              mutate(NumOfGPositionsCovered_Edited = TotalCoverageAtGPositions - NumOfGPositionsCovered) %>%
              mutate(NumOfTPositionsCovered_Edited = TotalCoverageAtTPositions - NumOfTPositionsCovered) %>%
              select(StrandDecidingMethod,Sample,NumOfAPositionsCovered_Edited,TotalCoverageAtAPositions,NumOfCPositionsCovered_Edited,TotalCoverageAtCPositions,NumOfGPositionsCovered_Edited,TotalCoverageAtGPositions,NumOfTPositionsCovered_Edited,TotalCoverageAtTPositions) %>% 
              pivot_longer(cols = NumOfAPositionsCovered_Edited:TotalCoverageAtTPositions,
                    names_to = "Coverage_type",
                    values_to = "Coverage") %>% 
              as.data.frame()

ggplot(data=sub_df, aes(x=factor(Sample, levels = ggplot_samples_x_order), y = Coverage, fill = factor(Coverage_type, levels = c("NumOfAPositionsCovered_Edited","TotalCoverageAtAPositions","NumOfCPositionsCovered_Edited","TotalCoverageAtCPositions","NumOfGPositionsCovered_Edited","TotalCoverageAtGPositions","NumOfTPositionsCovered_Edited","TotalCoverageAtTPositions")))) +
  geom_bar(stat="identity", position=position_dodge()) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Scores of edition") +
  xlab("") +
  guides(fill=guide_legend(title="")) + 
  scale_fill_brewer(palette="Paired")
ggsave(paste0(data_output, "Coverage_by_base_reference_by_sample.png"), width =10, height = 5)

