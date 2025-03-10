### script to summarise results from Editing analysis by SPRINT.

pdf(NULL)
library(optparse)
library(ggplot2)
library(stringr)
library(patchwork)
set.seed(1234)

## parameters
option_list <- list(
  ### Project
  make_option("--SPRINT_path", help="Path to the SPRINT analysis"),
  make_option("--samples_order_for_ggplot", help="Samples order for ggplot")
)
parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
ggplot_samples_x_order <- unlist(stringr::str_split(args$options$samples_order_for_ggplot, ","))
data_path <- args$options$SPRINT_path
data_output <- data_path

for(res_type in c("identified_all", "identified_hyper", "identified_regular")){
    
    ## get files names
    files <- paste0(data_path,list.files(path = data_path, pattern = paste0("SPRINT_",res_type,".res"), recursive = TRUE))
    nb_samples <- length(files)
    
    ## load and format data
    df <- data.frame()
    for(file in files){
      sample_name <- strsplit(file,"/",fixed=T)[[1]]
      sample_name <- sample_name[length(sample_name)-1]
      tmp_df <- read.table(file, sep="\t", header=FALSE)
      tmp_df$samples=sample_name
      df <- rbind(df,tmp_df)
    }
    colnames(df) <- c("chromonsome", "start", "end", "type", "supporting", "strand", "AD:DP", "samples")
    df$samples <- factor(df$samples, level = ggplot_samples_x_order)
    #change AC to A-to-C etc
    first_base <- str_split_fixed(df$type, "", 2)[,1]
    second_base <- str_split_fixed(df$type, "", 2)[,2]
    df$type <- paste0(first_base,"-to-",second_base)
    
    ##compute percentage of read supporting the edition
    df$total_reads <- str_split_fixed( df$`AD:DP`, ":", 2)[,2]
    df$percent_reads_supporting <- apply(df, 1, function(x){ round((as.numeric(x[5])*100)/as.numeric(x[9]),2) })
    sub_df <- table(df[,c("type","samples")])
    sub_df <- data.frame(edition_types = rownames(sub_df), data.frame(rbind(sub_df)))
    write.table(sub_df, file=paste0(data_output, "Summary_table_counts_SPRINT_",res_type,".tsv"), sep="\t", row.names = FALSE, quote=FALSE)
    sub_df <- table(df[,c("type","samples","strand")])
    write.table(sub_df, file=paste0(data_output, "Summary_table_counts_SPRINT_",res_type,"_by_strand.tsv"), sep="\t", row.names = FALSE, quote=FALSE)
    
    ## graphs
    #number of all edition by sample
    ggplot(data=df, aes(x=samples)) +
      geom_histogram(stat="count", position=position_dodge()) + 
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ylab("Number of edition") +
      xlab("")
    ggsave(paste0(data_output, "Number_of_all_edition_by_sample_",res_type,".png"), width = min(1*nb_samples,109), height = 5, limitsize = FALSE)
    
    #number of each Edition by sample
    ggplot(data=df, aes(x=samples, fill=type)) +
      geom_histogram(stat="count", position=position_dodge()) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ylab("Number of edition") +
      xlab("")
    ggsave(paste0(data_output, "Number_of_each_edition_by_sample_",res_type,".png"), width = min(2*nb_samples,109), height = 5, limitsize = FALSE)
    
    #number of each Edition by sample (log y scale)
    ggplot(data=df, aes(x=samples, fill=type)) +
      geom_histogram(stat="count", position=position_dodge()) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_y_continuous(trans="log10") + 
      ylab("Number of edition (log10 scale)") +
      xlab("")
    ggsave(paste0(data_output, "Number_of_each_edition_by_sample_log_y_scale_",res_type,".png"), width = min(2*nb_samples,109), height = 5, limitsize = FALSE)
    
    #number and percentage of read supporting the Edition by sample
    if(nb_samples == 1){
      plot_percent <- ggplot(data=df, aes(x=percent_reads_supporting, fill=type)) +
              geom_histogram() +
              theme_classic() +
              ggtitle(unique(sort(df$samples))) +
              xlab("Percentages of reads supporting the edition")
      plot_count <- ggplot(data=df, aes(x=supporting, fill=type)) +
              geom_histogram() +
              scale_x_log10(breaks=c(1,2,3,4,5,10,max(df$supporting))) +
              theme_classic() +
              ggtitle(unique(sort(df$samples))) +
              xlab("Number of reads supporting the edition")
      png(paste0(data_output, "Percentages_of_reads_supporting_the_edition_for_each_sample_",res_type,".png"), width = 500, height = 400)
      print(plot_percent)
      dev.off()
      png(paste0(data_output, "Number_of_reads_supporting_the_edition_for_each_sample_",res_type,".png"), width = 500, height = 400)
      print(plot_count)
      dev.off()
    }else{ #several samples
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
      png(paste0(data_output, "Percentages_of_reads_supporting_the_edition_for_each_sample_",res_type,".png"), width = 500, height = min(400*nb_samples,32000))
      print(wrap_plots(plot_percent, nrow=length(plot_percent)))
      dev.off()
      png(paste0(data_output, "Number_of_reads_supporting_the_edition_for_each_sample_",res_type,".png"), width = 500, height = min(400*nb_samples,32000))
      print(wrap_plots(plot_count, nrow=length(plot_count)))
      dev.off()
    }

}

Sys.sleep(10)
