setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")

list.of.packages <- c("easypackages", "ggplot2", "readODS", "colorspace")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# For the printed files
num_to_return <- 1
exe_num <- sample(1:as.numeric(10000), num_to_return)

thisFontSize <- 15

library("easypackages")
libraries(list.of.packages)

data_table <- read_ods("../results/pancancer_signature_results_recap.ods", sheet=5)
data_table$"cancer_type" <- gsub("_", " ", data_table$"cancer_type")  

data_table$"normMCC"  <- (data_table$"MCC" + 1) / 2

list_of_rates <- c("F1_score", "accuracy", "TPR", "TNR", "PPV", "NPV", "PR_AUC", "ROC_AUC", "normMCC")
# list_of_rates <- c("F1_score")
theseColors <- setNames(rainbow_hcl(nrow(data_table)), levels(data_table$"cancer_type")  )

data_table$"cancer_type_factor" <- factor(data_table$"cancer_type", levels = data_table$"cancer_type")
data_table$"cancer_type" <- data_table$"cancer_type_factor"

for(this_rate in list_of_rates){

    ylim_low <- 0
    ylim_upp <- 1
    
    this_rate_without_underscore <- gsub("_", " ", this_rate)
    
    data_table <- data_table[order(-data_table[[this_rate]]),]

   # p_this_rate_plot <- ggplot(data=data_table, mapping=aes(x=reorder(cancer_type, -.data[[this_rate]]), y=.data[[this_rate]], color=cancer_type))  + geom_bar(stat="identity") + scale_color_manual(values=theseColors)
    
    p_this_rate_plot <- ggplot(data_table, aes(x=reorder(cancer_type, -.data[[this_rate]]), y=.data[[this_rate]], fill=cancer_type))  + geom_bar(stat="identity")  + ylab("") + xlab("") + ggtitle(this_rate_without_underscore)  +  theme(plot.title = element_text(hjust = 0.5, size = thisFontSize), legend.title=element_text(size=thisFontSize), legend.text=element_text(size=thisFontSize), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
    
#     theme(plot.title = element_text(size = 12, face = "bold"),
#     legend.title=element_text(size=10), 
#     legend.text=element_text(size=9))
    
    p_this_rate_plot <- p_this_rate_plot + scale_fill_discrete(breaks=data_table$"cancer_type")
    
    p_this_rate_plot <- p_this_rate_plot + scale_y_continuous(limits=c(ylim_low, ylim_upp))
    
    folderName <- paste0("../results/barplots_", exe_num, "/")
    createFolderCommand <- paste0("mkdir -p ", folderName)
    system(createFolderCommand)
    cat("Executed command: ", createFolderCommand, "\n", sep="")
    
    pdfThisPlotFile <- paste0(folderName,"pancancer_", this_rate,"_barplot", exe_num,".pdf")
    cat("We're going to save the ", pdfThisPlotFile, " file\n", sep="")
    ggsave(pdfThisPlotFile)
}

# 
# # MCC
# 
# MCC_BARPLOT_FLAG <- FALSE
# if(MCC_BARPLOT_FLAG == TRUE) {
# 
#     ylim_mcc_low <- -1
#     ylim_mcc_upp <- 1
# 
#     data_table <- data_table[order(-data_table[c("MCC")]),]
#     data_table$"cancer_type_factor" <- factor(data_table$"cancer_type", levels = data_table$"cancer_type")
#     data_table$"cancer_type" <- data_table$"cancer_type_factor"
#     
#     p_mcc_plot <- ggplot(data_table, aes(x=reorder(cancer_type, -MCC), y=MCC, fill=cancer_type)) + geom_bar(stat="identity", color="black",  position=position_dodge())  + ylab("mean MCC") + xlab("") + ggtitle("survival binary prediction")  +
#     theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
#     
#     p_mcc_plot <- p_mcc_plot + scale_y_continuous(limits=c(ylim_mcc_low, ylim_mcc_upp)) #, breaks=c(1:10)) 
#     # p_mcc_plot <- p_mcc_plot + theme_classic()
# 
#     pdfMccPlotFile <- paste0("../results/pancancer_MCC_barplot", exe_num,".pdf")
#     ggsave(pdfMccPlotFile)
# }
# 
# # PR AUC
# 
# PR_AUC_BARPLOT_FLAG <- FALSE
# if(PR_AUC_BARPLOT_FLAG == TRUE) {
# 
#     ylim_low <- 0
#     ylim_upp <- 1
#     
#     data_table <- data_table[order(-data_table[c("PR_AUC")]),]
#     data_table$"cancer_type_factor" <- factor(data_table$"cancer_type", levels = data_table$"cancer_type")
#     data_table$"cancer_type" <- data_table$"cancer_type_factor"
#     
#     p_PR_AUC_plot <- ggplot(data_table, aes(x=reorder(cancer_type, -PR_AUC), y=PR_AUC, fill=cancer_type)) + geom_bar(stat="identity", color="black",  position=position_dodge())  + ylab("mean PR AUC") + xlab("") + ggtitle("survival binary prediction")  +
#     theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
#     
#     p_PR_AUC_plot <- p_PR_AUC_plot + scale_y_continuous(limits=c(ylim_low, ylim_upp)) #, breaks=c(1:10)) 
#     # p_PR_AUC_plot <- p_PR_AUC_plot + theme_classic()
# 
#     pdfPRAUCplotFile <- paste0("../results/pancancer_PR_AUC_barplot", exe_num,".pdf")
#     ggsave(pdfPRAUCplotFile)
# }
# 
# # ROC AUC
# 
# ROC_AUC_BARPLOT_FLAG <- FALSE
# if(ROC_AUC_BARPLOT_FLAG == TRUE) {
# 
#     ylim_low <- 0
#     ylim_upp <- 1
#     
#     data_table <- data_table[order(-data_table[c("ROC_AUC")]),]
#     data_table$"cancer_type_factor" <- factor(data_table$"cancer_type", levels = data_table$"cancer_type")
#     data_table$"cancer_type" <- data_table$"cancer_type_factor"
#     
#     p_ROC_AUC_plot <- ggplot(data_table, aes(x=reorder(cancer_type, -ROC_AUC), y=ROC_AUC, fill=cancer_type)) + geom_bar(stat="identity", color="black",  position=position_dodge())  + ylab("mean ROC AUC") + xlab("") + ggtitle("survival binary prediction")  +
#     theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
#     
#     p_ROC_AUC_plot <- p_ROC_AUC_plot + scale_y_continuous(limits=c(ylim_low, ylim_upp)) #, breaks=c(1:10)) 
#     # p_ROC_AUC_plot <- p_ROC_AUC_plot + theme_classic()
# 
#     pdfROCAUCplotFile <- paste0("../results/pancancer_ROC_AUC_barplot", exe_num,".pdf")
#     ggsave(pdfROCAUCplotFile)
# }
