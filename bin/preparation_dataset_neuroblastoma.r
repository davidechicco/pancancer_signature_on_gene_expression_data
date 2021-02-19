setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
# set.seed(11)

list.of.packages <- c("easypackages", "plyr") # other packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

source("utils.r")

# For the printed files
num_to_return <- 1
exe_num <- sample(1:as.numeric(10000), num_to_return)


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
listOfBiocPackages <- c("annotate", "GEOquery")

bioCpackagesNotInstalled <- which( !listOfBiocPackages %in% rownames(installed.packages()) )
cat("package missing listOfBiocPackages[", bioCpackagesNotInstalled, "]: ", listOfBiocPackages[bioCpackagesNotInstalled], "\n", sep="")

# check there's still something left to install
if( length(bioCpackagesNotInstalled) ) {
    BiocManager::install(listOfBiocPackages[bioCpackagesNotInstalled])
}

library("easypackages")
libraries(list.of.packages)
libraries(listOfBiocPackages)

GSE_code <- "GSE16237"
thisGEOplatform <- "GPL570"
datasetName <-  "Hiyama2010"


gset <- getGEO(GSE_code,  GSEMatrix =TRUE, getGPL=FALSE)

if (length(gset) > 1) idx <- grep(thisGEOplatform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
gene_expression <- as.data.frame(exprs(gset))

cat("str(gset@phenoData@data)\n")
print(str(gset@phenoData@data))
cat("str(gset@phenoData@data)\n")

LABEL_DETECTED <- FALSE
 
 if(LABEL_DETECTED == TRUE) {
 
    # # # # we add the labels
    library("plyr")
    label_list <- c()
    i <- 1
    for(thisTitle in gset@phenoData@data$"outcome of the patient:ch1") {
      
	if(grepl("Died of disease", thisTitle)) {
	      label_list[[i]] <-  0
	 } else if(grepl("Alive", thisTitle)) {
	    label_list[[i]] <- 1
	  }
	    i <- i + 1
        }
      
    cat("label_list:\n")
    print(label_list)
    
    targetName <- "survival"
    
    labels_df_temp <- as.data.frame(label_list)
    
    labels_df <- as.data.frame(t(labels_df_temp))
    colnames(labels_df) <- colnames(gene_expression)
    rownames(labels_df) <- targetName
    gene_expression_with_labels <- rbind(labels_df, gene_expression)
    
    gene_expression_with_labels$ID <- rownames(gene_expression_with_labels)

    
    ## selection based on signature probesets

    DCangelosi_probesets_signature <- c("200738_s_at", "17356_s_at", "206686_at", "226452_at", "223172_s_at",  "223193_x_at",  "224314_s_at", "230630_at", "202022_at")
    
    patients_data_filtered_ourSignature <- (gene_expression_with_labels[gene_expression_with_labels$ID %in% c(DCangelosi_probesets_signature, "survival"), ])
      
    print("dim(patients_data_filtered_ourSignature)\n")
    print(dim(patients_data_filtered_ourSignature))

    # patients_data_filtered_ourSignature <- rbind(gene_expression_with_labels[1,], patients_data_filtered_ourSignature)

    patients_data_filtered_ourSignature$"ID" <- NULL
    
    patients_data_filtered_ourSignature_t <- as.data.frame(t(patients_data_filtered_ourSignature))

    library("dplyr")
    patients_data_filtered_ourSignature_t <- patients_data_filtered_ourSignature_t %>% dplyr::select(-targetName,targetName)

    patients_data_filtered_ourSignature_t <-  slice(patients_data_filtered_ourSignature_t, 1:(n()-1))    
    patients_data_filtered_ourSignature_t[,targetName] <- as.numeric(patients_data_filtered_ourSignature_t[,targetName]) 

    tableOnlyOnesAndZeros <-  rbind(patients_data_filtered_ourSignature_t[patients_data_filtered_ourSignature_t$"survival"==0,], patients_data_filtered_ourSignature_t[patients_data_filtered_ourSignature_t$"survival"==1,])

    folderPath <- "../results/"
    
    outputFileName <- paste0(folderPath, datasetName, "_", GSE_code, "_probesets_dataset_", exe_num, ".csv")
    write.table(tableOnlyOnesAndZeros, file=outputFileName, row.names=TRUE, sep=",", col.names=NA)
    cat("saved file  ", outputFileName, "\n")
    
    
    ## selection based on signature gene symbols
    
    DCangelosi_gene_symbols_signature <- c("PGK1", "PDK1", "MTFP1", "FAM162A", "EGLN1", "AK4", "ALDOC")

    library("geneExpressionFromGEO")
    associateSymbolsToGenes <- TRUE
    verbose <- TRUE
    geneExpressionDF <- getGeneExpressionFromGEO(GSE_code,  associateSymbolsToGenes, verbose)
    
    labels_df$"GeneSymbol" <- "survival"
    gene_expression_with_labels_symbols <- rbind(labels_df, geneExpressionDF)
    
    patients_data_filtered_symbols_signature <- (gene_expression_with_labels_symbols[gene_expression_with_labels_symbols$"GeneSymbol" %in% c(DCangelosi_gene_symbols_signature, "survival"), ])
    
    rownames(patients_data_filtered_symbols_signature) <- make.unique(patients_data_filtered_symbols_signature$"GeneSymbol")
    patients_data_filtered_symbols_signature$"GeneSymbol" <- NULL
    
    patients_data_filtered_symbols_signature_t <- as.data.frame(t(patients_data_filtered_symbols_signature))

    library("dplyr")
    patients_data_filtered_symbols_signature_t <- patients_data_filtered_symbols_signature_t %>% dplyr::select(-targetName,targetName)
      
    print("dim(patients_data_filtered_symbols_signature)\n")
    print(dim(patients_data_filtered_symbols_signature))
    
    outputFileNameSymbols <- paste0(folderPath, datasetName, "_", GSE_code, "_gene_symbols_dataset_", exe_num, ".csv")
    write.table(patients_data_filtered_symbols_signature_t, file=outputFileNameSymbols, row.names=TRUE, sep=",", col.names=NA)
    cat("saved file  ", outputFileNameSymbols, "\n")

}


