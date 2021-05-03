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

GSE_code <- "GSE33118"
thisGEOplatform <- "GPL570"
datasetName <-  "Bilbault2017"


gset <- getGEO(GSE_code,  GSEMatrix =TRUE, getGPL=FALSE)

if (length(gset) > 1) idx <- grep(thisGEOplatform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
gene_expression <- as.data.frame(exprs(gset))

cat("str(gset@phenoData@data)\n")
print(str(gset@phenoData@data))
cat("str(gset@phenoData@data)\n")

LABEL_DETECTED <- TRUE
 
 if(LABEL_DETECTED == TRUE) {
 
    # # # # we add the labels
    library("plyr")
    label_list <- c()
    i <- 1
    for(thisTitle in gset@phenoData@data$"survival:ch1") {
      
	if(grepl("positive", thisTitle)) {
	      label_list[[i]] <-  1
	 } else if(grepl("negative", thisTitle)) {
	    label_list[[i]] <- 0
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
    
    gene_expression_with_labels$"ID" <- rownames(gene_expression_with_labels)

    ## selection based on signature probesets
    TingWang_probesets_signature <- c("1552263_at", "209857_s_at", "206437_at", "212046_x_at", "207163_s_at", "200602_at", "202615_at", "223454_at", "209906_at", "205118_at", "204533_at", "201040_at", "207008_at", "210772_at", "1567457_at", "203146_s_at", "206722_s_at", "1555814_a_at")

    
    patients_data_filtered_ourSignature <- (gene_expression_with_labels[gene_expression_with_labels$ID %in% c(TingWang_probesets_signature, "survival"), ])
      
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
    
    TingWang_gene_symbols_signature <- c("MAPK1", "SPHK2", "S1PR4", "MAPK3", "AKT1", "APP", "GNAQ", "CXCL16", "C3AR1", "FPR1", "CXCL10", "GNAI2", "CXCR2", "FPR2", "RAC1", "GABBR1", "LPAR2", "RHOA")

    library("geneExpressionFromGEO")
    associateSymbolsToGenes <- TRUE
    verbose <- TRUE
    geneExpressionDF <- getGeneExpressionFromGEO(GSE_code,  associateSymbolsToGenes, verbose)
    
    labels_df$"GeneSymbol" <- "survival"
    gene_expression_with_labels_symbols <- rbind(labels_df, geneExpressionDF)
    
    patients_data_filtered_symbols_signature <- (gene_expression_with_labels_symbols[gene_expression_with_labels_symbols$"GeneSymbol" %in% c(TingWang_gene_symbols_signature, "survival"), ])
    
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


