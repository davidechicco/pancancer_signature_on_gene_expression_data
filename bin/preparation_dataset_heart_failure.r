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

GSE_code <- "GSE76701"
thisGEOplatform <- "GPL570"
datasetName <-  "Kim2016"


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
    for(thisTitle in gset@phenoData@data$"source_name_ch1") {
      
	if(grepl("Non-failing adult human heart", thisTitle)) {
	      label_list[[i]] <-  1
	 } else if(grepl("Failing adult human heart", thisTitle)) {
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

    DavidLin_probesets_signature <- c("227228_s_at", "244218_at", "237001_at", "228512_at", "226764_at", "235036_at", "204021_s_at", "222876_s_at", "211783_s_at", "214149_s_at", "227682_at", "211300_s_at", "1564154_at", "1557487_at", "235380_at", "215743_at", "229695_at", "227900_at", "243341_at", "227410_at", "233341_s_at", "242197_x_at", "228026_at", "241343_at", "220985_s_at", "205456_at", "205976_at", "206420_at", "211900_x_at", "221860_at", "207521_s_at", "216246_at", "1564155_x_at", "205387_s_at", "205005_s_at", "244828_x_at", "208900_s_at", "215611_at", "218280_x_at", "219151_s_at", "1558827_a_at", "216609_at", "238439_at", "202905_x_at", "209972_s_at", "225793_at", "221011_s_at", "211893_x_at", "223502_s_at", "243810_at", "221058_s_at", "209092_s_at", "212599_at", "239673_at", "209871_s_at", "235543_at", "238121_at", "225746_at", "232262_at", "228520_s_at", "229235_at", "228466_at", "239346_at", "222820_at", "242051_at", "243982_at", "214163_at", "224518_s_at", "226158_at", "1560163_at", "225688_s_at", "207351_s_at", "242622_x_at", "235434_at", "1558809_s_at", "1558996_at", "213986_s_at", "224698_at", "", "238058_at", "234734_s_at", "220500_s_at", "235507_at", "213906_at", "218338_at", "221768_at", "37793_r_at", "243026_x_at", "206513_at", "230526_at", "229544_at", "204747_at", "224516_s_at", "218362_s_at", "1555830_s_at", "214316_x_at", "237346_at", "1555960_at", "243916_x_at", "224812_at", "206099_at")

    patients_data_filtered_ourSignature <- (gene_expression_with_labels[gene_expression_with_labels$ID %in% c(DavidLin_probesets_signature, "survival"), ])
      
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
    
    DavidLin_gene_symbols_signature <- c("CCDC88C", "PTCD3", "ZNF827", "LIX1L", "PURA", "ADAP2", "MTA1", "ATP6V0E1",  "TP53",    "NMT2",  "CBLB",  "FAM43A", "POLR1B", "CD36", "RP5- 1000E10.4", "RNASEH1", "RNF170", "CD3E", "FASTKD2", "IGSF6", "CD6", "HNRNPL", "ATP2A3",   "CGB", "CGB5”, “CGB7", "CGB8", "NMT2", "NAF1", "TOP1", "TCF12", "HIST2H2AA3", "HIST2H2AA4", "RABL2A", "RABL2B", "ZNF831", "TXN", "ANKRD22", "NBN", "JTV1", "LIX1L", "LBH", "CD6", "TNFSF13B",  "CKLF", "GLOD4", "AUTS2",  "APBA2",  "GK5", "RAB11FIP4", "PIGL", "APLP2",    "TNRC6C",  "KLHL28", "HSPB11", "ZNF559", "KLHL24", "LOC1001295 10", "PHLDB2", "SH2D2A", "PTEN",  "LOC284408", "FOXP1", "C19orf6", "FAM62B", "FLJ27365", "LOC150381", "TNRC6A", "RABL2A", "RABL2B", "PCMTD1", "MYBL1", "PHC1", "PHC1B", "SFPQ", "RAD51L3", "XIAP", "AIM2", "LOC100131096",  "IFIT3", "CXXC5", "DIS3", "FAM62B",  "TGDS", "HINT1", "UBLCP1", "HIBADH", "PRKCH")

    library("geneExpressionFromGEO")
    associateSymbolsToGenes <- TRUE
    verbose <- TRUE
    geneExpressionDF <- getGeneExpressionFromGEO(GSE_code,  associateSymbolsToGenes, verbose)
    
    labels_df$"GeneSymbol" <- "survival"
    gene_expression_with_labels_symbols <- rbind(labels_df, geneExpressionDF)
    
    patients_data_filtered_symbols_signature <- (gene_expression_with_labels_symbols[gene_expression_with_labels_symbols$"GeneSymbol" %in% c(DavidLin_gene_symbols_signature, "survival"), ])
    
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


