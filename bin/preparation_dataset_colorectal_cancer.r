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

GSE_code <- "GSE17536"
thisGEOplatform <- "GPL570"
datasetName <-  "Smith2010"


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
    for(thisTitle in gset@phenoData@data$"overall_event") {
      
	if(grepl("no death", thisTitle)) {
	      label_list[[i]] <-  1
	 } else if(grepl("death", thisTitle)) {
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

    VanLaar_probesets_signature <- c("1553954_at", "1554078_s_at", "1555832_s_at", "1555950_a_at", "1560089_at", "1560587_s_at", "1563796_s_at", "200006_at", "200632_s_at", "200665_s_at", "200827_at", "200838_at", "200839_s_at", "200931_s_at", "200983_x_at", "201012_at", "201141_at", "201170_s_at", "201185_at", "201261_x_at", "201289_at", "201323_at", "201422_at", "201426_s_at", "201578_at", "201590_x_at", "201666_at",
    "201925_s_at", "201926_s_at", "201939_at", "201951_at", "202068_s_at", "202237_at", "202238_s_at", "202419_at", "202457_s_at", "202478_at", "202839_s_at", "202887_s_at", "202904_s_at", "202939_at", "202949_s_at", "203072_at", "203083_at", "203382_s_at", "203476_at", "203895_at", "204264_at", "204472_at", "204620_s_at", "204679_at", "205677_s_at", "205963_s_at", "207543_s_at", "207574_s_at", "208891_at", "208892_s_at",
    "208893_s_at", "208918_s_at", "208961_s_at", "209043_at", "209101_at", "209184_s_at", "209185_s_at", "209193_at", "209345_s_at", "209386_at", "209387_s_at", "209457_at", "209545_s_at", "209624_s_at", "209711_at", "209875_s_at", "210095_s_at", "210275_s_at", "210427_x_at", "210495_x_at", "210512_s_at", "210517_s_at", "210592_s_at", "210652_s_at", "210845_s_at", "211074_at", "211719_x_at", "211924_s_at", "211928_at", "211988_at",
    "212013_at", "212143_s_at", "212171_x_at", "212463_at", "212464_s_at", "212501_at", "212632_at", "212884_x_at", "213274_s_at", "213503_x_at", "213905_x_at", "214581_x_at", "214620_x_at", "214866_at", "215033_at", "215034_s_at", "215792_s_at", "216392_s_at", "216442_x_at", "217762_s_at", "217773_s_at", "217996_at", "218213_s_at", "218698_at", "218856_at", "218902_at", "219038_at", "219206_x_at", "219539_at", "221419_s_at",
    "221479_s_at", "221563_at", "221648_s_at", "221656_s_at", "221730_at", "221731_x_at", "221745_at", "222421_at", "222994_at", "223003_at", "223122_s_at", "223163_s_at", "223312_at", "223454_at", "223455_at", "224602_at", "224606_at", "224657_at", "224777_s_at", "224806_at",
    "224890_s_at", "224911_s_at", "225010_at", "225011_at", "225337_at", "225494_at", "225670_at", "225750_at", "226041_at", "226594_at", "226648_at", "226727_at", "226987_at", "227143_s_at", "227338_at", "227735_s_at", "227736_at", "227961_at", "229676_at", "231576_at", "234983_at", "241355_at", "242648_at", "35156_at", "36711_at", "58780_s_at")

    
    patients_data_filtered_ourSignature <- (gene_expression_with_labels[gene_expression_with_labels$ID %in% c(VanLaar_probesets_signature, "survival"), ])
      
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
    
    VanLaar_gene_symbols_signature <- c("ALG14", "DNAJA3", "KLF6", "CD55", "LOC286208", "PRDX5", "EARS2", "PARK7", "NDRG1", "SPARC", "PLOD1", "CTSB", "CTSB", "VCL", "CD59", "ANXA1", "GPNMB", "BHLHE40", "HTRA1", "BGN", "CYR61", "EBNA1BP2", "IFI30", "VIM", "PODXL", "ANXA2", "TIMP1", "CD55", "CD55", "PLK2", "ALCAM", "LDLR", "NNMT", "NNMT", "KDSR", "PPP3CA", "TRIB2", "NDUFB7", "DDIT4", "LSM5", "ZMPSTE24", "FHL2", "MYO1E", "THBS2", "APOE", "TPBG", "PLCB4", "CPT2", "GEM", "VCAN", "KCNK1", "DLEU1", "DNAJA3", "P4HA1", "GADD45B", "DUSP6", "DUSP6", "DUSP6", "NADK", "KLF6", "PAPSS1", "CTGF", "IRS2", "IRS2", "PIM1", "PI4K2A", "TM4SF1", "TM4SF1", "DUSP5", "RIPK2", "MCCC2", "SLC35D1", "SPP1", "IGFBP3", "ZFAND5", "ANXA2", "FN1", "VEGFA", "AKAP12", "SAT1", "TTC39A", "PLAUR", "FOLR1", "FN1", "PLAUR", "DYNC1H1", "SMARCE1", "PXDN", "IGFBP3", "VEGFA", "CD59", "FN1", "CEBPB", "STX7", "APOE", "CTSB", "ANXA2", "BGN", "TNFRSF21", "PAM", "PLAUR", "TM4SF1", "TM4SF1", "DNAJC11", "SEC23IP", "FN1", "RAB31", "NDUFA4", "PHLDA1", "C11orf10", "APIP", "TNFRSF21", "NOTCH1", "MORC4", "TMBIM4", "GEMIN6", "BNIP3L", "DUSP10", "ARHGEF10L", "COL5A2", "VCAN", "DCAF7", "UBE2H", "PRDX5", "C19orf43", "SFRP2", "ZC3HC1", "C2orf7", "CXCL16", "TCHP", "C4orf3", "KLF6", "ERRFI1", "PAFAH1B2", "TRIM25", "C7orf59", "DCBLD2", "CCDC6", "PRKAR2A", "ABHD2", "DYNLL2", "FAM173B", "NAPEPLD", "HIF1AN", "CISD3", "RBM15B", "BID", "LOC440983", "C10orf99", "C10orf99", "CTSB", "MTPAP", "HR", "KLHL8", "R3HCC1", "MAFF", "FLJ10357")


    library("geneExpressionFromGEO")
    associateSymbolsToGenes <- TRUE
    verbose <- TRUE
    geneExpressionDF <- getGeneExpressionFromGEO(GSE_code,  associateSymbolsToGenes, verbose)
    
    labels_df$"GeneSymbol" <- "survival"
    gene_expression_with_labels_symbols <- rbind(labels_df, geneExpressionDF)
    
    patients_data_filtered_symbols_signature <- (gene_expression_with_labels_symbols[gene_expression_with_labels_symbols$"GeneSymbol" %in% c(VanLaar_gene_symbols_signature, "survival"), ])
    
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


