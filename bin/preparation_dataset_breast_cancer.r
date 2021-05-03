setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
# set.seed(11)

list.of.packages <- c("easypackages", "plyr") # other packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="https://utstat.toronto.edu/cran/")

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

GSE_code <- "GSE15852"
thisGEOplatform <- "GPL570"
datasetName <-  "Bong2010"


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
    for(thisTitle in gset@phenoData@data$"histopathological exam:ch1") {
      
	if(grepl("carcinoma", thisTitle)) {
	      label_list[[i]] <-  1
	 } else if(grepl("normal", thisTitle)) {
	    label_list[[i]] <- 0
	  }
	    i <- i + 1
        }
      
    cat("label_list:\n")
    print(label_list)
    
    targetName <- "cancer"
    
    labels_df_temp <- as.data.frame(label_list)
    
    labels_df <- as.data.frame(t(labels_df_temp))
    colnames(labels_df) <- colnames(gene_expression)
    rownames(labels_df) <- targetName
    gene_expression_with_labels <- rbind(labels_df, gene_expression)
    
    gene_expression_with_labels$"ID" <- rownames(gene_expression_with_labels)

    ## selection based on signature probesets

    
    Wen_probeset_signature <- c("212741_at", "213524_s_at", "219064_at", "210375_at", "219140_s_at", "210517_s_at", "221748_s_at", "210964_s_at", "201348_at", "207175_at", "209493_at", "200606_at", "205478_at", "203549_s_at", "216594_x_at", "219398_at", "221009_s_at", "222317_at", "218124_at",
"211696_x_at", "219295_s_at", "205220_at", "201839_s_at", "211653_x_at", "215695_s_at", "209116_x_at", "202988_s_at", "213071_at", "205610_at", "205382_s_at", "204151_x_at", "204388_s_at", "203323_at", "213068_at", "209612_s_at", "201540_at", "209763_at", "218237_s_at", "218087_s_at",
"266_s_at", "201650_at", "207092_at", "203548_s_at", "217232_x_at", "201690_s_at", "206243_at", "209699_x_at", "210832_x_at", "201596_x_at",
"203980_at", "204894_s_at", "205913_at", "204570_at", "221246_x_at", "208510_s_at", "216331_at", "217414_x_at", "43427_at", "219738_s_at", "210201_x_at", "212071_s_at", "209581_at", "211745_x_at", "204159_at", "207275_s_at", "213706_at", "208383_s_at", "209614_at", "214456_x_at", "209160_at", "203571_s_at", "206488_s_at", "216379_x_at", "210298_x_at", "219737_s_at", "214414_x_at", "209771_x_at", "218168_s_at",
"213693_s_at", "221928_at", "49452_at", "201131_s_at", "204134_at", "209298_s_at", "210715_s_at", "209555_s_at", "202931_x_at", "212281_s_at", "202942_at", "209496_at", "204997_at", "205824_at", "203851_at", "214093_s_at", "217028_at", "203065_s_at", "202286_s_at", "218186_at", "204719_at", "204653_at", "214063_s_at", "206068_s_at", "212636_at", "212886_at", "221747_at", "209458_x_at", "209183_s_at", "218677_at",
"208651_x_at", "201243_s_at", "209210_s_at", "209686_at", "220061_at", "214505_s_at", "209008_x_at", "213764_s_at", "201286_at", "201719_s_at", "203324_s_at", "203337_x_at", "212122_at", "209602_s_at", "210299_s_at", "212112_s_at", "211699_x_at", "213605_s_at", "204825_at", "219298_at", "212282_at", "204223_at", "222351_at", "209616_s_at", "213933_at", "218736_s_at", "58780_s_at", "213451_x_at", "212875_s_at", "209512_at",
"202283_at", "203215_s_at", "210740_s_at", "204018_x_at", "212653_s_at", "201417_at", "208789_at", "218723_s_at", "201450_s_at", "209138_x_at", "206093_x_at", "203336_s_at", "209016_s_at", "221295_at", "215039_at", "205932_s_at", "203799_at", "211946_s_at", "205609_at", "206030_at", "201689_s_at", "218175_at", "201539_s_at", "203296_s_at", "37117_at", "201664_at", "210869_s_at", "208029_s_at", "218031_s_at", "202768_at",
"212097_at", "204848_x_at", "206100_at", "219117_s_at", "200878_at", "203986_at", "213765_at", "205022_s_at", "208650_s_at", "222043_at", "216834_at", "208607_s_at", "204731_at", "219093_at", "201890_at", "205968_at", "208737_at", "209540_at", "214439_x_at", "203917_at", "201448_at", "215379_x_at", "216894_x_at", "215946_x_at", "212230_at", "210078_s_at", "209613_s_at", "202193_at", "202088_at", "202963_at",
"218625_at", "209513_s_at", "201963_at", "206385_s_at", "220327_at", "202112_at", "203680_at", "213913_s_at", "213518_at", "204326_x_at", "203693_s_at", "213364_s_at", "201761_at", "210046_s_at", "213150_at", "217771_at", "204798_at", "203397_s_at", "201202_at", "214677_x_at", "203358_s_at", "201242_s_at", "213793_s_at", "201340_s_at", "219288_at", "203400_s_at", "212036_s_at", "221024_s_at", "206201_s_at",
"201425_at", "213236_at", "211922_s_at", "202982_s_at", "201952_at")


    
    patients_data_filtered_ourSignature <- (gene_expression_with_labels[gene_expression_with_labels$ID %in% c(Wen_probeset_signature, "cancer"), ])
      
    print("dim(patients_data_filtered_ourSignature)\n")
    print(dim(patients_data_filtered_ourSignature))

    # patients_data_filtered_ourSignature <- rbind(gene_expression_with_labels[1,], patients_data_filtered_ourSignature)

    patients_data_filtered_ourSignature$"ID" <- NULL
    
    patients_data_filtered_ourSignature_t <- as.data.frame(t(patients_data_filtered_ourSignature))

    library("dplyr")
    patients_data_filtered_ourSignature_t <- patients_data_filtered_ourSignature_t %>% dplyr::select(-targetName,targetName)

    patients_data_filtered_ourSignature_t <-  slice(patients_data_filtered_ourSignature_t, 1:(n()-1))    
    patients_data_filtered_ourSignature_t[,targetName] <- as.numeric(patients_data_filtered_ourSignature_t[,targetName]) 

    tableOnlyOnesAndZeros <-  rbind(patients_data_filtered_ourSignature_t[patients_data_filtered_ourSignature_t$"cancer"==0,], patients_data_filtered_ourSignature_t[patients_data_filtered_ourSignature_t$"cancer"==1,])

    folderPath <- "../results/"
    
    outputFileName <- paste0(folderPath, datasetName, "_", GSE_code, "_probesets_dataset_", exe_num, ".csv")
    write.table(tableOnlyOnesAndZeros, file=outputFileName, row.names=TRUE, sep=",", col.names=NA)
    cat("saved file  ", outputFileName, "\n")
    
    
    ## selection based on signature gene symbols
    
    Wen_gene_symbols_signature <- c("MAOA", "G0S2", "ITIH5", "PTGER3", "RBP4", "AKAP12", "TNS1", "GYG2", "GPX3", "ADIPOQ", "PDZD2", "DSP", "PPP1R1A", "LPL", "AKR1C1", "CIDEC", "ANGPTL4", "PDE3B", "RETSAT", "HBB", "PCOLCE2", "HCAR3", "EPCAM", "AKR1C2", "LOC101060798", "GYG2", "HBB", "RGS1", "DPT", "MYOM1", "CFD", "AKR1C1", "MAOA", "CAV2", "DPT", "ADH1B", "FHL1", "CHRDL1", "SLC38A1", "SORBS1", "CD24", "KRT19", "LEP", "LPL", "HBB", "TPD52", "TIMP4", "AKR1C2", "LOC101060798", "PTGER3", "KRT18", "FABP4", "AOC3", "PLIN1", "COX7A1", "TNS1", "PPARG", "ITGA7", "HBA1", "HBA2", "ACACB", "PCDH9", "BIN1", "SPTBN1", "PLA2G16", "HBA1", "HBA2", "CDKN2C", "ACSL1", "GPD1", "PCK1", "ADH1B", "SAA1", "SAA2", "AKR1C3", "ADIRF", "CD36", "CD24", "FHL1", "PCDH9", "HBA1", "HBA2", "CD24", "ADCK3", "MUC1", "ACACB", "ACACB", "CDH1", "PDE2A", "ITSN1", "SPINT2", "CD36", "BIN1", "TMEM97", "ETFB", "RARRES2", "GPD1", "HSPB2", "IGFBP6", "FUBP1", "CXCR4", "CAV1", "TACSTD2", "RAB25", "ABCA8", "TFAP2A", "TF", "ACADL", "QKI", "CCDC69", "TNS1", "HBA1", "HBA2", "C10orf10", "S100A14", "CD24", "ATP1B1", "FERMT2", "S100B", "ACSM5", "FHL1", "KRT8", "MFAP5", "SDC1", "EPB41L2", "CAV2", "ITGB1BP1", "RHOQ", "GATA3", "FHL1", "STX12", "HBA1", "HBA2", "MELK", "ECHDC3", "TMEM97", "PRELP", "PPP2R1B", "CES1", "LOC100653057", "PTGER3", "PALMD", "ARHGEF40", "LOC101060681", "TNXA", "TNXB", "C2CD2", "HSDL2", "SERPINF1", "MYO6", "ITPK1", "HBA1", "HBA2", "EHBP1", "SOX4", "PTRF", "RGCC", "TIA1", "IGLC1", "LOC101060681", "TNXA", "TNXB", "ITGB1BP1", "KRT7", "CIDEA", "LOC339524", "MSX1", "CD302", "LY75-CD302", "PRRC2C", "ANGPT1", "ASPA", "TPD52", "CCDC92", "FHL1", "ATP1A2", "ARHGAP8", "PRR5-ARHGAP8", "SMC4", "MCAM", "LAPTM4B", "FOXN3", "FOSB", "CAV1", "HBG1", "HBG2", "CPM", "FKBP11", "EPAS1", "FAM47E", "FAM47E-STBD1", "STBD1", "MFAP5", "FOXN3", "CD24", "CLU", "RGS1", "SAA1", "SAA2", "SAA2-SAA4", "TGFBR3", "PID1", "RRM2", "KCNS3", "ATP6V1G1", "IGF1", "BIN1", "CXADR", "TIA1", "IGLV1-44", "CDKN1C", "IGLL3P", "PPAP2B", "KCNAB1", "ADH1B", "LIMK2", "SLC39A6", "RFX5", "NRN1", "HSDL2", "ACSL1", "ANK3", "VGLL3", "VWF", "PRKAR2B", "TBC1D30", "PRKCI", "MT1X", "E2F3", "SNX1", "MTHFD2", "IDH2", "HOXA10", "GOLM1", "MYB", "GALNT3", "PCNA", "IGLC1", "EZH2", "ATP1B1", "HOMER1", "ENC1", "C3orf14", "TF", "PNN", "SLC2A10", "MEOX2", "ALDH2", "SASH1", "CAT", "ACOT1", "ACOT2", "ALCAM")



    library("geneExpressionFromGEO")
    associateSymbolsToGenes <- TRUE
    verbose <- TRUE
    geneExpressionDF <- getGeneExpressionFromGEO(GSE_code,  associateSymbolsToGenes, verbose)
    
    labels_df$"GeneSymbol" <- "cancer"
    gene_expression_with_labels_symbols <- rbind(labels_df, geneExpressionDF)
    
    patients_data_filtered_symbols_signature <- (gene_expression_with_labels_symbols[gene_expression_with_labels_symbols$"GeneSymbol" %in% c(Wen_gene_symbols_signature, "cancer"), ])
    
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


