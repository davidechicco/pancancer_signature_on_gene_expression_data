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

outputFileName <- ""

GSE_code <- "GSE12417" # this line will change for each different dataset
thisGEOplatform <- "GPL96" # this line will change for each different dataset
datasetName <-  "Metzeler2018" # this line will change for each different dataset
cancer_type <- "leukemia" # this line will change for each different dataset

cat("\n\tGSE_code: ", GSE_code, "\n", sep="")
cat("\tthisGEOplatform: ", thisGEOplatform, "\n", sep="")
cat("\tdatasetName: ", datasetName, "\n", sep="")
cat("\tcancer_type: ", cancer_type, "\n", sep="")

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
    for(thisTitle in gset@phenoData@data$"characteristics_ch1") { # this line will change for each different dataset  
      
	if(grepl("(0=alive/1=dead): 1", thisTitle, fixed=TRUE)) { # this line will change for each different dataset
	      label_list[[i]] <-  0
          } else if(grepl("(0=alive/1=dead): 0", thisTitle, fixed=TRUE)) { # this line will change for each different dataset
	    label_list[[i]] <- 1
	  } else {
             label_list[[i]] <- -1
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

    source("cancer_prognosis_probesets_signatures.r")
    signature_name <- "our_pancancer_signature"
    our_merged_pancancer_prognostic_signature <- unique(our_merged_pancancer_prognostic_signature_original)
    
    patients_data_filtered_ourSignature <- (gene_expression_with_labels[gene_expression_with_labels$ID %in% c(our_merged_pancancer_prognostic_signature, "survival"), ])
      
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

    folderPath <- "../data/"
    mkdirCommand <- paste0("mkdir -p ", folderPath)
    system(mkdirCommand)
    cat("command invoked: ", mkdirCommand, "\n")
    
    outputFileName <- paste0(folderPath, datasetName, "_", GSE_code, "_", cancer_type,  "_", signature_name, "_probesets_dataset_", exe_num, ".csv")
    write.table(tableOnlyOnesAndZeros, file=outputFileName, row.names=FALSE, sep=",")
    cat("saved file  ", outputFileName, "\n")
    
}

