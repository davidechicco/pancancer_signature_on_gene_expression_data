setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")

# For the printed files
num_to_return <- 1
exe_num <- sample(1:as.numeric(10000), num_to_return)

list.of.packages <- c("easypackages", "plyr", "gprofiler2") # other packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
listOfBiocPackages <- c("annotate", "GEOquery", "biomaRt")

bioCpackagesNotInstalled <- which( !listOfBiocPackages %in% rownames(installed.packages()) )
cat("package missing listOfBiocPackages[", bioCpackagesNotInstalled, "]: ", listOfBiocPackages[bioCpackagesNotInstalled], "\n", sep="")

# check there's still something left to install
if( length(bioCpackagesNotInstalled) ) {
    BiocManager::install(listOfBiocPackages[bioCpackagesNotInstalled])
}

library("easypackages")
libraries(listOfBiocPackages)  

# Function that retrieves the annotations of a probeset list.
# Parameters: 
#  probeset_array  vector of affy_hg_u133a probeset ids
#  flank_len       gene flank length (default=25)
#  flank_stream    "downstream" (default) or "upstream"
#  flank_type      "gene" (default), "coding_gene", "transcript", or "coding_transcript"
retrieveAnnotations_GPL96 <- function(probeset_array, 
                                      flank_len = 25, 
                                      flank_stream = "downstream", 
                                      flank_type = "gene") {
  
  stopifnot(flank_stream %in% c("downstream", "upstream"))
  stopifnot(flank_type %in% c("gene", "coding_gene", "transcript", "coding_transcript"))
  flank_stream <- paste0(flank_stream, "_flank")
  flank_type <- paste0(flank_type, "_flank")
  
  cat("Attention: the function retrieveAnnotations_GPL96() works only for GPL96 datasets\n")
  
  # Gene list
  currSpecieMart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  dataset = 'hsapiens_gene_ensembl', 
                  host = "mar2016.archive.ensembl.org")
  
  thisAnnotLookup <- 
    getBM(mart = currSpecieMart, 
          attributes = c("affy_hg_u133a", "ensembl_gene_id", "gene_biotype", "external_gene_name", "chromosome_name", "start_position", "end_position", flank_type), 
          filter = c("affy_hg_u133a", flank_stream), 
          values = list(probeset_array, flank_len), 
          uniqueRows = TRUE, 
          checkFilters = FALSE, bmHeader = TRUE)
          
          
  
  return(thisAnnotLookup)
}
# 
# thisGEOplatform <- "GPL96"
# platform_ann <- readGEOAnn(GEOAccNum = thisGEOplatform)
# platform_ann_df <- as.data.frame(platform_ann, stringsAsFactors=FALSE)
# 
# allProbesetsIDs <- (platform_ann_df$"ID")

# MCM5_signature_probesets <- c("200785_s_at", "201856_s_at", "201926_s_at", "202085_at", "202170_s_at", "202922_at", "204511_at", "205364_at", "206472_s_at", "207730_x_at", "208664_s_at", "208892_s_at", "210404_x_at", "211741_x_at", "212110_at", "213778_x_at", "214261_s_at", "214730_s_at", "215066_at", "216735_x_at", "218380_at", "220389_at", "221610_s_at", "222086_s_at", "222168_at")

source("cancer_prognosis_probesets_signatures.r")
signature_name <- "our_pancancer_signature"
our_merged_pancancer_prognostic_signature <- unique(our_merged_pancancer_prognostic_signature_original)
 
annotations <- retrieveAnnotations_GPL96(our_merged_pancancer_prognostic_signature)

cat("annotations:\n")
print(annotations[,c("Affy HG U133A probeset", "Associated Gene Name")])

for(i in 1:nrow(annotations)) {

     this_gene_symbol <- annotations[i,6]
     cat(this_gene_symbol, "\", \"", sep="")
}
cat("\n")

# entrezGeneIDlist <- gconvert(annotations$"Ensembl Gene ID",  organism = "hsapiens", target="ENTREZGENE_ACC")

if(nrow(annotations) < length(our_merged_pancancer_prognostic_signature)) {
    cat("Attention: only ", nrow(annotations), " probesets out of ", length(our_merged_pancancer_prognostic_signature), " have annotations\n", sep="")
    probesWithoutAnnotations <-  length(our_merged_pancancer_prognostic_signature) - nrow(annotations)
    cat("and ",probesWithoutAnnotations, " probesets do not have annotations\n", sep="")
}


SAVE_FILE <- FALSE

if(SAVE_FILE) {

    outputFile <- paste0("../results/GPL96_probesets_sequences_symbols_MCM5_rand", exe_num,".csv")
    write.csv(annotations, file=outputFile, row.names=FALSE)
    cat("saved file ", outputFile, "\n")

    if(nrow(annotations) < length(MCM5_signature_probesets)) {
        cat("Attention: only ", nrow(annotations), " probes out of ", length(MCM5_signature_probesets), " have annotations\n", sep="")
        probesWithoutAnnotations <-  length(MCM5_signature_probesets) - nrow(annotations)
        cat("and ",probesWithoutAnnotations, " probes do not have annotations\n", sep="")
    }

}
