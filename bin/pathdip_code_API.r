setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")

# For the printed files
num_to_return <- 1
exe_num <- sample(1:as.numeric(10000), num_to_return)

list.of.packages <- c("easypackages", "httr") # other packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

library("easypackages")
libraries(list.of.packages)  


url <- "http://ophid.utoronto.ca/pathDIP/Http_API"

# Davide
gene_list <- c("TIMP1", "NDUFB7", "PQBP1", "LSM5", "STMN2", "TPBG", "CEBPB", "AGMAT", "GPR27", "CTGF", "PHLDA1", "PI4K2A", "TM4SF1", "GEMIN6", "HTRA1", "SLC35D1", "DUSP5", "TNFRSF21", "AP002518.1", "CH507-513H4.3", "CH507-513H4.5", "CH507-513H4.6", "HEXIM1", "RP11-197N18.8", "RP5-1028K7.3", "AL031777.2", "CH507-513H4.4", "ANP32C", "RP11-196I18.3", "RPL3P7", "DLEU1_2", "IRS2", "RPL3P4", "HIST1H2BF", "AKAP12", "HSPB8", "MYT1L", "IKBKAP", "GPNMB", "KCNK1", "NDUFA4", "BNIP3L", "CTSB", "FN1", "NDRG1", "GEM", "DUSP6", "PGK1", "PDK1", "HPS5", "ZMPSTE24", "EBNA1BP2", "CYR61", "PXDN", "RRAGD", "COL5A2", "THBS2", "SPARC", "TMEM258", "MLANA", "DYNC1H1", "PPP3CA", "GADD45B", "DUSP10", "FHL2", "ANXA1", "FOSL2", "DCAF7", "CPT2", "PARK7", "LDLR", "ZNF236", "PAM", "DDIT4", "PODXL", "DLEU1", "L1CAM", "TRIM14", "APOE", "TRIM22", "KDSR", "RAB31", "TRIB2", "ANXA2", "BHLHE40", "KLF6", "MB", "CD59", "PLOD1", "VCL", "SPP1", "MCCC2", "DSTN", "PAPSS1", "VCAN", "BGN", "MDM2", "RIPK2", "DNAJC11", "NAP1L1", "STX7", "APIP", "CD55", "NADK", "IGFBP3", "RBPMS", "PLCB4", "VEGFA", "RGS4", "R3HCC1", "PIM1", "ALDOC", "ALDH1A2", "VIM", "LCP1", "MYO1E", "SEC23IP", "PLK2", "PLAUR", "FAM64A", "SAT1", "NNMT", "TMBIM4", "RPL3", "ARHGEF10L", "PIAS1", "SMARCE1", "ATP1B1", "UMPS", "EDN3", "TDRD3", "ARHGEF40", "MAFF", "ALCAM", "MC2R", "ZFAND5", "DNAJA3", "MORC4", "TTC39A", "TMBIM4", "P4HA1", "NOTCH1")

# IDs = paste0(unique(cluster2$Gene.name), collapse=",")
IDs = paste0(unique(gene_list), collapse=",")

component <- "Literature curated (core) pathway memberships"

sources = "ACSN2,BioCarta,EHMN,HumanCyc,INOH,IPAVS,KEGG,NetPath,OntoCancro,Panther_Pathway,PharmGKB,PID,RB-Pathways,REACTOME,stke,systems-biology.org,SignaLink2.0,SIGNOR2.0,SMPDB,Spike,UniProt_Pathways,WikiPathways";

# Define search function
searchOnGenesymbols <- function(IDs, component, sources) {
  
  parameters <- list(
    typeChoice = "Gene Symbol",
    IDs = IDs,
    TableName = component,
    DataSet = sources
  )
  res <- POST(url, body = parameters, encode = "form", verbose())
}

# run search
res <- searchOnGenesymbols(IDs, component, sources)

#define Map function

makeMap <- function(res) {
  
  ENTRY_DEL = "\001"
  KEY_DEL = "\002"
  
  response = httr::content(res, "text")
  
  arr = unlist(strsplit(response, ENTRY_DEL, fixed = TRUE))
  
  list_map <- list("")
  vec_map_names <- c("");
  
  for (str in arr) {
    arrKeyValue = unlist(strsplit(str, KEY_DEL, fixed = TRUE));
    
    if (length(arrKeyValue) > 1) {
      list_map[length(list_map) + 1] <- arrKeyValue[2]
      vec_map_names[length(vec_map_names) + 1] <- arrKeyValue[1]
    }
  }
  names(list_map) <- vec_map_names
  list_map
}

# run map function
responseCode = status_code(res)
if (responseCode != 200) {
  
  cat("Error: Response Code : ", responseCode, "\r\n")
} else {
  
  list_map <- makeMap(res)
}

# Pathway table

cluster2_path = strsplit(list_map$Summary, split = "\r\n")
cluster2_path2 = lapply(cluster2_path, strsplit, split="\t")
cluster2_sign_path = t(as.data.frame(cluster2_path2))
rownames(cluster2_sign_path) <- c()
colnames(cluster2_sign_path) = cluster2_sign_path[1,]
cluster2_sign_path = cluster2_sign_path[-1,]
cluster2_sign_path = as.data.frame(cluster2_sign_path)
cluster2_sign_path[] = lapply(cluster2_sign_path, function(x) if(is.factor(x)) as.character(x)
                              else x)
cluster2_sign_path[,3:5] = lapply(cluster2_sign_path[,3:5], function(x) as.numeric(x))
cluster2_sign_path_b = subset(cluster2_sign_path,`q-value (FDR: BH-method)`<0.01)

write.table(cluster2_sign_path_b, file = "./pathDIP_files/cluster2vs1_SigPathLit001BH.txt", sep = "\t")

# # Gene annotation table 
# 
# cluster2_path_d = strsplit(list_map$Details, split = "\r\n")
# cluster2_path_d2 = lapply(cluster2_path_d, strsplit, split="\t")
# cluster2_path2gene = t(as.data.frame(cluster2_path_d2))
# rownames(cluster2_path2gene) <- c()
# colnames(cluster2_path2gene) = cluster2_path2gene[1,]
# cluster2_path2gene = cluster2_path2gene[-1,]
# cluster2_path2gene = as.data.frame(cluster2_path2gene)
# cluster2_path2gene[] = lapply(cluster2_path2gene, function(x) if(is.factor(x)) as.character(x)
#                               else x)
# # Filter annotation table to keep only significant
# cluster2_path2gene_b = cluster2_path2gene[cluster2_path2gene[,4]%in%cluster2_sign_path_b$`Pathway Source` & cluster2_path2gene[,5]%in%cluster2_sign_path_b$`Pathway Name`,]

write.table(cluster2_sign_path_b, file = "./pathDIP_files/cluster2vs1_SigPathLit001BH.txt", sep = "\t")
# write.table(cluster2_path2gene_b, file = "../pathDIP_files/cluster2vs1_PathMapLit001BH.txt", sep = "\t")
