setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)

execution_number <- 100

EXP_ARG_NUM <- 2



# fileName <- "/home/dave/my_projects/probesets_versus_gene_symbols/data/Hiyama2010_GSE16237_neuroblastoma_our_pancancer_signature_probesets_dataset_188.csv"

# fileName <- "/home/dave/my_projects/probesets_versus_gene_symbols/data/Pasini2011_GSE38749_stomach cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Shi2010_GSE24080_multiple_myeloma_our_pancancer_signature_probesets_dataset_1786.csv"

# source("preparation_dataset_multiple_myeloma_our_signature.r")
# fileName <- outputFileName

# fileName <- "../data/Shi2010_GSE24080_multiple_myeloma_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Chen2020_GSE161158_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Reister2012_GSE31684_bladder_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Heaton2011_GSE33371_adrenocortical_carcinoma_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Hatzis2009_GSE25066_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Iqbal2015_GSE58445_non_hodgkin_lymphoma_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Uehara2015_GSE65986_ovarian_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Leich2009_GSE16131_non_hodgkin_lymphoma_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/ZChen2020_GSE157011_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Sieber2010_GSE14333_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Marisa2013_GSE39582_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Shinto2020_GSE143985_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Gotoh2018_GSE92921_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Pintilie2013_GSE50081_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Ivshina2006_GSE4922_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Staub2009_GSE12945_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Heiskanen2015_GSE68465_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Rousseaux2013_GSE30219_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Sabatier2010_GSE21653_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Schmidt2008_GSE11121_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Desmedt2007_GSE7390_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Sinn2009_GSE124647_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Miller2013_GSE45255_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/MetzgerFilho2016_GSE88770_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Dedeurwaerder2011_GSE20711_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Huang2014_GSE48390_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Jezequel2015_GSE58812_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Karn2011_GSE31519_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Kim2020_GSE135565_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Lin2009_GSE19697_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Loi2008_GSE9195_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Wang2010_GSE19615_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Yenamandra2015_GSE61304_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Smith2009a_GSE17536_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Smith2009b_GSE17537_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/DelRoi2017_GSE72970_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Bild2005_GSE3141_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Botling2012_GSE37745_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Kohno2011_GSE31210_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Micke2011_GSE28571_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Philipsen2010_GSE19188_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Beauchamp2014_GSE38832_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <-"../data/Potti2006_GSE3593_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Son2007_GSE8894_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Tsao2010_GSE14814_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Xie2011_GSE29013_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Kawaguchi2012_GSE34771_non_Hodgkin_lymphoma_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Spivak2014_GSE47018_leukemia_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Herold2013_GSE37642_leukemia_our_pancancer_signature_probesets_dataset_1786.csv"

# fileName <- "../data/Herold2011_GSE22762_leukemia_our_pancancer_signature_probesets_dataset_1786.csv"

fileName <- "../data/Metzeler2018_GSE12417_leukemia_our_pancancer_signature_probesets_dataset_1786.csv"

TRAIN_SET_OVERSAMPLING_SYNTHETIC <- FALSE

targetName <- "survival"


cat("fileName: ", fileName, "\n", sep="")
cat("targetName: ", targetName, "\n", sep="")

list.of.packages <- c("easypackages", "PRROC", "e1071", "randomForest","class", "gmodels", "formula.tools", "dplyr", "pastecs", "ROSE")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="https://utstat.toronto.edu/cran/")

library("easypackages")
libraries(list.of.packages)


source("./confusion_matrix_rates.r")
source("./utils.r")

NUM_METRICS <- 9
confMatDataFrame <- matrix(ncol=NUM_METRICS, nrow=1)
colnames(confMatDataFrame) <- c("MCC", "F1 score", "accuracy", "TP_rate", "TN_rate", "PPV", "NPV", "PR_AUC", "ROC_AUC")

threshold <- 0.5

patients_data <- read.csv(file=fileName,head=TRUE,sep=",",stringsAsFactors=FALSE)
cat("fileName = ", fileName, "\n", sep="")

patients_data$"sampleID" <- NULL


# let's put the target label last on the right 
patients_data <- patients_data%>%select(-targetName,targetName)
target_index <- dim(patients_data)[2]    
patients_data_original <- patients_data

cat("original dataset:\n")
dataset_dim_retriever(patients_data_original)
imbalance_retriever(patients_data_original[,target_index])

patients_data[,target_index] <- as.factor(patients_data[,target_index])

# formula
allFeaturesFormula <- as.formula(paste(as.factor(colnames(patients_data)[target_index]), '.', sep=' ~ ' ))

# cycle of executions


cat("Number of executions = ", execution_number, "\n", sep="")
for(exe_i in 1:execution_number)
{
    cat("[Execution number ", exe_i, " out of ", execution_number, "]\n", sep="" )
    cat("[Randomizing the rows]\n")
    patients_data <- patients_data[sample(nrow(patients_data)),] # shuffle the rows

    totalElements <- dim(patients_data)[1]

    subsets_size <- totalElements

    target_label <- colnames(patients_data[target_index])
    cat("target_label = ", target_label, "\n", sep="")

    if (subsets_size != totalElements) {
        cat("ATTENTION: We are running the method on a subset of the original dataset, \n", sep="")
        cat(" containing only ", subsets_size, " elements \n", sep="")
        cat(" instead of ", totalElements, " elements \n", sep="")
    }

    patients_data <- patients_data[1:subsets_size, ]

    dataset_dim_retriever(patients_data)
    imbalance_retriever(patients_data[,target_index])

    training_set_perc <- 80
    INPUT_PERC_POS <- 50
    cat("[training set = ", training_set_perc,"%]\n", sep="")
    cat("[test set = ", (100-training_set_perc),"%]\n", sep="")

    artificialBalance <- FALSE
    balancedFlag <- FALSE # flag that sets everything to 50% 50% ratio

    if (artificialBalance == TRUE) {


        train_data_balancer_output <- train_data_balancer(patients_data, target_index, training_set_perc, INPUT_PERC_POS, balancedFlag)

        patients_data_train <- train_data_balancer_output[[1]]
        patients_data_test <- train_data_balancer_output[[2]]
        
        # Creating the subsets for the targets
        patients_data_train_labels <- patients_data_train[, target_index] # NEW
        patients_data_test_labels <- patients_data_test[, target_index]   # NEW

    } else {


        # the training set is the first 60% of the whole dataset
        training_set_first_index <- 1 # NEW
        training_set_last_index <- round(dim(patients_data)[1]*training_set_perc/100) # NEW

        # the test set is the last 40% of the whole dataset
        test_set_first_index <- training_set_last_index+1 # NEW
        test_set_last_index <- dim(patients_data)[1] # NEW

        cat("[Creating the training set and test set for the values]\n")
        patients_data_train <- patients_data[training_set_first_index:training_set_last_index, 1:(target_index)] # NEW
        patients_data_test <- patients_data[test_set_first_index:test_set_last_index, 1:(target_index)] # NEW
        
        # train_balanced_over <- ovun.sample(allFeaturesFormula, data = patients_data_train, method = "both",  p=0.5, N = (nrow(patients_data_train)*2))$data
        # patients_data_train <- train_balanced_over
         
         # https://www.analyticsvidhya.com/blog/2016/03/practical-guide-deal-imbalanced-classification-problems/

         
         if(TRAIN_SET_OVERSAMPLING_SYNTHETIC == TRUE)
         {
            thisP <- 0.5
         
            data.rose <- ROSE(allFeaturesFormula, data = patients_data_train, p=thisP, seed = 1)$data
            patients_data_train <- data.rose
         }
        
        cat("[training set dimensions: ", dim(patients_data_train)[1], " patients]\n")

        cat("[test set dimensions: ", dim(patients_data_test)[1], " patients]\n")

        cat("[Creating the training set and test set for the labels \"1\"-\"0\"]\n")
        patients_data_train_labels <- patients_data_train[, target_index] # NEW
        patients_data_test_labels <- patients_data[test_set_first_index:test_set_last_index, target_index]   # NEW

    }


    dataset_dim_retriever(patients_data_train)
    imbalance_retriever(patients_data_train[, targetName])    

    cat("\n[Training the random forest classifier on the training set]\n")

    rf_new <- NULL
    rf_new <- randomForest(allFeaturesFormula, data=patients_data_train, importance=FALSE, proximity=TRUE)
    
    cat("\n[Applying the trained random forest classifier on the test set]\n")
    patients_data_test_PRED <- predict(rf_new, patients_data_test, type="prob")[,"1"]

    thisConfMat <- confusion_matrix_rates(patients_data_test_labels, patients_data_test_PRED, "@@@ Test set @@@")
    
    if (exe_i == 1)  confMatDataFrame <-  thisConfMat
    else  confMatDataFrame <- rbind(confMatDataFrame, thisConfMat)
    
 }
 
cat("\n\n\n=== final results ===\n")
cat("Number of executions = ", execution_number, "\n", sep="")

# statistics on the dataframe of confusion matrices
statDescConfMatr <- stat.desc(confMatDataFrame)
meanSigmaRowResults <- (statDescConfMatr)[c("mean","std.dev"),]
print(statDescConfMatr)
cat("\n\n=== === === ===\n")
print(dec_three(meanSigmaRowResults))
cat("\n\n=== === === ===\n\n\n")

# printResultsLatex("Random forests", meanSigmaRowResults)

#   cat("\t", colnames(meanSigmaRowResults), "\\\\ \n", sep=" & ")
#     cat("mean ", as.character(dec_three((meanSigmaRowResults)["mean",])), sep=" & ")
#     cat("$\\sigma$", as.character(dec_three((meanSigmaRowResults)["std.dev",])), "\\\\ \n", sep=" & ")

cat("original dataset:\n")
dataset_dim_retriever(patients_data_original)
imbalance_retriever(patients_data_original[,target_index])

computeExecutionTime()