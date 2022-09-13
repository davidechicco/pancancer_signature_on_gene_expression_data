setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)

# fileName <- "../data/ZChen2020_GSE157011_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv"

args <- commandArgs(trailingOnly=TRUE)

# The dataset must be in a .csv file with patients on the rows and features on the column.
# The first row must contain the column feauture names.
# The last column on the right must contain the binary target
# fileName <- "/home/dave/my_projects/cardiovascular_heart_disease/data/dataset_edited_without_time.csv" # TO CHANGE
# targetName <- "death_event" #TO CHANGE

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  fileName <- args[1]
}


targetName <- "survival"

cat("fileName: ", fileName, "\n", sep="")
cat("targetName: ", targetName, "\n", sep="")

list.of.packages <- c("easypackages", "clusterSim", "PRROC", "e1071", "rpart",  "dplyr", "pastecs", "stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

library("easypackages")
libraries(list.of.packages)

source("./confusion_matrix_rates.r")
source("./utils.r")

threshold <- 0.5

# file reading
patients_data <- read.csv(fileName, header = TRUE, sep =",");
cat("Read data from file ", fileName, "\n", sep="")

NUM_METRICS <- 9
confMatDataFrame <- matrix(ncol=NUM_METRICS, nrow=1)
colnames(confMatDataFrame) <- c("MCC", "F1 score", "accuracy", "TP_rate", "TN_rate", "PPV", "NPV", "PR_AUC", "ROC_AUC")

# let's put the target label last on the right 
patients_data <- patients_data%>%select(-targetName,targetName)

target_index <- dim(patients_data)[2]
original_patients_data <- patients_data

execution_number <- 100
cat("Number of executions = ", execution_number, "\n", sep="")
for(exe_i in 1:execution_number)
{
    cat("[Execlution number ", exe_i, " out of ", execution_number, "]\n", sep="" )
    
    patients_data <- patients_data[sample(nrow(patients_data)),] # shuffle the rows
    totalElements <- dim(patients_data)[1]

    subsets_size <- totalElements

    if (subsets_size != totalElements) {
        cat("!!! ATTENTION: We are running the method on a subset of the original dataset, \n", sep="")
        cat("!!! containing only ", subsets_size, " elements \n", sep="")
        cat("!!! instead of ", totalElements, " elements \n", sep="")
    }

    patients_data <- patients_data[1:subsets_size, ]

    dataset_dim_retriever(patients_data)
    imbalance_retriever(patients_data$Metastasis)


    target_index <- dim(patients_data)[2]

    training_set_perce <- 60
    cat("training_set_perce = ", training_set_perce, "% \n", sep="")
    validation_set_perce <- 20
    cat("validation_set_perce = ", validation_set_perce, "% \n", sep="")
    test_set_perce <- 100 - training_set_perce - validation_set_perce
    cat("test_set_perce = ", test_set_perce, "% \n", sep="")

    # the training set is the first 60% of the whole dataset
    training_set_first_index <- 1 # NEW
    training_set_last_index <- round(dim(patients_data)[1]*training_set_perce/100) # NEW

    # the validation set is the following 20% of the whole dataset
    validation_set_first_index <- round(dim(patients_data)[1]*training_set_perce/100)+1 # NEW
    validation_set_last_index <- round(dim(patients_data)[1]*(training_set_perce+validation_set_perce)/100) # NEW

    # the test set is the last 20% of the whole dataset
    test_set_first_index <- round(dim(patients_data)[1]*(training_set_perce+validation_set_perce)/100)+1 # NEW
    test_set_last_index <- dim(patients_data)[1] # NEW

    cat("[Creating the subsets for the values]\n")
    patients_data_train <- patients_data[training_set_first_index:training_set_last_index, 1:(target_index-1)] # NEW
    patients_data_validation <- patients_data[validation_set_first_index:validation_set_last_index, 1:(target_index-1)] # NEW
    patients_data_test <- patients_data[test_set_first_index:test_set_last_index, 1:(target_index-1)] # NEW

    cat("[Creating the subsets for the labels \"1\"-\"0\"]\n")
    patients_data_train_labels <- patients_data[training_set_first_index:training_set_last_index, target_index] # NEW
    patients_data_validation_labels <- patients_data[validation_set_first_index:validation_set_last_index, target_index] # NEW
    patients_data_test_labels <- patients_data[test_set_first_index:test_set_last_index, target_index]   # NEW


    library(class)
    library(gmodels)

    # # The k value must be lower than the size of the trainingset
    maxK <- 10 #NEW

    mcc_array <- character(length(maxK))

    # NEW PART:

    cat("\n[Optimization of the hyper-parameter k start]\n")
    # optimizaion loop
    for(thisK in 2:maxK)
    {
    # apply k-NN with the current K value
    # train on the training set, evaluate in the validation set by computing the MCC
    # save the MCC corresponding to the current K value
    
            cat("\n[Training the kNN model (with k=",thisK,") on training set & applying the kNN model to validation set]\n", sep="")
            
            patients_data_validation_pred <- knn(train = patients_data_train, test = patients_data_validation, cl = patients_data_train_labels, k=thisK)
            patients_data_validation_pred_binary <- as.numeric (patients_data_validation_pred)-1
            
            mcc_outcome <- mcc(patients_data_validation_labels, patients_data_validation_pred_binary)
            cat("When k=",thisK,", the MCC value is ",mcc_outcome, "\t (worst possible: -1; best possible: +1)\n", sep="")
            
            mcc_array[thisK] <- mcc_outcome
    
    }

    # select the k corresponding to the highest MCC and call it k_best
    bestMCC <- max(mcc_array)
    bestK <- match(bestMCC, mcc_array)
    cat("\nThe best k value is ", bestK,", corresponding to MCC=", mcc_array[bestK],"\n", sep="")

    cat("[Optimization end]\n\n")

    cat("\n @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ @ \n")

    # apply k-NN with k_best to the test set

    cat("[Training the kNN model (with the OPTIMIZED hyper-parameter k=",bestK,") on training set & applying the kNN to the test set]\n", sep="")
    patients_data_test_pred <- knn(train = patients_data_train, test = patients_data_test, cl = patients_data_train_labels, k=bestK)

    patients_data_test_pred <- as.numeric(patients_data_test_pred)-1

    thisConfMat <- confusion_matrix_rates(patients_data_test_labels, patients_data_test_pred, "@@@ Test set @@@")


        if (exe_i == 1)  confMatDataFrame <-  thisConfMat
        else  confMatDataFrame <- rbind(confMatDataFrame, thisConfMat)
    
 }
 
cat("\n\n\n=== final results ===\n")
cat("Number of executions = ", execution_number, "\n", sep="")

# statistics on the dataframe of confusion matrices
statDescConfMatr <- stat.desc(confMatDataFrame)
meanAndSdRowResults <- (statDescConfMatr)[c("mean", "std.dev"),]
print(dec_three(statDescConfMatr))

datasetName <- (str_split(fileName, "/")[[1]][3] %>% str_split("_"))[[1]][1]
meanMCC <- meanAndSdRowResults[1,1]
cat("\ndataset  averageMCC\n")
cat(datasetName, " ", dec_three(meanMCC), "\n")
cat("\n\n=== === === ===\n")
