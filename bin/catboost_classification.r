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

cat("input fileName: ", fileName, "\n", sep="")
cat("targetName: ", targetName, "\n", sep="")

list.of.packages <- c("easypackages", "catboost", "formula.tools", "dplyr", "pastecs", "stringr")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

library("easypackages")
libraries(list.of.packages)


source("./confusion_matrix_rates.r")
source("./utils.r")

NUM_METRICS <- 9
confMatDataFrame <- matrix(ncol=NUM_METRICS, nrow=1)
colnames(confMatDataFrame) <- c("MCC", "F1 score", "accuracy", "TP rate", "TN rate", "PPV", "NPV", "PR AUC", "ROC AUC")

threshold <- 0.5

patients_data <- read.csv(file=fileName,head=TRUE,sep=",",stringsAsFactors=FALSE)
cat("fileName = ", fileName, "\n", sep="")

# let's put the target label last on the right 
patients_data <- patients_data%>%select(-targetName,targetName)
patients_data_original <- patients_data

# cycle of executions

execution_number <- 100
cat("Number of executions = ", execution_number, "\n", sep="")
for(exe_i in 1:execution_number)
{

    cat("[Execlution number ", exe_i, " out of ", execution_number, "]\n", sep="" )
    patients_data <- patients_data[sample(nrow(patients_data)),] # shuffle the rows

    target_index <- dim(patients_data)[2]
    cat("target_index = ", target_index, "\n", sep="") 

    training_set_perce <- 80

    cat("training_set_perce = ", training_set_perce, "%\n", sep="")

    # the training set is the first training_set_perce% of the whole dataset
    training_set_first_index <- 1 # NEW
    training_set_last_index <- round(dim(patients_data)[1]*training_set_perce/100) # NEW

    # the test set is the last 20% of the whole dataset
    test_set_first_index <- round(dim(patients_data)[1]*training_set_perce/100)+1 # NEW
    test_set_last_index <- dim(patients_data)[1] # NEW

    cat("[Creating the subsets for the values]\n")
    prc_data_train <- patients_data[training_set_first_index:training_set_last_index, 1:(target_index-1)] # NEW
    prc_data_test <- patients_data[test_set_first_index:test_set_last_index, 1:(target_index-1)] # NEW


    cat("[Creating the subsets for the labels \"1\"-\"0\"]\n")
    prc_data_train_labels <- patients_data[training_set_first_index:training_set_last_index, target_index] # NEW
    prc_data_test_labels <- patients_data[test_set_first_index:test_set_last_index, target_index]   # NEW

    # apply k-NN with k_best to the test set

    cat("\n[Training the catBoost model on training set & applying the catBoost to test set]\n", sep="")

    # lin_reg_model_new <- lm(prc_data_train_labels ~ ., data=prc_data_train)
    
    fit_params <- list(iterations = 100,
                   loss_function = 'Logloss',
                   logging_level = 'Silent')
                   
                   
    training_data_pool <- catboost.load_pool(prc_data_train, label = prc_data_train_labels)
    test_data_pool <- catboost.load_pool(prc_data_test, label = prc_data_test_labels)
    catboost_trained_model <- catboost.train(training_data_pool, params = fit_params)

    # prc_data_test_pred <- predict(lin_reg_model_new, prc_data_test)
    prc_data_test_pred <- catboost.predict(catboost_trained_model, test_data_pool,  prediction_type = "Class")

    prc_data_test_pred_bin <- as.numeric(prc_data_test_pred)
    prc_data_test_pred_bin[prc_data_test_pred_bin>=threshold]<-1
    prc_data_test_pred_bin[prc_data_test_pred_bin<threshold]<-0

    thisConfMat <- confusion_matrix_rates(prc_data_test_labels, prc_data_test_pred_bin, "@@@ Test set @@@")
    if (exe_i == 1)  confMatDataFrame <-  thisConfMat
    else  confMatDataFrame <- rbind(confMatDataFrame, thisConfMat)
    
 }
 
cat("\n\n\n=== final results ===\n")
cat("analyzed_dataset: ", fileName, "\n", sep="")
cat("Number of executions = ", execution_number, "\n", sep="")

# statistics on the dataframe of confusion matrices
statDescConfMatr <- stat.desc(confMatDataFrame)
meanAndSdRowResults <- (statDescConfMatr)[c("mean", "std.dev"),]
print(dec_three(meanAndSdRowResults))

datasetName <- (str_split(fileName, "/")[[1]][3] %>% str_split("_"))[[1]][1]
meanMCC <- meanAndSdRowResults[1,1]
cat("\ndataset  averageMCC\n")
cat(datasetName, " ", dec_three(meanMCC), "\n")
cat("\n\n=== === === ===\n")
