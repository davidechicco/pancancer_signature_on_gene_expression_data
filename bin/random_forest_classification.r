setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)

execution_number <- 100

EXP_ARG_NUM <- 2

TRAIN_SET_OVERSAMPLING_SYNTHETIC <- FALSE

# fileName <- "/home/dave/my_projects/probesets_versus_gene_symbols/results/Hiyama2010_GSE16237_gene_symbols_dataset_8454.csv"
# fileName <- "/home/dave/my_projects/probesets_versus_gene_symbols/results/Hiyama2010_GSE16237_probesets_dataset_8454.csv"

# fileName <- "/home/dave/my_projects/probesets_versus_gene_symbols/results/Smith2010_GSE17536_probesets_dataset_2181.csv"
# fileName <- "/home/dave/my_projects/probesets_versus_gene_symbols/results/Smith2010_GSE17536_gene_symbols_dataset_2181.csv"

# fileName <- "/home/dave/my_projects/probesets_versus_gene_symbols/results/Kim2016_GSE76701_gene_symbols_dataset_7113.csv"
# fileName <- "/home/dave/my_projects/probesets_versus_gene_symbols/results/Kim2016_GSE76701_probesets_dataset_7113.csv"


# fileName <- "/home/dave/my_projects/probesets_versus_gene_symbols/results/Bong2010_GSE15852_probesets_dataset_3875.csv"
# fileName <- "/home/dave/my_projects/probesets_versus_gene_symbols/results/Bong2010_GSE15852_gene_symbols_dataset_3875.csv"
# targetName <- "cancer"

# fileName <- "/home/dave/my_projects/probesets_versus_gene_symbols/results/Bilbault2017_GSE33118_probesets_dataset_4687.csv"
# fileName <- "/home/dave/my_projects/probesets_versus_gene_symbols/results/Bilbault2017_GSE33118_gene_symbols_dataset_4687.csv"

fileName <- "/home/dave/my_projects/probesets_versus_gene_symbols/results/Hiyama2010_GSE16237_our_pancancer_signature_probesets_dataset_188.csv"

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

patients_data[,target_index] <- as.factor(patients_data[,target_index])

# formula
allFeaturesFormula <- as.formula(paste(as.factor(colnames(patients_data)[target_index]), '.', sep=' ~ ' ))

# cycle of executions


cat("Number of executions = ", execution_number, "\n", sep="")
for(exe_i in 1:execution_number)
{
    cat("[Execlution number ", exe_i, " out of ", execution_number, "]\n", sep="" )
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

computeExecutionTime()