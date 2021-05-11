setwd(".")
options(stringsAsFactors = FALSE)
cat("\014")
set.seed(11)   

list.of.packages <- c("easypackages", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="https://utstat.toronto.edu/cran/")

library("easypackages")
libraries(list.of.packages)

# Function that reads in a vector made of binary values and prints the imbalance rates
imbalance_retriever_new <- function(thisVector)
{
  lun <- length(table(thisVector))
  if (lun != 2) {
  
    print("This vector is not binary. The imbalance_retriever() function will stop here");
    return(FALSE);
  
  }  

  cat("[Imbalance of this dataset]\n")
  number_of_elements_of_first_class <- unname(table(thisVector)[1])
  name_of_elements_of_first_class <- names(table(thisVector)[1])
  cat("[class: ",name_of_elements_of_first_class, "  #elements = ", number_of_elements_of_first_class, "]\t", sep="")
  cat(dec_three(unname(table(thisVector))[1]*100/length(thisVector)),"%\n", sep="")
  
  number_of_elements_of_second_class <-unname(table(thisVector)[2])
  name_of_elements_of_second_class <-names(table(thisVector)[2])
  cat("[class: ",name_of_elements_of_second_class, "  #elements = ", number_of_elements_of_second_class, "]\t", sep="")
  cat(dec_three(unname(table(thisVector))[2]*100/length(thisVector)),"%\n", sep="")
  
  cat("\n")

  return(TRUE);
}

# let's compute time
global_start_time <- Sys.time()

source("utils.r")

targetName <- "survival"

datasetFileNames <- list.files("../data", pattern="*.csv", full.names=TRUE)

for(thisFileName in datasetFileNames) {

        patients_data <- read.csv(file=thisFileName,head=TRUE,sep=",",stringsAsFactors=FALSE)
        cat("fileName = ", thisFileName, "\n", sep="")
        
        # let's put the target label last on the right 
        patients_data <- patients_data%>%select(-targetName,targetName)
        target_index <- ncol(patients_data)
        imbalance_retriever_new(patients_data[,target_index])

}