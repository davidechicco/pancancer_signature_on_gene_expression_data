#!/bin/bash
#
#$ -cwd
#$ -S /bin/bash
#
set -o nounset -o pipefail -o errexit
set -o xtrace

array_datasets_files=("../data/Beauchamp2014_GSE38832_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Bild2005_GSE3141_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Bogunovic2009_GSE19234_skin_cancer_our_pancancer_signature_probesets_dataset_5056.csv" "../data/Botling2012_GSE37745_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Chen2020_GSE161158_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Dedeurwaerder2011_GSE20711_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv" 
"../data/DelRoi2017_GSE72970_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Desmedt2007_GSE7390_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Gotoh2018_GSE92921_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv"  "../data/Hatzis2009_GSE25066_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Heaton2011_GSE33371_adrenocortical_carcinoma_our_pancancer_signature_probesets_dataset_1786.csv"  
"../data/Heiskanen2015_GSE68465_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv"  "../data/Herold2011_GSE22762_leukemia_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Herold2013_GSE37642_leukemia_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Hiyama2010_GSE16237_neuroblastoma_our_pancancer_signature_probesets_dataset_188.csv" "../data/Huang2014_GSE48390_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Iqbal2015_GSE58445_non_hodgkin_lymphoma_our_pancancer_signature_probesets_dataset_1786.csv" 
"../data/Ivshina2006_GSE4922_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Jezequel2015_GSE58812_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Karn2011_GSE31519_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Kawaguchi2012_GSE34771_non_Hodgkin_lymphoma_our_pancancer_signature_probesets_dataset_1786.csv" 
"../data/Kim2020_GSE135565_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Kohno2011_GSE31210_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Leich2009_GSE16131_non_hodgkin_lymphoma_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Lenz2008_GSE10846_non_hodgkin_lymphoma_our_pancancer_signature_probesets_dataset_1786.csv" 
"../data/Lin2009_GSE19697_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Loi2008_GSE9195_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Marisa2013_GSE39582_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Metzeler2018_GSE12417_leukemia_our_pancancer_signature_probesets_dataset_1786.csv"  "../data/MetzgerFilho2016_GSE88770_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv" 
"../data/Micke2011_GSE28571_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Miller2013_GSE45255_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Mulligan2007_GSE9782_multiple_myeloma_our_pancancer_signature_probesets_dataset_7648.csv" "../data/Pasini2011_GSE38749_stomach_cancer_our_pancancer_signature_probesets_dataset_1786.csv" 
"../data/Philipsen2010_GSE19188_lung_cancer_our_pancancer_signature_probesets_dataset_805.csv" "../data/Pintilie2013_GSE50081_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Potti2006_GSE3593_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Reister2012_GSE31684_bladder_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Rousseaux2013_GSE30219_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv"  
"../data/Sabatier2010_GSE21653_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Schmidt2008_GSE11121_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Shi2010_GSE24080_multiple_myeloma_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Shinto2020_GSE143985_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv" 
"../data/Sieber2010_GSE14333_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Sinn2009_GSE124647_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Smith2009a_GSE17536_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Smith2009b_GSE17537_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Son2007_GSE8894_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Spivak2014_GSE47018_leukemia_our_pancancer_signature_probesets_dataset_1786.csv" 
"../data/Staub2009_GSE12945_colorectal_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Tsao2010_GSE14814_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Uehara2015_GSE65986_ovarian_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/VanLoo2009_GSE7788_lymphoma_our_pancancer_signature_probesets_dataset_4134.csv"  "../data/Wang2010_GSE19615_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv" 
"../data/Xie2011_GSE29013_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/Yenamandra2015_GSE61304_breast_cancer_our_pancancer_signature_probesets_dataset_1786.csv" "../data/ZChen2020_GSE157011_lung_cancer_our_pancancer_signature_probesets_dataset_1786.csv")

random_number=$(shuf -i1-100000 -n1)
today=`date +%Y-%m-%d`

# catBoost
method="catBoost"
outputFile=""
outputFile="../results/"$today"_"$method"_random"$random_number".txt"

for dataset_name in ${array_datasets_files[@]}; do
  echo $dataset_name
  Rscript catboost_classification.r $dataset_name >> $outputFile 2>> $outputFile;
done

# lightGBM
method="lightGBM"
outputFile=""
outputFile="../results/"$today"_"$method"_random"$random_number".txt"

for dataset_name in ${array_datasets_files[@]}; do
  echo $dataset_name
  Rscript lightgbm_classification.r $dataset_name >> $outputFile 2>> $outputFile;
done


# DecisionTree
method="DecisionTree"
outputFile=""
outputFile="../results/"$today"_"$method"_random"$random_number".txt"

for dataset_name in ${array_datasets_files[@]}; do
  echo $dataset_name
  Rscript cart_classification.r $dataset_name >> $outputFile 2>> $outputFile;
done


# k-Nearest Neighbors
method="kNN"
outputFile=""
outputFile="../results/"$today"_"$method"_random"$random_number".txt"

for dataset_name in ${array_datasets_files[@]}; do
  echo $dataset_name
  Rscript knn_classification.r $dataset_name >> $outputFile 2>> $outputFile;
done
