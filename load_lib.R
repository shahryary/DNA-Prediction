# libraries source file 
list.pkg<- c("data.table", "abind","ggplot2","corrplot","stringr","plyr","caret","doParallel","tidyr","miscTools",
             "randomForest","Hmisc","RColorBrewer","formattable","ggthemes","reshape2","gtools","dplyr")

req_pkg<-function(packages){
  new.pkg <- list.pkg[!(list.pkg %in% installed.packages()[,"Package"])]
  if(length(new.pkg)) 
    install.packages(new.pkg)
  sapply(list.pkg, require, character.only = TRUE)
}

req_pkg(list.pkg)

