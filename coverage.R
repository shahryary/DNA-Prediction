# Find coverage per gene 

# loading libraries 
source("load_lib.R")
# get working directory 
current_directory<-getwd()

# --------- 
options(scipen=999)
# setting path to GSM
gsm="/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage/DNase_1Mbase/GSM1027306/"
# setting path to Stem 
stem="/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage/Stem15_1Mbase_Results/"
goodChrOrder <- paste("chr",c(1:22,"X","Y"),sep="")


    # Coverage Files Path
    files = list.files(path = gsm, pattern="*.coverage")
    # read just 1,2,3,7 columns 
    col.class   <- c(NA, NA, NA,"NULL","NULL","NULL",NA)
    # rename columns 
    col.names <- c("chr","start","finish","gsm")
    
    #  function to read all files from directory 
      df_gsm<-NULL
      for(i in 1:length(files)){ 
        # Reading files from source 
        tmp<-fread(paste(gsm,files[i],sep = ""),col.names = col.names,colClasses = col.class)
        tmp<- tmp[order(start),]
        df_gsm<-rbind(df_gsm,tmp)
      }
    #----------------- for Stem
    files = list.files(path = stem, pattern="*.coverage")
    # read just 1,2,3,7 columns 
    col.class   <- c(NA, NA, NA,"NULL","NULL","NULL",NA)
    #  function to read all files from directory 
    df_stem<-NULL
    col.names <- c("chr","start","finish","stem")
    for(i in 1:length(files)){ 
      # Reading files from source 
      tmp<-fread(paste(stem,files[i],sep = ""),col.names = col.names,colClasses = col.class)
      tmp<- tmp[order(start),]
      df_stem<-rbind(df_stem,tmp)
     
    }
    
# getting values column and creating matrix data frame to use in correlation 
df_cor<-df_gsm[,c("chr","gsm")]
df_cor$stem<-df_stem[,c("stem")]


