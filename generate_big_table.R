# loading libraries 
source("load_lib.R")
# get working directory 
current_directory<-getwd()

# --------- 
options(scipen=999)
goodChrOrder <- paste("chr",c(1:22,"X","Y"),sep="")

stem="/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage/Stem6_1Mbase_Results/"

# Find all Dnase Folders (folders starting with GSM in main directory) 
folders <- dir("/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage/DNase_1Mbase",pattern="GSM*")
#folders[1:2]
flag=TRUE
for (folder_name in 1:10){#length(folders)
  GSM_name=folders[folder_name]
  # Pass the folder name to GSM
  gsm=paste("/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage/DNase_1Mbase/",GSM_name,"/",sep = "")
  
  # Coverage Files Path
  files = list.files(path = gsm, pattern="*.coverage")
  # read just 1,2,3,7 columns 
  col.class   <- c(NA, NA, NA,NA,NA,"NULL",NA)
  # rename columns 
  col.names <- c("chr","start","end","cnt1","cnt2",GSM_name)
  
  #  function to read all files from directory 
  df_gsm<-NULL
  for(i in 1:length(files)){ 
    # Reading files from source 
    tmp<-fread(paste(gsm,files[i],sep = ""),col.names = col.names,colClasses = col.class)
    tmp<- tmp[order(start),]
    df_gsm<-rbind(df_gsm,tmp)
  }
  # flag==True it's mean first run (crating basic dataframe )
  if (flag == TRUE){
    tmp_big<-NULL
    tmp_big<-df_gsm[,c(1,2,3,6)]
    flag = FALSE
  } else {
  tmp_big$newcol1<-df_gsm[,6]
  colnames(tmp_big)[ncol(tmp_big)] <-paste(GSM_name,sep = "") 
  }
  #names(tmp_big)[-1]<-paste(GSM_name,"",sep = "")
  #---- just ignoring in first level---------- for Stem   
  # files_stem = list.files(path = stem, pattern="*.coverage")
  # # read just 1,2,3,7 columns 
  # col.class   <- c(NA, NA, NA,"NULL","NULL","NULL",NA)
  # #  function to read all files from directory 
  # df_stem<-NULL
  # col.names <- c("chr","start","end","stem")
  # for(i in 1:length(files_stem)){ 
  #   # Reading files from source 
  #   tmp_stem<-fread(paste(stem,files_stem[i],sep = ""),col.names = col.names,colClasses = col.class)
  #   tmp_stem<- tmp_stem[order(start),]
  #   df_stem<-rbind(df_stem,tmp_stem)
  #   
  # }
} 

  
