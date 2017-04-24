# loading libraries 
source("load_lib.R")
# get working directory 
current_directory<-getwd()

# --------- 
options(scipen=999)
goodChrOrder <- paste("chr",c(1:22,"X","Y"),sep="")
df_zero_ref<-fread("df_zero_ref.csv")

stem="/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage/Stem6_1Mbase_Results/"
# Find all Dnase Folders (folders starting with GSM in main directory) 
folders <- dir("/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage/DNase_1Mbase",pattern="GSM*")

flag=TRUE
for (folder_name in 1:length(folders)){#
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
    #----------- ignoring 0 positions in files according to the references from "df_zero_ref"
    tmp_char=tmp[1,1]
    aa<-df_zero_ref[which(df_zero_ref$chr==tmp_char$chr),]
    tmp<-tmp[!(tmp$start %in% aa$start),]
    
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
}

#------------------------------------------------------------------------------------
# for Stem folders


stem="/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage/"

# Find all Dnase Folders (folders starting with GSM in main directory) 
folders <- dir("/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage",pattern="Stem*")
#flag=TRUE
for (folder_name in 1:length(folders)){
  Stem_name=folders[folder_name]
  # Pass the folder name to GSM
  stem=paste("/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage/",Stem_name,"/",sep = "")
  main_name=gsub("_.*","",Stem_name)
  # Coverage Files Path
  files = list.files(path = stem, pattern="*.coverage")
  # read just 1,2,3,7 columns 
  col.class   <- c(NA, NA, NA,NA,NA,"NULL",NA)
  # rename columns 
  col.names <- c("chr","start","end","cnt1","cnt2",main_name)
  
  #  function to read all files from directory 
  df_gsm<-NULL
  for(i in 1:length(files)){ 
    # Reading files from source 
    tmp<-fread(paste(stem,files[i],sep = ""),col.names = col.names,colClasses = col.class)
    tmp<- tmp[order(start),]
    #----------- ignoring 0 positions in files
    tmp_char=tmp[1,1]
    aa<-df_zero_ref[which(df_zero_ref$chr==tmp_char$chr),]
    tmp<-tmp[!(tmp$start %in% aa$start),]
    df_gsm<-rbind(df_gsm,tmp)
  }
  # flag==True it's mean first run (crating basic dataframe )
  if (flag == TRUE){
    tmp_big<-NULL
    tmp_big<-df_gsm[,c(1,2,3,6)]
    flag = FALSE
  } else {
    tmp_big$newcol1<-df_gsm[,6]
    colnames(tmp_big)[ncol(tmp_big)] <-paste(main_name,sep = "") 
  }
}




#fwrite(df_zero_ref,file = "df_zero_ref.csv", sep = "\t",col.names = TRUE)
#fwrite(tmp_big,file = "finalTable.csv", sep = "\t",col.names = TRUE)
  
mytest<-fread("/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage_bysample/By_sample_finalTable.csv")
