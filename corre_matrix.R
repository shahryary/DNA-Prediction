# correct file V.1.1

# find correation values for all "sample" and save to the tabeles 

# loading libraries 
source("load_lib.R")
# get working directory 
current_directory<-getwd()
df_zero_ref<-fread("df_zero_ref.csv")

# --------- 
options(scipen=999)
goodChrOrder <- paste("chr",c(1:22,"X","Y"),sep="")
# setting path to Stem 
stem="/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage_stems/Stem6_1Mbase_Results/"

# Find all Dnase Folders (folders starting with GSM in main directory) 
folders <- dir("/Volumes/DISK_IN/BIGDATA_HSE/Master_These/Corr_matrix_all/BySample",pattern="GSM*")

ptm <- proc.time()
flag=TRUE
for (folder_name in 1:length(folders)){
  GSM_name=folders[folder_name]
  # Pass the folder name to GSM
  gsm=paste("/Volumes/DISK_IN/BIGDATA_HSE/Master_These/Corr_matrix_all/BySample/",GSM_name,"/",sep = "")
  
  # Coverage Files Path
  files = list.files(path = gsm, pattern="*.coverage")
  # read just 1,2,3,7 columns 
  col.class   <- c(NA, NA, NA,"NULL","NULL","NULL",NA)
  # rename columns 
  col.names <- c("chr","start","end","gsm")
  
  #  function to read all files from directory 
  df_gsm<-NULL
  for(i in 1:length(files)){ 
    # Reading files from source 
    tmp<-fread(paste(gsm,files[i],sep = ""),col.names = col.names,colClasses = col.class)
    tmp<- tmp[order(start),]
    df_gsm<-rbind(df_gsm,tmp)
  }
  #----------------- for Stem
  files_stem = list.files(path = stem, pattern="*.coverage")
  # read just 1,2,3,7 columns 
  col.class   <- c(NA, NA, NA,"NULL","NULL","NULL",NA)
  #  function to read all files from directory 
  df_stem<-NULL
  col.names <- c("chr","start","end","stem")
  for(i in 1:length(files_stem)){ 
    # Reading files from source 
    tmp_stem<-fread(paste(stem,files_stem[i],sep = ""),col.names = col.names,colClasses = col.class)
    tmp_stem<- tmp_stem[order(start),]
    df_stem<-rbind(df_stem,tmp_stem)
    
  }
  
  # getting values column and creating matrix data frame to use in correlation 
  df_cor<-df_gsm[,c("chr","start","end","gsm")]
  df_cor$stem<-df_stem[,c("stem")]
  ## ! correct load
  
  
  
  #-------------------- correlation
  df_correlation<-NULL
  for (i in 1:length(goodChrOrder)){
    ch=goodChrOrder[i]
    tmp_df<-df_cor[which(df_cor$chr==ch),]
    
    aa<-df_zero_ref[which(df_zero_ref$chr==ch),]
    tmp_df<-tmp_df[!(tmp_df$start %in% aa$start),]
    
    mat<-tmp_df[,c("gsm","stem")]
    cormatrix<-rcorr(as.matrix(mat),type=c("pearson")) # corr 1.1
    cordata = melt(cormatrix$r)
    te<-cordata$value[2]
    df_correlation<-rbind(df_correlation,c(ch,te))
  }
  
  df_correlation<-as.data.frame(df_correlation)
  names(df_correlation)<-c("chr", paste(GSM_name,sep = ""))
  
  # flag==True it's mean first run (crating basic dataframe )
  if (flag == TRUE){
    tmp_big<-NULL
    tmp_big<-df_correlation[,c(1,2)]
    flag = FALSE
  } else {
    tmp_big$newcol1<-df_correlation[,2]
    colnames(tmp_big)[ncol(tmp_big)] <-paste(GSM_name,sep = "") 
  }
  
}

fwrite(tmp_big,file = "correlation_by_sample_Stem6.csv", sep = "\t",col.names = TRUE)
proc.time() - ptm


