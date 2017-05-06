# correct file V.1.0
# 06-May

# find correation values for all "sample" -General correlation - (NOT 24 chromosome)

# loading libraries 
source("load_lib.R")
# get working directory 
current_directory<-getwd()
df_zero_ref<-fread("df_zero_ref.csv")

# --------- 
options(scipen=999)
goodChrOrder <- paste("chr",c(1:22,"X","Y"),sep="")
# setting path to Stem 
stem="/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage_stems/Stem15_1Mbase_Results/"

# Find all Dnase Folders (folders starting with GSM in main directory) 
folders <- dir("/Volumes/DISK_IN/BIGDATA_HSE/Master_These/Corr_matrix_all/BySample",pattern="GSM*")

ptm <- proc.time()
flag=TRUE
tmp_big<-NULL
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
    # excluding 0 points
    aa<-df_zero_ref[which(df_zero_ref$chr==tmp$chr[1]),]
    tmp<-tmp[!(tmp$start %in% aa$start),]
    # binding for large table
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
    # excluding 0 points
    aa<-df_zero_ref[which(df_zero_ref$chr==tmp_stem$chr[1]),]
    tmp_stem<-tmp_stem[!(tmp_stem$start %in% aa$start),]
    # binding for large table
    df_stem<-rbind(df_stem,tmp_stem)
    
  }
  
  # getting values column and creating matrix data frame to use in correlation 
  df_cor<-df_gsm[,c("chr","start","end","gsm")]
  df_cor$stem<-df_stem[,c("stem")]
  ## ! correct load
  
  
  
  #-------------------- correlation
  df_correlation<-NULL

    mat<-df_cor[,c("gsm","stem")]
    cormatrix<-rcorr(as.matrix(mat),type=c("pearson")) # corr 1.1
    cordata = melt(cormatrix$r)
    te<-cordata$value[2]
    df_correlation<-rbind(df_correlation,c(te))
  
  df_correlation<-as.data.frame(df_correlation)
  names(df_correlation)<-c( paste(GSM_name,sep = ""))
  

  tmp_big <-rbind(tmp_big,c(GSM_name,df_correlation[,1]))

  
}

tmp_big<-data.frame(tmp_big, stringsAsFactors=FALSE)
names(tmp_big)<-c("GSM","Stem15")
tmp_big[] <- lapply(tmp_big, as.character)
tmp_big$Stem15<-as.numeric(as.character(tmp_big$Stem15))


fwrite(tmp_big,file = "tmp_general_Stem15.csv", sep = "\t",col.names = TRUE)

xx<-fread(file = "tmp_general_Stem15.csv")

#---------------------- cleaning data and create a general table
# 2478 samples

Stem16 <- fread("tmp_general_Stem16.csv",col.names = c("GSM","Stem16") )
Stem15<-fread("tmp_general_Stem15.csv",col.names = c("GSM","Stem15") )
Stem6<-fread("tmp_general_Stem6.csv",col.names = c("GSM","Stem6") )

tmp_cv<-merge(Stem16, Stem15, by.x=c("GSM"), by.y=c("GSM"))

cv<-merge(tmp_cv, Stem6, by.x=c("GSM"), by.y=c("GSM"))

sDes <- fread("SampleDescription.csv",header = FALSE,col.names = c("tissue","experiment","gsm") )

final_result<-merge(sDes, cv, by.x=c("gsm"), by.y=c("GSM"))

fwrite(final_result,file = "general_corr_by_sample.csv", sep = "\t",col.names = TRUE)


# ------------------------  extract DNase from By_sample

Dnase_final_result<-final_result[which(final_result$experiment =="DNase_hypersensitivity"),]
fwrite(Dnase_final_result,file = "general_corr_DNase.csv", sep = "\t",col.names = TRUE)
xx<-fread(file = "general_corr_DNase.csv")
