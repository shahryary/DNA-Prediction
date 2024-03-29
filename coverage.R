# correct file Ver. 1.3 
# 02.May.2017
# plotting correlation for GSM 

source("load_lib.R")
# get working directory 
current_directory<-getwd()
df_zero_ref<-fread("df_zero_ref.csv")
tissue_descript<-fread("tissueTypes.csv",colClasses =c(NA, NA, NA,"NULL") )

# --------- 
options(scipen=999)
goodChrOrder <- paste("chr",c(1:22,"X","Y"),sep="")

stem="/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage_stems/Stem6_1Mbase_Results/"
stem_num="Stem6"

# Find all Dnase Folders (folders starting with GSM in main directory) 
folders <- dir("/Volumes/DISK_IN/BIGDATA_HSE/Master_These/Corr_matrix_all/BySample",pattern="GSM*")

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

  
    #-------------------- correlation
    df_correlation<-NULL
    for (i in 1:length(goodChrOrder)){
      ch=goodChrOrder[i]
      tmp_df<-df_cor[which(df_cor$chr==ch),]
      
      aa<-df_zero_ref[which(df_zero_ref$chr==ch),]
      tmp_df<-tmp_df[!(tmp_df$start %in% aa$start),]
      
      mat<-tmp_df[,c("gsm","stem")]
      cormatrix<-rcorr(as.matrix(mat),type=c("pearson"))
      cordata = melt(cormatrix$r)
      te<-cordata$value[2]
      df_correlation<-rbind(df_correlation,c(ch,te))
    }
    names(df_correlation)<-c("chr", "corelation")
    df_correlation<-as.data.frame(df_correlation)
    
    
    df.long<-melt(df_correlation)
    df.long$V2=as.numeric(levels(df.long$V2))[df.long$V2]
    # sort charachter vector in good format 
    df.long$V1 <- factor(df.long$V1,levels=goodChrOrder)
    
    # find description about the Tissue type
    
    tiss_name<-tissue_descript[which(tissue_descript$gsm==GSM_name),]
    #------------------------------------------------------------------
    
    # setting palate color 
    cols <- colorRampPalette(brewer.pal(9, "RdBu"))#PuBu
    
    if(min(df.long$V2) < 0) {
      const=-0.1
      break_point<-round(break_point<-seq(round((min(df.long$V2+const)),1), round(max((df.long$V2-const)),1), by=0.1),1)
      myPal <- cols(length(unique(df.long$V2<0.0)))
    }else{
      const=0.1
      break_point<-round(break_point<-seq(0, round(max((df.long$V2+const)),1), by=0.1),1)
      cols <- colorRampPalette(brewer.pal(9, "Blues"))#PuBu
      myPal <- cols(length(unique(df.long$V2<0.5)))
      
    }
    
    # loading theme and plottting 
    theme_set(theme_solarized())
    bplot <- ggplot(df.long,aes(x=V1,y=V2,fill=V2))
    bplot <- bplot + geom_bar(stat="identity")
    #bplot <- bplot + geom_tile(color="white",size=0.1)
    bplot <- bplot + geom_text(aes(label = formattable(V2,digit=2,format='f')), position=position_dodge(width=0.9), vjust=-0.25)
    bplot <- bplot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    bplot <- bplot + theme(axis.ticks=element_blank())
    bplot <- bplot + labs(x=NULL, y=NULL, title=paste("Pearson Correlation for: ",GSM_name,",From: ",tiss_name$tissue,tiss_name$type,",with StemLoops: ",stem_num,sep = "  "))
    bplot <- bplot + scale_y_continuous(limits=c(min(break_point), max(break_point)),breaks=break_point)
    bplot <- bplot + scale_fill_gradientn(colours = myPal)
    bplot <- bplot + theme(legend.position='none')
    bplot <- bplot + theme(panel.background = element_rect(), panel.grid.major.y = element_line( colour = "gray",linetype = "dashed"))
    bplot
# save to pngs file
png(paste("Image_Results_stem6/",GSM_name,".png",sep = ""),1200,800)
print(bplot)
dev.off()

}


