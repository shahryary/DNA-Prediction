# Find coverage per gene 

# loading libraries 
source("load_lib.R")
# get working directory 
current_directory<-getwd()

# --------- 
options(scipen=999)
goodChrOrder <- paste("chr",c(1:22,"X","Y"),sep="")
# setting path to Stem 
stem="/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage/Stem6_1Mbase_Results/"

# Find all Dnase Folders (folders starting with GSM in main directory) 
folders <- dir("/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage/DNase_1Mbase",pattern="GSM*")
#folders[1:2]
for (folder_name in 1:length(folders)){
    GSM_name=folders[folder_name]
    # Pass the folder name to GSM
    gsm=paste("/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage/DNase_1Mbase/",GSM_name,"/",sep = "")

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
  mat<-tmp_df[,c("gsm","stem")]
  cormatrix<-rcorr(as.matrix(mat),type=c("pearson")) # corr 1.1
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

# setting palate color 
cols <- colorRampPalette(brewer.pal(3, "OrRd"))
myPal <- cols(length(unique(df.long$V2>0.7)))

# loading theme and plottting 
theme_set(theme_solarized())
bplot <- ggplot(df.long,aes(x=V1,y=V2,fill=V2))
bplot <- bplot + geom_bar(stat="identity")
bplot <- bplot + geom_tile(color="white",size=0.1)
bplot <- bplot + geom_text(aes(label = formattable(V2,digit=2,format='f')), position=position_dodge(width=0.9), vjust=-0.25)
bplot <- bplot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
bplot <- bplot + theme(axis.ticks=element_blank())
bplot <- bplot + labs(x=NULL, y=NULL, title=paste("Spearman Correlation for: ",GSM_name))
bplot <- bplot + scale_y_continuous(limits=c(0, 1), breaks=c(0.0,0.2, 0.4, 0.6, 0.8,1.00))
bplot <- bplot + scale_fill_gradientn(colours = myPal)
bplot <- bplot + theme(legend.position='none')
bplot <- bplot + theme(panel.background = element_rect(), panel.grid.major.y = element_line( colour = "gray",linetype = "dashed"))
bplot
# save to pngs file
png(paste("Image_Results_stem6/",GSM_name,".png",sep = ""),1200,800)
print(bplot)
dev.off()

}


