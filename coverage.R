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
    files = list.files(path = stem, pattern="*.coverage")
    # read just 1,2,3,7 columns 
    col.class   <- c(NA, NA, NA,"NULL","NULL","NULL",NA)
    #  function to read all files from directory 
    df_stem<-NULL
    col.names <- c("chr","start","end","stem")
    for(i in 1:length(files)){ 
      # Reading files from source 
      tmp<-fread(paste(stem,files[i],sep = ""),col.names = col.names,colClasses = col.class)
      tmp<- tmp[order(start),]
      df_stem<-rbind(df_stem,tmp)
     
    }
    
# getting values column and creating matrix data frame to use in correlation 
df_cor<-df_gsm[,c("chr","start","finish","gsm")]
df_cor$stem<-df_stem[,c("stem")]

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


snpDensity<-ggplot(df_correlation)+ 
  geom_histogram(aes(x=start),binwidth=1e6,fill="#c0392b", alpha=0.75) +  # pick a binwidth that is not too small 
  fte_theme()+
  facet_wrap(~ chr,ncol = 2) + # seperate plots for each chr, x-scales can differ from chr to chr
  ggtitle("Density of chromosome  hg19") +
  xlab("Position in the genome") + 
  ylab("SNP density") 

plot(df_correlation)


test<-df_cor[which(df_cor$chr=="chr14"),]
mat<-test[,c("gsm","stem")]
cormatrix<-rcorr(as.matrix(mat),type=c("pearson")) # corr 1.1

cordata = melt(cormatrix$r)
te<-cordata$value[2]
hm.palette <- colorRampPalette(c('light green', 'dark green'))

gg <- ggplot(df_correlation, aes(x=Var1, y=Var2, fill=value))
gg <- gg + geom_tile(color="white",size=0.1)
#gg <- gg + geom_text(aes(label = formattable(value,digit=2,format='f')))
gg <- gg + theme(axis.ticks=element_blank())
gg <- gg + theme(axis.text=element_text(size=7))
gg <- gg + theme(legend.title=element_text(size=8))
gg <- gg + theme(legend.text=element_text(size=6))
gg <- gg + theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg <- gg + labs(x=NULL, y=NULL, title="Spearman Correlation Matrix")
gg <-  gg + scale_fill_gradientn(colours = hm.palette(100))

gg

cor(df_cor)
