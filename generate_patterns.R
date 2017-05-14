source("load_lib.R")

options(scipen=999)
goodChrOrder <- paste("chr",c(1:22,"X","Y"),sep="")
df_zero_ref<-fread("df_zero_ref.csv")

# root folder
# set this path accroding to the Stem folder located in your hard disk
stem_dir = "/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage_stems/Stem15_1Mbase_Results/"
# set stem number - Stem6 or Stem15 or Stem16
stem_num="Stem15"
# root folder
# you can set path for by_sample or Dnase
gsm_dir = "/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage_bysample/By_Sample_1Mbase/"
gsm_dir_name="skin_fetal"
gsm_ir_type="DNase_hypersensitivity"
gsm_number="GSM817158"

# detecting patterns 
# N= number of repetations to detect pattern  
N = 3
pattern<-function(gsm,stem){
  #gsm<-gsm[103:110]
  #stem<-stem[103:110]
  
  # first position of data frames
  startpos=0
  N=N-1
  # list 
  i<-1
  pos_start<-c()
  pos_stop<-c()
  # number of records
  lenset=NROW(stem)
  size=N
  while (startpos<lenset){
    if (size>lenset){
      size=lenset
      startpos=size # to finish 
    }
    tmp_gsm<-gsm[startpos:size]
    tmp_stem<-stem[startpos:size]

    # check if going to decress
  if(all(tmp_gsm$count == cummin(tmp_gsm$count))==TRUE){
    if(all(tmp_stem$count == cummin(tmp_stem$count))==TRUE){
      pos_start[[i]] <- head(tmp_gsm$start,1)
      pos_stop[[i]] <- tail(tmp_gsm$stop,1)
      i <- i + 1
      cat("added")
    }
  }
    # if increasing 
  if(all(tmp_gsm$count == cummax(tmp_gsm$count))==TRUE){
      if(all(tmp_stem$count == cummax(tmp_stem$count))==TRUE){
        pos_start[[i]] <- head(tmp_gsm$start,1)
        pos_stop[[i]] <- tail(tmp_gsm$stop,1)
        i <- i + 1
        cat("added")
        
      }
    }

    startpos=startpos+1
    size=size+1

  }
  # checking if there is no pattern for number N
  if(length(pos_start)==0){
    pos_start[[1]]<-0
    pos_stop[[1]]<-0
  }
  patt <- data.frame(xmin=pos_start, xmax=pos_stop, ymin=-Inf, ymax=Inf)
  return(patt)
}


for (i in 1:length(goodChrOrder)){
  chr=goodChrOrder[20]
  col.class   <- c(NA, NA, NA,"NULL","NULL","NULL",NA)
  col.names <- c("chr", "start", "stop", "count")

  stem <- fread(paste(stem_dir,chr,".",stem_num,".coverage",sep = ""),col.names = col.names,colClasses =col.class )
  stem<- stem[order(start),]
  gsm <- fread(paste(gsm_dir,gsm_dir_name,"/",gsm_ir_type,"/",gsm_number,"/",chr,".",gsm_number,".coverage",sep = ""),col.names = col.names,colClasses =col.class )
  #----------------- ignore 0 points
  tmp_char=chr
  aa<-df_zero_ref[which(df_zero_ref$chr==tmp_char),]
  stem<-stem[!(stem$start %in% aa$start),]
  gsm<-gsm[!(gsm$start %in% aa$start),]
  
  #-----------------
  df <- data.frame(position=stem$start,position_stop=stem$stop, stemloops=stem$count, Dnas=gsm$count)
  
  patt<-pattern(gsm,stem) # detecting patterns 
  gg <- ggplot(df, aes(x=position, y=signal, color=variable))
  gg <- gg + geom_line(aes(y=log(stemloops/mean(stem$count)),col="Stemloops")) 
  gg <- gg + geom_line(aes(y=log(Dnas/mean(gsm$count)), col="Dnas") )
  title<-paste(gsm_dir_name,gsm_ir_type,gsm_number,"Chromosome:",chr,"With Stemloops: ",stem_num,sep = "    ")
  gg <- gg + ggtitle(title)
  gg <- gg+ geom_rect(data=patt, aes(xmin=patt$xmin, xmax=patt$xmax, ymin=-Inf, ymax=Inf),
                      fill="orange",
                      alpha=0.2,
                      inherit.aes = FALSE)
  gg
  # save to pngs file
  # tmp is directory that you want to save files
  png(paste("tmp/",gsm_number,chr,".png",sep = ""),1200,800)
  print(gg)
  dev.off()
}

gsm1<-gsm[1:5,4]
all(gsm1 == cummin(gsm1))
stem1<-stem[1:5,4]
all(stem1 == cummin(stem1))

if(all(gsm1 == cummin(gsm1))==all(stem1 == cummin(stem1))){
  cat("True")
}else{
  cat("False")
}



rcorr(df$stemloops,df$Dnas)


tmp_df<-df_cor[which(df_cor$chr==ch),]
mat<-df[,c("stemloops","Dnas")]

cormatrix<-rcorr(as.matrix(mat/1000000),type=c("pearson")) # corr 1.1
cordata = melt(cormatrix$r)
te<-cordata$value[2]


