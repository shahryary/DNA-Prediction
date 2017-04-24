options(scipen=999)
goodChrOrder <- paste("chr",c(1:22,"X","Y"),sep="")

# root folder
# set this path accroding to the Stem folder located in your hard disk
stem_dir = "/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage_stems/Stem15_1Mbase_Results/"
# set stem number - Stem6 or Stem15 or Stem16
stem_num="Stem15"
# root folder
# you can set path for by_sample or Dnase
gsm_dir = "/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage_bysample/By_Sample_1Mbase/"
gsm_dir_name="pancreas"
gsm_ir_type="DNase_hypersensitivity"
gsm_number="GSM1027335"


for (i in 1:length(goodChrOrder)){
      chr=goodChrOrder[i]
      col.class   <- c(NA, NA, NA,"NULL",NA,"NULL","NULL")
      col.names <- c("chr", "start", "stop", "count")
      stem <- fread(paste(stem_dir,chr,".",stem_num,".coverage",sep = ""),col.names = col.names,colClasses =col.class )
      stem<- stem[order(start),]
      gsm <- fread(paste(gsm_dir,gsm_dir_name,"/",gsm_ir_type,"/",gsm_number,"/",chr,".",gsm_number,".coverage",sep = ""),col.names = col.names,colClasses =col.class )
      
      df <- data.frame(position=stem$start, stemloops=stem$count, Dnas=gsm$count)
      gg <- ggplot(df, aes(x=position, y=signal, color=variable))
      gg <- gg + geom_line(aes(y=log(stemloops/mean(stem$count)),col="Stemloops")) 
      gg <- gg + geom_line(aes(y=log(Dnas/mean(gsm$count)), col="Dnas") )
      title<-paste(gsm_dir_name,gsm_ir_type,gsm_number,"Chromosome:",chr,"With Stemloops: ",stem_num,sep = "    ")
      gg <- gg + ggtitle(title)
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
# move 5 by 5 records for each data set , gsm , stem 
# check if there are minimum or there are maximum 
# if there are same add position to list 
# hight light color in png 
rect <- data.frame(xmin=c(0,21000000), xmax=c(20000000,22000000), ymin=-Inf, ymax=Inf)
gg <- gg+ geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                   fill="yellow",
                   alpha=0.2,
                   inherit.aes = FALSE)
gg
