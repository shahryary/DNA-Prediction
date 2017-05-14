# for teacher and see chromosome by chromosome.

options(scipen=999)
goodChrOrder <- paste("chr",c(1:22,"X","Y"),sep="")
df_zero_ref<-fread("df_zero_ref.csv")

# root folder
# set this path accroding to the Stem folder located in your hard disk
stem_dir = "/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage_stems/Stem6_1Mbase_Results/"
# set stem number - Stem6 or Stem15 or Stem16
stem_num="Stem6"
# root folder
# you can set path for by_sample or Dnase
gsm_dir = "/Volumes/DISK_IN/BIGDATA_HSE/Master_These/coverage_bysample/By_Sample_1Mbase/"
gsm_dir_name="CD4_memory_primary_cells"
gsm_ir_type="MeDIP-Seq"
gsm_number="GSM669608"


for (i in 1:length(goodChrOrder)){
      chr=goodChrOrder[1]
      ch=chr
      col.class   <- c(NA, NA, NA,"NULL","NULL","NULL",NA)
      col.names <- c("chr", "start", "stop", "count")
      stem <- fread(paste(stem_dir,chr,".",stem_num,".coverage",sep = ""),col.names = col.names,colClasses =col.class )
      stem<- stem[order(start),]
      gsm <- fread(paste(gsm_dir,gsm_dir_name,"/",gsm_ir_type,"/",gsm_number,"/",chr,".",gsm_number,".coverage",sep = ""),col.names = col.names,colClasses =col.class )
      
      df <- data.frame(position=stem$start, stemloops=stem$count, Dnas=gsm$count)
      # removing 0 points
      aa<-df_zero_ref[which(df_zero_ref$chr==ch),]
      df<-df[!(df$position %in% aa$start),]
      
      gg <- ggplot(df, aes(x=position, y=signal, color=variable))
      gg <- gg + geom_line(aes(y=stemloops/mean(stemloops),col="Stemloops")) 
      gg <- gg + geom_line(aes(y=-Dnas/mean(Dnas)+2.1, col="Dnas") )
      title<-paste(gsm_dir_name,gsm_ir_type,gsm_number,"Chromosome:",chr,"With Stemloops: ",stem_num,sep = "    ")
      gg <- gg + ggtitle(title)
      gg
      # save to pngs file
      # tmp is directory that you want to save files
      png(paste("tmp/",gsm_number,chr,".png",sep = ""),1200,800)
      print(gg)
      dev.off()
}

newX<-df$stemloops/mean(df$stemloops)
newY<- df$Dnas/mean(df$Dnas)
cor(newY[1:5],newX[1:5])

rcorr((df$stemloops),(df$Dnas),type=c("pearson"))

#-------------------------------------------



