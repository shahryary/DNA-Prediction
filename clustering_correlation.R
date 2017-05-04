# last update: 02.May.2017
# clustering by k-mean method   

# clustering CORRELATION files!

source("load_lib.R")
options(scipen=999)

df <- fread("correlation_by_sample_Stem15.csv")
df[,2:NCOL(df)]<-df[,2:NCOL(df)] %>% mutate_each(funs(round(.,3)), starts_with("GSM"))


# replacing Rows/Columns
tmp_df<-data.frame(t(df))
colnames(tmp_df) = as.character(unlist(tmp_df[1,]))
tmp_df = tmp_df[-1, ]          # removing the first row.


# detect number of cluster 
mydata <- tmp_df
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 1:23) wss[i] <- sum(kmeans(mydata,
                                     centers=i)$withinss)
plot(1:23, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares",
     main="Assessing the Optimal Number of Clusters with the Elbow Method",
     pch=20, cex=2)
# 

k <- 13
kmeans_result <- kmeans(tmp_df, k)
#tmp_df$cluster<-kmeans_result$cluster

library(fpc)
dd <- dist(tmp_df, method ="euclidean")
# Statistics for k-means clustering
km_stats <- cluster.stats(dd,  kmeans_result$cluster)
# (k-means) within clusters sum of squares
km_stats$within.cluster.ss

# Enhanced hierarchical clustering, cut in 3 groups
library("factoextra")
res.hc <- eclust(tmp_df, "hclust", k = 23, graph = FALSE) 
# Visualize
fviz_dend(res.hc, rect = TRUE, show_labels = FALSE)

#install.packages("vegan")
library(vegan)
groups <- levels(factor(kmeans_result$cluster))
ordiplot(cmdscale(dist(tmp_df)), type = "n", display = 'sp') # might issue 
cols <- rainbow(nlevels(factor(kmeans_result$cluster)))
for(i in seq_along(groups)) {
  points(cmdscale(dist(tmp_df))[factor(kmeans_result$cluster) ==
                                  groups[i], ], col = cols[i], pch = 16)
}
ordispider(cmdscale(dist(tmp_df)), factor(kmeans_result$cluster), label = FALSE,cex=.8)
ordihull(cmdscale(dist(tmp_df)), factor(kmeans_result$cluster),lty = "dotted")

tmp_df$cluster<-kmeans_result$cluster
fwrite(tmp_df,file = "By_samples_cluster_Stem15.csv", sep = "\t",col.names = TRUE,row.names = TRUE )
x<-fread(file = "By_samples_cluster_Stem15.csv")


sDes <- fread("SampleDescription.csv",header = FALSE,col.names = c("tissue","type","gsm") )
tmp_clus<-x[,c(1,26)]
cv<-merge(sDes, tmp_clus, by.x=c("gsm"), by.y=c("V1")) 

fwrite(cv,file = "cv_15.csv", sep = "\t",col.names = TRUE )

