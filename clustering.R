# clustering by k-mean method
source("load_lib.R")
options(scipen=999)

df <- fread("By_sample_finalTable.csv")

tmp_df<-df[,4:2481]

# rounfding numbers 3 digit 
cols <- names(tmp_df)
tmp_df[,(cols) := round(.SD,3), .SDcols=cols]


# detect number of cluster 
mydata <- tmp_df
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:60) wss[i] <- sum(kmeans(mydata,
                                     centers=i)$withinss)
plot(1:60, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares",
     main="Assessing the Optimal Number of Clusters with the Elbow Method",
     pch=20, cex=2)


# ------ Kmeans
tmp_df<-df[,4:10]

k <- 5
kmeans_result <- kmeans(tmp_df, k)
kmeans_result
dd <- dist(tmp_df, method ="euclidean")
# Statistics for k-means clustering
km_stats <- cluster.stats(dd,  kmeans_result$cluster)
# (k-means) within clusters sum of squares
km_stats$within.cluster.ss


# ---- The result of this part in saved as  "cl_dn.png" - running code it will be take about 40 minute in server!
# Enhanced hierarchical clustering, cut in 3 groups
#library("factoextra")
#res.hc <- eclust(tmp_df, "hclust", k = 26, graph = FALSE) 
# Visualize
# ~ 20 minutes to plot - I exported this result
#fviz_dend(res.hc, rect = TRUE, show_labels = FALSE)

library(vegan)
groups <- levels(factor(kmeans_result$cluster))
ordiplot(cmdscale(dist(tmp_df)), type = "n", display = 'sp') # might issue 
cols <- rainbow(nlevels(factor(kmeans_result$cluster)))
for(i in seq_along(groups)) {
  points(cmdscale(dist(tmp_df))[factor(kmeans_result$cluster) ==
                                  groups[i], ], col = cols[i], pch = 16)
}
ordispider(cmdscale(dist(tmp_df)), factor(kmeans_result$cluster), label = TRUE,cex=.3)
ordihull(cmdscale(dist(tmp_df)), factor(kmeans_result$cluster),lty = "dotted")

