

#----------------------------------------------------
# Building Model using Xgboost 


library(xgboost)
library(readr)
library(stringr)
library(caret)
library(car)
setDT(train_data)
setDT(test_data)

labels <- train_data$Stem6
ts_label <- test_data$Stem6
new_tr <- model.matrix(~.+0,data = train_data[,-c("Stem6"),with=F])
new_ts <- model.matrix(~.+0,data = test_data[,-c("Stem6"),with=F])

#preparing matrix
dtrain <- xgb.DMatrix(data = new_tr,label = labels)
dtest <- xgb.DMatrix(data = new_ts,label=ts_label)

#default parameters
params <- list(booster = "gbtree", objective = "reg:linear", eta=0.3,
               gamma=0, max_depth=6, min_child_weight=1, subsample=1, colsample_bytree=1)

xgbcv <- xgb.cv( params = params, data = dtrain, nrounds = 100, nfold = 5, showsd = T, stratified = T, 
                 print_every_n = 10, early_stop_round = 20, maximize = F) ##best iteration = 79

min(xgbcv$evaluation_log$test_rmse_mean)


xgb1 <- xgb.train (params = params, data = dtrain, nrounds = 79, watchlist = list(val=dtest,train=dtrain), 
                   print_every_n = 10, early_stop_round = 10, maximize = F , eval_metric = "error")
#model prediction
xgbpred <- predict (xgb1,dtest)
xgbpred <- ifelse (xgbpred > 0.5,1,0)

#confusion matrix
library(caret)
confusionMatrix(xgbpred, ts_label )

mat <- xgb.importance (feature_names = colnames(new_tr),model = xgb1)
par(mar=c(.1,.1,.1,.1))

xgb.plot.importance(mat[1:10,])


dev.off()

#----------------------------------------------
lmMod <- lm(Stem6 ~. , data=train_data)
distPred <- predict(lmMod, test_data)  # predict distance
summary(lmMod)
actuals_preds <- data.frame(cbind(actuals=test_data$Stem6, predicteds=distPred))  # make actuals_predicteds dataframe.
correlation_accuracy <- cor(actuals_preds)  # 
# caret

library(DAAG)
cvResults <- suppressWarnings(CVlm(data=tmp_df[,4:151], form.lm=Stem6 ~ . , m=10, dots=FALSE, seed=29, 
                                   legend.pos="topleft",  printit=FALSE, 
                                   main="Small symbols are predicted values while bigger ones are actuals."));  # performs the CV
attr(cvResults, 'ms') 

#----------------------------------------------





rbind(data.frame(group = "train", train_data),
      data.frame(group = "test", test_data)) %>%
  gather(x, y, GSM1024605:GSM1027308) %>%
  ggplot(aes(x = y, color = group, fill = group)) +
  geom_density(alpha = 0.3) +
  facet_wrap( ~ x, scales = "free", ncol = 3)
# regression 



set.seed(42)
model_rf <- caret::train(Stem6 ~ .,
                         data = train_data[,4:301],
                         method = "rf",
                         importance=TRUE,
                         preProcess = c("scale", "center"),
                         trControl = trainControl(method = "repeatedcv",
                                                  number = 10,
                                                  repeats = 10,
                                                  savePredictions = TRUE,
                                                  verboseIter = FALSE))

model_rf$finalModel$importance
pr_rf<-predict(model_rf, test_data)



act_pred <- data.frame(cbind(actuals=test_data$Stem6, predicteds=pr_rf))  # make actuals_predicteds dataframe.

cor_acc <- cor(act_pred)  # 

min_max_accuracy <- mean(apply(act_pred, 1, min) / apply(act_pred, 1, max))  
mape <- mean(abs((act_pred$predicteds - act_pred$actuals))/act_pred$actuals) 

imp <- model_rf$finalModel$importance
imp[order(imp, decreasing = TRUE), ]


importance <- varImp(model_rf, scale = TRUE)
plot(importance)
#----------
df <- fread("By_sample_finalTable.csv")



tmp_df<-df[,2470:2482]

k <- 27
kmeans_result <- kmeans(tmp_df, k)
kmeans_result
#table(data.frame(tmp_df[,c(Stem15)], kmeans_result$cluster))

table(kmeans_result$cluster)


plot(tmp_df, col = kmeans_result$cluster, main= "(A) Plot with clusters")
dev.off()

distance <- dist(tmp_df, method="euclidean")
cluster <- hclust(distance, method="average")
plot(cluster, hang=-1, label=tmp_df$Stem15)



# ------------- detect clusters

mydata <- tmp_df
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:30) wss[i] <- sum(kmeans(mydata,
                                     centers=i)$withinss)
plot(1:30, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares",
     main="Assessing the Optimal Number of Clusters with the Elbow Method",
     pch=20, cex=2)






#------------------
df <- fread("By_sample_finalTable.csv")

tmp_df<-df[,4:2481]

cols <- names(tmp_df)
tmp_df[,(cols) := round(.SD,3), .SDcols=cols]
#setDF(tmp_df)

k <- 25
kmeans_result <- kmeans(tmp_df, k)
kmeans_result

# Enhanced hierarchical clustering, cut in 3 groups
library("factoextra")
res.hc <- eclust(tmp_df, "hclust", k = 26, graph = FALSE) 
# Visualize
fviz_dend(res.hc, rect = TRUE, show_labels = FALSE)


