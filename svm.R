# SVM method

# loading libraries 
source("load_lib.R")
current_directory<-getwd()
options(scipen=999)

df <- fread("By_experiment_finalTable.csv")

#--------------------
# SVM
#-------------------



#chr13-16 and chr20-22 and then chr1 and chr9

#ch="chr1"
tmp_df<-df[,1:352]
#tmp_df$Stem15<-df[,353]
#tmp_df$Stem16<-df[,354]
tmp_df$Stem6<-df[,355]

#tmp_df<-tmp_df[which(tmp_df$chr==ch),]

tmp_df<-tmp_df[,4:353]
#cols<-c(colnames(tmp_df))


# set.seed(42)
# index <- createDataPartition(tmp_df$Stem6, p = 0.7, list = FALSE)
# train_data <- tmp_df[index, ]
# test_data <- tmp_df[-index, ]
# mtry_def <- floor(sqrt(ncol(train_data))*.80)
# 
# model <- randomForest(Stem6 ~ . , data=train_data,importance=TRUE, ntree=1000,mtry=mtry_def)
# pred<-predict(model, test_data)
# 
# print(model)
# importance(model)
# plot(model)
# plot( importance(model), lty=2, pch=16)
# lines(importance(model))
# imp = importance(model)
# impvar = rownames(imp)[order(imp[, 1], decreasing=TRUE)]
# 
# cat("% Var explained: \n", 100* (1-sum((model$y-model$pred   )^2) /
#                                    sum((model$y-mean(model$y))^2)
# )
# )
# plot(model$y, model$predicted)
# 
# ##pearson correlation R²(pearson)
#cat("% Pearson cor: \n ", 100*cor(model$y,model$predicted)^2)
# ##spearman correlation R²(spearman)
#cat("% spearman cor: \n ", 100*cor(model$y,model$predicted,method="s")^2)
# 
# qqnorm((pred - test_data$Stem6)/sd(pred - test_data$Stem6))
# 
# qqline((pred-test_data$Stem6)/sd(pred - test_data$Stem6))
# 
# plot(test_data$Stem6, pred)
# # The root-mean-square deviation 
#RMSE <- (sum((pred-test_data$Stem6)^2)/length(test_data$Stem6))^(1/2)
#RMSE

#----------------------------------------------------
# using caret 

# configure multicore
cl <- makeCluster(detectCores())
registerDoParallel(cl)

set.seed(42)
index <- createDataPartition(tmp_df$Stem6, p = 0.8, list = FALSE)
train_data <- tmp_df[index, ]
test_data <- tmp_df[-index, ]

mtry_def <- floor(sqrt(ncol(train_data))*.80)

tune_grid <- expand.grid(mtry= c(mtry_def))

ptm <- proc.time()
set.seed(1234)
ml.svm <- train(Stem6 ~ .,
                data = train_data,
                preProc = c("center", "scale"),
                method = "svmLinear")

proc.time() - ptm
print(ml.svm)
ml.svm$results
summary(ml.svm)
#plot(model.rf$finalModel)
pred<-predict(ml.svm, test_data)

#------------- R square 
RMSE <- sqrt(sum((pred - test_data$Stem6)^2)/length(pred))
print(RMSE)
print(RMSE/mean(test_data$Stem6)) 
# 
postResample(pred = pred, obs = test_data$Stem6)



actual<-test_data$Stem6
result<-data.frame(actual=actual,predicted=pred)
#theme_set(theme_grey())

ggplot(result,aes(x=actual,y=predicted,color=predicted-actual),alpha=0.7)+
  geom_point()+
  geom_smooth(method="auto")+
  ggtitle('Predicted SVM Model')+
  labs(x="Actual",y="Predicted")


ImpMeasure<-varImp(ml.svm, scale = FALSE)
plot(ImpMeasure,top = 20 )



