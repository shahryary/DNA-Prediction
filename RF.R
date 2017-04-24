# loading libraries 
source("load_lib.R")
current_directory<-getwd()
options(scipen=999)

df <- fread("By_experiment_finalTable.csv")
#--------------------
# Random Forest 
#-------------------

# configure multicore
cl <- makeCluster(detectCores())
registerDoParallel(cl)

#chr13-16 and chr20-22 and then chr1 and chr9

ch="chr13"
tmp_df<-df[,1:352]
#tmp_df$Stem15<-df[,353]
#tmp_df$Stem16<-df[,354]
tmp_df$Stem6<-df[,355]

tmp_df<-tmp_df[which(tmp_df$chr==ch),]

tmp_df<-tmp_df[,4:353]
cols<-c(colnames(tmp_df))


set.seed(42)
index <- createDataPartition(tmp_df$Stem6, p = 0.7, list = FALSE)
train_data <- tmp_df[index, ]
test_data <- tmp_df[-index, ]
mtry_def <- floor(sqrt(ncol(train_data))*.70)

model <- randomForest(Stem6 ~ . , data=train_data,importance=TRUE, ntree=1000,mtry=mtry_def)
pred<-predict(model, test_data)

print(model)
importance(model)
plot(model)
plot( importance(model), lty=2, pch=16)
lines(importance(model))
imp = importance(model)
impvar = rownames(imp)[order(imp[, 1], decreasing=TRUE)]

cat("% Var explained: \n", 100* (1-sum((model$y-model$pred   )^2) /
                                   sum((model$y-mean(model$y))^2)
)
)
plot(model$y, model$predicted)

##pearson correlation R²(pearson)
cat("% Pearson cor: \n ", 100*cor(model$y,model$predicted)^2)
##spearman correlation R²(spearman)
cat("% spearman cor: \n ", 100*cor(model$y,model$predicted,method="s")^2)

qqnorm((pred - test_data$Stem6)/sd(pred - test_data$Stem6))

qqline((pred-test_data$Stem6)/sd(pred - test_data$Stem6))

plot(test_data$Stem6, pred)
# The root-mean-square deviation 
RMSE <- (sum((pred-test_data$Stem6)^2)/length(test_data$Stem6))^(1/2)
RMSE

#----------------------------------------------------
# using caret 
set.seed(42)
index <- createDataPartition(tmp_df$Stem6, p = 0.7, list = FALSE)
train_data <- tmp_df[index, ]
test_data <- tmp_df[-index, ]

mtry_def <- floor(sqrt(ncol(train_data))*.70)

tune_grid <- expand.grid(mtry= c(mtry_def))

ptm <- proc.time()
set.seed(1234)
model.rf <- train(Stem6 ~ .,
                  data = train_data,
                  method = "rf",
                  ntree = 4000, # How many trees to grow in total?
                  tuneGrid = tune_grid)#
proc.time() - ptm

print(model.rf)
plot(model.rf$finalModel)
pred<-predict(model.rf, test_data)
RMSE<- sqrt(sum((pred - test_data$Stem6)^2)/length(pred))
RMSE
print(RMSE/mean(test_data$Stem6))