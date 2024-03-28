#GOAL: to find the best classifier to classify correctly diabetics and non-diabetics
#Dataset is from Kaggle 

b <- read.csv("diabetes_data.csv", sep=","
              , dec = ".",  
              stringsAsFactors=TRUE, na.strings=c("NA","NaN", "", 'NULL'))


sapply(b, function(x)(sum(is.na(x))))
sapply(b, function(x)(length(unique(x))))
summary(b)
str(b)


library(funModeling)
status <- df_status(b, print_results = F)
status

?as.factor

b$Age <- as.factor(b$Age)
b$Sex <- as.factor(b$Sex)
b$HighChol <- as.factor(b$HighChol)
b$Smoker <- as.factor(b$Smoker)
b$CholCheck <- as.factor(b$CholCheck)
b$HeartDiseaseorAttack <- as.factor(b$HeartDiseaseorAttack)
b$PhysActivity <- as.factor(b$PhysActivity)
b$Fruits <- as.factor(b$Fruits)
b$Veggies <- as.factor(b$Veggies)
b$HvyAlcoholConsump <- as.factor(b$HvyAlcoholConsump)
b$GenHlth <- as.factor(b$GenHlth)
b$DiffWalk <- as.factor(b$DiffWalk)
b$Stroke <- as.factor(b$Stroke)
b$HighBP <- as.factor(b$HighBP)
b$Diabetes <- as.factor(b$Diabetes)

summary(b$BMI)
#ci sono BMI molto alti (molto oltre 40), il massimo è 98. Abbiamo visto che il record di BMI 
#è stato 146 --> è possibile dunque raggiungere un Bmi molto oltre 40 (obesità più grave).

Diabetes <- b$Diabetes
b$diab <- b$Diabetes
levels(b$diab)<-c('c0','c1')
b$Diabetes <- NULL


#vogliamo minimizzare il numero di previsti senza diabete, ma che sono in realtà malati (false negative rate)
#massimizziamo la sensitivity 
#non usiamo la precision perchè minimizzeremmo il numero di previsti malati, ma che non hanno il diabete. 
#sensitivity --> percentuale di riga (mentre precision --> percentuale di colonna)


#SPLIT IN SCORE DATASET AND TRAINING/TEST DATASET
library(caret)
set.seed(107)
scoreIndex <- createDataPartition(y = b$diab, p = .10, list = FALSE)
score <- b[scoreIndex,]
dataset <- b[-scoreIndex,]
#


#model selection 
library(Boruta)  
set.seed(123)
boruta<- Boruta(diab~., data = dataset, doTrace = 1)
plot(boruta, xlab = "features", xaxt = "n", ylab="MDI")
print(boruta)

boruta.metrics <- attStats(boruta)
head(boruta.metrics)
table(boruta.metrics$decision)

#dalla model selection con boruta non emergono variabili non importanti o tentative, ye!

set.seed(107)
Trainindex <- createDataPartition(y = dataset$diab, p = .70, list = FALSE)
train <- dataset[ Trainindex,]
test  <- dataset[-Trainindex,]



# STEP 1 ------------------------------------------------------------------
#Pre processing
#separation, zero variance, dati mancanti e collinearità
#zero variance l'abbiamo fatta dentro 

## LOGISTICO####
library(caret)
set.seed(1234)
metric <- "Spec"
control <- trainControl(method= "cv",number=10, summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE )
glm=train(diab~.,data=train , method = "glm", preProcess=c("corr", "nzv"),metric=metric,
          trControl = control, tuneLength=5, trace=FALSE)
glm #le 6 variabili tolte era per la zero variance
summary(glm)

logPred_train<-predict(glm, train)
confusionMatrix(logPred_train, train$diab)
#specificity->0.7671          

logPred_test <- predict(glm , test)# il target previsto
head(logPred_test)
confusionMatrix(logPred_test, test$diab)
#specificity->0.7679          
#non overfitta

## LASSO####
library(glmnet)
set.seed(1234)
ctrl =trainControl(method="cv", number = 10, classProbs = T,
                   summaryFunction=twoClassSummary, savePredictions = TRUE)
#grid = expand.grid(.alpha=1,.lambda=seq(0, 1, by = 0.01))
metric <- "Spec"
lasso=train(diab~.,
            data=train, method = "glmnet",
            trControl = ctrl, tuneLength=5, na.action = na.pass,metric=metric)
            
lasso

lassoPred_train= predict(lasso,train) #target previsto (no probabilità)
head(lassoPred_train)
confusionMatrix(lassoPred_train, train$diab)
#spec->0.7725

lassoPred_test= predict(lasso,test) #target previsto (no probabilità)
confusionMatrix(lassoPred_test, test$diab)
#spec->0.7728

#non overfitta

#PLS####
#influenza outliers
library(pls)
set.seed(1234)
Control=trainControl(method= "cv",number=10, classProbs=TRUE,
                     summaryFunction=twoClassSummary, savePredictions = TRUE)
metric<-"Spec"
pls=train(diab~. , data=train , method = "pls",metric=metric, 
          trControl = Control, tuneLength=5, preProcess= c("center","corr", "nzv"))
pls



plsPred_train= predict(pls,train) #target previsto (no probabilità)
confusionMatrix(plsPred_train, train$diab)
#spec->0.7473

plsPred_test= predict(pls,test) #target previsto (no probabilità)
confusionMatrix(plsPred_test, test$diab)
#spec->0.7412
#non overfitta

#Naive Bayes####
#scegliamo caret invece di klar perchè non sappiamo tunare la soglia
library(caret)
set.seed(1234)
ctrl =trainControl(method="cv", number = 10, classProbs = T,
                   summaryFunction=twoClassSummary,savePredictions = TRUE)
metric<-"Spec"
naivebayes =train(diab~.,
                 data=train,method = "naive_bayes",metric=metric,preProcess=c("corr"),
                 trControl = ctrl, tuneLength=5, na.action = na.pass) 
naivebayes

NBPred_train= predict(naivebayes,train) #target previsto (no probabilità)
confusionMatrix(NBPred_train, train$diab)
#spec->0.8257

NBPred_test= predict(naivebayes,test) #target previsto (no probabilità)
confusionMatrix(NBPred_test, test$diab)
#spec->0.8307

##Random forest####
library(caret)
set.seed(1234)
metric <- "Spec"   
control <- trainControl(method="cv", number=10, summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
random_forest <- train(diab~., data=train, method="rf", metric=metric, ntree=250, trControl=control)

random_forest
plot(random_forest)

RFPred_train= predict(random_forest,train) #target previsto (no probabilità)
confusionMatrix(RFPred_train, train$diab)
#spec -> 0.7980

RFPred_test= predict(random_forest,test) #target previsto (no probabilità)
confusionMatrix(RFPred_test, test$diab)
#spec --> 0.7868



#Vimportance <- varImp(random_forest)
#head(Vimportance)
#plot(Vimportance)


#Gradient Boost####

set.seed(1234) 
metric <- "Spec"   
control <- trainControl(method="cv", number=10, summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
gradient_boost <- train(diab ~ ., data = train, method = "gbm", trControl = control, metric = metric, verbose=FALSE) 

gradient_boost

GBPred_train= predict(gradient_boost,train) #target previsto (no probabilità)
confusionMatrix(GBPred_train, train$diab)
#spec -> 0.7769 

GBPred_test= predict(gradient_boost,test) #target previsto (no probabilità)
confusionMatrix(GBPred_test, test$diab)
#spec -->  0.7771

set.seed(1234) 
metric <- "Spec"   
control <- trainControl(method="boot", summaryFunction = twoClassSummary, classProbs = TRUE)
gradient_boost2 <- train(diab ~ ., data = train, method = "gbm", trControl = control, metric = metric, verbose=FALSE) 

gradient_boost2

GBPred_train2= predict(gradient_boost2,train) #target previsto (no probabilità)
confusionMatrix(GBPred_train2, train$diab)
#spec -> 0.7828 

GBPred_test2= predict(gradient_boost2,test) #target previsto (no probabilità)
confusionMatrix(GBPred_test2, test$diab)
# spec --> 0.7805

#Neural Netwrok####

#non proviamo con il KNN e QDA/LDA perchè le covariate di input sono in gran parte 
#categoriali

library(caret)
set.seed(1234)
metric <- "Spec"
ctrl = trainControl(method="cv", number=10, summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)  #, search = "grid")
nnet <- train(diab ~ ., data = train, 
              method = "nnet",
              preProcess = c('corr', 'nzv',"range"),
              metric=metric, trControl=ctrl,
              trace = TRUE, # use true to see convergence
              maxit = 1000)
nnet

NNETpred_train= predict(nnet,train) #target previsto (no probabilità)
confusionMatrix(NNETpred_train, train$diab)
#spec-> 0.8150

NNETpred_test= predict(nnet,test) #target previsto (no probabilità)
confusionMatrix(NNETpred_test, test$diab)
# spec -->  0.8135


#stacking####
library(caret)
library(caretEnsemble)
set.seed(1234)
modelli <- list("rf"=random_forest,"glm"=glm, "lasso"=lasso, 'pls' = pls, 'nb'=naivebayes, 'gb'=gradient_boost, 'nnet'=nnet) 
class(modelli) <- "caretList" 
control <- trainControl(method="cv", number=10, summaryFunction = twoClassSummary, classProbs = TRUE,savePredictions = TRUE)
stackglm<-caretStack(modelli,method='glm',metric="Spec", trControl=control) 

stackglm


stackPred_train= predict(stackglm,train) #target previsto (no probabilità)
confusionMatrix(stackPred_train, train$diab)
#spec-> 0.7756

stackPred_test= predict(stackglm,test) #target previsto (no probabilità)
confusionMatrix(stackPred_test, test$diab)
# spec --> 0.7813

test$target_stack <- stackPred_test


results <- resamples(list(glm=glm, lasso=lasso, naivebayes=naivebayes, random_forest=random_forest, pls=pls, gradient_boost=gradient_boost, nnet=nnet))
# Riassumiamo le distribuzioni delle metriche sulle folds
summary(results)
# boxplots dei risultati
bwplot(results)

#ROC####
test$glm=predict(glm,test, "prob")[,1]
test$lasso=predict(lasso,test, "prob")[,1]
test$naivebayes=predict(naivebayes,test, "prob")[,1]
test$pls=predict(pls,test, "prob")[,1]
test$random_forest=predict(random_forest,test, "prob")[,1]
test$gradient_boost=predict(gradient_boost,test, "prob")[,1]
test$nnet=predict(nnet,test, "prob")[,1]
test$stackglm <- predict(stackglm, newdata=test, type="prob")
test$stackglm

library(pROC)
roc.glm=roc(diab ~ glm, data = test)
roc.gradient_boost=roc(diab ~ gradient_boost, data = test)
roc.lasso=roc(diab ~ lasso, data = test)
roc.naivebayes=roc(diab ~ naivebayes, data = test)
roc.pls=roc(diab ~ pls, data = test)
roc.random_forest=roc(diab ~ random_forest, data = test)
roc.nnet=roc(diab ~ nnet, data = test)
roc.stackglm=roc(diab ~ stackglm, data = test)

class(stackglm)
class(pls)
train(stackglm)
roc.glm
roc.gradient_boost
roc.lasso
roc.naivebayes
roc.pls
roc.random_forest
roc.nnet
roc.stackglm

plot(roc.glm)
plot(roc.gradient_boost,add=T,col="red")
plot(roc.lasso,add=T,col="blue")
plot(roc.naivebayes,add=T,col="yellow")
plot(roc.pls,add=T,col="pink")
plot(roc.random_forest,add=T,col="green")
plot(roc.nnet,add=T,col="turquoise2")
plot(roc.stackglm,add=T,col="grey")
# for ROC smaller models has better ROC

#table(b$diab)/70692
#
rho1=0.5
rho0=0.5
#
true1 =0.096
true0 =0.904

copy = test
#
#posterior GLM####

test$glm_adj <- (test$glm*true0/rho0)/((test$glm*true0/rho0)+((1-test$glm)*true1/rho1))

library(funModeling)
gain_lift(data = test, score = 'glm_adj', target = 'diab')

#posterior pls####
test$pls_adj <- (test$pls*true0/rho0)/((test$pls*true0/rho0)+((1-test$pls)*true1/rho1))
gain_lift(data = test, score = 'pls_adj', target = 'diab')

#posterior lasso####
test$lasso_adj <- (test$lasso*true0/rho0)/((test$lasso*true0/rho0)+((1-test$lasso)*true1/rho1))
gain_lift(data = test, score = 'lasso_adj', target = 'diab')

#posterior random forest####
test$rf_adj <- (test$random_forest*true0/rho0)/((test$random_forest*true0/rho0)+((1-test$random_forest)*true1/rho1))
gain_lift(data = test, score = 'rf_adj', target = 'diab')

#posterior nnet####
test$nnet_adj <- (test$nnet*true0/rho0)/((test$nnet*true0/rho0)+((1-test$nnet)*true1/rho1))
gain_lift(data = test, score = 'nnet_adj', target = 'diab')

#posterior stacking#####

test$stacking_adj <- (test$stackglm*true0/rho0)/((test$stackglm*true0/rho0)+((1-test$stackglm)*true1/rho1))
gain_lift(data = test, score = 'stacking_adj', target = 'diab')

#posterior gb####
test$gb_adj <- (test$gradient_boost*true0/rho0)/((test$gradient_boost*true0/rho0)+((1-test$gradient_boost)*true1/rho1))
gain_lift(data = test, score = 'gb_adj', target = 'diab')

#soglia####

# save validation results: observed target an pred probs
df=test
df$Class=test$diab
head(df)
df$Class=ifelse(df$Class=="c0","M","R") # il nostro event c1 ? ora M
head(df)
df$ProbM=test$stacking_adj

library(dplyr)
# for each threshold, find tp, tn, fp, fn and the sens=prop_true_M, spec=prop_true_R, precision=tp/(tp+fp)

thresholds <- seq(from = 0, to = 1, by = 0.01)
prop_table <- data.frame(threshold = thresholds, prop_true_M = NA,  prop_true_R = NA, true_M = NA,  true_R = NA ,fn_M=NA)

for (threshold in thresholds) {
  pred <- ifelse(df$ProbM > threshold, "M", "R")  # be careful here!!!
  pred_t <- ifelse(pred == df$Class, TRUE, FALSE)
  
  group <- data.frame(df, "pred" = pred_t) %>%
    group_by(Class, pred) %>%
    dplyr::summarise(n = n())
  
  group_M <- filter(group, Class == "M")
  
  true_M=sum(filter(group_M, pred == TRUE)$n)
  prop_M <- sum(filter(group_M, pred == TRUE)$n) / sum(group_M$n)
  
  prop_table[prop_table$threshold == threshold, "prop_true_M"] <- prop_M
  prop_table[prop_table$threshold == threshold, "true_M"] <- true_M
  
  fn_M=sum(filter(group_M, pred == FALSE)$n)
  # true M predicted as R
  prop_table[prop_table$threshold == threshold, "fn_M"] <- fn_M
  
  
  group_R <- filter(group, Class == "R")
  
  true_R=sum(filter(group_R, pred == TRUE)$n)
  prop_R <- sum(filter(group_R, pred == TRUE)$n) / sum(group_R$n)
  
  prop_table[prop_table$threshold == threshold, "prop_true_R"] <- prop_R
  prop_table[prop_table$threshold == threshold, "true_R"] <- true_R
  
}

head(prop_table, n=10)


# prop_true_M = sensitivity
# prop_true_R = specificicy





# n of observations of the validation set    
prop_table$n=nrow(df)

# false positive (fp_M) by difference of   n and            tn,                 tp,         fn, 
prop_table$fp_M=nrow(df)-prop_table$true_R-prop_table$true_M-prop_table$fn_M

# find accuracy
prop_table$acc=(prop_table$true_R+prop_table$true_M)/nrow(df)

# find precision
prop_table$prec_M=prop_table$true_M/(prop_table$true_M+prop_table$fp_M)

# find F1 =2*(prec*sens)/(prec+sens)
# prop_true_M = sensitivity

prop_table$F1=2*(prop_table$prop_true_M*prop_table$prec_M)/(prop_table$prop_true_M+prop_table$prec_M)

# verify not having NA metrics at start or end of data 
tail(prop_table)
head(prop_table)
# we have typically some NA in the precision and F1 at the boundary..put,impute 1,0 respectively 

library(Hmisc)
#impute NA as 0, this occurs typically for precision
prop_table$prec_M=impute(prop_table$prec_M, 1)
prop_table$F1=impute(prop_table$F1, 0)
tail(prop_table)
colnames(prop_table)
colnames(prop_table) <- c("threshold", "sensitivity", "specificity", "true_M", "true_R", "fn_M", "n", 
                          "fp_M", "acc" ,        "precision" ,     "F1")
head(prop_table)

# drop counts, PLOT only metrics
prop_table2 = prop_table[,-c(4:8)] 
head(prop_table2)

# plot measures vs soglia
# before we must impile data vertically: one block for each measure
library(dplyr)
library(tidyr)
gathered=prop_table2 %>%
  gather(x, y, sensitivity:F1)
head(gathered)


# grafico con tutte le misure 
library(ggplot2)
gathered %>%
  ggplot(aes(x = threshold, y = y, color = x)) +
  geom_point() +
  geom_line() +
  scale_color_brewer(palette = "Set1") +
  labs(y = "measures",
       color = "M: event\nR: nonevent")



# follow sensitivity= prop true M...beccre i veri ricchi (soglie basse) 
# anche F1 ciferma soglie attorno a  0.20
prop_table2[prop_table2$specificity >= 0.85, ]
pred <- ifelse(test$stacking_adj > 0.94, "c0", "c1")
pred <- as.factor(pred)
table(actual=test$diab,pred)

confusionMatrix(pred, test$diab)



# step 4 ------------------------------------------------------------------

posterior_diab= predict(stackglm, score, "prob")
head(posterior_diab)
