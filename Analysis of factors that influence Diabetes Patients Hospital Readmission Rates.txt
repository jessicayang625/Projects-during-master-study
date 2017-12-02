diabeteshs<-read.csv("E:/diabetic_data.csv")
diabetesh<-diabeteshs[,-c(1,2,6,11:12,19:22,28,30,33,37:41,43:47)]
names(diabetesh)
dim(diabetesh)
##let readmitted be 0 and 1
diabetesh$readmitted<-ifelse(diabeteshs$readmitted=="<30",'1','0')
diabetesh[, 28]<- sapply(diabetesh[, 28], as.numeric)
summary(diabetesh)
##check whether it's numeric
class(diabetesh$readmitted) 
##fill the missing values in one column with NA
diabetesh$race[diabeteshs$race=='?'] <-NA
diabetems<-na.omit(diabetesh) ###delete missing values
dim(diabetems)
##random get 2000 values to do the analysis
set.seed(1)
diabete.sub<-diabetems[sample(nrow(diabetems), 2000),]
dim(diabete.sub)
summary (diabete.sub)

#random get 1500 observations of training data and 500 test data
set.seed(123)
train=sample(1:nrow(diabete.sub),nrow(diabete.sub)*0.75)
test<--train
diabete.train<-diabete.sub[train,]
diabete.test<-diabete.sub[test,]
View(diabete.train)


### logistic regression ###
library (ISLR)
library(leaps)
attach(diabete.train)
detach(diabete.train)
glm.fit<-glm(readmitted ~ .,data=diabete.train,family=binomial)
summary(glm.fit)
glm.step<-step(glm.fit)
summary(glm.step)
glm.fita1c<-glm(readmitted ~admission_type_id + discharge_disposition_id + 
                  time_in_hospital + num_procedures + number_emergency + number_inpatient + 
                  max_glu_serum + repaglinide + insulin,data=diabete.train,family=binomial)
summary(glm.fita1c)
Summary(diabete.train)
glm.probs <- predict(glm.fita1c,type ="response")
glm.probs[1:1500]
contrasts(factor(readmitted))
head(diabete.train)
glm.pred=rep ("0",1500)
glm.pred[glm.probs >.5]="1"
table(glm.pred,readmitted )
(1324+2) /1500
mean(glm.pred== readmitted )
#predict with test data
#glm.probsts <- predict (glm.fita1c,diabete.testrm, type="response")
glm.fita1cts<-glm(readmitted ~time_in_hospital + number_emergency + number_inpatient + metformin + nateglinide + glimepiride +A1Cresult,data=diabete.test,family =binomial)
glm.probsts<-predict (glm.fita1cts,diabete.test , type="response")
readmitted.ts= diabete.test[,28]
glm.predts<-rep ("0",500)
glm.predts[glm.probsts >.5]="1"
table(glm.predts,readmitted.ts)
mean(glm.predts== readmitted.ts)
mean(glm.predts!=readmitted.ts)


#### randomForest #####
diabete.test11<-diabete.test$readmitted
dim(diabete.test)
library (randomForest)
set.seed (111)
diabete.rf =randomForest(readmitted~.,mtry = 9,importance = TRUE ,data=diabete.train) #mtry=p/3
diabete.rf
yhat.rf = predict (diabete.rf,newdata =diabete.test)
mean((yhat.rf-diabete.test11)^2)
importance(diabete.rf)
varImpPlot(diabete.rf)


