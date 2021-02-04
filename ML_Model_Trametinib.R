library(caret)
library(tidyverse)
library(ggplot2)
library(ggpubr)

##################################################################

#Set Working Directory
setwd("C:/Users/Pedro/Documents/Work/work_from_home/Book Chapter") #CHANGE PATH
list.files()

#Load Files
Kinase_Activity <- read.csv("Kinase_activity_PDT.csv",h=T) #(loads kinase activity dataset)
Drug_Response<- read.csv("Drug_response.csv",h=T) #(loads drug response dataset)

#Centre and Scale Kinase Activity Dataset
row.names(Kinase_Activity)<-Kinase_Activity$ids 
Kinase_Activity$ids<-NULL
ML_test_no <- preProcess(Kinase_Activity, method = c("zv","center", "scale"))
Kinase_Activity_Processed <- predict(ML_test_no, Kinase_Activity)
Kinase_Activity_Processed$ids<-row.names(Kinase_Activity_Processed)

#Link Kinase Activity data to Drug Response to Trametinib
Trametinib_Response<- Drug_Response[,c(1,6)]%>%as.data.frame()
colnames(Trametinib_Response)<-(c("ids","Trametinib_1uM"))
Kinase_Activity_Trametinib<-merge(Kinase_Activity_Processed,Trametinib_Response)
row.names(Kinase_Activity_Trametinib)<-Kinase_Activity_Trametinib$ids
Kinase_Activity_Trametinib$ids <- NULL

#Create Partition of Train Set and Test Set
set.seed(41) #MEKi
index <- createDataPartition(Kinase_Activity_Trametinib$Trametinib_1uM, p=0.70, list=FALSE)
trainSet <- Kinase_Activity_Trametinib[ index,]
testSet <- Kinase_Activity_Trametinib[-index,]

#Generate PLS Model for Feature Selection
fit_control <- trainControl(
  method = "LOOCV",
  number = 10)
PLS_fit <- train(Trametinib_1uM ~ ., 
                data = trainSet, 
                metric = "RMSE",
                method = "pls",
                tuneLength = 7,
                maximize=TRUE,
                trControl = fit_control,
                linout = TRUE)

#Select Features with Importance > 50 in PLS Model
varImp_PLS <- varImp(PLS_fit)
Importance<-varImp_PLS[1]%>%as.data.frame()

  #Importance$kinase<-rownames(Importance)
  #Importance_1 <-Importance %>% top_n(15, Overall)
  #rownames(Importance_1)<-Importance_1$kinase

Importance_1 <- subset (Importance, Overall>50)
Importance_1$kinase<-rownames(Importance_1)

#Plot Features with Importance > 50 in PLS Model
Barplot_Importance_PLS<-ggplot(Importance_1, aes(x=reorder(kinase, -Overall), y=Overall))+
  geom_bar(stat="identity", color="black", position=position_dodge())
Barplot_Importance_PLS

#Fillter Features with Importance > 50 in PLS Model in Train and Test set for RF Model
Subset_K<-as.vector(rownames(Importance_1))
Subset_K<-c(Subset_K,"ids","Trametinib_1uM")
trainSet$ids<-rownames(trainSet)
trainSet<- subset(trainSet,select=Subset_K)
trainSet$ids<-NULL
testSet$ids<-rownames(testSet)
testSet<- subset(testSet,select=Subset_K)
testSet$ids<-NULL

#Generate RF Model
fit_control <- trainControl(## 10-fold CV
  method = "LOOCV",
  number = 20)
rf_fit <- train(Trametinib_1uM ~ ., 
                data = trainSet, 
                metric = "RMSE",
                method = "rf",
                tuneLength = 7,
                maximize=TRUE,
                trControl = fit_control,
                ntree=1000,
                linout = TRUE)

#Summary of the Generated RF Model
rf_fit
rf_fit$finalModel
rf_fit$results

#Generate Datasets of Predicted and Measured Values for Train and Test Sets
Predictors_train<- predict(rf_fit, trainSet)
Predictors_train<-as.data.frame(Predictors_train)
Predictors_test <- predict(rf_fit, testSet)
Predictors_test<-as.data.frame(Predictors_test)
Obserbed_train<-trainSet$Trametinib_1uM
Obserbed_train<-as.data.frame(Obserbed_train)
Obserbed_test<-testSet$Trametinib_1uM
Obserbed_test<-as.data.frame(Obserbed_test)
trainplot<- cbind(Obserbed_train, Predictors_train)
testplot<- cbind(Obserbed_test, Predictors_test)

#Calculate RMSE for Train and Test Sets
Train_RMSE<-RMSE(trainplot$Predictors_train, trainplot$Obserbed_train ,na.rm = T)
Test_RMSE<-RMSE(testplot$Predictors_test, testplot$Obserbed_test, na.rm = T)
Train_RMSE
Test_RMSE
 
#Plot Predicted vs Measured drug responses in train set
dotplot_train<- ggscatter(trainplot, x = "Obserbed_train", y = "Predictors_train",   fill = "red",  size= 5, 
                          add = "reg.line", conf.int = TRUE, ellipse = FALSE,
                          cor.coef = TRUE, cor.method = "spearman",add.params = list(color="grey", cor.coef.size=20),
                          xlab = "Measured", ylab = "Predicted")+
  geom_point(pch=16, colour="black", size=5,)+
  theme(legend.position = "none")
dotplot_train

#Plot Predicted vs Measured drug responses in test set
dotplot_test<- ggscatter(testplot, x = "Obserbed_test", y = "Predictors_test",   fill = "red",  size= 5, 
                         add = "reg.line", conf.int = TRUE, ellipse = FALSE,
                         cor.coef = TRUE, cor.method = "spearman",add.params = list(color="grey", cor.coef.size=20),
                         xlab = "Measured", ylab = "Predicted")+
  scale_color_manual(values=c("red",'brown','blue', "cyan", "orange","green", "purple",
                              "violet",  "yellow","green","limegreen","goldenrod","lawngreen"))+
  geom_point(pch=16, aes(color=row.names(testplot)), size=5,  )+
  theme(legend.position = "right")+
  theme(legend.title = element_blank())
dotplot_test

#Plot Feature Importance for Features with Importance >20 in RF Model
varImp_rf <- varImp(rf_fit)
Importance<-varImp_rf[1]%>%as.data.frame()
Importance_1 <- subset (Importance, Overall>20)
Importance_1$kinase<-rownames(Importance_1)
Bartplot_Importance_RF<-ggplot(Importance_1, aes(x=reorder(kinase, -Overall), y=Overall))+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()+
  theme(strip.text.x = element_text(size=18))+
  theme(strip.background = element_rect( fill="white", size=1))+
  theme(axis.title.x= element_text(size=18))+
  theme(axis.text.x=element_text(size=16, colour = "black", angle = 90, hjust=0.95,vjust=0.2))+
  theme(axis.text.y=element_text(size=16, colour = "black",  hjust=0.95,vjust=0.2))+
  theme(axis.line = element_line(color="black"))
Bartplot_Importance_RF

#Select Features with Importance >20 in RF Model
dfImportance_1<-Importance_1$kinase

#Generate Summary of Relevant Parameters
Lengh_sub<-length(Subset_K)
Spearvalue <- cor.test(testplot$Predictors_test, testplot$Obserbed_test, method ="spearman")
Spearvalue_B<-Spearvalue[3]
Spearvalue_B<-as.data.frame(Spearvalue_B)
Spearvalue_A<-Spearvalue[4]
Spearvalue_A<-as.data.frame(Spearvalue_A)
dfcolspearman <- cbind(Spearvalue_A,Spearvalue_B,Lengh_sub, Test_RMSE)
dfresults1<- dfcolspearman
A<-"replicate"
rownames(dfresults1)<-A



###################################################################################################################################


#Generate Replicates of the Workflow

for(i in 1:49)
   {i<-A    
   
  #Create Partition
    index <- createDataPartition(Kinase_Activity_Trametinib$Trametinib_1uM, p=0.70, list=FALSE)
    trainSet <- Kinase_Activity_Trametinib[ index,]
    testSet <- Kinase_Activity_Trametinib[-index,]

  #Generate PLS Model for Feature Selection
    fit_control <- trainControl(
    method = "LOOCV",
    number = 10)
    rf_fit <- train(Trametinib_1uM ~ ., 
                data = trainSet, 
                metric = "RMSE",
                method = "pls",
                tuneLength = 7,
                maximize=TRUE,
                trControl = fit_control,
                linout = TRUE)

  #Select Features with Importance > 50 in PLS Model    
    varImp_rf <- varImp(rf_fit)
    Importance<-varImp_rf[1]%>%as.data.frame()
    Importance_1 <- subset (Importance, Overall>50)
    Importance_1$kinase<-rownames(Importance_1)

    
  #Fillter Features with Importance > 50 in PLS Model in Train and Test set for RF Model
    Subset_K<-as.vector(rownames(Importance_1))
    Subset_K<-c(Subset_K,"ids","Trametinib_1uM")
    trainSet$ids<-rownames(trainSet)
    trainSet<- subset(trainSet,select=Subset_K)
    trainSet$ids<-NULL
    testSet$ids<-rownames(testSet)
    testSet<- subset(testSet,select=Subset_K)
    testSet$ids<-NULL

  #Generate RF Model
    fit_control <- trainControl(## 10-fold CV
    method = "LOOCV",
    number = 20)
    rf_fit <- train(Trametinib_1uM ~ ., 
                data = trainSet, 
                metric = "RMSE",
                method = "rf",
                tuneLength = 7,
                maximize=TRUE,
                trControl = fit_control,
                ntree=1000,
                linout = TRUE)
    
  #Generate Datasets of Predicted and Measured Values for Train and Test Sets
    Predictors_train<- predict(rf_fit, trainSet)
    Predictors_train<-as.data.frame(Predictors_train)
    Predictors_test <- predict(rf_fit, testSet)
    Predictors_test<-as.data.frame(Predictors_test)
    Obserbed_train<-trainSet$Trametinib_1uM
    Obserbed_train<-as.data.frame(Obserbed_train)
    Obserbed_test<-testSet$Trametinib_1uM
    Obserbed_test<-as.data.frame(Obserbed_test)
    trainplot<- cbind(Obserbed_train, Predictors_train)
    testplot<- cbind(Obserbed_test, Predictors_test)
    
  #Calculate RMSE for Train and Test Sets
    Train_RMSE<-RMSE(trainplot$Predictors_train, trainplot$Obserbed_train ,na.rm = T)
    Test_RMSE<-RMSE(testplot$Predictors_test, testplot$Obserbed_test, na.rm = T)

  #Select Features with Importance >20 in RF Model
    varImp_rf <- varImp(rf_fit)
    Importance<-varImp_rf[1]%>%as.data.frame()
    Importance_1 <- subset (Importance, Overall>20)
    Importance_1$kinase<-rownames(Importance_1)
    dfImportance<-Importance_1$kinase
    dfImportance
    dfImportance_1<-c(dfImportance_1, dfImportance) 

  #Generate Summary of Relevant Parameters
    Lengh_sub<-length(Subset_K)
    Spearvalue <- cor.test(testplot$Predictors_test, testplot$Obserbed_test, method ="spearman")
    Spearvalue_B<-Spearvalue[3]%>%as.data.frame()
    Spearvalue_A<-Spearvalue[4]%>%as.data.frame()
    dfcolspearman <- cbind(Spearvalue_A,Spearvalue_B,Lengh_sub,Test_RMSE)
    dfresults<- dfcolspearman
    rownames(dfresults)<-A
    dfresults1 <- rbind(dfresults1, dfresults)
  }


############################################################################################################################


#Generate a CSV File with the Summary of Relevant Parameters
dfresults_Tramentinib<-dfresults1
write.csv(dfresults_Tramentinib, "C:/Users/Pedro/Documents/Work/work_from_home/Book Chapter/dfresults_Tramentinib.csv") #CHANGE PATH

#Calculate Average RMSE across Replicates of the Worklow
PDT_average_RMSE<-mean(dfresults1$Test_RMSE)
PDT_average_RMSE

#Plot Average RMSE across Replicates of the Worflow
dfresults1$RMSE<-"RMSE"
Boxplot_Average_RMSE<-ggplot(dfresults1, aes(x=RMSE, y=Test_RMSE))+
  geom_boxplot(colour="black")+
  geom_point(colour="black", pch=16, size = 5)+
  ylim(0,0.5)+
  theme_minimal()+
  theme(strip.text.x = element_text(size=18))+
  theme(strip.background = element_rect( fill="white", size=1))+
  theme(axis.title.x= element_text(size=18))+
  theme(axis.text.x=element_text(size=16, colour = "black", hjust=0.95,vjust=0.2))+
  theme(axis.text.y=element_text(size=16, colour = "black",  hjust=0.95,vjust=0.2))+
  theme(axis.line = element_line(color="black"))
Boxplot_Average_RMSE

#Plot Frequency of Features with Importance in the RF Model >20 across Replicates of the Workflow 
dfImportance_2<-as.data.frame(dfImportance_1)

Barplot_Importance_Frequency<-ggplot(dfImportance_2,aes( x = reorder(dfImportance_1,dfImportance_1,function(x)-length(x))))+
  geom_bar( color="black", position=position_dodge())+
  theme_minimal()+
  theme(strip.text.x = element_text(size=18))+
  theme(strip.background = element_rect( fill="white", size=1))+
  theme(axis.title.x= element_text(size=18))+
  theme(axis.text.x=element_text(size=16, colour = "black", angle = 90, hjust=0.95,vjust=0.2))+
  theme(axis.text.y=element_text(size=16, colour = "black",  hjust=0.95,vjust=0.2))+
  theme(axis.line = element_line(color="black"))
Barplot_Importance_Frequency

