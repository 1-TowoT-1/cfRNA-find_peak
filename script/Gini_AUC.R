library(edgeR)
library(pROC)
library(stringr)
library(caret)
library(plyr)
argv <- commandArgs(T)


DEG_dat <- read.table(argv[1],sep="\t",header=T)
Treat <- argv[2]
Control <- argv[3]
output <- argv[4]



Treat <- unlist(str_split(pattern =",",Treat))
Control <- unlist(str_split(pattern =",",Control))



#分组的评均值
Treat_mean <- rowMeans(data.frame(DEG_dat[,Treat]))
Control_mean <- rowMeans(data.frame(DEG_dat[,Control]))



#Gini系数评估
Control_gini <- gini(t(DEG_dat[,Control]))
Treat_gini <- gini(t(DEG_dat[,Treat]))



#AUC评估
model_pre <- DEG_dat[,c("Gene_ID",Treat,Control)]
rownames(model_pre) <- model_pre$Gene_ID
model_pre <- as.data.frame(t(model_pre[,-1]))
type <- c(rep(0,length(Control)),rep(1,length(Treat)))
model_pre$Type <- type

train <- model_pre
test <- model_pre

Single_Factor <- function(x){
  #第一步:训练数据模型评估
  FML <- as.formula(paste0("Type~",x))
  #构建回归模型
  GLogit <- glm(FML, data = train, family=binomial('logit'))
  #提取主要信息
  GSum <- summary(GLogit)
  #提取HR
  #HR <- round(GSum$coefficients[,2],2)
  #提取P值
  PValue <- (GSum$coefficients[,4])[2]
  #提取期望值
  Estimate=paste0(round(GSum$coefficients[,'Estimate'],5), collapse = '+')
  #提取CI值
  CI <- paste0(round(confint(GLogit)[2,],5),collapse = '~')
  
  #第二步：测试数据测试是否符合模型标准
  test_rdata <- test[,c(x,"Type")]
  #predict函数来评估我们的test数据是否符合GLogit模型
  probability <- predict(object = GLogit, newdata = test, type = 'response')
  pre_test <- cbind(test_rdata,probability)
  pre_test <- transform(pre_test, predict = ifelse(probability <= 0.5, 0, 1))
  ##总结出混淆矩阵
  pre_table <- table(pre_test$Type,pre_test$predict)
  pre_colname <- colnames(pre_table)
  
  ###计算该gene（因素）构造模型对测试数据的NPV，SN,SP,accuracy
  #NPV：阴性预测值，指的是在所有被预测为阴性的样本中，实际为阴性的比例
  #SN：敏感性，也称为真阳性率或查全率
  #SP：特异度，评估预测负例的准确率
  #accuracy：评估模型的正确分类比例
  if (length(pre_colname) == 2) {
    NPV <- pre_table[1,1]/(pre_table[1,1]+pre_table[2,1])
    SN <- pre_table[2,2]/(pre_table[2,1]+pre_table[2,2])
    SP <- pre_table[1,1]/(pre_table[1,1]+pre_table[1,2])
    accuracy <- sum(diag(pre_table))/sum(pre_table)
  }else {
    if (pre_colname[1]==0) {
      NPV <- pre_table[1,1]/(pre_table[1,1]+pre_table[2,1])
      SN <- 0
      SP <- pre_table[1,1]/pre_table[1,1]
      accuracy <- pre_table[1,1]/(pre_table[1,1]+pre_table[2,1])
    }else {
      if (pre_colname[1]==1) {
        NPV <- NA
        SN <- pre_table[2,1]/pre_table[2,1]
        SP <- 0
        accuracy <- pre_table[2,1]/(pre_table[1,1]+pre_table[2,1])
      }else {
        NPV <- NA
        SN <- NA
        SP <- NA
        accuracy <- NA
      }
    }
  }
  
  pre_sum <- sum(pre_table)
  data_sum <- dim(train)[1]
  
  #模型预测的ROC曲线,AUC来评估模型的好坏
  roc_curve <- roc(pre_test$Type ~ probability)
  roc_x <- 1 - roc_curve$specificities
  roc_y <- roc_curve$sensitivities
  # plot(roc_x,roc_y)
  AUC = roc_curve$auc[1]
  
  UniLogit <- data.frame('Characteristics' = x, 'Train_Num' = data_sum, 'Predict_Num' = pre_sum, accuracy, NPV, SP, SN, AUC, 'Estimate' = Estimate, 'CI95' = CI, 'P Value' = PValue)
  return(UniLogit)
}

VarNames <- colnames(train)[1:(ncol(train)-1)]
Univar <- lapply(VarNames,Single_Factor)
#使用ldply将数据转化为 data.frame 格式
Univar <- ldply(Univar,data.frame)
Univar_mid<- Univar[,c("Characteristics","AUC")]
colnames(Univar_mid) <- c("Gene_ID","AUC")

new_DEG_dat <- DEG_dat[,c(1,3,6,7)]
new_DEG_dat <- cbind(new_DEG_dat,"Control_mean"=Control_mean,
                     "Control_gini"=Control_gini,
                     "Treat_mean"=Treat_mean,
                     "Treat_gini"=Treat_gini)
new_DEG_dat<- merge(new_DEG_dat,Univar_mid,by="Gene_ID")


write.table(new_DEG_dat,output,sep="\t",quote=F,row.names=F)

