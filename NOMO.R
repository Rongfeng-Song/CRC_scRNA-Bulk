library(rms)
library(regplot)
library(tidyverse)
library(survival)
library(regplot)

load("OKModel_risk.rdata")

entire_risk
TCGA_risk
train_risk = TCGA_risk
train_risk = entire_risk
#修改数据类型
train_risk[1:4,1:4]
dim(train_risk)
str(train_risk[,c(1:10,1224,1225)])
train_risk$Age = as.numeric(train_risk$Age)
train_risk = na.omit(train_risk)
table(train_risk$Gender)
train_risk$Gender = ifelse(train_risk$Gender == "female","Female",
                           ifelse(train_risk$Gender == "male","Male",ifelse(train_risk$Gender == "Male","Male","Female")))
train_risk$Gender =as.factor(train_risk$Gender)
train_risk$Clinical_stage = as.character(train_risk$Clinical_stage)
train_risk$risk = as.character(train_risk$risk)
#构建多因素模型
multiCox <- coxph(Surv(OS_time, Status) ~  Age + Gender +  Clinical_stage + 
                    risk, data = train_risk)
risk = train_risk
pdf("entire_nomo.pdf",width = 6,height = 6)
regplot(multiCox,
        observation=risk[4,], #用哪行观测
        obscol = "#326db1",
        failtime = c(1,3,5), 
        plots = c("bars","boxes"),
        droplines = T, # 是否画竖线
        points = T,
        title = "TCGA COAD cohort", # 更换标题
        # odds = T, # 是否显示OR值
        showP = T, # 是否显示变量的显著性标记（默认：T）
        rank = "sd", # 根据sd给变量排序
        # interval="confidence", # 展示可信区间
        clickable = F, # 是否可以交互（默认：F）
        prfail = F)
dev.off()
#矫正曲线
risk[4,1:9]
#1-year
if(T){
  
  cox1 <- cph(Surv(OS_time, Status) ~ Age + Gender + Clinical_stage +risk,
              surv=T,x=T, y=T,time.inc = 1,data=risk)
  cal <- calibrate(cox1, cmethod="KM", method="boot", u=1,
                   m= 100, B=1000)
  #3-year
  cox3 <- cph(Surv(OS_time, Status) ~ Age + Gender + Clinical_stage +risk,
              surv=T,x=T, y=T,time.inc = 3,data=risk)
  ca3 <- calibrate(cox3, cmethod="KM", method="boot", u=3,
                   m= 100, B=1000)
  #5-year
  cox5 <- cph(Surv(OS_time, Status) ~ Age + Gender + Clinical_stage +risk,
              surv=T,x=T, y=T,time.inc = 5,data=risk)
  ca5 <- calibrate(cox5, cmethod="KM", method="boot", u=5,
                   m= 100, B=1000)
}


#绘制曲线
if(T){
  pdf("TCGA Correction curve.pdf",width = 5,height = 5)
  plot(cal,lwd=2,lty=1,errbar.col="black",
       xlim = c(0.5,1) ,
       ylim = c(0.4,1), 
       xlab ="Nomogram-Predicted Probability of 1,3,5Year Survival",
       ylab="Actual 1,3,5Year Survival",col="blue",sub=F)
  plot(ca3,
       add = T ,
       lwd=2,lty=1,errbar.col="black",col="red",sub=F)
  plot(ca5,
       add = T ,
       lwd=2,lty=1,errbar.col="black",col="green",sub=F)
  legend("bottomright", legend = c("1-year", "3-year", "5-year"), 
         col = c("blue", "red", "green"), lty = 1, lwd = 2, 
         text.col = "black", cex = 0.8, box.col = "transparent")
  
  dev.off()
}

