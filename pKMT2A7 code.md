```R
library(rms)
library(survival)
library(glmnet)
library(dplyr)
library(AnnoProbe)

#Univariable Cox regression using EFS as endpoint
training <- read.csv("~/matrix_training.csv")
validation <- read.csv("~/matrix_validation.csv")
variables <- colnames(training)[c(6:18821)]
variables
uni_cox <- function(variables){
  formula <- as.formula(paste0("Surv(EFS,EFS_Status)~",variables))
  surv_uni_cox <- summary(coxph(formula,data=training))
  uni_cox_rep <- data.frame("Factor"=variables,
                            "beta"=round(surv_uni_cox$coefficients[,1],2),
                            "HR"=round(exp(surv_uni_cox$coefficients[,1]),2),
                            "p_value"=round(surv_uni_cox$coefficients[,5],3),
                            "CI5" <-round(surv_uni_cox$conf.int[,3],2),
                            "CI95" <-round(surv_uni_cox$conf.int[,4],2),
                            "Wald_pval"=round(as.numeric(surv_uni_cox$waldtest[3]),1))
  return(uni_cox_rep)
}
Uni_cox <-lapply(variables,uni_cox)
Uni_cox <- ldply(Uni_cox,data.frame)
Uni_cox <- Uni_cox[c(Uni_cox$p_value<0.05),]

#Identify protein coding genes
gene_id_anno <- annoGene(Uni_cox$Factor,ID_type="ENSEMBL" )
gene_pc <- subset(gene_id_anno,gene_id_anno$biotypes=="protein_coding")
cox_var <- gene_pc$ENSEMBL
index <- match(cox_var,colnames(training))
training <- training[,c(1:5,index)]

index <- match(colnames(training),colnames(validation))
validation <- validation[,index]

#LASSO
set.seed(34)
v1 <- data.matrix(training[,c(6:2050)])
v2 <- data.matrix(Surv(training$EFS,training$EFS_Status))
fit <- glmnet(v1,v2,family = "cox")
plot(fit, xvar = "lambda",label = F)
fit.cv <- cv.glmnet(v1,v2,family = "cox",nfolds = 10)
coe <- coef(fit,s=fit.cv$lambda.min)
act_index <- which(coe!=0)
row.names(coe)[act_index]
coefficient <- c(coe[act_index])
coefficient

#Calculate risk score for each patient in training set
risk_score <- predict(fit.cv,newx=data.matrix(v1),s=fit.cv$lambda.min,type="response")
training <- as.data.frame(cbind(training[1:5],risk_score))
training$risk_level <- ifelse(training$risk_score<1.11,"G1","G2")

#Calculate risk score for each patient in validation set
v1_val <- data.matrix(validation_cox_efs[,c(6:2050)])
risk_score_val <- predict(fit.cv,newx=data.matrix(v1_val),s=fit.cv$lambda.min,type="response")
validation <- as.data.frame(cbind(validation[,1:5],risk_score_val))
colnames(validation)[6] <- "risk_score"
validation$risk_level <- ifelse(validation$risk_score<1.11,"G1","G2")
```

