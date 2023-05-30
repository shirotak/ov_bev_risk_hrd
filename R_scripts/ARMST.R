# R script for adjusted RMST
# Ref; https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fpst.2151&file=pst2151-sup-0001-Supinfo.pdf

library("dplyr")
library("survival")
library("survRM2")
library("RISCA")
source('./ARMST_functions.R')

# All time
in_df='input/ICON7_all.tsv'
out_prefix='R_results/ARSMT_alltime'
duration='PFS'
data=read.delim(in_df)
cat='Bev'
covs=c('Stage34','Optimal','Age_high','Serous')

# unadjusted
rmst_result_sub1_unadj =rmst2(data[,duration],data$Rec,data[,cat])
add1=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[1,] ))
row.names(add1)=cat
dfw1=add1
add0=t( data.frame(rmst_result_sub1_unadj$RMST.arm1$rmst, rmst_result_sub1_unadj$RMST.arm0$rmst))
row.names(add0)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw0=add0
add4=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[2,] ))
row.names(add4)=cat
dfw4=add4

# Adjusted
Pr1 <- glm( Bev ~ Stage34 + Optimal + Age_high+ Serous,
           data=data,family=binomial(link = "logit"))$fitted.values
W <- ( data[,cat]==1) * (1/Pr1) +( data[,cat]==0) * (1)/(1-Pr1)
res.akm <- ipw.survival(times=data[,duration],failures=data$Rec,
                        variable= data[,cat], weights=W)$table.surv
data$W <- W
survdat1 <- res.akm[res.akm$variable==1,]
wtdat1 <- data[data[,cat]==1,]
wtdat1$time =data[data[,cat]==1,duration]
wtdat1$status = data[data[,cat]==1,]$Rec
survdat0 <- res.akm[res.akm$variable==0,]
wtdat0 <- data[data[,cat]==0,]
wtdat0$time = data[data[,cat]==0,duration]
wtdat0$status <-data[data[,cat]==0,]$Rec
tau= min(c(max(wtdat1$time),max(wtdat0$time)))
rmst1akme_1 <- rmst1_AKME(survdata =survdat1,wtdata =wtdat1,tau=tau)
rmst1akme_0 <- rmst1_AKME(survdata =survdat0,wtdata =wtdat0,tau=tau)
rmst1akme_1$rmst[1:2]
rmst1akme_0$rmst[1:2]
armst <- ARMST_diff(rmst_trt = rmst1akme_1$rmst[1],
                    rmst_ctrl = rmst1akme_0$rmst[1],
                    rmst_var_trt = rmst1akme_1$rmst.var,
                    rmst_var_ctrl = rmst1akme_0$rmst.var)
# write out
add2=armst$rmst.diff.result
row.names(add2)=cat
dfw2= add2
add3=armst$rmst.ratio.result
row.names(add3)=cat
dfw3= add3
add5=t( data.frame(rmst1akme_1$rmst, rmst1akme_0$rmst))
row.names(add5)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw5=add5


## repeat
cat='Stage34'
covs=c('Bev','Optimal','Age_high','Serous')

# unadjusted
rmst_result_sub1_unadj =rmst2(data[,duration],data$Rec,data[,cat])
add1=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[1,] ))
row.names(add1)=cat
dfw1= rbind(dfw1,add1)
add0=t( data.frame(rmst_result_sub1_unadj$RMST.arm1$rmst, rmst_result_sub1_unadj$RMST.arm0$rmst))
row.names(add0)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw0= rbind(dfw0,add0)
add4=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[2,] ))
row.names(add4)=cat
dfw4= rbind(dfw4,add4)

# adjusted
Pr1 <- glm( Stage34 ~ Bev + Optimal + Age_high+ Serous,
            data=data,family=binomial(link = "logit"))$fitted.values
W <- ( data[,cat]==1) * (1/Pr1) +( data[,cat]==0) * (1)/(1-Pr1)
res.akm <- ipw.survival(times=data[,duration],failures=data$Rec,
                        variable= data[,cat], weights=W)$table.surv
data$W <- W
survdat1 <- res.akm[res.akm$variable==1,]
wtdat1 <- data[data[,cat]==1,]
wtdat1$time =data[data[,cat]==1,duration]
wtdat1$status = data[data[,cat]==1,]$Rec
survdat0 <- res.akm[res.akm$variable==0,]
wtdat0 <- data[data[,cat]==0,]
wtdat0$time = data[data[,cat]==0,duration]
wtdat0$status <-data[data[,cat]==0,]$Rec
tau= min(c(max(wtdat1$time),max(wtdat0$time)))
rmst1akme_1 <- rmst1_AKME(survdata =survdat1,wtdata =wtdat1,tau=tau)
rmst1akme_0 <- rmst1_AKME(survdata =survdat0,wtdata =wtdat0,tau=tau)
armst <- ARMST_diff(rmst_trt = rmst1akme_1$rmst[1],
                    rmst_ctrl = rmst1akme_0$rmst[1],
                    rmst_var_trt = rmst1akme_1$rmst.var,
                    rmst_var_ctrl = rmst1akme_0$rmst.var)
# write out
add2=armst$rmst.diff.result
row.names(add2)=cat
dfw2= rbind(dfw2,add2)
add3=armst$rmst.ratio.result
row.names(add3)=cat
dfw3= rbind(dfw3,add3)
add5=t( data.frame(rmst1akme_1$rmst, rmst1akme_0$rmst))
row.names(add5)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw5= rbind(dfw5,add5)

## repeat
cat='Optimal'
covs=c('Bev','Stage34','Age_high','Serous')
# unadjusted
rmst_result_sub1_unadj =rmst2(data[,duration],data$Rec,data[,cat])
add1=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[1,] ))
row.names(add1)=cat
dfw1= rbind(dfw1,add1)
add0=t( data.frame(rmst_result_sub1_unadj$RMST.arm1$rmst, rmst_result_sub1_unadj$RMST.arm0$rmst))
row.names(add0)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw0= rbind(dfw0,add0)
add4=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[2,] ))
row.names(add4)=cat
dfw4= rbind(dfw4,add4)
# adjusted
Pr1 <- glm( Optimal ~ Bev + Stage34 + Age_high + Serous,
            data=data,family=binomial(link = "logit"))$fitted.values
W <- ( data[,cat]==1) * (1/Pr1) +( data[,cat]==0) * (1)/(1-Pr1)
res.akm <- ipw.survival(times=data[,duration],failures=data$Rec,
                        variable= data[,cat], weights=W)$table.surv
data$W <- W
survdat1 <- res.akm[res.akm$variable==1,]
wtdat1 <- data[data[,cat]==1,]
wtdat1$time =data[data[,cat]==1,duration]
wtdat1$status = data[data[,cat]==1,]$Rec
survdat0 <- res.akm[res.akm$variable==0,]
wtdat0 <- data[data[,cat]==0,]
wtdat0$time = data[data[,cat]==0,duration]
wtdat0$status <-data[data[,cat]==0,]$Rec
tau= min(c(max(wtdat1$time),max(wtdat0$time)))
rmst1akme_1 <- rmst1_AKME(survdata =survdat1,wtdata =wtdat1,tau=tau)
rmst1akme_0 <- rmst1_AKME(survdata =survdat0,wtdata =wtdat0,tau=tau)
armst <- ARMST_diff(rmst_trt = rmst1akme_1$rmst[1],
                    rmst_ctrl = rmst1akme_0$rmst[1],
                    rmst_var_trt = rmst1akme_1$rmst.var,
                    rmst_var_ctrl = rmst1akme_0$rmst.var)
# write out
add2=armst$rmst.diff.result
row.names(add2)=cat
dfw2= rbind(dfw2,add2)
add3=armst$rmst.ratio.result
row.names(add3)=cat
dfw3= rbind(dfw3,add3)
add5=t( data.frame(rmst1akme_1$rmst, rmst1akme_0$rmst))
row.names(add5)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw5= rbind(dfw5,add5)

## repeat
cat='Age_high'
covs=c('Bev','Stage34','Optimal','Serous')
# unadjusted
rmst_result_sub1_unadj =rmst2(data[,duration],data$Rec,data[,cat])
add1=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[1,] ))
row.names(add1)=cat
dfw1= rbind(dfw1,add1)
add0=t( data.frame(rmst_result_sub1_unadj$RMST.arm1$rmst, rmst_result_sub1_unadj$RMST.arm0$rmst))
row.names(add0)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw0= rbind(dfw0,add0)
add4=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[2,] ))
row.names(add4)=cat
dfw4= rbind(dfw4,add4)
# adjusted
Pr1 <- glm( Age_high ~ Bev + Stage34 + Optimal + Serous,
            data=data,family=binomial(link = "logit"))$fitted.values
W <- ( data[,cat]==1) * (1/Pr1) +( data[,cat]==0) * (1)/(1-Pr1)
res.akm <- ipw.survival(times=data[,duration],failures=data$Rec,
                        variable= data[,cat], weights=W)$table.surv
data$W <- W
survdat1 <- res.akm[res.akm$variable==1,]
wtdat1 <- data[data[,cat]==1,]
wtdat1$time =data[data[,cat]==1,duration]
wtdat1$status = data[data[,cat]==1,]$Rec
survdat0 <- res.akm[res.akm$variable==0,]
wtdat0 <- data[data[,cat]==0,]
wtdat0$time = data[data[,cat]==0,duration]
wtdat0$status <-data[data[,cat]==0,]$Rec
tau= min(c(max(wtdat1$time),max(wtdat0$time)))
rmst1akme_1 <- rmst1_AKME(survdata =survdat1,wtdata =wtdat1,tau=tau)
rmst1akme_0 <- rmst1_AKME(survdata =survdat0,wtdata =wtdat0,tau=tau)
armst <- ARMST_diff(rmst_trt = rmst1akme_1$rmst[1],
                    rmst_ctrl = rmst1akme_0$rmst[1],
                    rmst_var_trt = rmst1akme_1$rmst.var,
                    rmst_var_ctrl = rmst1akme_0$rmst.var)
# write out
add2=armst$rmst.diff.result
row.names(add2)=cat
dfw2= rbind(dfw2,add2)
add3=armst$rmst.ratio.result
row.names(add3)=cat
dfw3= rbind(dfw3,add3)
add5=t( data.frame(rmst1akme_1$rmst, rmst1akme_0$rmst))
row.names(add5)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw5= rbind(dfw5,add5)

## repeat
cat='Serous'
covs=c('Bev','Stage34','Optimal','Age_high')
# unadjusted
rmst_result_sub1_unadj =rmst2(data[,duration],data$Rec,data[,cat])
add1=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[1,] ))
row.names(add1)=cat
dfw1= rbind(dfw1,add1)
add0=t( data.frame(rmst_result_sub1_unadj$RMST.arm1$rmst, rmst_result_sub1_unadj$RMST.arm0$rmst))
row.names(add0)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw0= rbind(dfw0,add0)
add4=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[2,] ))
row.names(add4)=cat
dfw4= rbind(dfw4,add4)
# adjusted
Pr1 <- glm( Serous ~ Bev + Stage34 + Optimal + Age_high,
           data=data,family=binomial(link = "logit"))$fitted.values
W <- ( data[,cat]==1) * (1/Pr1) +( data[,cat]==0) * (1)/(1-Pr1)
res.akm <- ipw.survival(times=data[,duration],failures=data$Rec,
                       variable= data[,cat], weights=W)$table.surv
data$W <- W
survdat1 <- res.akm[res.akm$variable==1,]
wtdat1 <- data[data[,cat]==1,]
wtdat1$time =data[data[,cat]==1,duration]
wtdat1$status = data[data[,cat]==1,]$Rec
survdat0 <- res.akm[res.akm$variable==0,]
wtdat0 <- data[data[,cat]==0,]
wtdat0$time = data[data[,cat]==0,duration]
wtdat0$status <-data[data[,cat]==0,]$Rec
tau= min(c(max(wtdat1$time),max(wtdat0$time)))
rmst1akme_1 <- rmst1_AKME(survdata =survdat1,wtdata =wtdat1,tau=tau)
rmst1akme_0 <- rmst1_AKME(survdata =survdat0,wtdata =wtdat0,tau=tau)
armst <- ARMST_diff(rmst_trt = rmst1akme_1$rmst[1],
                    rmst_ctrl = rmst1akme_0$rmst[1],
                    rmst_var_trt = rmst1akme_1$rmst.var,
                    rmst_var_ctrl = rmst1akme_0$rmst.var)
# write out
add2=armst$rmst.diff.result
row.names(add2)=cat
dfw2= rbind(dfw2,add2)
add3=armst$rmst.ratio.result
row.names(add3)=cat
dfw3= rbind(dfw3,add3)
add5=t( data.frame(rmst1akme_1$rmst, rmst1akme_0$rmst))
row.names(add5)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw5= rbind(dfw5,add5)

# Result Check
dfw1
dfw4
dfw2
dfw3
dfw0
dfw5
out_prefix

# Write out
write.table(dfw0,paste0(out_prefix,'_unadjusted_raw.tsv'),sep='\t',quote = F,
            row.names = T,col.names = NA)
write.table(dfw1,paste0(out_prefix,'_unadjusted_diff.tsv'),sep='\t',quote = F,
            row.names = T,col.names = NA)
write.table(dfw4,paste0(out_prefix,'_unadjusted_ratio.tsv'),sep='\t',quote = F,
            row.names = T,col.names = NA)
write.table(dfw5,paste0(out_prefix,'_adjusted_raw.tsv'),sep='\t',quote = F,
            row.names = T,col.names = NA)
write.table(dfw2,paste0(out_prefix,'_adjusted_diff.tsv'),sep='\t',quote = F,
            row.names = T,col.names = NA)
write.table(dfw3,paste0(out_prefix,'_adjusted_ratio.tsv'),sep='\t',quote = F,
            row.names = T,col.names = NA)

##################################### Before ###################################
in_df='input/ICON7_all_before.tsv'
out_prefix='./R_results/ARSMT_before'
duration='PFS'
data=read.delim(in_df)

cat='Bev'
covs=c('Stage34','Optimal','Age_high','Serous')

# unadjusted
rmst_result_sub1_unadj =rmst2(data[,duration],data$Rec,data[,cat])
add1=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[1,] ))
row.names(add1)=cat
dfw1=add1
add0=t( data.frame(rmst_result_sub1_unadj$RMST.arm1$rmst, rmst_result_sub1_unadj$RMST.arm0$rmst))
row.names(add0)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw0=add0
add4=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[2,] ))
row.names(add4)=cat
dfw4=add4
# Adjusted
Pr1 <- glm( Bev ~ Stage34 + Optimal + Age_high+ Serous,
            data=data,family=binomial(link = "logit"))$fitted.values
W <- ( data[,cat]==1) * (1/Pr1) +( data[,cat]==0) * (1)/(1-Pr1)
res.akm <- ipw.survival(times=data[,duration],failures=data$Rec,
                        variable= data[,cat], weights=W)$table.surv
data$W <- W
survdat1 <- res.akm[res.akm$variable==1,]
wtdat1 <- data[data[,cat]==1,]
wtdat1$time =data[data[,cat]==1,duration]
wtdat1$status = data[data[,cat]==1,]$Rec
survdat0 <- res.akm[res.akm$variable==0,]
wtdat0 <- data[data[,cat]==0,]
wtdat0$time = data[data[,cat]==0,duration]
wtdat0$status <-data[data[,cat]==0,]$Rec
tau= min(c(max(wtdat1$time),max(wtdat0$time)))
rmst1akme_1 <- rmst1_AKME(survdata =survdat1,wtdata =wtdat1,tau=tau)
rmst1akme_0 <- rmst1_AKME(survdata =survdat0,wtdata =wtdat0,tau=tau)
rmst1akme_1$rmst[1:2]
rmst1akme_0$rmst[1:2]
armst <- ARMST_diff(rmst_trt = rmst1akme_1$rmst[1],
                    rmst_ctrl = rmst1akme_0$rmst[1],
                    rmst_var_trt = rmst1akme_1$rmst.var,
                    rmst_var_ctrl = rmst1akme_0$rmst.var)
add2=armst$rmst.diff.result
row.names(add2)=cat
dfw2= add2
add3=armst$rmst.ratio.result
row.names(add3)=cat
dfw3= add3
add5=t( data.frame(rmst1akme_1$rmst, rmst1akme_0$rmst))
row.names(add5)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw5=add5

## repeat
cat='Stage34'
covs=c('Bev','Optimal','Age_high','Serous')

# unadjusted
rmst_result_sub1_unadj =rmst2(data[,duration],data$Rec,data[,cat])
add1=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[1,] ))
row.names(add1)=cat
dfw1= rbind(dfw1,add1)
add0=t( data.frame(rmst_result_sub1_unadj$RMST.arm1$rmst, rmst_result_sub1_unadj$RMST.arm0$rmst))
row.names(add0)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw0= rbind(dfw0,add0)
add4=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[2,] ))
row.names(add4)=cat
dfw4= rbind(dfw4,add4)
# adjusted
Pr1 <- glm( Stage34 ~ Bev + Optimal + Age_high+ Serous,
            data=data,family=binomial(link = "logit"))$fitted.values
W <- ( data[,cat]==1) * (1/Pr1) +( data[,cat]==0) * (1)/(1-Pr1)
res.akm <- ipw.survival(times=data[,duration],failures=data$Rec,
                        variable= data[,cat], weights=W)$table.surv
data$W <- W
survdat1 <- res.akm[res.akm$variable==1,]
wtdat1 <- data[data[,cat]==1,]
wtdat1$time =data[data[,cat]==1,duration]
wtdat1$status = data[data[,cat]==1,]$Rec
survdat0 <- res.akm[res.akm$variable==0,]
wtdat0 <- data[data[,cat]==0,]
wtdat0$time = data[data[,cat]==0,duration]
wtdat0$status <-data[data[,cat]==0,]$Rec
tau= min(c(max(wtdat1$time),max(wtdat0$time)))
rmst1akme_1 <- rmst1_AKME(survdata =survdat1,wtdata =wtdat1,tau=tau)
rmst1akme_0 <- rmst1_AKME(survdata =survdat0,wtdata =wtdat0,tau=tau)
armst <- ARMST_diff(rmst_trt = rmst1akme_1$rmst[1],
                    rmst_ctrl = rmst1akme_0$rmst[1],
                    rmst_var_trt = rmst1akme_1$rmst.var,
                    rmst_var_ctrl = rmst1akme_0$rmst.var)
add2=armst$rmst.diff.result
row.names(add2)=cat
dfw2= rbind(dfw2,add2)
add3=armst$rmst.ratio.result
row.names(add3)=cat
dfw3= rbind(dfw3,add3)
add5=t( data.frame(rmst1akme_1$rmst, rmst1akme_0$rmst))
row.names(add5)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw5= rbind(dfw5,add5)

## repeat
cat='Optimal'
covs=c('Bev','Stage34','Age_high','Serous')
# unadjusted
rmst_result_sub1_unadj =rmst2(data[,duration],data$Rec,data[,cat])
add1=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[1,] ))
row.names(add1)=cat
dfw1= rbind(dfw1,add1)
add0=t( data.frame(rmst_result_sub1_unadj$RMST.arm1$rmst, rmst_result_sub1_unadj$RMST.arm0$rmst))
row.names(add0)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw0= rbind(dfw0,add0)
add4=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[2,] ))
row.names(add4)=cat
dfw4= rbind(dfw4,add4)
# adjusted
Pr1 <- glm( Optimal ~ Bev + Stage34 + Age_high + Serous,
            data=data,family=binomial(link = "logit"))$fitted.values
W <- ( data[,cat]==1) * (1/Pr1) +( data[,cat]==0) * (1)/(1-Pr1)
res.akm <- ipw.survival(times=data[,duration],failures=data$Rec,
                        variable= data[,cat], weights=W)$table.surv
data$W <- W
survdat1 <- res.akm[res.akm$variable==1,]
wtdat1 <- data[data[,cat]==1,]
wtdat1$time =data[data[,cat]==1,duration]
wtdat1$status = data[data[,cat]==1,]$Rec
survdat0 <- res.akm[res.akm$variable==0,]
wtdat0 <- data[data[,cat]==0,]
wtdat0$time = data[data[,cat]==0,duration]
wtdat0$status <-data[data[,cat]==0,]$Rec
tau= min(c(max(wtdat1$time),max(wtdat0$time)))
rmst1akme_1 <- rmst1_AKME(survdata =survdat1,wtdata =wtdat1,tau=tau)
rmst1akme_0 <- rmst1_AKME(survdata =survdat0,wtdata =wtdat0,tau=tau)
armst <- ARMST_diff(rmst_trt = rmst1akme_1$rmst[1],
                    rmst_ctrl = rmst1akme_0$rmst[1],
                    rmst_var_trt = rmst1akme_1$rmst.var,
                    rmst_var_ctrl = rmst1akme_0$rmst.var)
add2=armst$rmst.diff.result
row.names(add2)=cat
dfw2= rbind(dfw2,add2)
add3=armst$rmst.ratio.result
row.names(add3)=cat
dfw3= rbind(dfw3,add3)
add5=t( data.frame(rmst1akme_1$rmst, rmst1akme_0$rmst))
row.names(add5)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw5= rbind(dfw5,add5)

## repeat
cat='Age_high'
covs=c('Bev','Stage34','Optimal','Serous')
# unadjusted
rmst_result_sub1_unadj =rmst2(data[,duration],data$Rec,data[,cat])
add1=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[1,] ))
row.names(add1)=cat
dfw1= rbind(dfw1,add1)
add0=t( data.frame(rmst_result_sub1_unadj$RMST.arm1$rmst, rmst_result_sub1_unadj$RMST.arm0$rmst))
row.names(add0)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw0= rbind(dfw0,add0)
add4=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[2,] ))
row.names(add4)=cat
dfw4= rbind(dfw4,add4)
# adjusted
Pr1 <- glm( Age_high ~ Bev + Stage34 + Optimal + Serous,
            data=data,family=binomial(link = "logit"))$fitted.values
W <- ( data[,cat]==1) * (1/Pr1) +( data[,cat]==0) * (1)/(1-Pr1)
res.akm <- ipw.survival(times=data[,duration],failures=data$Rec,
                        variable= data[,cat], weights=W)$table.surv
data$W <- W
survdat1 <- res.akm[res.akm$variable==1,]
wtdat1 <- data[data[,cat]==1,]
wtdat1$time =data[data[,cat]==1,duration]
wtdat1$status = data[data[,cat]==1,]$Rec
survdat0 <- res.akm[res.akm$variable==0,]
wtdat0 <- data[data[,cat]==0,]
wtdat0$time = data[data[,cat]==0,duration]
wtdat0$status <-data[data[,cat]==0,]$Rec
tau= min(c(max(wtdat1$time),max(wtdat0$time)))
rmst1akme_1 <- rmst1_AKME(survdata =survdat1,wtdata =wtdat1,tau=tau)
rmst1akme_0 <- rmst1_AKME(survdata =survdat0,wtdata =wtdat0,tau=tau)
armst <- ARMST_diff(rmst_trt = rmst1akme_1$rmst[1],
                    rmst_ctrl = rmst1akme_0$rmst[1],
                    rmst_var_trt = rmst1akme_1$rmst.var,
                    rmst_var_ctrl = rmst1akme_0$rmst.var)
add2=armst$rmst.diff.result
row.names(add2)=cat
dfw2= rbind(dfw2,add2)
add3=armst$rmst.ratio.result
row.names(add3)=cat
dfw3= rbind(dfw3,add3)
add5=t( data.frame(rmst1akme_1$rmst, rmst1akme_0$rmst))
row.names(add5)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw5= rbind(dfw5,add5)

## repeat
cat='Serous'
covs=c('Bev','Stage34','Optimal','Age_high')
# unadjusted
rmst_result_sub1_unadj =rmst2(data[,duration],data$Rec,data[,cat])
add1=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[1,] ))
row.names(add1)=cat
dfw1= rbind(dfw1,add1)
add0=t( data.frame(rmst_result_sub1_unadj$RMST.arm1$rmst, rmst_result_sub1_unadj$RMST.arm0$rmst))
row.names(add0)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw0= rbind(dfw0,add0)
add4=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[2,] ))
row.names(add4)=cat
dfw4= rbind(dfw4,add4)
# adjusted
Pr1 <- glm( Serous ~ Bev + Stage34 + Optimal + Age_high,
            data=data,family=binomial(link = "logit"))$fitted.values
W <- ( data[,cat]==1) * (1/Pr1) +( data[,cat]==0) * (1)/(1-Pr1)
res.akm <- ipw.survival(times=data[,duration],failures=data$Rec,
                        variable= data[,cat], weights=W)$table.surv
data$W <- W
survdat1 <- res.akm[res.akm$variable==1,]
wtdat1 <- data[data[,cat]==1,]
wtdat1$time =data[data[,cat]==1,duration]
wtdat1$status = data[data[,cat]==1,]$Rec
survdat0 <- res.akm[res.akm$variable==0,]
wtdat0 <- data[data[,cat]==0,]
wtdat0$time = data[data[,cat]==0,duration]
wtdat0$status <-data[data[,cat]==0,]$Rec
tau= min(c(max(wtdat1$time),max(wtdat0$time)))
rmst1akme_1 <- rmst1_AKME(survdata =survdat1,wtdata =wtdat1,tau=tau)
rmst1akme_0 <- rmst1_AKME(survdata =survdat0,wtdata =wtdat0,tau=tau)
armst <- ARMST_diff(rmst_trt = rmst1akme_1$rmst[1],
                    rmst_ctrl = rmst1akme_0$rmst[1],
                    rmst_var_trt = rmst1akme_1$rmst.var,
                    rmst_var_ctrl = rmst1akme_0$rmst.var)
add2=armst$rmst.diff.result
row.names(add2)=cat
dfw2= rbind(dfw2,add2)
add3=armst$rmst.ratio.result
row.names(add3)=cat
dfw3= rbind(dfw3,add3)
add5=t( data.frame(rmst1akme_1$rmst, rmst1akme_0$rmst))
row.names(add5)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw5= rbind(dfw5,add5)

# Check
dfw1
dfw4
dfw2
dfw3
dfw0
dfw5
out_prefix
# Write out
write.table(dfw0,paste0(out_prefix,'_unadjusted_raw.tsv'),sep='\t',quote = F,
            row.names = T,col.names = NA)
write.table(dfw1,paste0(out_prefix,'_unadjusted_diff.tsv'),sep='\t',quote = F,
            row.names = T,col.names = NA)
write.table(dfw4,paste0(out_prefix,'_unadjusted_ratio.tsv'),sep='\t',quote = F,
            row.names = T,col.names = NA)
write.table(dfw5,paste0(out_prefix,'_adjusted_raw.tsv'),sep='\t',quote = F,
            row.names = T,col.names = NA)
write.table(dfw2,paste0(out_prefix,'_adjusted_diff.tsv'),sep='\t',quote = F,
            row.names = T,col.names = NA)
write.table(dfw3,paste0(out_prefix,'_adjusted_ratio.tsv'),sep='\t',quote = F,
            row.names = T,col.names = NA)

##################################### After ###################################
in_df='input/ICON7_all_after.tsv'
out_prefix='./R_results/ARSMT_after'
duration='PFS_r'
data=read.delim(in_df)

cat='Bev'
covs=c('Stage34','Optimal','Age_high','Serous')

# unadjusted
rmst_result_sub1_unadj =rmst2(data[,duration],data$Rec,data[,cat])
add1=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[1,] ))
row.names(add1)=cat
dfw1=add1
add0=t( data.frame(rmst_result_sub1_unadj$RMST.arm1$rmst, rmst_result_sub1_unadj$RMST.arm0$rmst))
row.names(add0)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw0=add0
add4=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[2,] ))
row.names(add4)=cat
dfw4=add4
# Adjusted
Pr1 <- glm( Bev ~ Stage34 + Optimal + Age_high+ Serous,
            data=data,family=binomial(link = "logit"))$fitted.values
W <- ( data[,cat]==1) * (1/Pr1) +( data[,cat]==0) * (1)/(1-Pr1)
res.akm <- ipw.survival(times=data[,duration],failures=data$Rec,
                        variable= data[,cat], weights=W)$table.surv
data$W <- W
survdat1 <- res.akm[res.akm$variable==1,]
wtdat1 <- data[data[,cat]==1,]
wtdat1$time =data[data[,cat]==1,duration]
wtdat1$status = data[data[,cat]==1,]$Rec
survdat0 <- res.akm[res.akm$variable==0,]
wtdat0 <- data[data[,cat]==0,]
wtdat0$time = data[data[,cat]==0,duration]
wtdat0$status <-data[data[,cat]==0,]$Rec
tau= min(c(max(wtdat1$time),max(wtdat0$time)))
rmst1akme_1 <- rmst1_AKME(survdata =survdat1,wtdata =wtdat1,tau=tau)
rmst1akme_0 <- rmst1_AKME(survdata =survdat0,wtdata =wtdat0,tau=tau)
rmst1akme_1$rmst[1:2]
rmst1akme_0$rmst[1:2]
armst <- ARMST_diff(rmst_trt = rmst1akme_1$rmst[1],
                    rmst_ctrl = rmst1akme_0$rmst[1],
                    rmst_var_trt = rmst1akme_1$rmst.var,
                    rmst_var_ctrl = rmst1akme_0$rmst.var)
add2=armst$rmst.diff.result
row.names(add2)=cat
dfw2= add2
add3=armst$rmst.ratio.result
row.names(add3)=cat
dfw3= add3
add5=t( data.frame(rmst1akme_1$rmst, rmst1akme_0$rmst))
row.names(add5)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw5=add5

## repeat
cat='Stage34'
covs=c('Bev','Optimal','Age_high','Serous')

# unadjusted
rmst_result_sub1_unadj =rmst2(data[,duration],data$Rec,data[,cat])
add1=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[1,] ))
row.names(add1)=cat
dfw1= rbind(dfw1,add1)
add0=t( data.frame(rmst_result_sub1_unadj$RMST.arm1$rmst, rmst_result_sub1_unadj$RMST.arm0$rmst))
row.names(add0)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw0= rbind(dfw0,add0)
add4=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[2,] ))
row.names(add4)=cat
dfw4= rbind(dfw4,add4)
# adjusted
Pr1 <- glm( Stage34 ~ Bev + Optimal + Age_high+ Serous,
            data=data,family=binomial(link = "logit"))$fitted.values
W <- ( data[,cat]==1) * (1/Pr1) +( data[,cat]==0) * (1)/(1-Pr1)
res.akm <- ipw.survival(times=data[,duration],failures=data$Rec,
                        variable= data[,cat], weights=W)$table.surv
data$W <- W
survdat1 <- res.akm[res.akm$variable==1,]
wtdat1 <- data[data[,cat]==1,]
wtdat1$time =data[data[,cat]==1,duration]
wtdat1$status = data[data[,cat]==1,]$Rec
survdat0 <- res.akm[res.akm$variable==0,]
wtdat0 <- data[data[,cat]==0,]
wtdat0$time = data[data[,cat]==0,duration]
wtdat0$status <-data[data[,cat]==0,]$Rec
tau= min(c(max(wtdat1$time),max(wtdat0$time)))
rmst1akme_1 <- rmst1_AKME(survdata =survdat1,wtdata =wtdat1,tau=tau)
rmst1akme_0 <- rmst1_AKME(survdata =survdat0,wtdata =wtdat0,tau=tau)
armst <- ARMST_diff(rmst_trt = rmst1akme_1$rmst[1],
                    rmst_ctrl = rmst1akme_0$rmst[1],
                    rmst_var_trt = rmst1akme_1$rmst.var,
                    rmst_var_ctrl = rmst1akme_0$rmst.var)
add2=armst$rmst.diff.result
row.names(add2)=cat
dfw2= rbind(dfw2,add2)
add3=armst$rmst.ratio.result
row.names(add3)=cat
dfw3= rbind(dfw3,add3)
add5=t( data.frame(rmst1akme_1$rmst, rmst1akme_0$rmst))
row.names(add5)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw5= rbind(dfw5,add5)

## repeat
cat='Optimal'
covs=c('Bev','Stage34','Age_high','Serous')
# unadjusted
rmst_result_sub1_unadj =rmst2(data[,duration],data$Rec,data[,cat])
add1=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[1,] ))
row.names(add1)=cat
dfw1= rbind(dfw1,add1)
add0=t( data.frame(rmst_result_sub1_unadj$RMST.arm1$rmst, rmst_result_sub1_unadj$RMST.arm0$rmst))
row.names(add0)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw0= rbind(dfw0,add0)
add4=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[2,] ))
row.names(add4)=cat
dfw4= rbind(dfw4,add4)
# adjusted
Pr1 <- glm( Optimal ~ Bev + Stage34 + Age_high + Serous,
            data=data,family=binomial(link = "logit"))$fitted.values
W <- ( data[,cat]==1) * (1/Pr1) +( data[,cat]==0) * (1)/(1-Pr1)
res.akm <- ipw.survival(times=data[,duration],failures=data$Rec,
                        variable= data[,cat], weights=W)$table.surv
data$W <- W
survdat1 <- res.akm[res.akm$variable==1,]
wtdat1 <- data[data[,cat]==1,]
wtdat1$time =data[data[,cat]==1,duration]
wtdat1$status = data[data[,cat]==1,]$Rec
survdat0 <- res.akm[res.akm$variable==0,]
wtdat0 <- data[data[,cat]==0,]
wtdat0$time = data[data[,cat]==0,duration]
wtdat0$status <-data[data[,cat]==0,]$Rec
tau= min(c(max(wtdat1$time),max(wtdat0$time)))
rmst1akme_1 <- rmst1_AKME(survdata =survdat1,wtdata =wtdat1,tau=tau)
rmst1akme_0 <- rmst1_AKME(survdata =survdat0,wtdata =wtdat0,tau=tau)
armst <- ARMST_diff(rmst_trt = rmst1akme_1$rmst[1],
                    rmst_ctrl = rmst1akme_0$rmst[1],
                    rmst_var_trt = rmst1akme_1$rmst.var,
                    rmst_var_ctrl = rmst1akme_0$rmst.var)
add2=armst$rmst.diff.result
row.names(add2)=cat
dfw2= rbind(dfw2,add2)
add3=armst$rmst.ratio.result
row.names(add3)=cat
dfw3= rbind(dfw3,add3)
add5=t( data.frame(rmst1akme_1$rmst, rmst1akme_0$rmst))
row.names(add5)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw5= rbind(dfw5,add5)

## repeat
cat='Age_high'
covs=c('Bev','Stage34','Optimal','Serous')
# unadjusted
rmst_result_sub1_unadj =rmst2(data[,duration],data$Rec,data[,cat])
add1=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[1,] ))
row.names(add1)=cat
dfw1= rbind(dfw1,add1)
add0=t( data.frame(rmst_result_sub1_unadj$RMST.arm1$rmst, rmst_result_sub1_unadj$RMST.arm0$rmst))
row.names(add0)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw0= rbind(dfw0,add0)
add4=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[2,] ))
row.names(add4)=cat
dfw4= rbind(dfw4,add4)
# adjusted
Pr1 <- glm( Age_high ~ Bev + Stage34 + Optimal + Serous,
            data=data,family=binomial(link = "logit"))$fitted.values
W <- ( data[,cat]==1) * (1/Pr1) +( data[,cat]==0) * (1)/(1-Pr1)
res.akm <- ipw.survival(times=data[,duration],failures=data$Rec,
                        variable= data[,cat], weights=W)$table.surv
data$W <- W
survdat1 <- res.akm[res.akm$variable==1,]
wtdat1 <- data[data[,cat]==1,]
wtdat1$time =data[data[,cat]==1,duration]
wtdat1$status = data[data[,cat]==1,]$Rec
survdat0 <- res.akm[res.akm$variable==0,]
wtdat0 <- data[data[,cat]==0,]
wtdat0$time = data[data[,cat]==0,duration]
wtdat0$status <-data[data[,cat]==0,]$Rec
tau= min(c(max(wtdat1$time),max(wtdat0$time)))
rmst1akme_1 <- rmst1_AKME(survdata =survdat1,wtdata =wtdat1,tau=tau)
rmst1akme_0 <- rmst1_AKME(survdata =survdat0,wtdata =wtdat0,tau=tau)
armst <- ARMST_diff(rmst_trt = rmst1akme_1$rmst[1],
                    rmst_ctrl = rmst1akme_0$rmst[1],
                    rmst_var_trt = rmst1akme_1$rmst.var,
                    rmst_var_ctrl = rmst1akme_0$rmst.var)
add2=armst$rmst.diff.result
row.names(add2)=cat
dfw2= rbind(dfw2,add2)
add3=armst$rmst.ratio.result
row.names(add3)=cat
dfw3= rbind(dfw3,add3)
add5=t( data.frame(rmst1akme_1$rmst, rmst1akme_0$rmst))
row.names(add5)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw5= rbind(dfw5,add5)

## repeat
cat='Serous'
covs=c('Bev','Stage34','Optimal','Age_high')
# unadjusted
rmst_result_sub1_unadj =rmst2(data[,duration],data$Rec,data[,cat])
add1=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[1,] ))
row.names(add1)=cat
dfw1= rbind(dfw1,add1)
add0=t( data.frame(rmst_result_sub1_unadj$RMST.arm1$rmst, rmst_result_sub1_unadj$RMST.arm0$rmst))
row.names(add0)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw0= rbind(dfw0,add0)
add4=t(as.data.frame( rmst_result_sub1_unadj$unadjusted.result[2,] ))
row.names(add4)=cat
dfw4= rbind(dfw4,add4)
# adjusted
Pr1 <- glm( Serous ~ Bev + Stage34 + Optimal + Age_high,
            data=data,family=binomial(link = "logit"))$fitted.values
W <- ( data[,cat]==1) * (1/Pr1) +( data[,cat]==0) * (1)/(1-Pr1)
res.akm <- ipw.survival(times=data[,duration],failures=data$Rec,
                        variable= data[,cat], weights=W)$table.surv
data$W <- W
survdat1 <- res.akm[res.akm$variable==1,]
wtdat1 <- data[data[,cat]==1,]
wtdat1$time =data[data[,cat]==1,duration]
wtdat1$status = data[data[,cat]==1,]$Rec
survdat0 <- res.akm[res.akm$variable==0,]
wtdat0 <- data[data[,cat]==0,]
wtdat0$time = data[data[,cat]==0,duration]
wtdat0$status <-data[data[,cat]==0,]$Rec
tau= min(c(max(wtdat1$time),max(wtdat0$time)))
rmst1akme_1 <- rmst1_AKME(survdata =survdat1,wtdata =wtdat1,tau=tau)
rmst1akme_0 <- rmst1_AKME(survdata =survdat0,wtdata =wtdat0,tau=tau)
armst <- ARMST_diff(rmst_trt = rmst1akme_1$rmst[1],
                    rmst_ctrl = rmst1akme_0$rmst[1],
                    rmst_var_trt = rmst1akme_1$rmst.var,
                    rmst_var_ctrl = rmst1akme_0$rmst.var)
add2=armst$rmst.diff.result
row.names(add2)=cat
dfw2= rbind(dfw2,add2)
add3=armst$rmst.ratio.result
row.names(add3)=cat
dfw3= rbind(dfw3,add3)
add5=t( data.frame(rmst1akme_1$rmst, rmst1akme_0$rmst))
row.names(add5)=c( paste0(cat,'_arm1'),paste0(cat,'_arm0'))
dfw5= rbind(dfw5,add5)

# Check
dfw1
dfw4
dfw2
dfw3
dfw0
dfw5
out_prefix
# Write out
write.table(dfw0,paste0(out_prefix,'_unadjusted_raw.tsv'),sep='\t',quote = F,
            row.names = T,col.names = NA)
write.table(dfw1,paste0(out_prefix,'_unadjusted_diff.tsv'),sep='\t',quote = F,
            row.names = T,col.names = NA)
write.table(dfw4,paste0(out_prefix,'_unadjusted_ratio.tsv'),sep='\t',quote = F,
            row.names = T,col.names = NA)
write.table(dfw5,paste0(out_prefix,'_adjusted_raw.tsv'),sep='\t',quote = F,
            row.names = T,col.names = NA)
write.table(dfw2,paste0(out_prefix,'_adjusted_diff.tsv'),sep='\t',quote = F,
            row.names = T,col.names = NA)
write.table(dfw3,paste0(out_prefix,'_adjusted_ratio.tsv'),sep='\t',quote = F,
            row.names = T,col.names = NA)

