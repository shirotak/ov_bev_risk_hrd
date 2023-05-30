# R script for RMST analysis
# Ref:https://cran.r-project.org/web/packages/survRM2/vignettes/survRM2-vignette3-2.html

library(survival)
library(survRM2)

# All
D=read.delim('./input/ICON7_all.tsv')
time   = D$PFS
status = D$Rec
arm    = D$Bev
tau=min(max(D$PFS[D$Bev==1]),max(D$PFS[D$Bev==0]))
time=pmin(time,tau)
obj = rmst2(time/30, status = status, arm = arm)
res=obj$unadjusted.result
res
pdf(file = './R_results/All.pdf',height = 1.7,width = 3.2)
par(mar=c(3, 2, 0.5, 0.1))
plot(obj, xlab='', ylab='',
     col.RMTL = "skyblue",col = 'black',col.main='white',
     cex.sub=1.1,cex.axis=0.8,cex.main=0.1,
     mgp=c(1,0.5,0))
dev.off()
write.table(res,'./R_results/All_RMST.tsv',sep='\t',quote = F,
            col.names = NA,row.names = T)

## Before
D=read.delim('./input/ICON7_all_before.tsv')
time   = D$PFS
time=pmin(time,54*7)
status = D$Rec
arm    = D$Bev
obj = rmst2(time/30, status = status, arm = arm)
res=obj$unadjusted.result
res
pdf(file = './R_results/All_before.pdf',height = 1.7,width = 3.2)
par(mar=c(3, 2, 0.5, 0.1))
plot(obj, xlab='', ylab='',
     col.RMTL = "skyblue",col = 'black',col.main='white',
     cex.sub=1.1,cex.axis=0.8,cex.main=0.1,
     mgp=c(1,0.5,0))
dev.off()
write.table(res,'./R_results/All_before_RMST.tsv',sep='\t',quote = F,
            col.names = NA,row.names = T)

## After
D=read.delim('./input/ICON7_all_after.tsv')
time   = D$PFS_r
status = D$Rec
arm    = D$Bev
tau=min(max(D$PFS_r[D$Bev==1]),max(D$PFS_r[D$Bev==0]))
time=pmin(time,tau)
obj = rmst2(time/30, status, arm)
res=obj$unadjusted.result
res
pdf(file = './R_results/All_after.pdf',height = 1.7,width = 3.2)
par(mar=c(3, 2, 0.5, 0.1))
plot(obj, xlab='', ylab='',
     col.RMTL = "skyblue",col = 'black',col.main='white',
     cex.sub=1.1,cex.axis=0.8,cex.main=0.1,
     mgp=c(1,0.5,0))
dev.off()
write.table(res,'./R_results/All_after_RMST.tsv',sep='\t',quote = F,
            col.names = NA,row.names = T)

# Serous
D=read.delim('./input/ICON7_serous.tsv')
time   = D$PFS
status = D$Rec
arm    = D$Bev
obj = rmst2(time/30, status = status, arm = arm)
res=obj$unadjusted.result
res
pdf(file = './R_results/Serous_all.pdf',height = 1.7,width = 3.2)
par(mar=c(3, 2, 0.5, 0.1))
plot(obj, xlab='', ylab='',
     col.RMTL = "skyblue",col = 'black',col.main='white',
     cex.sub=1.1,cex.axis=0.8,cex.main=0.1,
     mgp=c(1,0.5,0))
dev.off()
write.table(res,'./R_results/Serous_all_RMST.tsv',sep='\t',quote = F,
            col.names = NA,row.names = T)

D=read.delim('./input/ICON7_serous_before.tsv')
time   = D$PFS
time=pmin(time,54*7)
status = D$Rec
arm    = D$Bev
tau    = 54*7
obj = rmst2(time/30, status = status, arm = arm,tau = tau/30)
res=obj$unadjusted.result
res
pdf(file = './R_results/Serous_before.pdf',height = 1.7,width = 3.2)
par(mar=c(3, 2, 0.5, 0.1))
plot(obj, xlab='', ylab='',
     col.RMTL = "skyblue",col = 'black',col.main='white',
     cex.sub=1.1,cex.axis=0.8,cex.main=0.1,
     mgp=c(1,0.5,0))
dev.off()
write.table(res,'./R_results/Serous_before_RMST.tsv',sep='\t',quote = F,
            col.names = NA,row.names = T)

D=read.delim('./input/ICON7_serous_after.tsv')
time   = D$PFS_r
status = D$Rec
arm    = D$Bev
tau=min(max(D$PFS_r[D$Bev==1]),max(D$PFS_r[D$Bev==0]))
time=pmin(time,tau)
obj = rmst2(time/30, status, arm,tau = tau/30)
res=obj$unadjusted.result
res
pdf(file = './R_results/Serous_after.pdf',height = 1.7,width = 3.2)
par(mar=c(3, 2, 0.5, 0.1))
plot(obj, xlab='', ylab='',
     col.RMTL = "skyblue",col = 'black',col.main='white',
     cex.sub=1.1,cex.axis=0.8,cex.main=0.1,
     mgp=c(1,0.5,0))
dev.off()
write.table(res,'./R_results/Serou_after_RMST.tsv',sep='\t',quote = F,
            col.names = NA,row.names = T)

# NonSerous
D=read.delim('./input/ICON7_nonserous.tsv')
time   = D$PFS
status = D$Rec
arm    = D$Bev
obj = rmst2(time/30, status = status, arm = arm)
res=obj$unadjusted.result
res
pdf(file = './R_results/NonSerous_all.pdf',height = 1.7,width = 3.2)
par(mar=c(3, 2, 0.5, 0.1))
plot(obj, xlab='', ylab='',
     col.RMTL = "skyblue",col = 'black',col.main='white',
     cex.sub=1.1,cex.axis=0.8,cex.main=0.1,
     mgp=c(1,0.5,0))
dev.off()
write.table(res,'./R_results/NonSerous_all_RMST.tsv',sep='\t',quote = F,
            col.names = NA,row.names = T)

D=read.delim('./input/ICON7_nonserous.tsv')
time   = D$PFS
time=pmin(time,54*7)
status = D$Rec
arm    = D$Bev
tau    = 54*7
obj = rmst2(time/30, status = status, arm = arm,tau = tau/30)
res=obj$unadjusted.result
res
pdf(file = './R_results/NonSerous_before.pdf',height = 1.7,width = 3.2)
par(mar=c(3, 2, 0.5, 0.1))
plot(obj, xlab='', ylab='',
     col.RMTL = "skyblue",col = 'black',col.main='white',
     cex.sub=1.1,cex.axis=0.8,cex.main=0.1,
     mgp=c(1,0.5,0))
dev.off()
write.table(res,'./R_results/NonSerous_before_RMST.tsv',sep='\t',quote = F,
            col.names = NA,row.names = T)

D=read.delim('./input/ICON7_nonserous_after.tsv')
time   = D$PFS_r
status = D$Rec
arm    = D$Bev
tau=min(max(D$PFS_r[D$Bev==1]),max(D$PFS_r[D$Bev==0]))
time=pmin(time,tau)
obj = rmst2(time/30, status, arm,tau = tau/30)
res=obj$unadjusted.result
res
pdf(file = './R_results/NonSerous_after.pdf',height = 1.7,width = 3.2)
par(mar=c(3, 2, 0.5, 0.1))
plot(obj, xlab='', ylab='',
     col.RMTL = "skyblue",col = 'black',col.main='white',
     cex.sub=1.1,cex.axis=0.8,cex.main=0.1,
     mgp=c(1,0.5,0))
dev.off()
write.table(res,'./R_results/NonSerous_after_RMST.tsv',sep='\t',quote = F,
            col.names = NA,row.names = T)

# serous HRD
D=read.delim('./input/ICON7_HRD.tsv')
time   = D$PFS
status = D$Rec
arm    = D$Bev
obj = rmst2(time/30, status = status, arm = arm)
res=obj$unadjusted.result
res
pdf(file = './R_results/HRD_all.pdf',height = 1.7,width = 3.2)
par(mar=c(3, 2, 0.5, 0.1))
plot(obj, xlab='', ylab='',
     col.RMTL = "skyblue",col = 'black',col.main='white',
     cex.sub=1.1,cex.axis=0.8,cex.main=0.1,
     mgp=c(1,0.5,0))
dev.off()
write.table(res,'./R_results/HRD_all_RMST.tsv',sep='\t',quote = F,
            col.names = NA,row.names = T)

D=read.delim('./input/ICON7_HRD_before.tsv')
time   = D$PFS
time=pmin(time,54*7)
status = D$Rec
arm    = D$Bev
tau    = 54*7
obj = rmst2(time/30, status = status, arm = arm,tau = tau/30)
res=obj$unadjusted.result
res
pdf(file = './R_results/HRD_before.pdf',height = 1.7,width = 3.2)
par(mar=c(3, 2, 0.5, 0.1))
plot(obj, xlab='', ylab='',
     col.RMTL = "skyblue",col = 'black',col.main='white',
     cex.sub=1.1,cex.axis=0.8,cex.main=0.1,
     mgp=c(1,0.5,0))
dev.off()
write.table(res,'./R_results/HRD_before_RMST.tsv',sep='\t',quote = F,
            col.names = NA,row.names = T)

D=read.delim('./input/ICON7_HRD_after.tsv')
time   = D$PFS_r
status = D$Rec
arm    = D$Bev
tau=min(max(D$PFS_r[D$Bev==1]),max(D$PFS_r[D$Bev==0]))
time=pmin(time,tau)
obj = rmst2(time/30, status, arm,tau = tau/30)
res=obj$unadjusted.result
res
pdf(file = './R_results/HRD_after.pdf',height = 1.7,width = 3.2)
par(mar=c(3, 2, 0.5, 0.1))
plot(obj, xlab='', ylab='',
     col.RMTL = "skyblue",col = 'black',col.main='white',
     cex.sub=1.1,cex.axis=0.8,cex.main=0.1,
     mgp=c(1,0.5,0))
dev.off()
write.table(res,'./R_results/HRD_after_RMST.tsv',sep='\t',quote = F,
            col.names = NA,row.names = T)

# nonHRD
D=read.delim('./input/ICON7_nonHRD.tsv')
time   = D$PFS
status = D$Rec
arm    = D$Bev
obj = rmst2(time/30, status = status, arm = arm)
res=obj$unadjusted.result
res
pdf(file = './R_results/NonHRD_all.pdf',height = 1.7,width = 3.2)
par(mar=c(3, 2, 0.5, 0.1))
plot(obj, xlab='', ylab='',
     col.RMTL = "skyblue",col = 'black',col.main='white',
     cex.sub=1.1,cex.axis=0.8,cex.main=0.1,
     mgp=c(1,0.5,0))
dev.off()
write.table(res,'./R_results/NonHRD_all_RMST.tsv',sep='\t',quote = F,
            col.names = NA,row.names = T)

D=read.delim('./input/ICON7_nonHRD_before.tsv')
time   = D$PFS
time=pmin(time,54*7)
status = D$Rec
arm    = D$Bev
tau    = 54*7
obj = rmst2(time/30, status = status, arm = arm,tau = tau/30)
res=obj$unadjusted.result
res
pdf(file = './R_results/NonHRD_before.pdf',height = 1.7,width = 3.2)
par(mar=c(3, 2, 0.5, 0.1))
plot(obj, xlab='', ylab='',
     col.RMTL = "skyblue",col = 'black',col.main='white',
     cex.sub=1.1,cex.axis=0.8,cex.main=0.1,
     mgp=c(1,0.5,0))
dev.off()
write.table(res,'./R_results/NonHRD_before_RMST.tsv',sep='\t',quote = F,
            col.names = NA,row.names = T)

D=read.delim('./input/ICON7_nonHRD_after.tsv')
time   = D$PFS_r
status = D$Rec
arm    = D$Bev
tau=min(max(D$PFS_r[D$Bev==1]),max(D$PFS_r[D$Bev==0]))
time=pmin(time,tau)
obj = rmst2(time/30, status, arm,tau = tau/30)
res=obj$unadjusted.result
res
pdf(file = './R_results/NonHRD_after.pdf',height = 1.7,width = 3.2)
par(mar=c(3, 2, 0.5, 0.1))
plot(obj, xlab='', ylab='',
     col.RMTL = "skyblue",col = 'black',col.main='white',
     cex.sub=1.1,cex.axis=0.8,cex.main=0.1,
     mgp=c(1,0.5,0))
dev.off()
write.table(res,'./R_results/NonHRD_after_RMST.tsv',sep='\t',quote = F,
            col.names = NA,row.names = T)
