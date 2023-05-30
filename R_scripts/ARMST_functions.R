## Function from #https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fpst.2151&file=pst2151-sup-0001-Supinfo.pdf
######## ARMST difference and ratio
ARMST_diff <- function(rmst_trt=rmst_trt,
                       rmst_ctrl=rmst_ctrl,
                       rmst_var_trt=rmst_var_trt,
                       rmst_var_ctrl=rmst_var_ctrl,
                       alpha=0.05){
  rmst.diff.10 = rmst_trt - rmst_ctrl
  rmst.diff.10.se = sqrt(rmst_var_ctrl + rmst_var_trt)
  rmst.diff.10.low = rmst.diff.10 - qnorm(1 - alpha/2) * rmst.diff.10.se
  rmst.diff.10.upp = rmst.diff.10 + qnorm(1 - alpha/2) * rmst.diff.10.se
  rmst.diff.pval = pnorm(-abs(rmst.diff.10)/rmst.diff.10.se) * 2
  rmst.diff.result=cbind(round(rmst.diff.10,6),round(rmst.diff.10.se,6),
                         round(rmst.diff.10.low,6), round(rmst.diff.10.upp,6),
                         round(rmst.diff.pval,6))
  rmst.log.ratio.10 = log(rmst_trt) - log(rmst_ctrl)
  rmst.log.ratio.10.se=
    sqrt(rmst_var_trt/rmst_trt/rmst_trt + rmst_var_ctrl/rmst_ctrl/rmst_ctrl)
  rmst.log.ratio.10.low =
    rmst.log.ratio.10 - qnorm(1-alpha/2)*rmst.log.ratio.10.se
  rmst.log.ratio.10.upp =
    rmst.log.ratio.10 + qnorm(1-alpha/2)*rmst.log.ratio.10.se
  rmst.log.ratio.pval=pnorm(-abs(rmst.log.ratio.10)/rmst.log.ratio.10.se)*2
  
  rmst.ratio.result = 
    cbind(round(exp(rmst.log.ratio.10),6),round((rmst.log.ratio.10.se),6),
          round(exp(rmst.log.ratio.10.low),6),round(exp(rmst.log.ratio.10.upp),6),
          round(rmst.log.ratio.pval,6))
  
  colnames(rmst.ratio.result)=c("est","se",'lcl','ucl','p-value')
  colnames(rmst.diff.result) <- c('est','se','lcl','ucl','p-value')
  
  return_list <- list(rmst.diff.result = rmst.diff.result,
                      rmst.ratio.result = rmst.ratio.result)
  return(return_list)
}

#ARMST function
#This is a modification of the R function ‘rmst1’ in ‘survRM2’
rmst1_AKME=function(survdata, wtdata, tau, alpha=0.05){
  #-- time
  #-- statuts
  #-- tau -- truncation time
  #-- alpha -- gives (1-alpha) confidence interval
  #-- W weighting from AKME
  wtdata<-wtdata[order(wtdata$time),]
  idx=survdata$time<=tau
  wk.time=sort(c(survdata$time[idx],tau))
  w1 <-sapply(wk.time,function(x){sum(wtdata$W[wtdata$time >= x])})
  w2 <-sapply(wk.time,function(x){sum((wtdata$W[wtdata$time >= x])^2)})
  wk.surv=survdata$survival[idx]
  wk.n.risk =survdata$n.risk[idx]
  wk.n.event=survdata$n.event[idx]
  time.diff <- diff(c(0, wk.time))
  areas <- time.diff * c(1, wk.surv)
  rmst = sum(areas)
  rmst
  wk.var <- ifelse((wk.n.risk-wk.n.event)==0, 0,
                   wk.n.event /((wk.n.risk - wk.n.event)))
  wk.var =c(wk.var,0)
  rmst.var =
    sum(((cumsum(rev(areas[-1]))^2*rev(wk.var)[-1])*rev(w2)[-1])/rev(w1)[-1]^2)
  rmst.se = sqrt(rmst.var)
  #--- check ---
  # print(ft, rmean=tau)
  #--- output ---
  out=matrix(0,2,4)
  out[1,]=c(rmst, rmst.se, rmst-qnorm(1-alpha/2)*rmst.se,
            rmst+qnorm(1-alpha/2)*rmst.se)
  out[2,]=c(tau-out[1,1], rmst.se, tau-out[1,4], tau-out[1,3])
  rownames(out)=c("RMST","RMTL")
  colnames(out)=c("Est.","se", paste('lower .', round((1-alpha)*100,digits=0), sep = ''),
                  paste("upper .",round((1-alpha)*100, digits=0), sep=""))
  Z=list()
  Z$result=out
  Z$rmst = out[1,]
  Z$rmtl = out[2,]
  Z$tau=tau
  Z$rmst.var = rmst.var
  return(Z)
}

