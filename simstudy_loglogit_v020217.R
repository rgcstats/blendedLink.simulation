###########################################################################
# Simulation Study for Clark and Barr (2017):
#    A Blended Link Approach to Relative Risk Regression
# to appear in Statistical Methods in Medical Research
# This code was originally run using Rstudio server with R version 3.3.2
# on a 60 core high performance computer with operating system
# Fedora Linux 25.
# Tables for the main paper and supplementary materials are output
# as Latex files.
###########################################################################

###########################################################
# install and load packages
###########################################################

library(devtools)
# download my blendedLink package. Version 1.0 is on CRAN,
#  but the latest (1.1) isn't yet, so get it from Github.
# Version 1.1 will be on CRAN eventually
devtools::install_github("rgcstats/blendedLink_1.1")
# install other packages
install_version("logbin","2.0.1")
install_version("glm2","1.1.2")
install_version("doMC","1.3.4")
install_version("foreach","1.4.3")
install_version("reshape2","1.4.2")
install_version("dplyr","0.5.0")
library(blendedLink)
library(glm2)
library(logbin)
library(doMC)
library(foreach)
library(reshape2)
library(dplyr)

###########################################################
# create an objective function for use in calculating beta to 
#   give a specified relative risk - see paper for details
###########################################################

objfn.for.beta <- function(beta,target.RR,alpha,truelinkinv,upper.value=0.5,lower.value=-0.5){
  truelinkinv(alpha+beta*upper.value) / truelinkinv(alpha+beta*lower.value) - target.RR
}

###########################################################
# function to create a link object from a given model (specified as text)
###########################################################

choose.link <- function(modelname){
  switch(modelname,
         "loglogit0.8"=blendedLink("log","logit",0.8),
         "logit"=make.link("logit") ,
         "log"=make.link("log") ,
         "poisson"=make.link("log") ,
         "probit"=make.link("probit") ,
         "cloglog"=make.link("cloglog") )
}

###########################################################
# create a data frame of the scenarios to be simulated 
###########################################################

scenarios <- expand.grid(RR=c(1.2,1.5,3),median.prob=c(0.3,0.5),sampsize=c(100,500,1000),
                         dsn.x=c("binary","normal","t4","uniform"),
                         numcovariates=2,
                         truemodel=c("logit"),
                         method=c("logit","log","poisson","loglogit0.8","adaptive"),
                         stringsAsFactors=F)

###########################################################
# run the simulations in parallel
###########################################################

registerDoMC(cores=40)
t1 <- proc.time()[3]
R <- 1000

allsimresults <- foreach(scenario=1:nrow(scenarios),.combine=rbind) %dopar% {
  RR <- scenarios[scenario,"RR"]
  median.prob <- scenarios[scenario,"median.prob"]
  sampsize <- scenarios[scenario,"sampsize"]
  truemodel <- scenarios[scenario,"truemodel"]
  method <- scenarios[scenario,"method"]
  dsn.x <- scenarios[scenario,"dsn.x"]
  numcovariates <- scenarios[scenario,"numcovariates"]
  fullsimresults <- NULL
  # set up true and assumed link
  # calculate alpha and beta corresponding to supplied median prob and relative risk
  alpha <- choose.link(truemodel)$linkfun(median.prob)
  if(dsn.x!="binary")  beta <- uniroot(f=objfn.for.beta,interval=c(0,20),alpha=alpha,target.RR=RR,
                  truelinkinv=choose.link(truemodel)$linkinv)$root
  if(dsn.x=="binary")  beta <- uniroot(f=objfn.for.beta,interval=c(0,20),alpha=alpha,target.RR=RR,
                                       truelinkinv=choose.link(truemodel)$linkinv,upper.value=1,
                                       lower.value=-1)$root
  beta2 <- uniroot(f=objfn.for.beta,interval=c(0,20),alpha=alpha,target.RR=1.2,
                   truelinkinv=choose.link(truemodel)$linkinv)$root
  RR1a <- choose.link(truemodel)$linkinv(beta*0.5+alpha) / choose.link(truemodel)$linkinv(-0.5*beta+alpha) # x1 goes from -0.5 to 0.5, should be equal to RR
  RR1b <- choose.link(truemodel)$linkinv(beta*2+alpha) / choose.link(truemodel)$linkinv(beta*1+alpha) # x1 goes from 1 to 2
  RR1c <- choose.link(truemodel)$linkinv(beta*(-1)+alpha) / choose.link(truemodel)$linkinv(beta*(-2)+alpha) # x1 goes from -2 to -1
  RR1d <- choose.link(truemodel)$linkinv(beta*(1)+alpha) / choose.link(truemodel)$linkinv(beta*(-1)+alpha) # x1 goes from -1 to 1
  simresults <- data.frame(ahat=rep(NA,R),bhat1=NA,bhat2=NA,se.bhat1=NA,se.bhat2=NA,comptime=NA,
                           RRhat1a=NA,RRhat1b=NA,RRhat1c=NA,RRhat1d=NA,converged=NA,boundary=NA,
                           pval.cutover=NA)
  set.seed(5646)
  for(r in c(1:R)){
    x1 <- x2 <- NA
    x1 <- switch( dsn.x , "uniform"=runif(sampsize,min=-sqrt(3),max=sqrt(3)), "normal"=rnorm(sampsize,mean=0,sd=1) ,
                  "t4"=rt(sampsize,df=4)*sqrt((4-2)/4) ,
                  "binary"=sample(c(-1,1),sampsize,replace=T))
    if(numcovariates==2) x2 <- rnorm(sampsize,mean=0,sd=1)
    u <- runif(sampsize)
    if(numcovariates==1) lp <- alpha + beta*x1
    if(numcovariates==2) lp <- beta*x1 + beta2*x2 + alpha
    if(numcovariates==1) modformula <- y~x1
    if(numcovariates==2) modformula <- y~x1+x2
    p <- choose.link(truemodel)$linkinv(lp)
    y <- 1*(u<=p)
    glmfit <- NA
    t1.glmfit <- proc.time()[3]
    if(method=="logit"){
      glmfit <- glm(modformula,family=binomial("logit"))
      estlinkinv <- make.link("logit")$linkinv
    }
    if(method=="cloglog"){
      glmfit <- glm2(modformula,family=binomial("cloglog"))
      estlinkinv <- make.link("cloglog")$linkinv
    }
    if(method=="log"){
      glmfit <- logbin(modformula)
      estlinkinv <- make.link("log")$linkinv
    }
    if(method=="poisson"){
      glmfit <- glm(modformula,family=poisson("log"))
      estlinkinv <- make.link("log")$linkinv
    }
    if(method=="loglogit0.8"){
      glmfit <- glm(modformula,family=binomial(blendedLink("log","logit",0.8)))
      if(numcovariates==1) simresults[r,c("pval.cutover")] <- pvalue.cutover(data.frame(y=y,x1=x1),y~x1,"log","logit",0.8,eps=0.025)$pval
      if(numcovariates==2) simresults[r,c("pval.cutover")] <- pvalue.cutover(data.frame(y=y,x1=x1,x2=x2),y~x1+x2,"log","logit",0.8,eps=0.025)$pval
      estlinkinv <- blendedLink("log","logit",0.8)$linkinv
    }
    if(method=="adaptive"){
      if(numcovariates==1) simresults[r,c("pval.cutover")] <- pvalue.cutover(data.frame(y=y,x1=x1),y~x1,"log","logit",0.8,eps=0.025)$pval
      if(numcovariates==2) simresults[r,c("pval.cutover")] <- pvalue.cutover(data.frame(y=y,x1=x1,x2=x2),y~x1+x2,"log","logit",0.8,eps=0.025)$pval
      if(simresults$pval.cutover[r]<0.05){
        glmfit <- glm(modformula,family=binomial(link="logit"))
        estlinkinv <- make.link("logit")$linkinv
      } else{
        glmfit <- glm2(modformula,family=binomial(blendedLink("log","logit",0.8)))
        estlinkinv <- blendedLink("log","logit",0.8)$linkinv
      }
    }
    t2.glmfit <- proc.time()[3]
    simresults[r,"comptime"] <- t2.glmfit-t1.glmfit
    simresults[r,"converged"] <- glmfit$converged
    simresults[r,"boundary"] <- glmfit$boundary
    simresults[r,c("ahat","bhat1")] <- glmfit$coef[1:2]
    simresults[r,c("se.bhat1")] <- sqrt(vcov(glmfit)[2,2]) 
    simresults[r,"RRhat1a"] <- estlinkinv(simresults$ahat[r]+simresults$bhat1[r]*0.5) /
      estlinkinv(simresults$ahat[r]-simresults$bhat1[r]*0.5)
    simresults[r,"RRhat1b"] <- estlinkinv(simresults$ahat[r]+simresults$bhat1[r]*2) /
      estlinkinv(simresults$ahat[r]+simresults$bhat1[r]*1)
    simresults[r,"RRhat1c"] <- estlinkinv(simresults$ahat[r]+simresults$bhat1[r]*(-1)) /
      estlinkinv(simresults$ahat[r]+simresults$bhat1[r]*(-2))
    simresults[r,"RRhat1d"] <- estlinkinv(simresults$ahat[r]+simresults$bhat1[r]*(1)) /
      estlinkinv(simresults$ahat[r]+simresults$bhat1[r]*(-1))
  }
  cat("scenario ",scenario , "\n" )
  fullsimresults <- rbind( fullsimresults , 
                           data.frame( R=R , RR=RR , RR1a=RR1a , RR1b=RR1b , RR1c=RR1c , RR1d=RR1d ,
                                       sampsize=sampsize , median.prob=median.prob ,
                                       truemodel=truemodel , dsn.x=dsn.x,
                                       method=method , e.ahat=mean(simresults$ahat) , beta=beta,
                                       e.bhat1=mean(simresults$bhat1) , #e.bhat2=mean(simresults$bhat2) ,
                                       mse.bhat1=mean((simresults$bhat1-beta)^2) , 
                                       mse.RRhat1a=mean((simresults$RRhat1a-RR1a)^2) ,
                                       mse.RRhat1b=mean((simresults$RRhat1b-RR1b)^2) ,
                                       mse.RRhat1c=mean((simresults$RRhat1c-RR1c)^2) ,
                                       mse.RRhat1d=mean((simresults$RRhat1d-RR1d)^2) ,
                                       e.RRhat1a=mean(simresults$RRhat1a),
                                       e.RRhat1b=mean(simresults$RRhat1b),
                                       e.RRhat1c=mean(simresults$RRhat1c),
                                       e.RRhat1d=mean(simresults$RRhat1d),
                                       sd.RRhat1a=sqrt(var(simresults$RRhat1a)),
                                       sd.RRhat1b=sqrt(var(simresults$RRhat1b)),
                                       sd.RRhat1c=sqrt(var(simresults$RRhat1c)),
                                       sd.RRhat1d=sqrt(var(simresults$RRhat1d)),
                                       p.boundary=mean(simresults$boundary),
                                       p.converged=mean(simresults$converged),
                                       mean.comptime=mean(simresults$comptime),
                                       bhat1.noncover.pct=100 * mean(abs(simresults$bhat1-beta)/simresults$se.bhat1>1.96) ,
                                       type1.pct.cutover=100 * mean(simresults$pval.cutover<=0.05) ) )
  fullsimresults
}
allsimresults$bias.RRhat1a <- allsimresults$e.RRhat1a-allsimresults$RR1a
allsimresults$bias.RRhat1b <- allsimresults$e.RRhat1b-allsimresults$RR1b
allsimresults$bias.RRhat1c <- allsimresults$e.RRhat1c-allsimresults$RR1c
allsimresults$bias.RRhat1d <- allsimresults$e.RRhat1d-allsimresults$RR1d
allsimresults$rmse.RRhat1a <- sqrt(allsimresults$mse.RRhat1a)
allsimresults$rmse.RRhat1b <- sqrt(allsimresults$mse.RRhat1b)
allsimresults$rmse.RRhat1c <- sqrt(allsimresults$mse.RRhat1c)
allsimresults$rmse.RRhat1d <- sqrt(allsimresults$mse.RRhat1d)
allsimresults$bias.RRhat1c <- allsimresults$e.RRhat1c-allsimresults$RR1c
allsimresults$bias.RRhat1d <- allsimresults$e.RRhat1d-allsimresults$RR1d
t2 <- proc.time()[3]
t2-t1

save.image("simresults_010217.rdata")
#load("simresults_010217.rdata")

#####################################################
# Function to Produce a Latex Table
#####################################################

table.producer <- function(allsimresults,n,dsn.x,signal,whichRR){
  simresults <- allsimresults[(allsimresults$truemodel=="logit")&(allsimresults$sampsize==n)&
                                (allsimresults$dsn.x==dsn.x),]
  which.measure <- paste0(signal,".RRhat",whichRR)
  tableformula <- as.formula(paste0("median.prob+RR",whichRR,"~method"))
  simres.spread <- dcast(simresults,tableformula,value.var=which.measure)
  simres.spread[,-1] <- format(round(simres.spread[,-1],digits=3),nsmall=3)
  prob.reject <- dcast(simresults,tableformula,value.var="type1.pct.cutover")
  simres.spread <- cbind(simres.spread,p.reject=format(round(prob.reject$adaptive,digits=1),nsmall=1))
  simres.spread[,ncol(simres.spread)] <- paste0( simres.spread[,ncol(simres.spread)] , " \\\\ " )
  if(signal=="bias") longsignal <- "Biases"
  if(signal=="sd") longsignal <- "Standard deviations (SDs)"
  if(signal=="rmse") longsignal <- "Root mean squared errors (RMSEs)"
  if(dsn.x=="normal") long.dsn.x <- "$N(0,1)$"
  if(dsn.x=="uniform") long.dsn.x <- "$U(-\\sqrt{3},-\\sqrt{3})$"
  if(dsn.x=="t4") long.dsn.x <- "$t_4/\\sqrt{2}$"
  if(dsn.x=="binary") long.dsn.x <- "binary (values $\\pm 1$ with equal probability)"
  filename <- paste0("simtab_n",n,"_x",dsn.x,"_",signal,whichRR,".tex")
  if(whichRR=="1a"){xval1<- 0.5; xval0 <- -0.5}
  if(whichRR=="1b"){xval1<- 2; xval0 <- 1}
  if(whichRR=="1c"){xval1<- -1; xval0 <- -2}
  if(whichRR=="1d"){xval1<- 1; xval0 <- -1}
  cat("\\begin{table}[H] \n",file=filename)
  cat("\\small\\sf\\centering \n",file=filename,append=T)
  cat("\\caption{",longsignal," of relative risk (RR) estimators for $x_1=",xval1,"$ vs $x_1=",xval0,
      "$ when $\\mbox{logit}(P[Y=1])=\\alpha+\\beta_1 x_1 + \\beta_2 x_2$. Sample size is 1000. ",
      "$x_1 \\sim $",long.dsn.x,", $x_2 \\sim N(0,1)$. {\\it Median prob} refers to the probability that ",
      "$Y=1$ when both covariates are at their median value. {\\it Logit, log, poisson and log-logit} ",
      "refer to binary regression with logit link (logistic regression) and log link, ",
      "Poisson regression with log link and binary regression with blended log-logit link, respectively. ",
      "{\\it Prob.reject} is the simulation probability (\\%) of rejecting the null that the cutover ",
      "probability is $0.8$ in the blended log-logit model.} \n" , sep="" , file=filename,append=T) 
  cat("\\begin{tabular}{llrrrrr} \n",file=filename,append=T)
  cat("\\toprule \n",file=filename,append=T)
  cat(paste0("median & true RR & \\multicolumn{4}{c}{",signal," of estimated RR} & prob.reject \\\\ \n"),
      file=filename,append=T)  
  cat("prob. & & logit & log & poisson & log-logit  & (\\%) \\\\ \\midrule \n",file=filename,append=T)
  write.table(simres.spread[,c("median.prob",paste0("RR",whichRR),"logit","log","poisson","loglogit0.8","p.reject")],file=filename,sep=" & ",eol=" \n" , row.names=FALSE ,
              col.names=FALSE,quote=F,append=T)
  cat("\\bottomrule \n",file=filename,append=T)
  cat("\\end{tabular} \n",file=filename,append=T)
  cat("\\end{table} \n",file=filename,append=T)
  cat(filename,"\n")
}



###########################################################
# CREATE ALL SIMULATION TABLES AND SAVE AS LATEX FILES
###########################################################

# test table production function
table.producer(allsimresults,1000,"normal","bias","1a")
writeLines(readLines("simtab_n1000_xnormal_bias.RRhat1a.tex"))
writeLines(readLines("simtab_n1000_xnormal_bias.RRhat1a.tex"),con=file("test.tex"))


# produce all tables
for(n in c(100,500,1000)){
  for(dsn.x in c("normal","binary","t4","uniform")){
    for(signal in c("bias","sd","rmse")){
      for(whichRR in c("1a","1b","1c","1d")){
        table.producer(allsimresults,n=n,dsn.x=dsn.x,signal=signal,whichRR=whichRR)
      }
    }
  }
}

# create latex files for supplementary tables, by going through each n and signal
all.supp.filenames <- NULL
for(n in c(100,500,1000)){
  for(signal in c("bias","sd","rmse")){
    text <- c(readLines(paste0("simtab_n",n,"_xbinary_",signal,"1d.tex")),
              readLines(paste0("simtab_n",n,"_xnormal_",signal,"1a.tex")),
              readLines(paste0("simtab_n",n,"_xnormal_",signal,"1b.tex")),
              readLines(paste0("simtab_n",n,"_xnormal_",signal,"1c.tex")),
              readLines(paste0("simtab_n",n,"_xt4_",signal,"1a.tex")),
              readLines(paste0("simtab_n",n,"_xt4_",signal,"1b.tex")),
              readLines(paste0("simtab_n",n,"_xt4_",signal,"1c.tex")),
              readLines(paste0("simtab_n",n,"_xuniform_",signal,"1a.tex")),
              readLines(paste0("simtab_n",n,"_xuniform_",signal,"1b.tex")),
              readLines(paste0("simtab_n",n,"_xuniform_",signal,"1c.tex")) )
    if(signal=="bias") longsignal <- "Biases"
    if(signal=="sd") longsignal <- "Standard Deviations"
    if(signal=="rmse") longsignal <- "Root Mean Squared Errors"
    docheading <- paste0(longsignal," of Relative Risk Estimators in Simulation Study when $n=",n,"$ }")
    filename <- paste0("supp_tables_",signal,"_n",n,".tex")
    all.supp.filenames <- c(all.supp.filenames,filename)
    writeLines(c(readLines("preamble_supp_tables1.tex"),
                 docheading,
                 readLines("preamble_supp_tables2.tex"),
                 text,
                 "\\end{document}"),
               con=file(filename))
    close(con=file(filename))
  }
}

cat(rep(paste0("pdflatex ",all.supp.filenames," \n"),each=4),
    file="latex_calls.txt")


#############################################
# Plot of inverse link function for paper
#############################################

postscript("blendedlink.eps",onefile=T,width=6,height=5,horizontal=F)
curve(blendedLink("log","logit",0.8)$linkinv(eta),from=-3,to=3,xname="eta",xlab=expression(paste("linear predictor (",eta,")")),
      ylab=expression(paste("E[Y]=",mu=f(eta))),ylim=c(-0.01,1))
dev.off()

