# Simulation.function.r : a script to simulate assembly process (reflected in species abundances) in a community
# Author: ENAT (Adapted from Ning et al., 2020)
# Date: June 3, 2025

# NOTE.
# In this script: 
# FH represent condition with butyrate; FL: acetate
# Nk = number of individual in the community
# op = species pool
# spnumi = species number
# meta.ab = metacommunity (species w/ simulated abundance); in this study: following zero-sum multinominal distribution

# Function to simulate selection; 
# This function follows Gaussian distribution as a "fitness kernel".
# The peak of the bell curve represent the optimal trait/condition. 
# Thus, species closest to the mean are most fit (selected); Species in the tails are less likely to be selected.

# Notes:
# op = optimum niche, a value (0-1) assigned based on phylogeny tree that could be used to associate fitness trait to a species. 
# (e.g., species with op around 0.8 could survive high temperature)
# Nk = number of individual in each sample
# FL and FH = the "trait" that will be selected. In our scenario, species with op ~ 0.95 will be selected in butyrate condition.
# sig2. L and sig2.H = spread (variance) of the Gaussian function — how strongly selective the environment is. 

sim_select<-function(op,Nk,sampname,code.wd,
                     FL=0.05,FH=0.95,sig2.L=0.015,sig2.H=0.015)
{
  # 2.1 # totally deterministic, no variation among samples
  fit.L=exp(-((op-FL)^2)/(2*sig2.L))
  fit.H=exp(-((op-FH)^2)/(2*sig2.H))
  com.fit=rbind(matrix(fit.L,nr=10,nc=length(op),byrow = TRUE),
                matrix(fit.H,nr=10,nc=length(op),byrow = TRUE))
  source(paste0(code.wd,"/int.round.r"))
  comm=t(sapply(1:nrow(com.fit),function(i){int.round(com.fit[i,],sum.exp = Nk)}))
  dim(comm)
  rownames(comm)=sampname
  colnames(comm)=names(op)
  comm=comm[,colSums(comm)>0]
  #rowSums(comm);rowSums(comm>0)
  comm
}

# Function to simulate dispersal
dispersal.ZSM.sloan<-function(op,Nk,sampname,spnumi=100,meta.ab,
                              m1=0.01,m2=0.99,distinct=FALSE,fix.rich=TRUE)
{
  # 3.3 # Sloan's dispersal model, meta community per MZSM  
  pi1=meta.ab/sum(meta.ab)
  
  disp.sloan<-function(n,Nk,m1,m2=NULL,pi1,pi2=NULL,spnumi,fix.rich=TRUE)
  {
    if(is.null(m2))
    {
      yis=sapply(1:length(pi1),function(i){rgamma(n,Nk*m1*pi1[i])})
    }else{
      yis=sapply(1:length(pi1),function(i){rgamma(n,Nk*(m1*pi1[i]+m2*pi2[i]))})
    }
    yis=matrix(yis,nr=n,nc=length(pi1))
    xis=yis/rowSums(yis)
    t(sapply(1:nrow(xis),
             function(i)
             {
               outi=rep(0,ncol(xis))
               if(fix.rich)
               {
                 idi=sample(1:ncol(xis),spnumi,replace = FALSE,prob = xis[i,])
                 outi[idi]=1
                 rdi=sample(idi,Nk-spnumi,replace = TRUE,prob = xis[i,idi])
               }else{
                 rdi=sample(1:ncol(xis),Nk,replace = TRUE, prob = xis[i,])
               }
               tabi=table(rdi)
               outi[as.numeric(names(tabi))]=outi[as.numeric(names(tabi))]+as.vector(tabi)
               outi
             }))
  }
  if(distinct)
  {
    pool.A=disp.sloan(n=1,Nk=Nk,m1=m1,pi1=pi1,spnumi = spnumi,fix.rich=fix.rich)[1,]
    pi1b=pi1;pi1b[which(pool.A>0)]=0
    pool.B=disp.sloan(n=1,Nk=Nk,m1=m1,pi1=pi1b,spnumi = spnumi,fix.rich=fix.rich)[1,]
  }else{
    pool.AB=disp.sloan(n=2,Nk=Nk,m1=m1,pi1=pi1,spnumi = spnumi,fix.rich=fix.rich)
    pool.A=pool.AB[1,]
    pool.B=pool.AB[2,]
  }
  
  pi2.A=pool.A/sum(pool.A)
  pi2.B=pool.B/sum(pool.B)
  
  comm=matrix(0,nr=length(sampname),nc=length(op))
  rownames(comm)=sampname;colnames(comm)=names(op)
  dim(comm)
  
  comm["Almere_butyrate1",]=pool.A
  comm["Bath_butyrate1",]=pool.B
  comm[grep("Almere_butyrate[2-5]",rownames(comm)),]=disp.sloan(n=4,Nk=Nk,m1=m1,m2=m2,pi1 = pi1,pi2 = pi2.A,spnumi = spnumi,fix.rich=fix.rich)
  comm[grep("Bath_butyrate1[2-5]",rownames(comm)),]=disp.sloan(n=4,Nk=Nk,m1=m1,m2=m2,pi1 = pi1,pi2 = pi2.B,spnumi = spnumi,fix.rich=fix.rich)
  comm[grep("Almere_acetate",rownames(comm)),]=disp.sloan(n=5,Nk=Nk,m1=m1,m2=m2,pi1 = pi1,pi2 = pi2.A,spnumi = spnumi,fix.rich=fix.rich)
  comm[grep("Bath_acetate",rownames(comm)),]=disp.sloan(n=5,Nk=Nk,m1=m1,m2=m2,pi1 = pi1,pi2 = pi2.B,spnumi = spnumi,fix.rich=fix.rich)
  rowSums(comm)
  dim(comm)
  comm=comm[,colSums(comm)>0]
  dim(comm)
  comm
}

# Function to simulate drift
drift.ZSM.sloan<-function(op,Nk,sampname,spnumi=100,meta.ab,
                          m1=0.5,m2=0.01,FH=0.95,FL=0.05,sig2.W=4)
{
  # 4.3 # Sloan's dispersal model: moderate dispersal, meta community per MZSM.
  pi1=meta.ab/sum(meta.ab)
  
  disp.sloan<-function(n,Nk,m1,m2=NULL,pi1,pi2=NULL,spnumi)
  {
    if(is.null(m2))
    {
      yis=sapply(1:length(pi1),function(i){rgamma(n,Nk*m1*pi1[i])})
    }else{
      yis=sapply(1:length(pi1),function(i){rgamma(n,Nk*(m1*pi1[i]+m2*pi2[i]))})
    }
    xis=matrix(yis/rowSums(matrix(yis,nr=n,nc=length(pi1))),nr=n,nc=length(pi1))
    t(sapply(1:nrow(xis),
             function(i)
             {
               outi=rep(0,ncol(xis))
               idi=sample(1:ncol(xis),spnumi,replace = FALSE,prob = xis[i,])
               outi[idi]=1
               rdi=sample(idi,Nk-spnumi,replace = TRUE,prob = xis[i,idi])
               tabi=table(rdi)
               outi[as.numeric(names(tabi))]=outi[as.numeric(names(tabi))]+as.vector(tabi)
               outi
             }))
  }
  
  pool.A=disp.sloan(n=1,Nk=Nk,m1=m1,pi1=pi1,spnumi = spnumi)[1,]
  
  pi2=pool.A/sum(pool.A)
  
  comm=matrix(0,nr=length(sampname),nc=length(op))
  rownames(comm)=sampname;colnames(comm)=names(op)
  dim(comm)
  
  comm[grep("butyrate",rownames(comm)),]=disp.sloan(n=10,Nk=Nk,m1=m1,m2=m2,pi1 = pi1,pi2 = pi2,spnumi = spnumi)
  comm[grep("acetate",rownames(comm)),]=disp.sloan(n=10,Nk=Nk,m1=m1,m2=m2,pi1 = pi1,pi2 = pi2,spnumi = spnumi)
  rowSums(comm)
  dim(comm)
  comm=comm[,colSums(comm)>0]
  dim(comm)
  comm
}



