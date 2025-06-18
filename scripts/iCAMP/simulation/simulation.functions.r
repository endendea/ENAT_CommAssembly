##########################################
# 1.2 # functions
##########################################
compet<-function(op,Nk,sampname,code.wd,sp.bin,
                FL=0.05,FH=0.95,sig2.L=0.015,sig2.H=0.015,
                topn=3,minder=0.3)
{
  fit.L=exp(-((op-FL)^2)/(2*sig2.L))
  fit.H=exp(-((op-FH)^2)/(2*sig2.H))
  bsop<-function(opi,fsp,dopi)
  {
    which.maxn<-function(v)
    {
      out=which(v==max(v))
      if(length(out)>1){out=sample(out,1)}
      out
    }
    fsp2=rownames(dopi)[which.maxn(dopi[,fsp])]
    outsp=c(fsp,fsp2)
    absp=0.5^(seq(from=0,to=length(outsp)-1))
    rmsp=names(opi)[which(!(names(opi) %in% outsp))]
    while(length(rmsp)>0)
    {
      dopim=dopi[outsp,rmsp,drop=FALSE]
      outsp=c(outsp,rmsp[which.maxn(colSums((dopim^2)*absp))])
      absp=0.5^(seq(from=0,to=length(outsp)-1))
      rmsp=names(opi)[which(!(names(opi) %in% outsp))]
    }
    names(absp)=outsp
    absp[match(names(opi),outsp)]
  }
  topi<-function(fiti,...)
  {
    outi=sort(fiti,decreasing = TRUE)[1:topn]
    outi[which(outi>=(outi[1]-minder))]
  }
  
  # to each bin
  bin.lev=unique(sp.bin[,1])
  abl=lapply(1:length(bin.lev),
             function(i)
             {
               #message("bin i=",i," in ",length(bin.lev),". ")
               spi=rownames(sp.bin)[sp.bin[,1]==bin.lev[i]]
               opi=op[match(spi,names(op))]
               fit.Li=fit.L[match(spi,names(fit.L))]
               fit.Hi=fit.H[match(spi,names(fit.H))]
               # pick up a focal species
               top.Li=topi(fit.Li)
               top.Hi=topi(fit.Hi)
               fit.fspi=c(sample(top.Li,size=12,replace = TRUE),
                          sample(top.Hi,size=12,replace = TRUE))
               fspis=names(fit.fspi)
               fspi.lev=unique(fspis)
               # assign break-stick abundance according to OP difference
               dopi=as.matrix(dist(opi))
               abim=sapply(fspi.lev,bsop,opi = opi,dopi = dopi)
               abi=(t(abim[,match(fspis,fspi.lev)]))*fit.fspi
               rownames(abi)=sampname
               abi
             })
  # assemble bins
  abm=Reduce(cbind,abl)
  source(paste0(code.wd,"/int.round.r"))
  comm=t(sapply(1:nrow(abm),function(i){int.round(abm[i,],sum.exp = Nk)}))
  rownames(comm)=sampname
  colnames(comm)=colnames(abm)
  comm=comm[,colSums(comm)>0]
  #rowSums(comm);rowSums(comm>0)
  comm
}

sim_select<-function(op,Nk,sampname,code.wd,
                 FL=0.05,FH=0.95,sig2.L=0.015,sig2.H=0.015)
{
  # 2.1 # totally deterministic, no variation among samples
  fit.L=exp(-((op-FL)^2)/(2*sig2.L))
  fit.H=exp(-((op-FH)^2)/(2*sig2.H))
  com.fit=rbind(matrix(fit.L,nr=5,nc=length(op),byrow = TRUE),
                matrix(fit.H,nr=5,nc=length(op),byrow = TRUE))
  source(paste0(code.wd,"/int.round.r"))
  comm=t(sapply(1:nrow(com.fit),function(i){int.round(com.fit[i,],sum.exp = Nk)}))
  rownames(comm)=sampname
  colnames(comm)=names(op)
  comm=comm[,colSums(comm)>0]
  #rowSums(comm);rowSums(comm>0)
  comm
}

com.sel<-function(sp.bin.sel, compet.sel, op, Nk.sel,sampname,code.wd)
{
  Nk.cmp=compet.sel*Nk.sel
  if(Nk.cmp>0)
  {
    bins=unique(sp.bin.sel[,1])
    bin.cmp=sample(bins,size =max(floor(length(bins)*compet.sel),1))
    
    sp.bin.cmp=sp.bin.sel[which(sp.bin.sel[,1] %in% bin.cmp),,drop=FALSE]
    op.cmp=op[rownames(sp.bin.cmp)]
    
    comm.cmp=compet(op = op.cmp, Nk = Nk.cmp, sampname = sampname, code.wd = code.wd,
                    sp.bin = sp.bin.cmp, FL = 0.05, FH = 0.95, sig2.L = 0.015,
                    sig2.H = 0.015,topn=3,minder=0.3)
  }else{
    comm.cmp=NULL
    bin.cmp=NULL
  }
  
  Nk.flt=Nk.sel-Nk.cmp
  if(Nk.flt>0)
  {
    sp.bin.flt=sp.bin.sel[which(!(sp.bin.sel[,1] %in% bin.cmp)),,drop=FALSE]
    op.flt=op[rownames(sp.bin.flt)]
    comm.flt=select(op = op.flt,Nk = Nk.flt,sampname = sampname,code.wd = code.wd,
                    FL=0.05,FH=0.95,sig2.L=0.015,sig2.H=0.015)
    comm.sel=cbind(comm.flt,comm.cmp)
  }else{
    comm.sel=comm.cmp
    comm.flt=NULL
  }
  list(comm.sel=comm.sel,comm.flt=comm.flt,comm.cmp=comm.cmp)
}





dispersal.stegen<-function(op,Nk,spnumi=100,sampname,
                           FH=0.95,FL=0.05,sig2.W=4,disp.a=0.05,disp.h=1.5)
{
  # 3.1 # Stegen's model
  fit.WH=exp(-((op-FH)^2)/(2*sig2.W))
  fit.WL=exp(-((op-FL)^2)/(2*sig2.W))
  
  com.rand<-function(sampn,spnumi,Nk,fits)
  {
    t(sapply(1:sampn,
             function(i)
             {
               out=rep(0,length(fits))
               idi=sample(1:length(fits),spnumi,replace = FALSE,prob = fits)
               out[idi]=1
               randi=sample(idi,Nk-spnumi,replace=TRUE,prob=fits[idi])
               tabi=table(randi)
               out[as.numeric(names(tabi))]=out[as.numeric(names(tabi))]+as.vector(tabi)
               out
             }))
  }
  pool.A=com.rand(sampn=1,spnumi=spnumi,Nk=Nk,fits=fit.WL)
  fit.WLB=fit.WL
  pool.B=com.rand(sampn=1,spnumi=spnumi,Nk=Nk,fits=fit.WLB)
  
  fit.LAs=fit.WL+disp.a*(pool.A[1,]^disp.h)
  fit.LBs=fit.WL+disp.a*(pool.B[1,]^disp.h)
  fit.HAs=fit.WH+disp.a*(pool.A[1,]^disp.h)
  fit.HBs=fit.WH+disp.a*(pool.B[1,]^disp.h)
  
  comm=matrix(0,nr=length(sampname),nc=length(op))
  rownames(comm)=sampname;colnames(comm)=names(op)
  dim(comm)
  
  comm["Almere1",]=pool.A
  comm["Bath1",]=pool.B
  comm[grep("Almere[2-5]",rownames(comm)),]=com.rand(sampn=5,spnumi = spnumi,Nk=Nk,fits=fit.LAs)
  comm[grep("Bath[2-5]",rownames(comm)),]=com.rand(sampn=5,spnumi = spnumi,Nk=Nk,fits=fit.LBs)
  comm[grep("HA",rownames(comm)),]=com.rand(sampn=6,spnumi = spnumi,Nk=Nk,fits=fit.HAs)
  comm[grep("HB",rownames(comm)),]=com.rand(sampn=6,spnumi = spnumi,Nk=Nk,fits=fit.HBs)
  rowSums(comm)
  dim(comm)
  comm=comm[,colSums(comm)>0]
  dim(comm)
  comm
}

dispersal.sloan<-function(op,Nk,sampname,spnumi=100,
                          m1=0.01,m2=0.99,FH=0.95,FL=0.05,sig2.W=4)
{
  # 3.2 # Sloan's dispersal model
  fit.WH=exp(-((op-FH)^2)/(2*sig2.W))
  fit.WL=exp(-((op-FL)^2)/(2*sig2.W))
  pi1.WL=fit.WL/sum(fit.WL)
  pi1.WH=fit.WH/sum(fit.WH)
  
  disp.sloan<-function(n,Nk,m1,m2=NULL,pi1,pi2=NULL,spnumi)
  {
    if(is.null(m2))
    {
      yis=sapply(1:length(pi1),function(i){rgamma(n,Nk*m1*pi1[i])})
    }else{
      yis=sapply(1:length(pi1),function(i){rgamma(n,Nk*(m1*pi1[i]+m2*pi2[i]))})
    }
    xis=yis/rowSums(yis)
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
  
  pool.AB=disp.sloan(n=2,Nk=Nk,m1=m1,pi1=pi1.WL,spnumi = spnumi)
  pool.A=pool.AB[1,]
  pool.B=pool.AB[2,]
  
  pi2.A=pool.A/sum(pool.A)
  pi2.B=pool.B/sum(pool.B)
  
  comm=matrix(0,nr=length(sampname),nc=length(op))
  rownames(comm)=sampname;colnames(comm)=names(op)
  dim(comm)
  
  comm["Almere1",]=pool.A
  comm["Bath1",]=pool.B
  comm[grep("Almere[2-5]",rownames(comm)),]=disp.sloan(n=5,Nk=Nk,m1=m1,m2=m2,pi1 = pi1.WL,pi2 = pi2.A,spnumi = spnumi)
  comm[grep("Bath[2-5]",rownames(comm)),]=disp.sloan(n=5,Nk=Nk,m1=m1,m2=m2,pi1 = pi1.WH,pi2 = pi2.B,spnumi = spnumi)
  rowSums(comm)
  dim(comm)
  comm=comm[,colSums(comm)>0]
  dim(comm)
  comm
}
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
  
  comm["Almere1",]=pool.A
  comm["Bath1",]=pool.B
  comm[grep("Almere[2-5]",rownames(comm)),]=disp.sloan(n=4,Nk=Nk,m1=m1,m2=m2,pi1 = pi1,pi2 = pi2.A,spnumi = spnumi,fix.rich=fix.rich)
  comm[grep("Bath[2-5]",rownames(comm)),]=disp.sloan(n=4,Nk=Nk,m1=m1,m2=m2,pi1 = pi1,pi2 = pi2.B,spnumi = spnumi,fix.rich=fix.rich)
  rowSums(comm)
  dim(comm)
  comm=comm[,colSums(comm)>0]
  dim(comm)
  comm
}


drift.stegen<-function(op,Nk,sampname,spnumi=100,
                       FH=0.95,FL=0.05,sig2.W=4,disp.a=0.05,disp.l=0.2)
{
  # 4.1 # Stegen's model: weak selection with moderate dispersal
  fit.WH=exp(-((op-FH)^2)/(2*sig2.W))
  fit.WL=exp(-((op-FL)^2)/(2*sig2.W))
  
  com.rand<-function(sampn,spnumi,Nk,fits)
  {
    t(sapply(1:sampn,
             function(i)
             {
               out=rep(0,length(fits))
               idi=sample(1:length(fits),spnumi,replace = FALSE,prob = fits)
               out[idi]=1
               randi=sample(idi,Nk-spnumi,replace=TRUE,prob=fits[idi])
               tabi=table(randi)
               out[as.numeric(names(tabi))]=out[as.numeric(names(tabi))]+as.vector(tabi)
               out
             }))
  }
  pool.A=com.rand(sampn=1,spnumi=spnumi,Nk=Nk,fits=fit.WL)
  
  fit.Ls=fit.WL+disp.a*(pool.A[1,]^disp.l)
  fit.Hs=fit.WH+disp.a*(pool.A[1,]^disp.l)
  
  comm=matrix(0,nr=length(sampname),nc=length(op))
  rownames(comm)=sampname;colnames(comm)=names(op)
  dim(comm)
  
  comm["Almere1",]=pool.A
  comm[c(grep("Almere[2-5]",rownames(comm)),grep("LB",rownames(comm))),]=com.rand(sampn=11,spnumi = spnumi,Nk=Nk,fits=fit.Ls)
  comm[grep("H",rownames(comm)),]=com.rand(sampn=12,spnumi = spnumi,Nk=Nk,fits=fit.Hs)
  
  rowSums(comm)
  dim(comm)
  comm=comm[,colSums(comm)>0]
  dim(comm)
  comm
}

drift.sloan<-function(op,Nk,sampname,spnumi=100,
                      m1=0.01,m2=0.01,FH=0.95,FL=0.05,sig2.W=4)
{
  # 4.2 # Sloan's dispersal model: moderate dispersal
  fit.WH=exp(-((op-FH)^2)/(2*sig2.W))
  fit.WL=exp(-((op-FL)^2)/(2*sig2.W))
  pi1.WL=fit.WL/sum(fit.WL)
  pi1.WH=fit.WH/sum(fit.WH)
  
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
  
  pool.A=disp.sloan(n=1,Nk=Nk,m1=m1,pi1=pi1.WL,spnumi = spnumi)[1,]
  
  pi2=pool.A/sum(pool.A)
  
  comm=matrix(0,nr=length(sampname),nc=length(op))
  rownames(comm)=sampname;colnames(comm)=names(op)
  dim(comm)
  
  comm[grep("L",rownames(comm)),]=disp.sloan(n=12,Nk=Nk,m1=m1,m2=m2,pi1 = pi1.WL,pi2 = pi2,spnumi = spnumi)
  comm[grep("H",rownames(comm)),]=disp.sloan(n=12,Nk=Nk,m1=m1,m2=m2,pi1 = pi1.WH,pi2 = pi2,spnumi = spnumi)
  
  rowSums(comm)
  dim(comm)
  comm=comm[,colSums(comm)>0]
  dim(comm)
  comm
}
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
  
  comm[grep("Almere",rownames(comm)),]=disp.sloan(n=5,Nk=Nk,m1=m1,m2=m2,pi1 = pi1,pi2 = pi2,spnumi = spnumi)
  comm[grep("Bath",rownames(comm)),]=disp.sloan(n=5,Nk=Nk,m1=m1,m2=m2,pi1 = pi1,pi2 = pi2,spnumi = spnumi)
  rowSums(comm)
  dim(comm)
  comm=comm[,colSums(comm)>0]
  dim(comm)
  comm
}


simmu<-function(prefix,op,bin.size,...)
{
  save.wd=paste0(save.wd,"/",prefix)
  if(!dir.exists(save.wd)){dir.create(save.wd)}
  
  ######################
  if(file.exists(paste0("MZSM.meta.J",J/1000/1000,"M.theta",theta,".rda")))
  {
    MZSM.meta=lazyopen(paste0("MZSM.meta.J",J/1000/1000,"M.theta",theta,".rda"))
  }else{
    MZSM.meta=sads::rmzsm(n=length(op),J=Nk*1000,theta=theta)
    save(MZSM.meta,file=paste0("MZSM.meta.J",J/1000/1000,"M.theta",theta,".rda"))
  }
}


save.sim<-function(prefixi,comm,pd,tree,save.wd,ABP=NULL)
{
  library(iCAMP)
  if(!is.null(ABP))
  {
    spc=iCAMP::match.name(cn.list=list(comm=comm),rn.list = list(ABP=ABP),
                          both.list = list(pd=pd),tree.list=list(tree=tree),
                          silent = TRUE)
    tree.use=spc$tree
    pd.use=spc$pd
    comm.use=spc$comm
    ABP.use=spc$ABP
    sim.data=list(comm=comm.use,pd=pd.use,tree=tree.use,ABP=ABP.use)
    filen=paste0(save.wd,"/",prefixi,".comm.tree.pd.ABP.rda")
  }else{
    spc=iCAMP::match.name(cn.list=list(comm=comm),both.list = list(pd=pd),tree.list=list(tree=tree),silent = TRUE)
    tree.use=spc$tree
    pd.use=spc$pd
    comm.use=spc$comm
    sim.data=list(comm=comm.use,pd=pd.use,tree=tree.use)
    filen=paste0(save.wd,"/",prefixi,".comm.tree.pd.rda")
  }
  save(sim.data,file = filen)
  filen
}

ABPcount<-function(comm,sel.spname=NULL,disp.spname=NULL,drift.spname=NULL,sp.bin=NULL,binid.col=3)
{
  ABP=data.frame(matrix(NA,nrow = ncol(comm),ncol = 2),stringsAsFactors = FALSE)
  rownames(ABP)=colnames(comm)
  colnames(ABP)=c("ABin","AProcess")
  ABP[match(sel.spname,colnames(comm)),2]="Selection"
  ABP[match(disp.spname,colnames(comm)),2]="Dispersal"
  ABP[match(drift.spname,colnames(comm)),2]="Drift"
  if(sum(is.na(ABP[,2]))>0){warning("Wrong! ABP has unexpected NA.")}
  if(!is.null(sp.bin))
  {
    ABP[match(sel.spname,colnames(comm)),1]=sp.bin[match(sel.spname,rownames(sp.bin)),binid.col]
    ABP[match(disp.spname,colnames(comm)),1]=sp.bin[match(disp.spname,rownames(sp.bin)),binid.col]
    ABP[match(drift.spname,colnames(comm)),1]=sp.bin[match(drift.spname,rownames(sp.bin)),binid.col]
    ABP[,1]=paste0("A",ABP[,1])
    if(sum(is.na(ABP))>0){warning("Wrong! ABP has unexpected NA.")}
  }
  ABP
}