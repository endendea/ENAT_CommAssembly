## Script from iCAMP publication 1 journal Ning et al., 2020

rm(list=ls())


###########################
# 0 # file/folder paths and package loading
# you may need to change them to the paths on your computer before testing the code.

library(here)

wd = here("data", "processed", "SimulatedData", "SpeciesPool")
code.wd= here("scripts", "iCAMP", "simulation")# the folder save simulation functions

source(paste0(code.wd,"/tools.r"))
library(ape)
library(phytools)
library(vegan)
library(iCAMP)
library(glue)
library(ape)


###########################################################
# 1 # 
# 1.1 # tree
tree.file = glue(wd, "/YXIN_lowN_2094ASV.nwk") # phylogenetic tree from YXIN lowN dataset (2094ASVs)
tree <- read.tree(tree.file)
#tree$tip.label=paste0("OTU",tree$tip.label) # set species IDs as OTU1 ... OTU1139

# Check if the tree is binary (no polytomies)
if (!is.binary(tree)) {
  tree <- ape::multi2di(tree)  # randomly resolve polytomies
}

# Check for zero branch lengths
if (any(tree$edge.length == 0)) {
  # Add a small jitter to zero-length branches
  tree$edge.length[tree$edge.length == 0] <- 1e-8
}

#Verify no zero-length branches remain
any(tree$edge.length == 0)

drt=iCAMP::tree.droot(tree,nworker = 4)
(drt.max=max(drt$distRoot)) # check max distance to root
tree$edge.length=tree$edge.length/drt.max # correct the distance to root, to range within 1.
tree
pd=iCAMP::pdist.big(tree = tree, wd = wd, output = TRUE, nworker = 4)
tree.save=list(tree=tree,pd=pd)
save(tree.save,file=paste0(wd,"/YXIN.lowN.tree.pd.rda"))
  
# 1.2 # OPEN, optimum environment, the key trait
# 1.2.1 # method 1: Stegen 2015, low phylogenetic signal across tree
opens=list()
#enop.old=lazyopen("All_pops_sim_10002.csv") # data file from Stegen et al 2015
#enopv.old=enop.old[,2]
#names(enopv.old)=paste0("OTU",rownames(enop.old))
#hist(enop.old[,2],breaks=100)
#(Kvalue.old=phytools::phylosig(tree, x=enopv.old, method="K", test=FALSE, nsim=1000))
#opens$JS=list(openv=enopv.old,K=Kvalue.old)

# 1.2.2 # method 2: Brownian, medium phylogenetic signal across tree
source(paste0(code.wd,"/BM.ACDC.r"))
BMt=BM.ACDC(tree = tree,mean.trait = 0.5,sig2 = 0.25^2,
            bounds = c(0,1),g=NULL,rand.num = 1,
            nworker = 4,code.wd = code.wd)

enop=BMt$trait
enopv=enop[,1];names(enopv)=rownames(enop)
treen=BMt$tree.new
hist(enop[,1],breaks = 50)
(Kvalue=phytools::phylosig(tree, x=enopv, method="K", test=FALSE, nsim=1000))
opens$BM=list(openv=enopv,K=Kvalue)

# 1.2.3 # method 3: ACDC model, high phylogenetic signal across tree
source(paste0(code.wd,"/BM.ACDC.r"))
g=2000
BMt=BM.ACDC(tree = tree,mean.trait = 0.5,sig2 = 0.25^2,
            bounds = c(0,1),g=g,rand.num = 1,
            nworker = 4,code.wd = code.wd)

enop=BMt$trait
enopv=enop[,1];names(enopv)=rownames(enop)
treen=BMt$tree.new
hist(enop[,1],breaks = 100)
(Kvalue=phytools::phylosig(tree, x=enopv, method="K", test=FALSE, nsim=1000))
opens$ACDC=list(openv=enopv,K=Kvalue)

save(opens,file=paste0(wd,"/YXIN.lowN.opens.rda"))

# End #
