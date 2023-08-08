######## Run BayesPrism
##### @author: Cristiane Esteves
######### 09/2022

.libPaths("~/R/4.5")
library(BayesPrism)
library(data.table)
library(Seurat)

h1 <- Sys.time()

setwd("/Path/")


breast = readRDS("myPrism_breast.rds")

bp.res.initial <- run.prism(prism = breast, n.cores=70, update.gibbs=FALSE)
bp.res.update <- update.theta(bp = bp.res.initial)


save(bp.res.update, file="Prism_breast.RData")

rm(breast, bp.res.initial,bp.res.update)
########

coad = readRDS("myPrism_coad.rds")

bp.res.initial <- run.prism(prism = coad, n.cores=70, update.gibbs=FALSE)
bp.res.update <- update.theta(bp = bp.res.initial)

save(bp.res.update, file="Prism_coad.RData")
#######

path0 = ""
path1=""

liver = readRDS(paste0(path0,"myPrism_liver.rds"))

bp.res.initial <- run.prism(prism = liver, n.cores=70, update.gibbs=FALSE)
bp.res.update <- update.theta(bp = bp.res.initial)

save(bp.res.update, file=paste0(path1,"Prism_liver.RData"))
#######

luad = readRDS("myPrism_luad.rds")

bp.res.initial <- run.prism(prism = luad, n.cores=70, update.gibbs=FALSE)
bp.res.update <- update.theta(bp = bp.res.initial)

save(bp.res.update, file = "Prism_luad.RData")
#######

lusc = readRDS("myPrism_lusc.rds")

bp.res.initial <- run.prism(prism = lusc, n.cores=70, update.gibbs=FALSE)
bp.res.update <- update.theta(bp = bp.res.initial)

save(bp.res.update, file="Prism_lusc.RData")
#######

melanoma = readRDS("myPrism_melanoma.rds")

bp.res.initial <- run.prism(prism = melanoma, n.cores=70, update.gibbs=FALSE)
bp.res.update <- update.theta(bp = bp.res.initial)

save(bp.res.update, file="Prism_melanoma.RData")
#######
ovary = readRDS("myPrism_ovary.rds")

bp.res.initial <- run.prism(prism = ovary, n.cores=70, update.gibbs=FALSE)
bp.res.update <- update.theta(bp = bp.res.initial)

save(bp.res.update, file="Prism_ovary.RData")
#######
#falta esse p baixo
read = readRDS("myPrism_read.rds")

bp.res.initial <- run.prism(prism = read, n.cores=70, update.gibbs=FALSE)
bp.res.update <- update.theta(bp = bp.res.initial)

save(bp.res.update, file="Prism_read.RData")

#######

uveal = readRDS("myPrism_uveal.rds")

bp.res.initial <- run.prism(prism = uveal, n.cores=70, update.gibbs=FALSE)
bp.res.update <- update.theta(bp = bp.res.initial)

save(bp.res.update, file="Prism_uveal.RData")
#######

h2 <- Sys.time()
print(h2 - h1)
