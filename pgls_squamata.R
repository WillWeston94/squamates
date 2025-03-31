install.packages('nlme')
install.packages('rms')
install.packages('phytools')
install.packages('geiger')
install.packages('caper')
install.packages('tidyverse')
install.packages('cowplot')
install.packages('ggrepel')
install.packages('ggsci')
install.packages('patchwork')
install.packages('poolr')
install.packages('rr2')

install.packages('cpp11')
install.packages('iterators')
install.packages('fastmatch')

library(nlme)
#library(rms)
library(phytools)
library(geiger)
library(caper)
library(tidyverse)
library(cowplot)
library(ggrepel)
library(ggsci)
library(patchwork)
library(poolr)
library(rr2)
library(dplyr)

Data <- read.csv("GitHub/squamates/min20-2022.05.16.csv")
View(Data)


#Cut the data
squamata_data <- Data %>%
  filter(grepl("Squamata", Orders, ignore.case = TRUE))

write.csv(squamata_data, "squamata_data.csv", row.names = FALSE)


#pgls sey base (just run all of this)
modPgls.SEy = function (model, data, corClass = corBrownian, tree, se = NULL, 
                        method = c("REML", "ML"), interval = c(0, 1000), corClassValue=1, sig2e=NULL, ...) 
{
  Call <- match.call()
  corfunc <- corClass
  spp <- rownames(data)
  data <- cbind(data, spp)
  if (is.null(se)) 
    se <- setNames(rep(0, Ntip(tree)), tree$tip.label)[spp]
  else se <- se[spp]
  
  lk <- function(sig2e, data, tree, model, ve, corfunc, spp) {
    tree$edge.length <- tree$edge.length * sig2e
    ii <- sapply(1:Ntip(tree), function(x, e) which(e == 
                                                      x), e = tree$edge[, 2])
    tree$edge.length[ii] <- tree$edge.length[ii] + ve[tree$tip.label]
    vf <- diag(vcv(tree))[spp]
    w <- varFixed(~vf)
    COR <- corfunc(corClassValue, tree, form = ~spp, ...)
    fit <- gls(model, data = cbind(data, vf), correlation = COR, 
               method = method, weights = w)
    -logLik(fit)
  }
  
  if (is.null(sig2e)) {
    fit <- optimize(lk, interval = interval, data = data, tree = tree, 
                    model = model, ve = se^2, corfunc = corfunc, spp = spp)
    sig2e=fit$minimum
  }
  
  tree$edge.length <- tree$edge.length * sig2e
  ii <- sapply(1:Ntip(tree), function(x, e) which(e == x), 
               e = tree$edge[, 2])
  tree$edge.length[ii] <- tree$edge.length[ii] + se[tree$tip.label]^2
  vf <- diag(vcv(tree))[spp]
  w <- varFixed(~vf)
  obj <- gls(model, data = cbind(data, vf), correlation = corfunc(corClassValue, 
                                                                  tree, form = ~spp, ...), weights = w, method = method)
  obj$call <- Call
  obj$sig2e <- sig2e
  obj
}


#Internal function
pglsSEyPagelToOptimizeLambda=function(lambda,model,data,tree,...) {
  -logLik(modPgls.SEy(model=model,data=data,tree=tree,corClassValue=lambda,corClass=corPagel,fixed=T,...)) #Returns -logLikelihood of the pgls.SEy model with lambda fixed to the value of the lambda argument. sig2e will be optimized within modPgls.SEy unless given as an argument here
}

#Function intended for users
pglsSEyPagel=function(model, data, tree, lambdaInterval=c(0,1),...){
  optimizedModel=optimize(pglsSEyPagelToOptimizeLambda,lambdaInterval,model=model,data=data,tree=tree,...) #Optimizes lambda in the lambdaInterval using the pglsSEyPagelToOptimizeLambda function
  return(modPgls.SEy(model=model,data=data,tree=tree,corClass=corPagel,fixed=T,corClassValue=optimizedModel$minimum,...)) #Returns the final model fit
}


#read data
Squamata <- read.csv("squamata_data.csv")
View(Squamata)

#adult weight models
#adult weight neo

cutData <- Squamata[,c(5,9,10,11,13,38,42),drop=FALSE] 
cutData[cutData$adult_weight == -1, ] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("GitHub/squamates/min20Fixed516.nwk")

cutData$Species <- gsub(" ", "_", cutData$Species) 
cutData$common_name<-gsub("_", "", cutData$common_name)
includedSpecies<-cutData$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
cutData <- cutData[!(cutData$Keep==FALSE),]
rownames(cutData)<-cutData$Species
SE<-setNames(cutData$SE_simple,cutData$Species)[rownames(cutData)]
View(cutData)


#pgls model
adult.weight.neo<-pglsSEyPagel(NeoplasiaPrevalence~log10(adult_weight.g.),data=cutData,tree=pruned.tree,se=SE,method = "ML")

summary(adult.weight.neo) 




#adult weight mal
cutData <- Squamata[,c(5,9,10,11,17,38,42),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("GitHub/squamates/min20Fixed516.nwk")

cutData$Species <- gsub(" ", "_", cutData$Species) 
cutData$common_name<-gsub("_", "", cutData$common_name)
includedSpecies <- cutData$Species

pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
cutData <- cutData[!(cutData$Keep==FALSE),]
rownames(cutData)<-cutData$Species
SE<-setNames(cutData$SE_simple,cutData$Species)[rownames(cutData)]
View(cutData)

#pgls model
adult.weight.mal<-pglsSEyPagel(MalignancyPrevalence~log10(adult_weight.g.),data=cutData,
                               tree=pruned.tree,se=SE,method="ML")
summary(adult.weight.mal)
