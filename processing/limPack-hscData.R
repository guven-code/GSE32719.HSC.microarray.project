source("http://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("gcrma")
biocLite("limma")
biocLite("R.utils")
biocLite("affycoretools")
biocLite("Biobase")
biocLite("GEOquery")
library(affy)
library(gcrma)
library(limma)
library(R.utils)
library(affycoretools)
library(Biobase)
library(GEOquery)
geoNumber = "GSE32719"
gse<- getGEO(geoNumber, destdir=".",GSEMatrix=TRUE)
annotation(gse[[1]])
gpl=getGEO(annotation(gse[[1]]), destdir=".")
id2Name=Table(gpl)[,c(1,10)] #Here you can check which column contains the information you want, here I used 10
id2NameMatrix=matrix(NA, nrow(id2Name), 2)
for(i in 1:nrow(id2Name)){
  id2NameMatrix[i,1] = as.character(id2Name[i,1])
  split = unlist(strsplit(as.character(id2Name[i,2]), split="//")[[1]])
  if(length(split) == 1){
    id2NameMatrix[i,2] = as.character(split[1]);
  }
  else{
    id2NameMatrix[i,2] = as.character(split[2]);
  }
}
colnames(id2NameMatrix) = c("ID", "Name")
id2NameMatrix = as.data.frame(id2NameMatrix)
id2Name = id2NameMatrix$Name
names(id2Name) = id2NameMatrix$ID
pData(gse[[1]])
Group1 =c(rep("Control",17 ), rep("Case", 10)) #Here you have to define the grouping, you can read the pData of gse to get the information e.g. pData(gse[[1]])
Group2=Group1 =c(rep("Control",20 ), rep("Case", 7))
design1<- model.matrix(~Group1)
design2<-model.matrix(~Group2)
fit1 <-lmFit(exprs(gse[[1]]),design1)
fit2 <-lmFit(exprs(gse[[1]]),design2)

fit1<-eBayes(fit1)
fit2<-eBayes(fit2)
fitgenesSymbol = id2Name
top.all1 <- topTable(fit1,n=nrow(exprs(gse[[1]])),adjust="BH",coef=2) #This will give you all the pvalue when comparing the samples (for 2 condition)
top.all2 <- topTable(fit2,n=nrow(exprs(gse[[1]])),adjust="BH",coef=2) #This will give you all the pvalue when comparing the samples (for 2 condition)

summary(top.all1)


summary(top.all2)


dim(fit1)
colnames(fit1)
rownames(fit1)[1:10]
names(fit1)

# Fold-change thresholding
fitFold <- treat(fit1,lfc=0.1)
topTreat(fit1,coef=2)

# Volcano plot
volcanoplot(fit1,coef=2,highlight=2)

#volcanoplot(fit_contrast)
# Generate a list of top 100 differentially expressed genes
top_genes <- topTable(fit_contrast, number = 100, adjust = "BH")
# Summary of results (number of differentially expressed genes)
result <- decideTests(fit1)
summary(result)

# Mean-difference plot
plotMD(fit1,column=2)

# Q-Q plot of moderated t-statistics
qqt(fit1$t[,2],df=fit1$df.residual+fit1$df.prior)
abline(0,1)

# Various ways of writing results to file
## Not run: write.fit(fit,file="exampleresults.txt")
## Not run: write.table(fit,file="exampleresults2.txt")

# Fit with correlated arrays
# Suppose each pair of arrays is a block
#block <- c(1,1,2,2,3,3)
#dupcor <- duplicateCorrelation(z,design1,block=block)
#dupcor$consensus.correlation
#fitArrays <- lmFit(z,design1,block=block,correlation=dupcor$consensus)

# Fit with duplicate probes
# Suppose two side-by-side duplicates of each gene
#rownames(z) <- paste("Gene",rep(1:27339,each=2))
#dupcor <- duplicateCorrelation(z,design1,ndups=2)
#dupcor$consensus.correlation
#fitProbes <- lmFit(z,design1,ndups=2,correlation=dupcor$consensus)
#dim(fitProbes)
#fitProbes <- eBayes(fitProbes)
#topTable(fitProbes,coef=2)





