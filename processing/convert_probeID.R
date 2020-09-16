source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
biocLite("limma")
biocLite("Biobase")
biocLite("affy")

library(GEOquery)
library(biomaRt)


gset <- getGEO("GSE32719", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

class(gset)

length(gset)

slotNames(gset)

dim(gset)

str(gset)

pData(phenoData(gset))

dim(pData(phenoData(gset)))

colnames(pData(phenoData(gset)))

pData(phenoData(gset))[ , c(12,13)]

fvarLabels(gset) <- make.names(fvarLabels(gset))

head(exprs(gset))

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

par(mfrow=c(1,2))
# Before log2 - since it's already been log2 trasnform 
# we use 2^ to compute back the original value:
hist(2^exprs(gset), breaks=25)
#After log2
hist(exprs(gset), breaks=25)


indicesLookup <- match(rownames(gset), data_annotLookup$affy_hg_u133_plus_2)

rownames(gset) <- paste(data_annotLookup[indicesLookup, "external_gene_name"], 
                        c(1:length(indicesLookup)), sep="_")


data=exprs(gset)

