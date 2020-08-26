

#GSE32719
#http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32719


rm(list=ls())
#setwd("C:/Users/lenovo/Dropbox/analiz")
library(Biobase)
library(GEOquery)
library(limma)
library(foreach)

#install.packages("doMC", repos="http://R-Forge.R-project.org")
library(doMC)

#load platform data from bioconductor
#http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32719
#GPL570	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
#https://bioconductor.org/packages/release/data/annotation/html/hgu133plus2.db.html

FLAG1=TRUE
FLAG2=FALSE
FLAG3=TRUE

#Load series and platform data from GEO
gset <- getGEO("GSE32719", GSEMatrix =TRUE)
class(gset)
str(gset)


 gset <- gset$GSE32719_series_matrix.txt.gz
# make proper column names to match toptable 
# fvarLabels(gset) <- make.names(fvarLabels(gset))


##Exression seviyelerini kontrol eder ve normalize eder.
ex= exprs(gset) ##gset içinden expres verisini çekip deðiþken ile tanýmladýk.
apply( ex, 2, median) 
means = apply( ex, 2, mean) ##Data içindeki kiþilerin tanýmlayýcý istatistiklerini kutu grafikte çizdirir. 
means
#pdf(paste("C:/Users/lenovo/Dropbox/analiz" , "boxplot.pdf", sep=''),width=5,height=5)
if (FLAG1==TRUE){
boxplot( ex )
}
#dev.off()


#Bu normalizasyon datanýn toplam yoðunluðuna göre gerçekleþir.
#Scale ile bir ölçekleme saðlanmaktadýr.
scale = max(means)
for(i in 1:length(ex[1,])){
  ex[,i] = ex[,i] * scale/ means[i] ##for döngüsü için scale ile datadaki max. ortalama deðeri ile 
}                                   ###tanýmladýk.
apply( ex, 2, mean) / scale #Data içindeki tüm ort. deðerlerini / scale(max.ort)
apply( ex, 2, median) /scale #Data içindeki ortalama deðerleri scale ile bölerek bir normalizasyon yapýlýr.
boxplot( ex )


str(gset) ##data yý analizden önce inceler
gset@phenoData@data ##phenodata nýn içindeki datayý çeker.(title, last update date, type, prtocol vb. bilgsini verir)
experimental_design = gset@phenoData@data
experimental_design[1:3,] #Yukarýda belirttiðimiz bilgileri datadan sadece ilk 3 kiþiye iliþkin verir.
# gset@phenoData@varMetadata
# experimental_design[, "source_name_ch1"][1:10]
experimental_design[1:3,  c("title", "source_name_ch1")] ##Matrixdeki sütunlara title ekler.
unique( experimental_design$source_name_ch1 ) 


gpl <- annotation(gset) ##Biobase'de bulunan annotation(gset) gpl deðiþkenine atayýp getGEO'dan verileri aldýk.
platf <- getGEO(gpl, AnnotGPL=TRUE)
ncbifd <- data.frame(attr(dataTable(platf), "table")) ##Yukarýda tanýmladýðýmýz platf içinden table[54675:21] matrixi
                                                      ###burada çekip ncbifd deðiþkenine atadýk.Bu tablo bize Gene ID,
                                                      ####Gene title,Gene symbol gibi baþlýklarda çýktýlar sunar.

##Yaþa baðlý regresyon analizine baþlýyoruz.
sample.names = unique(experimental_design$source_name_ch1) ###Örneklem üzerindeki bireleri 3 grupta inceledik.(Y,M,O)
mylevels = c(1,2,3)
names(mylevels) = sample.names ##Yukarýda tanýmsal kodlarý eþitleyerek 3 grubu eþleþtirdik.
gset.levels = mylevels[ experimental_design$source_name_ch1 ]
my.pvalues = numeric(length=length(ex[,1])) ##P deðerlerinin ex datasetinden çekilip numeric deðerler içerdiðini
                                            ###ifade ettik.

# x2 = foreach(i=1:(3*2), .combine='cbind') %dopar% sqrt(i)

registerDoMC(cores=4) ##doMC package çaðýrýyoruz.
#my.pvalues = foreach( i = 1:100, .combine='rbind') %dopar% {
my.pvalues = foreach( i = 1:length(ex[,1]), .combine='rbind') %dopar% {
m = lm( ex[i,] ~ gset.levels ) ##lm komutu (lineer model) regresyon analizinde kullanýlýr. i'deðeri p deðerlerini 
                               ###sýrasýyla döndürecektir. ex gset(exprs) datalarýný ifade ederek bizim baðýmlý 
                               ####deðiþkenimizdir. Yukarýda tanýmladýðýmýz gset.levels göre reg. analizini uygular.
sm = summary(m)
  pf(sm$fstatistic[1], sm$fstatistic[2], sm$fstatistic[3], lower.tail = FALSE)
}
row.names(my.pvalues) = row.names(gset@assayData$exprs)

#pdf(paste("grafikler/" ,"figure-pvalues.pdf",sep=''),width=5,height=5)
hist(my.pvalues)
#dev.off()
summary(my.pvalues)

my.pvalues.BH = p.adjust(my.pvalues, "BH")
names(my.pvalues.BH) = row.names(gset@assayData$exprs)
#pdf(paste("grafikler/", "figure-pvaluesBH.pdf",sep=''),width=5,height=5)
hist(my.pvalues.BH)
#dev.off()
summary(my.pvalues.BH)

sig= my.pvalues.BH[my.pvalues.BH < 0.05]
sig= data.frame(sig)
sig$ID = row.names(sig)
#ncbifd[ match(names(sig), as.character( ncbifd$ID) ),  ]
sig2 = merge(sig, ncbifd, by="ID")

#GO analysis
sig.genes = unique(sig2$Gene.symbol)
#sig.geneIDs = unique(sig2$Gene.ID)
write.table(sig.genes, "__sig.genes.tsv", sep="\t", quote=F, row.names=F, col.names=F  )

background.genes = unique( ncbifd$Gene.symbol)
write.table(background.genes, "__background.genes.tsv", sep="\t", quote=F, row.names=F, col.names=F  )

#library(topGO)
#do GOrilla analysis


# source("https://bioconductor.org/biocLite.R")
# biocLite("hgu133plus2.db")
# library("hgu133plus2.db")
# help(package="hgu133plus2.db")
# ls("package:hgu133plus2.db")
# x <- hgu133plus2GENENAME
# x <- hgu133plus2SYMBOL
# mapped_probes <- mappedkeys(x)
# xx <- as.list(x[mapped_probes])
# xx[1:5]

View(ex)
View(experimental_design)
plot(experimental_design$source_name_ch1)
plot(experimental_design$characteristics_ch1.2)

sml <- c("Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y",
         "M","M","M","M","M","O","O","O","O","O","O","O","O")

fl <- as.factor(sml)
labels <- c("MIDDLE","OLD","YOUNG")

palette(c("#dfeaf4","#f4dfdf","#f2cb98"))
title <- paste ("GSE32719", '/', gpl, " selected samples", sep ='')
#pdf(paste("grafikler/", "colboxplot.pdf",sep=''),width=5,height=5)
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")
#dev.off()

hist(ex, col = "gray", main="GSE32719 - Histogram")
log = log2(ex)
hist(log, col = "gray", main="GSE32719 (log2) - Histogram")

hc = hclust(as.dist(1-cor(ex)))
plot(hc, main="GSE32719 - Hierarchical Clustering")


gds3942 <- getGEO('GDS3942', destdir=".")
#gds3942 <- getGEO(filename='GDS3942.soft.gz')

Meta(gds3942)$channel_count
Meta(gds3942)$description

colnames(Table(gds3942))
Table(gds3942)[1:10,1:7]

eset <- GDS2eSet(gds3942, do.log2=TRUE)

eset@featureData@data$ID[1:10]
eset@experimentData@other$sample_id

eset["1007_s_at","GSM812988"]
exprs(eset["1007_s_at","GSM812988"]) ##örn:gen deðeri
##makalede ilgilenilen genleri kullanarak analizler devam ettirilebilir.
#pdf(paste("grafikler/" , "figure-samplegeneanalysis.pdf", sep=''),width=5,height=5)
barplot(exprs(eset["215446_s_at"])) ##Örn:tümör oluþmuna inhibitör etki ifade eden gen deðeri
#dev.off()

Meta(gds3942)$platform

gpl570 <- getGEO('GPL570', destdir=".")

Meta(gpl570)$title
colnames(Table(gpl570))
Table(gpl570)[1:10,1:5] #or
Table(gpl570)[1:10,c("ID","Gene Title","Gene Symbol","ENTREZ_GENE_ID")]
  
attr(dataTable(platf), "table")[100,]
attr(dataTable(platf), "table")[1:100,c("ID", "Gene symbol")]
IDs <- attr(dataTable(platf), "table")[,c("ID", "Gene symbol")]
IDs

rows <- rownames(ex)
IDs[which(IDs[,1] == rows[1000]), 2]

selected  <- my.pvalues.BH < 0.05
esetSel <- eset [selected, ]
heatmap(exprs(esetSel))
#pdf(paste("grafikler/" ,"figure-heatmap.pdf",sep=''),width=5,height=5)
heatmap(exprs(esetSel), col=topo.colors(100))
#dev.off()
