# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Thu Aug 29 10:50:11 EDT 2019


Server: www.ncbi.nlm.nih.gov
Query: acc=GSE32719&platform=GPL570&type=txt&groups=&colors=&selection=XXXXXXXXXXXXXXXXXXXXXXXXXXX&padj=fdr&logtransform=auto&columns=ID&columns=adj.P.Val&columns=P.Value&columns=F&columns=Gene+symbol&columns=Gene+title&num=250&annot=ncbi

# Unable to generate script analyzing differential expression.
#      Invalid input: at least two groups of samples should be selected.

################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)
library(doMC)

# load series and platform data from GEO

gset <- getGEO("GSE32719", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# set parameters and draw the plot

dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE32719", '/', annotation(gset), " selected samples", sep ='')
boxplot(exprs(gset), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)

##Exression seviyelerini kontrol eder ve normalize eder.
ex= exprs(gset) ##gset i?inden expres verisini ?ekip de?i?ken ile tan?mlad?k.
class(ex)
dim(ex)
colnames(ex)
head(ex)
minimalSet <- ExpressionSet(assayData=ex)
#pDataFile <- file.path("/Users/user/AMH Dropbox/current-EG/hsc-jbcb/processing", "pData.txt")
#pData <- read.table(pDataFile,
 #                   row.names=1, header=TRUE, sep="\t")

apply( ex, 2, median) 
means = apply( ex, 2, mean) ##Data i?indeki ki?ilerin tan?mlay?c? istatistiklerini kutu grafikte ?izdirir. 
means
#pdf(paste("C:/Users/lenovo/Dropbox/analiz" , "boxplot.pdf", sep=''),width=5,height=5)
#if (FLAG1==TRUE){
  boxplot( ex )
#}
#dev.off()


#Bu normalizasyon datan?n toplam yo?unlu?una g?re ger?ekle?ir.
#Scale ile bir ?l?ekleme sa?lanmaktad?r.
scale = max(means)
for(i in 1:length(ex[1,])){
  ex[,i] = ex[,i] * scale/ means[i] ##for d?ng?s? i?in scale ile datadaki max. ortalama de?eri ile 
}                                   ###tan?mlad?k.
apply( ex, 2, mean) / scale #Data i?indeki t?m ort. de?erlerini / scale(max.ort)
apply( ex, 2, median) /scale #Data i?indeki ortalama de?erleri scale ile b?lerek bir normalizasyon yap?l?r.
boxplot( ex )


str(gset) ##data y? analizden ?nce inceler
gset@phenoData@data ##phenodata n?n i?indeki datay? ?eker.(title, last update date, type, prtocol vb. bilgsini verir)
experimental_design = gset@phenoData@data
experimental_design[1:3,] #Yukar?da belirtti?imiz bilgileri datadan sadece ilk 3 ki?iye ili?kin verir.
# gset@phenoData@varMetadata
# experimental_design[, "source_name_ch1"][1:10]
experimental_design[1:3,  c("title", "source_name_ch1")] ##Matrixdeki s?tunlara title ekler.
unique( experimental_design$source_name_ch1 ) 


gpl <- annotation(gset) ##Biobase'de bulunan annotation(gset) gpl de?i?kenine atay?p getGEO'dan verileri ald?k.
platf <- getGEO(gpl, AnnotGPL=TRUE)
ncbifd <- data.frame(attr(dataTable(platf), "table")) ##Yukar?da tan?mlad???m?z platf i?inden table[54675:21] matrixi
###burada ?ekip ncbifd de?i?kenine atad?k.Bu tablo bize Gene ID,
####Gene title,Gene symbol gibi ba?l?klarda ??kt?lar sunar.

##Ya?a ba?l? regresyon analizine ba?l?yoruz.
sample.names = unique(experimental_design$source_name_ch1) ###?rneklem ?zerindeki bireleri 3 grupta inceledik.(Y,M,O)
mylevels = c(1,2,3)
names(mylevels) = sample.names ##Yukar?da tan?msal kodlar? e?itleyerek 3 grubu e?le?tirdik.
gset.levels = mylevels[ experimental_design$source_name_ch1 ]
my.pvalues = numeric(length=length(ex[,1])) ##P de?erlerinin ex datasetinden ?ekilip numeric de?erler i?erdi?ini
###ifade ettik.



#limma package
fit <- lmFit(ex[1,], )

# x2 = foreach(i=1:(3*2), .combine='cbind') %dopar% sqrt(i)

registerDoMC(cores=4) ##doMC package ?a??r?yoruz.
#my.pvalues = foreach( i = 1:100, .combine='rbind') %dopar% {
my.pvalues = foreach( i = 1:length(ex[,1]), .combine='rbind') %dopar% {
  m = lm( ex[i,] ~ gset.levels ) ##lm komutu (lineer model) regresyon analizinde kullan?l?r. i'de?eri p de?erlerini 
  ###s?ras?yla d?nd?recektir. ex gset(exprs) datalar?n? ifade ederek bizim ba??ml? 
  ####de?i?kenimizdir. Yukar?da tan?mlad???m?z gset.levels g?re reg. analizini uygular.
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


#source("https://bioconductor.org/biocLite.R")
#biocLite("hgu133plus2.db")
#library("hgu133plus2.db")
#help(package="hgu133plus2.db")
 #ls("package:hgu133plus2.db")
 #x <- hgu133plus2GENENAME
 #x <- hgu133plus2SYMBOL
 #mapped_probes <- mappedkeys(x)
 #xx <- as.list(x[mapped_probes])
 #xx[1:5]

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
pdf(paste("grafikler/", "colboxplot.pdf",sep=''),width=5,height=5)
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)

par(cex.lab=0.8)#or y-axis
par(cex.axis=0.8)#for x-axis
legend("topleft", labels, fill=palette(), bty="n")
dev.off()

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
exprs(eset["1007_s_at","GSM812988"]) ##?rn:gen de?eri
##makalede ilgilenilen genleri kullanarak analizler devam ettirilebilir.
#pdf(paste("grafikler/" , "figure-samplegeneanalysis.pdf", sep=''),width=5,height=5)
barplot(exprs(eset["215446_s_at"])) ##?rn:t?m?r olu?muna inhibit?r etki ifade eden gen de?eri
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
