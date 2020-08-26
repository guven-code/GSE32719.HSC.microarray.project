#############################################################################
# Gene expression analysis of expression data from human bone
# marrow hematopoietic stem cells
# Organism	Homo sapiens
# Dataset: GSE32719 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32719)
# Pang WW, Price EA, Sahoo D, Beerman I et al. Human bone marrow hematopoietic stem cells 
# are increased in frequency and myeloid-biased with age.
# Proc Natl Acad Sci U S A 2011 Dec 13;108(50):20012-7.
# PMID: 22123971 (http://www.ncbi.nlm.nih.gov/pubmed/22123971)
# R code: Emine Guven
##############################################################################


GSEdata=read.csv("GSE32719.expression.matrix.csv",stringsAsFactors = FALSE)

geneIDs=GSEdata$X
conditions=colnames(GSEdata[,-1])

# Check the loaded dataset
dim(GSEdata) # Dimension of the dataset
head(GSEdata) # First few rows
tail(GSEdata) # Last few rows

raw_data=data.matrix(GSEdata,rownames.force = NA)
rownames(raw_data)=geneIDs

#write.csv(raw_data,"raw_data.csv")

data <- raw_data[,-1]

###################
# Exploratory plots
###################

# Check the behavior of the data
pdf("plots/histogramRawData.pdf")
hist(data, col = "gray", main="GSE32719 - Histogram")
dev.off()

# Log2 transformation (why?)
data2 = log2(data)

# Check the behavior of the data after log-transformation
pdf("plots/Log2histogramRawData.pdf")
hist(data2, col = "gray", main="GSE32719 (log2) - Histogram")
dev.off()
# Boxplot
pdf("plots\boxPlotRawData.pdf")
boxplot(data2, col=c("darkgreen", "darkgreen", "darkgreen","darkgreen","darkgreen",
                     "darkgreen", "darkgreen", "darkgreen","darkgreen","darkgreen",
                     "darkgreen", "darkgreen","darkgreen","darkred",
                     "darkred", "darkred", "darkred","darkred",
                     "blue", "blue", "blue","blue", "blue", "blue","blue", "blue"),
        main="GSE32719 - boxplots", las=2,cex.axis=0.6)
dev.off()
# Hierarchical clustering of the "samples" based on
# the correlation coefficients of the expression values
hc = hclust(as.dist(1-cor(data2)))
pdf("plots/HclustRawData.pdf")
plot(hc, main="GSE32719 - Hierarchical Clustering")
dev.off()
#######################################
# Differential expression (DE) analysis
#######################################

# Separate each of the conditions into three smaller data frames
young = data2[,1:13]
middle = data2[,14:18]
old = data2[,19:26]

# Compute the means of the samples of each condition
young.mean = apply(young, 1, mean)
middle.mean = apply(middle, 1, mean)
old.mean = apply(old,1, mean)

head(young.mean)

head(middle.mean)

head(old.mean)

# Just get the maximum of all the means
limit = max(young.mean, middle.mean,old.mean)

# Scatter plots of young vs old and young vs middle aged
pdf("plots/comparisonConditions.pdf")
par(mfrow=c(2,3))
plot(young.mean ~ old.mean, xlab = "young", ylab = "old",
     main = "GSE32719 - Scatter", xlim = c(0, limit), ylim = c(0, limit))
# Diagonal line
abline(0, 1, col = "red")

######compare young and middle##############
plot(young.mean ~ middle.mean, xlab = "young", ylab = "middle",
     #main = "GSE32719 - Scatter", 
     xlim = c(0, limit), ylim = c(0, limit))
# Diagonal line
abline(0, 1, col = "blue")
##########compare middle and old aged patients
plot(middle.mean ~ old.mean, xlab = "middle", ylab = "old",
     #main = "GSE32719 - Scatter", 
     xlim = c(0, limit), ylim = c(0, limit))
# Diagonal line
abline(0, 1, col = "green")
# Compute fold-change (biological significance)
# Difference between the means of the conditions
#par(mfrow=c(1,3))
foldYO = old.mean - young.mean

# Histogram of the fold differences
hist(foldYO, col = "red")

## between M vs Y ####
foldYM = middle.mean - young.mean

# Histogram of the fold differences
hist(foldYM, col = "blue")

## between O vs M ####
foldMO = old.mean - middle.mean

# Histogram of the fold differences
hist(foldMO, col = "green")
dev.off()


# Compute statistical significance (using t-test)
pvalueYO = NULL # Empty list for the p-values
tstatYO = NULL # Empty list of the t test statistics
pvalueYM = NULL # Empty list for the p-values
tstatYM = NULL # Empty list of the t test statistics
pvalueMO = NULL # Empty list for the p-values
tstatMO = NULL # Empty list of the t test statistics

for(i in 1 : nrow(data)) { # For each gene : 
  Y = young[i,] # WT of gene number i
  M = middle[i,] # KO of gene number i
  O = old[i,]
  # Compute t-test between the bi- conditions
  tYO = t.test(Y, O)
  
  tYM = t.test(Y, M)
  
  tMO = t.test(M, O)
  
  # Put the current p-value in the pvalues list
  pvalueYO[i] = tYO$p.value
  # Put the current t-statistic in the tstats list
  tstatYO[i] = tYO$statistic
  
  # Put the current p-value in the pvalues list
  pvalueYM[i] = tYM$p.value
  # Put the current t-statistic in the tstats list
  tstatYM[i] = tYM$statistic
  
  # Put the current p-value in the pvalues list
  pvalueMO[i] = tMO$p.value
  # Put the current t-statistic in the tstats list
  tstatMO[i] = tMO$statistic
}

head(pvalueYO)
head(pvalueYM)
head(pvalueMO)

pdf("plots/pvalueByConditionHist.pdf")
par(mfrow=c(3,1))
# Histogram of p-values (-log10)
hist(-log10(pvalueYO), col = "red")
hist(-log10(pvalueYM), col = "blue")
hist(-log10(pvalueMO), col = "green")
dev.off()

# Volcano: put the biological significance (fold-change)
# and statistical significance (p-value) in one plot
pdf("plots/foldYO_volcano_I.pdf")
plot(foldYO, -log10(pvalueYO), main = "GSE32719 YO conditions - Volcano")
dev.off()

fold_cutoff_YO = 0.075
pvalue_cutoff_YO = 0.001
abline(v = fold_cutoff_YO, col = "blue", lwd = 3)
abline(v = -fold_cutoff_YO, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff_YO), col = "green", lwd = 3)

# Screen for the genes that satisfy the filtering criteria

# Fold-change filter for "biological" significance
filter_by_fold_YO = abs(foldYO) >= fold_cutoff_YO
dim(data2[filter_by_fold_YO, ])

# P-value filter for "statistical" significance
filter_by_pvalue_YO = pvalueYO <= pvalue_cutoff_YO
dim(data2[filter_by_pvalue_YO, ])

# Combined filter (both biological and statistical)
filter_combined_YO = filter_by_fold_YO & filter_by_pvalue_YO

filtered_YO = data2[filter_combined_YO,]
dim(filtered_YO)

head(filtered_YO)


# Let's generate the volcano plot again,
# highlighting the significantly differential expressed genes
pdf("plots/foldYO_volcanoYO)_II.pdf")
plot(foldYO, -log10(pvalueYO), main = "GSE32719 YO conditions - Volcano #2")
points (foldYO[filter_combined_YO], -log10(pvalueYO[filter_combined_YO]),
        pch = 16, col = "red")
dev.off()
# Highlighting up-regulated in red and down-regulated in blue
pdf("plots/foldYO_volcanoYO_III.pdf")
plot(foldYO, -log10(pvalueYO), main = "GSE32719 YO conditions - Volcano #3")
points (foldYO[filter_combined_YO & foldYO < 0],
        -log10(pvalueYO[filter_combined_YO & foldYO < 0]),
        pch = 16, col = "blue")
points (foldYO[filter_combined_YO & foldYO > 0],
        -log10(pvalueYO[filter_combined_YO & foldYO > 0]),
        pch = 16, col = "red")
legend(1, 95, legend=c("Up", "Down"),
       col=c("red", "blue"), lty=1:2, cex=0.8)

dev.off()

upRegulated_DEGs_YO = foldYO[filter_combined_YO & foldYO < 0]
length(upRegulated_DEGs_YO)

downRegulated_DEGs_YO =  foldYO[filter_combined_YO & foldYO > 0]
length(downRegulated_DEGs_YO)

# Cluster the rows (genes) & columns (samples) by correlation
rowvYO = as.dendrogram(hclust(as.dist(1-cor(t(filtered_YO)))))
colvYO = as.dendrogram(hclust(as.dist(1-cor(filtered_YO))))

# Generate a heatmap
pdf("plots/heatmap_filtered_YO.pdf")
heatmap(filtered_YO, Rowv=rowvYO, Colv=colvYO, cexCol=0.7)
dev.off()
# install.packages("gplots")		# Uncomment if not already installed
# install.packages("RColorBrewer")	# Uncomment if not already installed

library(gplots)

# Enhanced heatmap
pdf("plots/heatmap2_filtered_YO.pdf")
heatmap.2(filtered_YO, Rowv=rowvYO, Colv=colvYO, cexCol=0.7,
          col = rev(redblue(256)), scale = "row")
dev.off()
# Save the heatmap to a PDF file
#pdf ("GSE32719_YO_DE_Heatmap.pdf")
heatmap.2(filtered_YO, Rowv=rowvYO, Colv=colvYO, cexCol=0.7,
          col = rev(redblue(256)), scale = "row")
#dev.off()

saveRDS(filtered_YO,"youngOld.RDS")

# Save the DE genes to a text file
write.csv(filtered_YO, "GSE32719_YO_DE.csv") #sep = "\t",
             #quote = FALSE)

#############################################################################

n = nrow(filtered_YO)

cor.table = NULL
x = NULL
y = NULL
cor.val = NULL
cor.sig = NULL

for (i in 1 : (n-1)) {
  x_name = rownames(filtered_YO)[i]
  x_exps = filtered_YO[i, ]	
  
  for (j in (i+1) : n) {
    y_name = rownames(filtered_YO)[j]
    y_exps = filtered_YO[j, ]
    
    output = cor.test(x_exps,y_exps)
    
    x = c(x, x_name)
    y = c(y, y_name)
    cor.val = c(cor.val, output$estimate)
    cor.sig = c(cor.sig, output$p.value)
  }
}

cor.table = data.frame (x, y, cor.val, cor.sig)

dim(cor.table)
head(cor.table)

sig_cutoff = 0.001

cor.filtered = subset (cor.table, cor.sig < sig_cutoff)

dim(cor.filtered)
head(cor.filtered)


########## YM conditions ################

# Volcano: put the biological significance (fold-change)
# and statistical significance (p-value) in one plot
pdf("plots/foldYM_volcanoYM_I.pdf")
plot(foldYM, -log10(pvalueYM), main = "GSE32719 YM conditions - Volcano")
dev.off()


fold_cutoff_YM = 0.05
pvalue_cutoff_YM = 0.001
abline(v = fold_cutoff_YM, col = "blue", lwd = 3)
abline(v = -fold_cutoff_YM, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff_YM), col = "green", lwd = 3)

# Screen for the genes that satisfy the filtering criteria

# Fold-change filter for "biological" significance
filter_by_fold_YM = abs(foldYM) >= fold_cutoff_YM
dim(data2[filter_by_fold_YM, ])

# P-value filter for "statistical" significance
filter_by_pvalue_YM = pvalueYM <= pvalue_cutoff_YM
dim(data2[filter_by_pvalue_YM, ])

# Combined filter (both biological and statistical)
filter_combined_YM = filter_by_fold_YM & filter_by_pvalue_YM

filtered_YM = data2[filter_combined_YM,]
dim(filtered_YM)

head(filtered_YM)


# Let's generate the volcano plot again,
# highlighting the significantly differential expressed genes
pdf("plots/foldYM_volcanoYM_II.pdf")
plot(foldYM, -log10(pvalueYM), main = "GSE32719 YM conditions - Volcano #2")
points (foldYM[filter_combined_YM], -log10(pvalueYM[filter_combined_YM]),
        pch = 16, col = "red")
dev.off()
# Highlighting up-regulated in red and down-regulated in blue
pdf("plots/foldYM_volcanoYM_III.pdf")
plot(foldYM, -log10(pvalueYM), main = "GSE32719 YM conditions - Volcano #3")
points (foldYM[filter_combined_YM & foldYM < 0],
        -log10(pvalueYM[filter_combined_YM & foldYM < 0]),
        pch = 16, col = "blue")
points (foldYM[filter_combined_YM & foldYM > 0],
        -log10(pvalueYM[filter_combined_YM & foldYM > 0]),
        pch = 16, col = "red")
legend(1, 95, legend=c("up", "down"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
dev.off()

upRegulated_DEGs_YM  = foldYO[filter_combined_YM & foldYM < 0]
length(upRegulated_DEGs_YM)

downRegulated_DEGs_YM =  foldYO[filter_combined_YM & foldYM > 0]
length(downRegulated_DEGs_YM)

# Cluster the rows (genes) & columns (samples) by correlation
rowvYM = as.dendrogram(hclust(as.dist(1-cor(t(filtered_YM)))))
colvYM = as.dendrogram(hclust(as.dist(1-cor(filtered_YM))))

# Generate a heatmap
pdf("plots/heatmap_filtered_YM.pdf")
heatmap(filtered_YM, Rowv=rowvYM, Colv=colvYM, cexCol=0.7)
dev.off()
# install.packages("gplots")		# Uncomment if not already installed
# install.packages("RColorBrewer")	# Uncomment if not already installed

library(gplots)

# Enhanced heatmap
pdf("plots/heatmap2_filtered_YM.pdf")
heatmap.2(filtered_YM, Rowv=rowvYM, Colv=colvYM, cexCol=0.7,
          col = rev(redblue(256)), scale = "row")
dev.off()
# Save the heatmap to a PDF file
#pdf ("GSE32719_YM_DE_Heatmap.pdf")
heatmap.2(filtered_YM, Rowv=rowvYM, Colv=colvYM, cexCol=0.7,
          col = rev(redblue(256)), scale = "row")
#dev.off()

saveRDS(filtered_YM,"youngMiddle.RDS")
# Save the DE genes to a text file
write.csv (filtered_YM, "GSE32719_YM_DE.csv")#, sep = "\t",
             #quote = FALSE)

#############################################################################

n = nrow(filtered_YM)

cor.table = NULL
x = NULL
y = NULL
cor.val = NULL
cor.sig = NULL

for (i in 1 : (n-1)) {
  x_name = rownames(filtered_YM)[i]
  x_exps = filtered_YM[i, ]	
  
  for (j in (i+1) : n) {
    y_name = rownames(filtered_YM)[j]
    y_exps = filtered_YM[j, ]
    
    output = cor.test(x_exps,y_exps)
    
    x = c(x, x_name)
    y = c(y, y_name)
    cor.val = c(cor.val, output$estimate)
    cor.sig = c(cor.sig, output$p.value)
  }
}

cor.table = data.frame (x, y, cor.val, cor.sig)

dim(cor.table)
head(cor.table)

sig_cutoff = 0.001

cor.filtered = subset (cor.table, cor.sig < sig_cutoff)

dim(cor.filtered)
head(cor.filtered)


########## MO conditions ################

# Volcano: put the biological significance (fold-change)
# and statistical significance (p-value) in one plot
pdf("plots/foldMO_volcanoMO_I.pdf")
plot(foldMO, -log10(pvalueMO), main = "GSE32719 MO conditions - Volcano")
dev.off()

fold_cutoff_MO = 0.005
pvalue_cutoff_MO = 0.001
abline(v = fold_cutoff_MO, col = "blue", lwd = 3)
abline(v = -fold_cutoff_MO, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff_MO), col = "green", lwd = 3)

# Screen for the genes that satisfy the filtering criteria

# Fold-change filter for "biological" significance
filter_by_fold_MO = abs(foldMO) >= fold_cutoff_MO
dim(data2[filter_by_fold_MO, ])

# P-value filter for "statistical" significance
filter_by_pvalue_MO = pvalueMO <= pvalue_cutoff_MO
dim(data2[filter_by_pvalue_MO, ])

# Combined filter (both biological and statistical)
filter_combined_MO = filter_by_fold_MO & filter_by_pvalue_MO

filtered_MO = data2[filter_combined_MO,]
dim(filtered_MO)

head(filtered_MO)


# Let's generate the volcano plot again,
# highlighting the significantly differential expressed genes
pdf("plots/foldMO_volcanoMO_II.pdf")
plot(foldMO, -log10(pvalueMO), main = "GSE32719 MO conditions - Volcano #2")
points (foldMO[filter_combined_MO], -log10(pvalueMO[filter_combined_MO]),
        pch = 16, col = "red")
dev.off()
# Highlighting up-regulated in red and down-regulated in blue
pdf("plots/foldMO_volcanoMO_III.pdf")
plot(foldMO, -log10(pvalueMO), main = "GSE32719 MO conditions - Volcano #3")
points (foldMO[filter_combined_MO & foldMO < 0],
        -log10(pvalueMO[filter_combined_MO & foldMO < 0]),
        pch = 16, col = "red")
points (foldMO[filter_combined_MO & foldMO > 0],
        -log10(pvalueMO[filter_combined_MO & foldMO > 0]),
        pch = 16, col = "blue")

dev.off()

upRegulated_DEGs_MO = foldYO[filter_combined_MO & foldMO < 0]
length(upRegulated_DEGs_MO)


downRegulated_DEGs_MO =  foldMO[filter_combined_MO & foldMO > 0]
length(downRegulated_DEGs_MO)


write.csv(upRegulated_DEGs_YO,"up_YO.csv")
write.csv(downRegulated_DEGs_YO,"down_YO.csv")
write.csv(upRegulated_DEGs_YM,"up_YM.csv")
write.csv(downRegulated_DEGs_YM,"down_YM.csv")
write.csv(upRegulated_DEGs_MO,"up_MO.csv")
write.csv(downRegulated_DEGs_MO,"down_MO.csv")


# Cluster the rows (genes) & columns (samples) by correlation
rowvMO = as.dendrogram(hclust(as.dist(1-cor(t(filtered_MO)))))
colvMO = as.dendrogram(hclust(as.dist(1-cor(filtered_MO))))

# Generate a heatmap
pdf("plots/heatmap_filtered_MO.pdf")
heatmap(filtered_MO, Rowv=rowvMO, Colv=colvMO, cexCol=0.7)
dev.off()

# install.packages("gplots")		# Uncomment if not already installed
# install.packages("RColorBrewer")	# Uncomment if not already installed

library(gplots)

# Enhanced heatmap
pdf("plots/heatmap_filtered_MO.pdf")
heatmap.2(filtered_MO, Rowv=rowvMO, Colv=colvMO, cexCol=0.7,
          col = rev(redblue(256)), scale = "row")
dev.off()

# Save the heatmap to a PDF file
#pdf ("GSE32719_MO_DE_Heatmap.pdf")
heatmap.2(filtered_MO, Rowv=rowvMO, Colv=colvMO, cexCol=0.7,
          col = rev(redblue(256)), scale = "row")
#dev.off()

saveRDS(filtered_MO,"middleOld.RDS")

# Save the DE genes to a text file
write.csv (filtered_MO, "GSE32719_MO_DE.csv")#, sep = "\t",
             #quote = FALSE)

#############################################################################

n = nrow(filtered_MO)

cor.table = NULL
x = NULL
y = NULL
cor.val = NULL
cor.sig = NULL

for (i in 1 : (n-1)) {
  x_name = rownames(filtered_MO)[i]
  x_exps = filtered_MO[i, ]	
  
  for (j in (i+1) : n) {
    y_name = rownames(filtered_MO)[j]
    y_exps = filtered_MO[j, ]
    
    output = cor.test(x_exps,y_exps)
    
    x = c(x, x_name)
    y = c(y, y_name)
    cor.val = c(cor.val, output$estimate)
    cor.sig = c(cor.sig, output$p.value)
  }
}

cor.table = data.frame (x, y, cor.val, cor.sig)

dim(cor.table)
head(cor.table)

sig_cutoff = 0.001

cor.filtered = subset (cor.table, cor.sig < sig_cutoff)

dim(cor.filtered)
head(cor.filtered)



