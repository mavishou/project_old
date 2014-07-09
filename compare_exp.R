fCompare <- function(x, y, m="pearson"){
  for(i in 2:5){
    print(cor(x[, i], y[, i], method=m))
  }
}

setwd('/lustre/user/houm/projects/AnnoLnc/expression')
# cuff nornal gene count
gcNorm <- read.table('cuff/cuffnorm/genes.count_table', header=T, stringsAsFactors=F)
# cuff normal gene fpkm
gfNorm <- read.table('cuff/cuffnorm/genes.fpkm_table', header=T, stringsAsFactors=F)
# gcNorm[, c(2:4)] <- apply(gcNorm[, c(2:4)], 2, as.numeric)
# gfNorm[, c(2:4)] <- apply(gfNorm[, c(2:4)], 2, as.numeric)

# ----------------比较自己的count和fpkm--------------------------
plot(gcNorm[, 2], gfNorm[, 2], pch=20)
for(i in 2:5){
  print(cor(gcNorm[, i], gfNorm[, i], method='spearman'))
}
# 说明同一个组内的count和fpkm其实差不多，也就是说count是对reads长度进行了normalzie的
# 也就是说cuffnorm的作用其实只是在组间进行normalzie
# [1] 0.9772484
# [1] 0.9849432
# [1] 0.9907596
# [1] 0.9897649

# --------------比较normal和mask---------------------------------
# 为了避免组间normalize带来的影响，这里用count来比较
gcMask <- read.table('cuff/mask_rRNA/cuffnorm/genes.count_table', header=T, stringsAsFactors=F)
# check whether the order of genes are the same
all(gcMask[, 1] == gcNorm[, 1])
for(i in 2:5){
  print (cor(gcNorm[, i], gcMask[, i]))
}
plot(gcNorm[, 2], gcMask[, 2], pch=20)

# 基本没有区别，说明mask基本上没有影响
# [1] 0.9988739
# [1] 0.9949093
# [1] 0.9999801
# [1] 0.9948179

# Let's check the transcript level
tcNorm <- read.table('cuff/cuffnorm/isoforms.count_table', header=T, stringsAsFactors=F)
tcMask <- read.table('cuff/mask_rRNA/cuffnorm/isoforms.count_table', header=T, stringsAsFactors=F)
# check the order
all(tcNorm[, 1]==tcMask[, 1])
for(i in 2:5){
  print(cor(tcNorm[, i], tcMask[, i]))
}
# transcript level上也是基本上没有区别，以后可以不用mask了
# [1] 0.9983689
# [1] 0.9930249
# [1] 0.9999702
# [1] 0.9917543

# --------------比较normal和correct-------------------------------
### gene count
gcCorrect <- read.table('cuff/correct/cuffnorm/genes.count_table', header=T, stringsAsFactors=F)
fCompare(gcNorm, gcCorrect)
# 还是一样的
# [1] 0.9988739
# [1] 0.9949093
# [1] 0.9999801
# [1] 0.9948179

### gene fpkm
gfCorrect <- gcCorrect <- read.table('cuff/correct/cuffnorm/genes.fpkm_table', header=T, stringsAsFactors=F)
fCompare(gfNorm, gfCorrect)
# [1] 0.9982951
# [1] 0.9562654
# [1] 0.9999876
# [1] 0.9528909

### transcript level
tcCorrect <- read.table('cuff/correct/cuffnorm/isoforms.count_table', header=T, stringsAsFactors=F)
fCompare(tcNorm, tcCorrect)
# [1] 0.9983689
# [1] 0.9930249
# [1] 0.9999702
# [1] 0.9917543

# --------------比较normal和no_length-------------------------------
### gene count
gcNl <- read.table('cuff/no_length/cuffnorm/genes.count_table', header=T, stringsAsFactors=F)
fCompare(gcNorm, gcNl)
# [1] 0.9988739
# [1] 0.9949093
# [1] 0.9999801
# [1] 0.9948179
plot(gcNorm[, 2], gcNl[, 2], pch=20)

### gene fpkm
gfNl <- read.table('cuff/no_length/cuffnorm/genes.fpkm_table', header=T, stringsAsFactors=F)
# check the order
all(gfNl[, 1] == gfNorm[, 1])
fCompare(gfNorm, gfNl)
# [1] 0.1080504
# [1] 0.4147376
# [1] 0.1051083
# [1] 0.6595122
fCompare(gfNorm, gfNl, 'spearman')
# [1] 0.9767756
# [1] 0.9843196
# [1] 0.9904253
# [1] 0.9894595

plot(gfNorm[, 2], gfNl[, 2], pch=20)

# ---------------比较norm和HTSeq-----------------------------
### gene count
gcHtseq <- read.table('HTSeq/gene.txt', sep='\t', stringsAsFactors=F)
# check ID
gCuff <- gcNorm[, 1]
gHtseq <- gcHtseq[, 1]
all(gCuff %in% gHtseq)
# setdiff(gHtseq, gCuff)
gcHtseq <- gcHtseq[gcHtseq[, 1] %in% gCuff, ]
# check order
all(gcHtseq[, 1] == gcNorm[, 1])

fCompare(gcNorm, gcHtseq)
# [1] 0.6106263
# [1] 0.3612847
# [1] 0.3582761
# [1] 0.9885178
fCompare(gcNorm, gcHtseq, 'spearman')
# [1] 0.9195368
# [1] 0.9267516
# [1] 0.9134145
# [1] 0.9461353

plot(gcNorm[, 2], gcHtseq[, 2], pch=20, xlim=c(0, 8e5), ylim=c(0, 8e5), 
     xlab="Cuff", ylab="HTSeq", main="ERR030882 genes")
text(6e5, 6e5, 'Pearson cor: 0.61', cex=1.2)
text(6e5, 7e5, 'Spearman cor: 0.92', cex=1.2)

plot(gcNorm[, 5], gcHtseq[, 5], pch=20, 
     xlim=c(0, 4e5), ylim=c(0, 4e5), 
     xlab="Cuff", ylab="HTSeq", main="GSE24399_GSM601407 genes")
text(2e5, 2e5, 'Pearson cor: 0.99', cex=1.2)
text(2e5, 3e5, 'Spearman cor: 0.95', cex=1.2)

### transcript count
tcHtseq <- read.table('HTSeq/transcript.txt', sep='\t', stringsAsFactors=F)
tcHtseq <- tcHtseq[tcHtseq[, 1] %in% tcNorm[, 1], ]
# check ordr
all(tcNorm[, 1] == tcHtseq[, 1])

fCompare(tcNorm, tcHtseq)
# [1] 0.1980898
# [1] 0.1412311
# [1] 0.09203862
# [1] 0.6937974

fCompare(tcNorm, tcHtseq, 'spearman')
# [1] 0.5125536
# [1] 0.5286943
# [1] 0.5552942
# [1] 0.4875728

plot(tcNorm[, 2], tcHtseq[, 2], pch=20, 
     xlim=c(0, 8e5), ylim=c(0, 8e5), 
     xlab="Cuff", ylab="HTSeq", main="ERR030882 transcripts")
text(4e5, 4e5, 'Pearson cor: 0.20', cex=1.2)
text(4e5, 5e5, 'Spearman cor: 0.51', cex=1.2)

# ---------------------------Norm2--------------------------------
gcNorm2 <- read.table('cuff/cuffnorm2/genes.count_table', header=T, stringsAsFactors=F)
gfNorm2 <- read.table('cuff/cuffnorm2/genes.fpkm_table', header=T, stringsAsFactors=F)

# -----------------------featureCount---------------------------
gcFc <- read.table('featureCounts/all_gene_counts.txt', sep='\t', stringsAsFactors=F)
# genes of featurecount
gFc <- gcFc[, 1]
all(gCuff %in% gFc)
rownames(gcFc) <- gFc
gcFc <- gcFc[, ]
gcFc <- gcFc[gCuff, ]
# check order
all(gcFc[, 1] == gcNorm[, 1])

# compare featureCount with HTSeq
fCompare(gcHtseq, gcFc)
> fCompare(gcHtseq, gcFc)
# [1] 0.9930448
# [1] 0.9989145
# [1] 0.996702
# [1] 1

# calculate FPKM using reads count
# total mass
totalMass <- c(4.92436e+07, 5.08075e+07, 1.15251e+08, 1.20944e+07)
# gene length
geneLength <- read.table('featureCounts/gene_length.txt', stringsAsFactors=F)
geneLength <- geneLength[, -1]
names(geneLength) <- gFc
geneLength <- geneLength[gCuff]

gfMy <- gfNorm
gfMy[, 2:5] <- 0
for(i in 2:5){
  gfMy[, i] <- gcFc[, i] * 10^9 / (totalMass[i-1] * geneLength)
}

fCompare(gfMy, gfNorm2)
# [1] 0.02718042
# [1] 0.1237496
# [1] 0.01356684
# [1] 0.7207268

fCompare(gfMy, gfNorm2, 'spearman')
# [1] 0.8766704
# [1] 0.8983167
# [1] 0.8718634
# [1] 0.9215698

plot(gfNorm2[, 2], gfMy[, 2], pch=20, xlim=c(0, 1000), ylim=c(0, 1000))


