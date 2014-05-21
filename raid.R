# protein-RNA interaction
pr <- read.table("RNA-Protein-download.txt", header = T, stringsAsFactors=F, quote = "", sep = '\t')

prInteractionTypes <- unique(pr[, c(5, 9)])

# lncRNA protein interaction
lncp <- unique(pr[pr[, 5] == 'lncRNA', c(3,7,12)])
nrow(lncp)
# the number of lncRNA related protein interactions
# 333

# the number of lncRNAs
length(unique(lncp[, 1]))
# 77

# the number of proteins
length(unique(lncp[, 2]))
# 154

# the number of literatures
length(unique(lncp[, 3]))
# 139

# rna rna interaction
rr <- read.table("RNA-RNA-download.txt", header = T, stringsAsFactors=F, quote = "", sep = '\t')

rrInteractionTypes <- unique(rr[, c(5, 9)])

#lncRNA rna interaction
lncr <- unique(rr[rr[, 5] == 'lncRNA' | rr[, 9] == 'lncRNA', c(2,5,7,9,12)])
# lncRNA rna interaction types
lncrTypes <- apply(lncr[, c(2, 4)], 1, function(x) {
  y <- sort(x)
  y <- paste(sort(x), collapse="; ")
  return(y)
})

# the number of lncRNA related

length(unique(lncr[, 4]))

sort(table(lncrTypes), decreasing=T)
# lncrTypes
# lncRNA; mRNA   lncRNA; miRNA  lncRNA; lncRNA   lncRNA; other 
# 165              54               8               1 
# lncRNA; Protein 
# 1 

length(unique(lncr[, 1]))


