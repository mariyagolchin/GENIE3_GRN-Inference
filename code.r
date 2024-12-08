
exprMatr <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)
rownames(exprMatr) <- paste("Gene", 1:20, sep="")
colnames(exprMatr) <- paste("Sample", 1:5, sep="")
head(exprMatr)

# > head(exprMatr)
#       Sample1 Sample2 Sample3 Sample4 Sample5
# Gene1       7       3       9       5       1
# Gene2       6       8       7       6       9
# Gene3       9       6       7       8       2
# Gene4       8       7       3       1       8
# Gene5       9       6       6       3       6
# Gene6       7       8       2      10       1
# > 

# =================================================
# Run GENEI3 on gene expression data exprMatr
# =================================================
library(GENIE3)
set.seed(123) # For reproducibility of results
weightMat <- GENIE3(exprMatr)
dim(weightMat)
# > weightMat[1:5,1:5]
#             Gene1      Gene2      Gene3      Gene4      Gene5
# Gene1 0.000000000 0.13093249 0.10065088 0.04939223 0.02431176
# Gene2 0.081126576 0.00000000 0.11123464 0.03868114 0.02565248
# Gene3 0.087517620 0.14730375 0.00000000 0.05471999 0.06792992
# Gene4 0.022539722 0.02635689 0.01698693 0.00000000 0.09233107
# Gene5 0.009497612 0.01267937 0.01272033 0.11474264 0.00000000
# > 

The algorithm outputs a matrix containing the weights of the putative regulatory links,
with higher weights corresponding to more likely regulatory links. 
weightMat[i,j] is the weight of the link directed from the i-th gene to j-th gene.

Restrict the candidate regulators to a subset of genes
By default, all the genes in exprMatr are used as candidate regulators.
The list of candidate regulators can however be restricted to a subset of genes. 
This can be useful when you know which genes are transcription factors.

# Genes that are used as candidate regulators
regulators <- c(2, 4, 7)
# Or alternatively:
regulators <- c("Gene2", "Gene4", "Gene7")
weightMat <- GENIE3(exprMatr, regulators=regulators)
> weightMat[1:3,1:5]
          Gene1     Gene2      Gene3    Gene4     Gene5
Gene2 0.3894916 0.0000000 0.55763560 0.424013 0.1651934
Gene4 0.1221831 0.3561855 0.07950583 0.000000 0.5393023
Gene7 0.4883253 0.6438145 0.36285858 0.575987 0.2955042
> 

regulatorsList <- list("Gene1"=rownames(exprMatr)[1:10],
                       "Gene2"=rownames(exprMatr)[10:20],
                       "Gene20"=rownames(exprMatr)[15:20])
                    
#    This setup suggests that these genes are considered as potential regulators
#     for different groups within the network inference process.

set.seed(123)
weightList <- GENIE3(exprMatr, nCores=1, targets=names(regulatorsList), regulators=regulatorsList, returnMatrix=FALSE)


# You can obtain the list of all the regulatory links (from most likely to least likely) with this command:

linkList <- getLinkList(weightMat)
dim(linkList)
## [1] 57  3
head(linkList)

#   regulatoryGene targetGene    weight
# 1          Gene2     Gene19 0.7975999
# 2          Gene7      Gene2 0.6438145
# 3          Gene7     Gene17 0.6314861
# 4          Gene7     Gene14 0.6198380
# 5          Gene4     Gene16 0.6177857
# 6          Gene2      Gene7 0.6138833
# > 