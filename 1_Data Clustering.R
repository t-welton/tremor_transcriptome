# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
workingDir = "~"
setwd(workingDir);
# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#
exprMat = read.csv("expression_matrix.csv");
# Take a quick look at what is in the data set:
dim(exprMat);
names(exprMat);
datExpr0 = as.data.frame(t(exprMat[, -1]));
names(datExpr0) = exprMat$gene;
rownames(datExpr0) = names(exprMat)[-1];
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
#abline(h = 15, col = "red");
# Determine cluster under the line
#clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
#table(clust)
# clust 1 contains the samples we want to keep.
#keepSamples = (clust==1)
datExpr = datExpr0#[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
##read clin data
allTraits = read.csv("binary_vectors_Brain_Cortex.csv");
##dim(traitData)
##names(traitData)
# remove columns that hold information we do not need.
##allTraits = traitData[, -c(31, 16)];
##allTraits = allTraits[, c(2, 11:36) ];
dim(allTraits)
names(allTraits)
# Form a data frame analogous to expression data that will hold the clinical traits.
#phenotype = rownames(datExpr);
#traitRows = match(phenotype, allTraits$phenotype);
datTraits = allTraits[, -1];
rownames(datTraits) = allTraits[, 1];
collectGarbage();
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = TRUE); # change signed=T for our analysis
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
save(datExpr, datTraits, file = "BrainCortex-01a-dataInput.RData")
