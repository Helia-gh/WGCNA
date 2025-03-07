##### WGCNA  #####
# By Helia Ghazinejad 

library(WGCNA)
library(ggplot2)

###### Filtering & data prep #####

# Read TCGA files
BRCA_Normal<-read.table("TCGA_Normal.txt",row.names = 1,header=T)
BRCA_Tumor<-read.table("TCGA_Tumor.txt",row.names = 1,header=T)

# Merge the two datasets
Merged_Expression_Matrix<-cbind(BRCA_Normal,BRCA_Tumor)


# Count number of Tumor/Normal samples
Number_of_BRCA_Tumor<-length(BRCA_Tumor)
Number_of_Controls<-length(BRCA_Normal)

# Create phenotype (0 for control, 1 for case)
Phenotype<-c(rep(0,Number_of_Controls),rep(1,Number_of_BRCA_Tumor))

# Measure variance using MAD
gene_variance<-apply(Merged_Expression_Matrix,FUN = mad,1)

# Plot a histogram of variance
pdf("Variance Histogram.pdf")
hist(gene_variance,breaks = 100)
dev.off()

# Filter by median variance
Filtered_Expression <- Merged_Expression_Matrix[which(gene_variance > median(gene_variance)),]
filter_var <- apply(Filtered_Expression, FUN = mad, 1)
pdf("Filtered Histogram.pdf")
hist(filter_var, breaks=100)
dev.off()
    
  ###### WGCNA #####

# WGCNA expects a matrix or data frame in which columns are genes and rows ar samples. Hence, transpose!
WGCNA_Expression_Matrix<-t(Filtered_Expression)

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(WGCNA_Expression_Matrix, powerVector = powers, verbose = 5, networkType = "signed")

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
pdf("Similarity Matrix.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# This line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# We now calculate the adjacencies (similarity), using the soft thresholding power we chose from above
Power = 12
  
adjacency = adjacency(WGCNA_Expression_Matrix, power = Power,type = "signed");

# Turn adjacency into topological overlap
### Warning: This step can run for a while! Have patience ### 
TOM = TOMsimilarity(adjacency,TOMType = "signed");
dissTOM = 1-TOM
  
# Run hierarchical clustering on TOM matrix
geneTree = hclust(as.dist(dissTOM), method = "average");

# Module distinction using dynamic tree cut algorithm:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 30);

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
pdf("Gene Dendrogram.pdf")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(WGCNA_Expression_Matrix, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the resulting eigengene dendrogram
sizeGrWindow(7, 6)
pdf("Eigengene Dendrogram.pdf")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()

  ### Relating modules to trait ###

# Define numbers of genes and samples
nGenes = ncol(WGCNA_Expression_Matrix);
nSamples = nrow(WGCNA_Expression_Matrix);

# Run pearson correlation between module eigengene and trait of interest
moduleTraitCor = cor(MEs, Phenotype, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# Create correlation plot between trait and modules
sizeGrWindow(10,6)
# Display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
figureData <- data.frame(Modules = rownames(moduleTraitCor),
                         Tissue = matrix(nrow = nrow(moduleTraitCor), ncol = 1, data = "Tumor"),
                         Correlation = moduleTraitCor,
                         textMatrix = textMatrix)

pdf("Module Correlation.pdf")
ggplot(data = figureData, aes(x = Tissue, y = Modules, fill = Correlation))+
  scale_fill_gradient2(low = "blue", high = "red", midpoint = 0, limit = c(-1,1), space = "Lab", name = "Pearson\nCorrelation")+
  theme_minimal()+
  geom_tile()+
  geom_text(aes(Tissue, Modules, label = textMatrix), color = "black", size = 2.5)
dev.off()

# Identify hub-gene
hub_genes<-as.data.frame(chooseTopHubInEachModule(datExpr = WGCNA_Expression_Matrix,colorh = dynamicColors,type="signed"))



  ###### Summary output of network analysis results #####
# Define variable weight containing the weight column of datTrait

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(WGCNA_Expression_Matrix, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(WGCNA_Expression_Matrix, Phenotype, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", "TumorStatus", sep="");
names(GSPvalue) = paste("p.GS.","TumorStatus", sep="");


geneNames = colnames(WGCNA_Expression_Matrix)

# Create the starting data frame
geneInfo0 = data.frame(geneSymbol = geneNames,
                       moduleColor = dynamicColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, Phenotype, use = "p")));

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.TumorStatus));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo.csv")




  ################# Export to Cytoscape ####################


# Select modules
targetmodules = "salmon"

# Select module probes
inModule = is.finite(match(dynamicColors, targetmodules));
modGenes = geneNames[inModule];

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modGenes, modGenes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(targetmodules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.005
)

csvfile <- read.csv("geneInfo.csv")
