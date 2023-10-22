library(WGCNA)
#Network dendrogram----
datExpr1_combine <- read.table("datExpr1_combine_test.txt",header = T)
dim(datExpr1_combine)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()
allowWGCNAThreads()
enableWGCNAThreads(nThreads=6)
#We first check for genes and samples with too many missing values
gsg = goodSamplesGenes(datExpr1_combine, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  # Optioncombiney, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0)
    a <-
      printFlush(paste("Removing genes:", paste(names(datExpr1_combine)[!gsg$goodGenes], collapse = ", ")))
  
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr1_combine)[!gsg$goodSamples], collapse = ", ")))
  
  # Remove the offending genes and samples from the data:
  datExpr1_combine = datExpr1_combine[gsg$goodSamples, gsg$goodGenes]
}
#Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers.
sampleTree = hclust(dist(datExpr1_combine), method = "median")
dim(datExpr1_combine)
#[1]    107 14516
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too smcombine.
sizeGrWindow(12, 9)
pdf(file = "sample_clustering.pdf", width = 12, height = 9);
par(mar = c(0, 4, 2, 0), cex = 0.6)
plot(
  sampleTree,
  main = "Sample clustering to detect outliers",
  sub = "",
  xlab = "",
  cex.lab = 1.5,
  cex.axis = 1.5,
  cex.main = 2
)
abline(h = 60, col = "red");
dev.off()

#Choosing the soft-thresholding power---------------------------------------------------

# Choose a set of soft-thresholding powers
powers = c(1:30)
# Call the network topology analysis function
sft = pickSoftThreshold(
  datExpr1_combine,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file = "fit_index.pdf", width = 12, height = 9)
plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit,signed R^2",
  type = "n",
  main = paste("Scale independence")
)

text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red"
)
dev.off()

# Mean connectivity as a function of the soft-thresholding power
pdf(file = "mean connectivity.pdf", width = 12, height = 9)
plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = paste("Mean connectivity")
)
text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1,
  col = "red"
)
dev.off()

power = sft$powerEstimate
power
#One-step network construction and module detection----
net_combine_signed = blockwiseModules(
  datExpr1_combine,
  power = power,
  maxBlockSize = 20000,
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.2,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "TOM",
  verbose = 3
)


# Convert labels to colors for plotting
moduleLabels_combine_signed = net_combine_signed$colors
moduleColors_combine_signed = labels2colors(moduleLabels_combine_signed)
length(table(moduleColors_combine_signed))
#[1] 13
# Convert labels to colors for plotting
mergedColors = labels2colors(net_combine_signed$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net_combine_signed$dendrograms[[1]], mergedColors[net_combine_signed$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#Quantifying module–trait associations---------------------------------------
# Define numbers of genes and samples
nGenes = ncol(datExpr1_combine)
nSamples = nrow(datExpr1_combine)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr1_combine,
                        moduleColors_combine_signed,
                        nPC = 1,
                        excludeGrey = T)$eigengenes
MEs = orderMEs(MEs0)

# module eigengene
MEs = net_combine_signed$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# Correlation diagram between modules obtained by clustering according to expression levels between genes
# marDendro/marHeatmap Set bottom, left, top and right margins
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
#Gene relationship to trait and important modules----
modNames = substring(names(MEs_col), 3)
geneModuleMembership = as.data.frame(cor(datExpr1_combine, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")
#genes with the highest correlation with the MEs----
module = "black"
probes = colnames(datExpr1_combine)
inModule = (moduleColors_combine_signed == module)
modProbes_black = probes[inModule]
MEs_black <-
  geneModuleMembership[which(rownames(geneModuleMembership) %in% modProbes_black), ]
MEs_black <-
  data.frame(gene = rownames(MEs_black), kME = MEs_black$MMblack)
MEs_black_high <- MEs_black[which(abs(MEs_black$kME) >= 0.7), ]
MEs_black_high <-
  MEs_black_high[order(MEs_black_high$kME, decreasing = T), ]

module = "blue"
probes = colnames(datExpr1_combine)
inModule = (moduleColors_combine_signed == module)
modProbes_blue = probes[inModule]
MEs_blue <-
  geneModuleMembership[which(rownames(geneModuleMembership) %in% modProbes_blue), ]
MEs_blue <-
  data.frame(gene = rownames(MEs_blue), kME = MEs_blue$MMblue)
MEs_blue_high <- MEs_blue[which(abs(MEs_blue$kME) >= 0.7), ]
MEs_blue_high <-
  MEs_blue_high[order(MEs_blue_high$kME, decreasing = T), ]


module = "brown"
probes = colnames(datExpr1_combine)
inModule = (moduleColors_combine_signed == module)
modProbes_brown = probes[inModule]
MEs_brown <-
  geneModuleMembership[which(rownames(geneModuleMembership) %in% modProbes_brown), ]
MEs_brown <-
  data.frame(gene = rownames(MEs_brown), kME = MEs_brown$MMbrown)
MEs_brown_high <- MEs_brown[which(abs(MEs_brown$kME) >= 0.7), ]
MEs_brown_high <-
  MEs_brown_high[order(MEs_brown_high$kME, decreasing = T), ]

module = "green"
probes = colnames(datExpr1_combine)
inModule = (moduleColors_combine_signed == module)
modProbes_green = probes[inModule]
MEs_green <-
  geneModuleMembership[which(rownames(geneModuleMembership) %in% modProbes_green), ]
MEs_green <-
  data.frame(gene = rownames(MEs_green), kME = MEs_green$MMgreen)
MEs_green_high <- MEs_green[which(abs(MEs_green$kME) >= 0.7), ]
MEs_green_high <-
  MEs_green_high[order(MEs_green_high$kME, decreasing = T), ]

module = "greenyellow"
probes = colnames(datExpr1_combine)
inModule = (moduleColors_combine_signed == module)
modProbes_greenyellow = probes[inModule]
MEs_greenyellow <-
  geneModuleMembership[which(rownames(geneModuleMembership) %in% modProbes_greenyellow), ]
MEs_greenyellow <-
  data.frame(gene = rownames(MEs_greenyellow), kME = MEs_greenyellow$MMgreenyellow)
MEs_greenyellow_high <- MEs_greenyellow[which(abs(MEs_greenyellow$kME) >= 0.7), ]
MEs_greenyellow_high <-
  MEs_greenyellow_high[order(MEs_greenyellow_high$kME, decreasing = T), ]

module = "magenta"
probes = colnames(datExpr1_combine)
inModule = (moduleColors_combine_signed == module)
modProbes_magenta = probes[inModule]
MEs_magenta <-
  geneModuleMembership[which(rownames(geneModuleMembership) %in% modProbes_magenta), ]
MEs_magenta <-
  data.frame(gene = rownames(MEs_magenta), kME = MEs_magenta$MMmagenta)
MEs_magenta_high <- MEs_magenta[which(abs(MEs_magenta$kME) >= 0.7), ]
MEs_magenta_high <-
  MEs_magenta_high[order(MEs_magenta_high$kME, decreasing = T), ]

module = "pink"
probes = colnames(datExpr1_combine)
inModule = (moduleColors_combine_signed == module)
modProbes_pink = probes[inModule]
MEs_pink <-
  geneModuleMembership[which(rownames(geneModuleMembership) %in% modProbes_pink), ]
MEs_pink <-
  data.frame(gene = rownames(MEs_pink), kME = MEs_pink$MMpink)
MEs_pink_high <- MEs_pink[which(abs(MEs_pink$kME) >= 0.7), ]
MEs_pink_high <-
  MEs_pink_high[order(MEs_pink_high$kME, decreasing = T), ]

module = "purple"
probes = colnames(datExpr1_combine)
inModule = (moduleColors_combine_signed == module)
modProbes_purple = probes[inModule]
MEs_purple <-
  geneModuleMembership[which(rownames(geneModuleMembership) %in% modProbes_purple), ]
MEs_purple <-
  data.frame(gene = rownames(MEs_purple), kME = MEs_purple$MMpurple)
MEs_purple_high <- MEs_purple[which(abs(MEs_purple$kME) >= 0.7), ]
MEs_purple_high <-
  MEs_purple_high[order(MEs_purple_high$kME, decreasing = T), ]

module = "red"
probes = colnames(datExpr1_combine)
inModule = (moduleColors_combine_signed == module)
modProbes_red = probes[inModule]
MEs_red <-
  geneModuleMembership[which(rownames(geneModuleMembership) %in% modProbes_red), ]
MEs_red <-
  data.frame(gene = rownames(MEs_red), kME = MEs_red$MMred)
MEs_red_high <- MEs_red[which(abs(MEs_red$kME) >= 0.7), ]
MEs_red_high <-
  MEs_red_high[order(MEs_red_high$kME, decreasing = T), ]

module = "tan"
probes = colnames(datExpr1_combine)
inModule = (moduleColors_combine_signed == module)
modProbes_tan = probes[inModule]
MEs_tan <-
  geneModuleMembership[which(rownames(geneModuleMembership) %in% modProbes_tan), ]
MEs_tan <-
  data.frame(gene = rownames(MEs_tan), kME = MEs_tan$MMtan)
MEs_tan_high <- MEs_tan[which(abs(MEs_tan$kME) >= 0.7), ]
MEs_tan_high <-
  MEs_tan_high[order(MEs_tan_high$kME, decreasing = T), ]

module = "turquoise"
probes = colnames(datExpr1_combine)
inModule = (moduleColors_combine_signed == module)
modProbes_turquoise = probes[inModule]
MEs_turquoise <-
  geneModuleMembership[which(rownames(geneModuleMembership) %in% modProbes_turquoise), ]
MEs_turquoise <-
  data.frame(gene = rownames(MEs_turquoise), kME = MEs_turquoise$MMturquoise)
MEs_turquoise_high <- MEs_turquoise[which(abs(MEs_turquoise$kME) >= 0.7), ]
MEs_turquoise_high <-
  MEs_turquoise_high[order(MEs_turquoise_high$kME, decreasing = T), ]

module = "yellow"
probes = colnames(datExpr1_combine)
inModule = (moduleColors_combine_signed == module)
modProbes_yellow = probes[inModule]
MEs_yellow <-
  geneModuleMembership[which(rownames(geneModuleMembership) %in% modProbes_yellow), ]
MEs_yellow <-
  data.frame(gene = rownames(MEs_yellow), kME = MEs_yellow$MMyellow)
MEs_yellow_high <- MEs_yellow[which(abs(MEs_yellow$kME) >= 0.7), ]
MEs_yellow_high <-
  MEs_yellow_high[order(MEs_yellow_high$kME, decreasing = T), ]

#coexpression enrichment----
allgene <- colnames(datExpr1_combine)
top50 <- read.table("top50.txt")
candidate1 <- top50$V1
candidate1 <- rare_mac3_saige_lof_dmis_all1_top50[which(rare_mac3_saige_lof_dmis_all1_top50$odds_ratio>=1),"gene"]
candidate1 <- rare_mac3_saige_lof_dmis_all1_top50$gene

MEs_brown <- read.table("MEs_brown.txt")
colnames(MEs_brown) <- c("gene")
dim(MEs_brown)
MEs_brown_high <- read.table("MEs_brown_high.txt")
colnames(MEs_brown_high) <- c("gene")
dim(MEs_brown_high)
intersect(MEs_brown$gene,candidate1)
intersect(MEs_brown_high$gene,candidate1)

gene_overlap <- intersect(MEs_brown_high$gene,candidate1)
gene_overlap_other <- intersect(setdiff(allgene,MEs_brown$gene),candidate1)
compare <-
  matrix(c(
    length(gene_overlap),
    length(MEs_brown_high$gene),
    length(gene_overlap_other),
    length(setdiff(allgene,MEs_brown_high$gene))),
    nrow = 2)
result <- fisher.test(compare,conf.level = 0.95,alternative = "g")
result_brown <- cbind(length(gene_overlap),length(MEs_brown_high$gene),length(gene_overlap_other),
                      length(setdiff(allgene,MEs_brown_high$gene)),result$p.value,result$estimate)
#HM
{
  MEs_black <- read.table("MEs_black.txt")
  colnames(MEs_black) <- c("gene")
  MEs_black_high <- read.table("MEs_black_high.txt")
  colnames(MEs_black_high) <- c("gene")
  intersect(MEs_black$gene,top50)
  gene_overlap <- intersect(MEs_black_high$gene,top50)
  gene_overlap_other <- intersect(setdiff(allgene,MEs_black$gene),top50)
  compare <-
    matrix(c(
      length(gene_overlap),
      length(MEs_black_high$gene),
      length(gene_overlap_other),
      length(setdiff(allgene,MEs_black_high$gene))),
      nrow = 2)
  result <- fisher.test(compare ,conf.level = 0.95,alternative = "g")
  result_black <- cbind(length(gene_overlap),length(MEs_black_high$gene),length(gene_overlap_other),
                        length(setdiff(allgene,MEs_black_high$gene)),result$p.value,result$estimate)
  
  MEs_blue <- read.table("MEs_blue.txt")
  colnames(MEs_blue) <- c("gene")
  MEs_blue_high <- read.table("MEs_blue_high.txt")
  colnames(MEs_blue_high) <- c("gene")
  intersect(MEs_blue$gene,top50)
  gene_overlap <- intersect(MEs_blue_high$gene,top50)
  gene_overlap_other <- intersect(setdiff(allgene,MEs_blue$gene),top50)
  compare <-
    matrix(c(
      length(gene_overlap),
      length(MEs_blue_high$gene),
      length(gene_overlap_other),
      length(setdiff(allgene,MEs_blue_high$gene))),
      nrow = 2)
  result <- fisher.test(compare ,conf.level = 0.95,alternative = "g")
  result_blue <- cbind(length(gene_overlap),length(MEs_blue_high$gene),length(gene_overlap_other),
                       length(setdiff(allgene,MEs_blue_high$gene)),result$p.value,result$estimate)
  
  MEs_brown <- read.table("MEs_brown.txt")
  colnames(MEs_brown) <- c("gene")
  MEs_brown_high <- read.table("MEs_brown_high.txt")
  colnames(MEs_brown_high) <- c("gene")
  intersect(MEs_brown$gene,top50)
  gene_overlap <- intersect(MEs_brown_high$gene,top50)
  gene_overlap_other <- intersect(setdiff(allgene,MEs_brown$gene),top50)
  compare <-
    matrix(c(
      length(gene_overlap),
      length(MEs_brown_high$gene),
      length(gene_overlap_other),
      length(setdiff(allgene,MEs_brown_high$gene))),
      nrow = 2)
  result <- fisher.test(compare ,conf.level = 0.95,alternative = "g")
  result_brown <- cbind(length(gene_overlap),length(MEs_brown_high$gene),length(gene_overlap_other),
                        length(setdiff(allgene,MEs_brown_high$gene)),result$p.value,result$estimate)
  
  MEs_green <- read.table("MEs_green.txt")
  colnames(MEs_green) <- c("gene")
  MEs_green_high <- read.table("MEs_green_high.txt")
  colnames(MEs_green_high) <- c("gene")
  intersect(MEs_green$gene,top50)
  gene_overlap <- intersect(MEs_green_high$gene,top50)
  gene_overlap_other <- intersect(setdiff(allgene,MEs_green$gene),top50)
  compare <-
    matrix(c(
      length(gene_overlap),
      length(MEs_green_high$gene),
      length(gene_overlap_other),
      length(setdiff(allgene,MEs_green_high$gene))),
      nrow = 2)
  result <- fisher.test(compare ,conf.level = 0.95,alternative = "g")
  result_green <- cbind(length(gene_overlap),length(MEs_green_high$gene),length(gene_overlap_other),
                        length(setdiff(allgene,MEs_green_high$gene)),result$p.value,result$estimate)
  
  MEs_greenyellow <- read.table("MEs_greenyellow.txt")
  colnames(MEs_greenyellow) <- c("gene")
  MEs_greenyellow_high <- read.table("MEs_greenyellow_high.txt")
  colnames(MEs_greenyellow_high) <- c("gene")
  intersect(MEs_greenyellow$gene,top50)
  gene_overlap <- intersect(MEs_greenyellow_high$gene,top50)
  gene_overlap_other <- intersect(setdiff(allgene,MEs_greenyellow$gene),top50)
  compare <-
    matrix(c(
      length(gene_overlap),
      length(MEs_greenyellow_high$gene),
      length(gene_overlap_other),
      length(setdiff(allgene,MEs_greenyellow_high$gene))),
      nrow = 2)
  result <- fisher.test(compare ,conf.level = 0.95,alternative = "g")
  result_greenyellow <- cbind(length(gene_overlap),length(MEs_greenyellow_high$gene),length(gene_overlap_other),
                              length(setdiff(allgene,MEs_greenyellow_high$gene)),result$p.value,result$estimate)
  
  MEs_magenta <- read.table("MEs_magenta.txt")
  colnames(MEs_magenta) <- c("gene")
  MEs_magenta_high <- read.table("MEs_magenta_high.txt")
  colnames(MEs_magenta_high) <- c("gene")
  intersect(MEs_magenta$gene,top50)
  gene_overlap <- intersect(MEs_magenta_high$gene,top50)
  gene_overlap_other <- intersect(setdiff(allgene,MEs_magenta$gene),top50)
  compare <-
    matrix(c(
      length(gene_overlap),
      length(MEs_magenta_high$gene),
      length(gene_overlap_other),
      length(setdiff(allgene,MEs_magenta_high$gene))),
      nrow = 2)
  result <- fisher.test(compare ,conf.level = 0.95,alternative = "g")
  result_magenta <- cbind(length(gene_overlap),length(MEs_magenta_high$gene),length(gene_overlap_other),
                          length(setdiff(allgene,MEs_magenta_high$gene)),result$p.value,result$estimate)
  
  MEs_pink <- read.table("MEs_pink.txt")
  colnames(MEs_pink) <- c("gene")
  MEs_pink_high <- read.table("MEs_pink_high.txt")
  colnames(MEs_pink_high) <- c("gene")
  intersect(MEs_pink$gene,top50)
  gene_overlap <- intersect(MEs_pink_high$gene,top50)
  gene_overlap_other <- intersect(setdiff(allgene,MEs_pink$gene),top50)
  compare <-
    matrix(c(
      length(gene_overlap),
      length(MEs_pink_high$gene),
      length(gene_overlap_other),
      length(setdiff(allgene,MEs_pink_high$gene))),
      nrow = 2)
  result <- fisher.test(compare ,conf.level = 0.95,alternative = "g")
  result_pink <- cbind(length(gene_overlap),length(MEs_pink_high$gene),length(gene_overlap_other),
                       length(setdiff(allgene,MEs_pink_high$gene)),result$p.value,result$estimate)
  
  MEs_purple <- read.table("MEs_purple.txt")
  colnames(MEs_purple) <- c("gene")
  MEs_purple_high <- read.table("MEs_purple_high.txt")
  colnames(MEs_purple_high) <- c("gene")
  intersect(MEs_purple$gene,top50)
  gene_overlap <- intersect(MEs_purple_high$gene,top50)
  gene_overlap_other <- intersect(setdiff(allgene,MEs_purple$gene),top50)
  compare <-
    matrix(c(
      length(gene_overlap),
      length(MEs_purple_high$gene),
      length(gene_overlap_other),
      length(setdiff(allgene,MEs_purple_high$gene))),
      nrow = 2)
  result <- fisher.test(compare ,conf.level = 0.95,alternative = "g")
  result_purple <- cbind(length(gene_overlap),length(MEs_purple_high$gene),length(gene_overlap_other),
                         length(setdiff(allgene,MEs_purple_high$gene)),result$p.value,result$estimate)
  
  MEs_red <- read.table("MEs_red.txt")
  colnames(MEs_red) <- c("gene")
  MEs_red_high <- read.table("MEs_red_high.txt")
  colnames(MEs_red_high) <- c("gene")
  intersect(MEs_red$gene,top50)
  gene_overlap <- intersect(MEs_red_high$gene,top50)
  gene_overlap_other <- intersect(setdiff(allgene,MEs_red$gene),top50)
  compare <-
    matrix(c(
      length(gene_overlap),
      length(MEs_red_high$gene),
      length(gene_overlap_other),
      length(setdiff(allgene,MEs_red_high$gene))),
      nrow = 2)
  result <- fisher.test(compare ,conf.level = 0.95,alternative = "g")
  result_red <- cbind(length(gene_overlap),length(MEs_red_high$gene),length(gene_overlap_other),
                      length(setdiff(allgene,MEs_red_high$gene)),result$p.value,result$estimate)
  
  MEs_tan <- read.table("MEs_tan.txt")
  colnames(MEs_tan) <- c("gene")
  MEs_tan_high <- read.table("MEs_tan_high.txt")
  colnames(MEs_tan_high) <- c("gene")
  intersect(MEs_tan$gene,top50)
  gene_overlap <- intersect(MEs_tan_high$gene,top50)
  gene_overlap_other <- intersect(setdiff(allgene,MEs_tan$gene),top50)
  compare <-
    matrix(c(
      length(gene_overlap),
      length(MEs_tan_high$gene),
      length(gene_overlap_other),
      length(setdiff(allgene,MEs_tan_high$gene))),
      nrow = 2)
  result <- fisher.test(compare ,conf.level = 0.95,alternative = "g")
  result_tan <- cbind(length(gene_overlap),length(MEs_tan_high$gene),length(gene_overlap_other),
                      length(setdiff(allgene,MEs_tan_high$gene)),result$p.value,result$estimate)
  
  MEs_turquoise <- read.table("MEs_turquoise.txt")
  colnames(MEs_turquoise) <- c("gene")
  MEs_turquoise_high <- read.table("MEs_turquoise_high.txt")
  colnames(MEs_turquoise_high) <- c("gene")
  intersect(MEs_turquoise$gene,top50)
  gene_overlap <- intersect(MEs_turquoise_high$gene,top50)
  gene_overlap_other <- intersect(setdiff(allgene,MEs_turquoise$gene),top50)
  compare <-
    matrix(c(
      length(gene_overlap),
      length(MEs_turquoise_high$gene),
      length(gene_overlap_other),
      length(setdiff(allgene,MEs_turquoise_high$gene))),
      nrow = 2)
  result <- fisher.test(compare ,conf.level = 0.95,alternative = "g")
  result_turquoise <- cbind(length(gene_overlap),length(MEs_turquoise_high$gene),length(gene_overlap_other),
                            length(setdiff(allgene,MEs_turquoise_high$gene)),result$p.value,result$estimate)
  
  MEs_yellow <- read.table("MEs_yellow.txt")
  colnames(MEs_yellow) <- c("gene")
  MEs_yellow_high <- read.table("MEs_yellow_high.txt")
  colnames(MEs_yellow_high) <- c("gene")
  intersect(MEs_yellow$gene,top50)
  gene_overlap <- intersect(MEs_yellow_high$gene,top50)
  gene_overlap_other <- intersect(setdiff(allgene,MEs_yellow$gene),top50)
  compare <-
    matrix(c(
      length(gene_overlap),
      length(MEs_yellow_high$gene),
      length(gene_overlap_other),
      length(setdiff(allgene,MEs_yellow_high$gene))),
      nrow = 2)
  result <- fisher.test(compare ,conf.level = 0.95,alternative = "g")
  result_yellow <- cbind(length(gene_overlap),length(MEs_yellow_high$gene),length(gene_overlap_other),
                         length(setdiff(allgene,MEs_yellow_high$gene)),result$p.value,result$estimate)
}
result_HM <-
  rbind(
    result_black,
    result_blue,
    result_brown,
    result_green,
    result_greenyellow,
    result_magenta,
    result_pink,
    result_purple,
    result_red,
    result_tan,
    result_turquoise,
    result_yellow
  )
result_HM <- as.data.frame(result_HM)
colnames(result_HM) <- c("c1", "c2", "c3", "c4", "p", "or")
rownames(result_HM) <-
  c(
    "black",
    "blue",
    "brown",
    "green",
    "greenyellow",
    "magenta",
    "pink",
    "purple",
    "red",
    "tan",
    "turquoise",
    "yellow"
  )
result_HM$FDR <- p.adjust(result_HM$p, method = "fdr")
result_HM[3,]
