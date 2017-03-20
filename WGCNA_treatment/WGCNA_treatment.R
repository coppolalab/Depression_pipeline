#For WGCNA
library(WGCNA)
library(flashClust)
enableWGCNAThreads()

#For baseline processing
library(limma)
library(sva)
library(R.utils)
library(Biobase)

#Functional programming
library(magrittr)
library(purrr)

#Data arrangement
library(dplyr)
library(tidyr)
library(broom)

#String operations
library(stringr)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)
library(igraph)
library(TeachingDemos)

#Reading and writing tables
library(readr)
library(openxlsx)
library(parallel)

source('../../FRDA project/common_functions.R')
source('../../code/GO/enrichr.R')

gen.pcaplot <- function(filename, dataset, facet.bool, size.height, size.width)
{
    colnames(dataset)[2] <- "Module"
    dataset$Module %<>% str_replace("ME", "") 
    p <- ggplot(dataset, aes(x = as.numeric(x), y = value, fill = Module, color = Module)) + geom_point()
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank()) + ylab("First Principal Component")
    p <- p + theme(axis.title.x = element_blank()) + scale_x_continuous(as.numeric(unique(dataset$x)))
    p <- p + scale_color_manual(values = sort(unique(dataset$Module)))
    if (facet.bool == TRUE)
    {
        p <- p + facet_wrap(~ Module)
        p <- p + theme(legend.position = "none")
    } 
    CairoPDF(filename, height = size.height, width = size.width)
    print(p)
    dev.off()
}

gen.heatmap <- function(dataset, ME.genes)
{
    color <- as.character(unique(dataset$module.colors))
    dataset %<>% select(-module.colors) %>% scale
    max.dataset <- max(abs(dataset))
    print(dim(dataset))
    CairoPDF(paste("./modules/", color, sep = ""), width = 21, height = 12)
    par(mar = c(3.5,3,2,3))
    par(oma = c(4,0,2,0))
    plotMat(dataset, zlim = c(-max.dataset, max.dataset), main = paste(color, " (", nrow(dataset), ")", sep = ""))

    ME.genes.plot <- select(ME.genes, Sample.ID, matches(color))
    p <- ggplot(ME.genes.plot, aes_string(x = "Sample.ID", y = color))
    p <- p + geom_bar(stat = "identity") + xlab("Eigengene Expression")#+ ylim(c(-6, 16)) 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.text.x = element_text(angle = 90, size = 2))  
    p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
    print(p)
    dev.off()
}

enrichr.submit <- function(index, full.df, enrichr.terms, use.weights)
{
    dataset <- filter(full.df, module.colors == index)
    dir.create(file.path("./enrichr", index), showWarnings = FALSE, recursive = TRUE)
    enrichr.data <- map(enrichr.terms, get.enrichrdata, dataset, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]
    names(enrichr.data) <- enrichr.names
    trap1 <- map(names(enrichr.data), enrichr.wkbk, enrichr.data, index)
}

enrichr.wkbk <- function(subindex, full.df, index)
{
    dataset <- full.df[[subindex]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)

    filename = paste("./enrichr/", index, "/", index, "_", subindex, ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

plot.eigencor <- function(module.traits.pval, status.col, status.vector)
{
    sig.col <- paste(status.col, ".p.value", sep = "")
    cor.status.labeled <- data.frame(Color = rownames(module.traits.pval), select_(data.frame(module.traits.pval), sig.col))
    filter.cond <- paste(sig.col, "< 0.05")
    status.sig <- filter_(cor.status.labeled, filter.cond)
    me.genes.status <- select(ME.genes, one_of(as.character(status.sig$Color)))
    me.genes.status$Status <- status.vector
    split.cols <- str_split(status.col, "\\.")[[1]]
    me.status.melt <- melt(me.genes.status, id.vars = "Status") %>% filter(Status == split.cols[1] | Status == split.cols[2])
    colnames(me.status.melt)[2] <- "Module"

    sum.fun <- function(data.vector){ data.frame(ymin = min(data.vector), ymax = max(data.vector), y = mean(data.vector)) }
    me.status.melt$Module %<>% as.character
    p <- ggplot(me.status.melt, aes(x = factor(Status), y = value, col = Module)) + geom_point(position = "jitter")
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank()) + ylab("Eigengene")
    p <- p + theme(axis.title.x = element_blank()) + stat_summary(aes(group = 1), fun.y = mean, geom = "line", col = "black", position = position_dodge(width = 0.9))
    p <- p + scale_color_manual(values = sort(unique(me.status.melt$Module)))
    p <- p + facet_wrap(~ Module, scales = "free_y")
    p <- p + theme(legend.position = "none")

    filename <- paste(split.cols[1], "_", split.cols[2], "_eigengenes_05", sep = "")
    CairoPDF(filename, height = 13, width = 20)
    print(p)
    dev.off()
}

gen.enrichrplot <- function(enrichr.df, filename, plot.height = 5, plot.width = 8)
{
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map_int(Compose(unlist, length))
    enrichr.df$Log.pvalue <- -(log10(enrichr.df$P.value))
    enrichr.df$Term %<>% str_replace_all("\\ \\(.*$", "") %>% str_replace_all("\\_Homo.*$", "") %>% tolower #Remove any thing after the left parenthesis and convert to all lower case
    enrichr.df$Format.Name <- paste(enrichr.df$Database, ": ", enrichr.df$Term, " (", enrichr.df$Gene.Count, ")", sep = "")
    enrichr.df %<>% arrange(Log.pvalue)
    enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)
    enrichr.df.plot <- select(enrichr.df, Format.Name, Log.pvalue) %>% melt(id.vars = "Format.Name") 

    p <- ggplot(enrichr.df.plot, aes(Format.Name, value, fill = variable)) + geom_bar(stat = "identity") + geom_text(label = enrichr.df$Format.Name, hjust = "left", aes(y = 0.1)) + coord_flip() + theme_bw() + theme(legend.position = "none") 
    p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste('-', Log[10], ' P-value')))
    CairoPDF(filename, height = plot.height, width = plot.width)
    print(p)
    dev.off()
}

module.iterate <- function(color, modules.out) {
    module.overlaps <- map(gene.lists, get.intersect, color, modules.out) %>% reduce(rbind) %>% data.frame
    colnames(module.overlaps) <- c(str_c("Overlap.", color), str_c("P.value.", color))
    return(module.overlaps)
}

get.intersect <- function(gene.list, color, modules.out) {
    gene.vector <- gene.list[[1]]
    shared.symbols <- intersect(modules.out$Symbol, gene.vector)
    module.object <- filter(modules.out, Module == color)    
    overlap <- intersect(shared.symbols, module.object$Symbol)
    overlap.pval <- phyper(length(overlap), nrow(module.object), (nrow(modules.out) - nrow(module.object)), length(shared.symbols), lower.tail = FALSE)
    return(c(length(overlap), overlap.pval))
}

lumi.ect <- readRDS.gz("../treatment/save/ect.collapse.rda")
lumi.ect <- lumi.ect[,lumi.ect$Time.Point != "TP4"]
lumi.ect$Time.Point %<>% droplevels

model.ect <- model.matrix( ~ Age + Gender + Ethnicity + RIN, data = pData(lumi.ect))[,-1]
ect.design <- model.matrix( ~ Time.Point, data = pData(lumi.ect) )
rmcov.ect.expr <- removeBatchEffect(exprs(lumi.ect), covariates = model.ect, design = ect.design)
rmcov.ect <- lumi.ect #Make a copy of lumi object
exprs(rmcov.ect) <- rmcov.ect.expr #Transfer cleaned expression values into new lumi object

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

expr.collapse <- exprs(rmcov.ect) %>% t
sft <- pickSoftThreshold(expr.collapse, powerVector = powers, verbose = 5, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), networkType = "signed")
sft.df <- sft$fitIndices
saveRDS.gz(sft, file = "./save/sft.rda")

#Plot scale indendence and mean connectivity as functions of power
sft.df$multiplied <- sft.df$SFT.R.sq * -sign(sft.df$slope)
p <- ggplot(sft.df, aes(x = Power,  y = multiplied, label = Power)) + geom_point() + geom_text(vjust = -0.6, size = 4, col = "red")
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_hline(aes(yintercept = 0.9))
p <- p + xlab("Soft Threshold") + ylab("Scale Free Topology Model Fit, signed R^2") + ggtitle("Scale Independence")
CairoPDF(file = "./scaleindependence", width = 6, height = 6)
print(p)
dev.off()

p <- ggplot(sft.df, aes(x = Power,  y = mean.k.)) + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Soft Threshold") + ylab("Mean Connectivity") + ggtitle("Mean Connectivity")
CairoPDF(file = "./meanconnectivity", width = 6, height = 6)
print(p)
dev.off()

softPower <- 6
adjacency.expr <- adjacency(expr.collapse, power = softPower, type = "signed", corFnc = "bicor", corOptions = "maxPOutliers = 0.05")

TOM <- TOMsimilarity(adjacency.expr, verbose = 5)
dissimilarity.TOM <- 1 - TOM

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")
saveRDS.gz(geneTree, file = "./save/gene.tree.rda")

CairoPDF(file = "./genecluster", height = 10, width = 15)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

min.module.size <- 20

#Identify modules using dynamic tree cutting with hybrid clustering
dynamic.modules <- cutreeDynamic(dendro = geneTree, method = "hybrid", distM = dissimilarity.TOM, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
dynamic.colors <- labels2colors(dynamic.modules)
saveRDS.gz(dynamic.colors, file = "./save/dynamic.colors.rda")

CairoPDF(file = "./gene_dendrogram_and_module_colors", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Calculate module eigengenes
ME.list <- moduleEigengenes(expr.collapse, colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - bicor(ME.genes, maxPOutliers = 0.05)
METree <- flashClust(as.dist(MEDiss), method = "average")
saveRDS.gz(METree, file = "./save/me.tree.rda")

CairoPDF(file = "./module_eigengene_clustering_tree", height = 10, width = 15)
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.2
merge.all <- mergeCloseModules(expr.collapse, dynamic.colors, cutHeight = ME.dissimilarity.threshold, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), verbose = 3) 
merged.colors <- merge.all$colors
merged.genes <- merge.all$newMEs

CairoPDF("module_eigengene_clustering", height = 10, width = 15)
plotDendroAndColors(geneTree, cbind(dynamic.colors, merged.colors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Use merged eigengenes 
module.colors <- merged.colors
saveRDS.gz(module.colors, file = "./save/module.colors.rda")
color.order <- c("grey", standardColors(50))
modules.labels <- match(module.colors, color.order)
saveRDS.gz(modules.labels, file = "./save/modules.labels.rda")
ME.genes <- merged.genes
saveRDS.gz(ME.genes, file = "./save/me.genes.rda")

CairoPDF("eigengenes", height = 6, width = 8)
plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), marHeatmap = c(3,5,1,2), plotPreservation = "standard")
dev.off()

modules.out <- data.frame(Symbol = colnames(expr.collapse), Module = module.colors)
write.xlsx(modules.out, "modules_out.xlsx")

EigengeneModel <- function(ME.vector, trait.vector) {
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.aov <- aov(ME ~ Trait, trait.df) %>% TukeyHSD
    trait.tukey <- data.frame(trait.aov$Trait)
    colnames(trait.tukey) %<>% str_replace(" ", ".")
    return(select(trait.tukey, diff, p.adj))
}

EigengeneAOV <- function(ME.vector, pheno.data) {
    trait.df <- data.frame(pheno.data, ME = ME.vector)
    trait.tukey <- aov(ME ~ Time.Point * responder_to_treatment, trait.df) %>% TukeyHSD #%>% tidy
    trait.tukey
}

EigengeneAnova <- function(ME.vector, pheno.data) {
    trait.df <- data.frame(pheno.data, ME = ME.vector)
    trait.anova <- lm(ME ~ Ethnicity + Gender + Age + RIN + Time.Point + responder_to_treatment + Time.Point * responder_to_treatment, trait.df) %>% anova %>% tidy
    trait.anova
}

EigengeneBayes <- function(ME.vector, pheno.data) {
    trait.df <- data.frame(pheno.data, ME = ME.vector)
    trait.anova <- lmBF(ME ~ Ethnicity + Gender + Age + RIN + Time.Point + responder_to_treatment + Time.Point * responder_to_treatment, trait.df) 
    trait.anova
}

pdata <- pData(rmcov.ect)
#targets.final.known$ %<>% factor(levels = c("CONT", "ECT"))
anova.time.response <- map(ME.genes, EigengeneAnova, pdata) #%>% map(slice, 7) %>% reduce(rbind)#%>% map(extract, 5)
bayes.time.response <- map(ME.genes, EigengeneBayes, pdata) 
bf.time.response <- map(bayes.time.response, extractBF)
posterior.bf <- map(bayes.time.response, posterior, iterations = 10000) 
posterior.quantiles <- map(posterior.bf, summary) %>% map(extract2, "quantiles") %>% map(data.frame) %>% map(select, X2.5., X97.5.)

cor.time.response <- map(ME.genes, EigengeneAOV, pdata)
status.diff <- map(cor.status, select, diff) %>% map(t) %>% reduce(rbind) %>% data.frame
rownames(status.diff) <- names(cor.status)
status.pval <- map(cor.status, select, p.adj) %>% map(t) %>% reduce(rbind) %>% data.frame
rownames(status.pval) <- names(cor.status)

pval.adjust <- map(status.pval, p.adjust, method = "fdr", n = nrow(status.pval)) %>% reduce(cbind) %>% data.frame
rownames(pval.adjust) <- rownames(status.pval)
colnames(pval.adjust) <- paste(colnames(status.pval), ".pval")

#status.pval$Module <- rownames(status.pval)
text.matrix.traits <- paste(signif(as.matrix(status.diff), 2), '\n(', signif(as.matrix(pval.adjust), 1), ')', sep = '')
dim(text.matrix.traits) = dim(status.diff)

heatmap.range <- c(min(as.matrix(status.diff)) * 1.1, max(as.matrix(status.diff)) * 1.1)
gen.text.heatmap(as.matrix(status.diff), text.matrix.traits, colnames(status.diff), colnames(ME.genes), "", "module-trait relationships", heatmap.range)

all.degrees <- intramodularConnectivity(adjacency.expr, module.colors)
gene.info <- data.frame(Symbol = rownames(all.degrees), module.color = module.colors, all.degrees)

write_csv(data.frame(table(module.colors)), path = "./final_eigengenes.csv") 
gene.info$kscaled <- by(gene.info, gene.info$module.color, select, kWithin) %>% llply(function(x) { x / max (x) }) %>% reduce(c)
saveRDS.gz(gene.info, file = "./save/gene.info.rda")

gene.module.membership <- as.data.frame(bicor(expr.collapse, ME.genes, maxPOutliers = 0.05))
module.membership.pvalue <- as.data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(expr.collapse)))
names(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")

module.membership <- cbind(select(gene.info, Symbol:module.color), gene.module.membership, module.membership.pvalue)
write_csv(module.membership, "module_membership.csv")

colnames(gene.module.membership) %<>% str_replace("MM.", "")
colnames(module.membership.pvalue) %<>% str_replace("MM.pvalue.", "")
gene.module.membership$Symbol <- rownames(gene.module.membership)
module.membership.pvalue$Symbol <- rownames(module.membership.pvalue)

color.values <- unique(module.colors)
color.key <- paste(head(color.values, n = 1), tail(color.values, n = 1), sep = ":")
gene.module.membership.long <- gather_(data = gene.module.membership, "module.comparison", "correlation", color.values)
module.membership.pvalue.long <- gather_(data = module.membership.pvalue, "module.comparison", "p.value", color.values)
membership.join <- join(gene.module.membership.long, module.membership.pvalue.long)
eigengene.connectivity <- join(membership.join, gene.info) %>% select(Symbol, module.color:kscaled, module.comparison:p.value)
write_csv(eigengene.connectivity, "eigengene_connectivity.csv")

all.smooth <- apply(ME.genes, 2, smooth.spline, spar = 0.4) %>% llply(`[`, "y")
smooth.df <- data.frame(all.smooth)
colnames(smooth.df) <- names(all.smooth)
smooth.df$x <- as.factor(1:nrow(smooth.df))
smooth.plot <- melt(smooth.df, id.vars = "x")

gen.pcaplot("all_principal_components", smooth.plot, FALSE, 10, 15)
gen.pcaplot("facet_principal_components", smooth.plot, TRUE, 13, 25)

sample.ids <- factor(rownames(expr.collapse), levels = rownames(expr.collapse))
colnames(ME.genes) %<>% str_replace("ME", "")
ME.genes.plot <- mutate(data.frame(ME.genes), Sample.ID = sample.ids)
expr.data.plot <- data.frame(t(expr.collapse), module.colors)

cluster <- makeForkCluster(8)
split(expr.data.plot, expr.data.plot$module.colors) %>% map(gen.heatmap, ME.genes.plot)
by(expr.data.plot, expr.data.plot$module.colors, gen.heatmap, ME.genes.plot)

#modules.out <- select(gene.info, Symbol, module.color)
#targets.final.known$Sample.Name %<>% str_replace(" ", "")

enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "WikiPathways_2016", "Reactome_2016", "BioCarta_2016", "PPI_Hub_Proteins", "HumanCyc_2016", "NCI-Nature_2016", "Panther_2016") 
color.names <- unique(module.colors) %>% sort
trap1 <- map(color.names, enrichr.submit, modules.out, enrichr.terms, FALSE)

rownames(ME.genes) <- rownames(expr.collapse)

ME.df <- data.frame(Time.Point = pdata$Time.Point, Responder = pdata$responder_to_treatment, ME.genes)
#p <- ggplot(ME.df, aes(factor(Status), green) ) + geom_boxplot(col = "green")
#p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank(), legend.position = "none")
#p <- p + ylab("Eigengene") + ggtitle(paste("P < ", anova.status["MEgreen"]))
#CairoPDF("green.eigengene", width = 4, height = 4)
#print(p)
#dev.off()

#p <- ggplot(ME.df, aes(factor(Status), pink) ) + geom_boxplot(col = "pink" )
#p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank(), legend.position = "none")
#p <- p + ylab("Eigengene") + ggtitle(paste("P < ", anova.status["MEpink"]))
#CairoPDF("pink.eigengene", width = 4, height = 4)
#print(p)
#dev.off()

p <- ggplot(ME.df, aes(Responder, MEgreen, fill = Responder) ) + geom_violin() + facet_wrap( ~ Time.Point, ncol = 4 )
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank(), legend.position = "none")
p <- p + ylab("Eigengene") 
CairoPDF("green.eigengene", width = 12, height = 4)
print(p)
dev.off()

EstimatePlot <- function(eigengene, posterior.list){
    time.response <- select(data.frame(posterior.list[[eigengene]]), matches("^Time.Point.responder.*")) %>% gather(Group, Estimate)
    time.response$Time.Point <- str_split_fixed(time.response$Group, "\\.", 6)[,4]
    time.response$Responder <- str_split_fixed(time.response$Group, "\\.", 6)[,6] %>% str_replace("\\.", "")

    p <- ggplot(time.response, aes(Responder, Estimate, fill = Responder) ) + geom_violin() + facet_wrap( ~ Time.Point, ncol = 4 )
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank(), legend.position = "none")
    p <- p + ylab("Estimate") + ylim(c(-0.1, 0.1))
    CairoPDF(str_c(eigengene, ".estimate"), width = 12, height = 4)
    print(p)
    dev.off()
}

map(colnames(ME.genes), EstimatePlot, posterior.bf)

#plot.eigencor(module.traits.pval, "Control.Carrier", lumi.import$Status)
#plot.eigencor(module.traits.pval, "Control.Patient", lumi.import$Status)
#plot.eigencor(module.traits.pval, "Carrier.Patient", lumi.import$Status)

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

#Overlap
list.names <- list.files("../MDD_ReferenceGeneSets/GENE_LISTS/", full.names = TRUE)
gene.lists <- map(list.names, read_delim, delim = "\n")

treatment.overlaps.list <- map(color.values, module.iterate, modules.out) %>% reduce(cbind) %>% data.frame
treatment.overlaps <- data.frame(Gene.List = map_chr(gene.lists, names), Num.Genes = map_int(gene.lists, nrow), treatment.overlaps.list)
write.xlsx(treatment.overlaps, "./treatment.overlaps.xlsx")

treatment.pvals <- select(treatment.overlaps.list, contains("P.value")) %>% as.matrix
treatment.logpvals <- -log10(treatment.pvals)
treatment.count <- select(treatment.overlaps.list, contains("Overlap")) %>% as.matrix
text.matrix.overlaps <- str_c(treatment.count, '\n(', signif(treatment.pvals, 2), ')')
dim(text.matrix.overlaps) <- dim(treatment.count)
heatmap.labels <- str_c(treatment.overlaps$Gene.List, " (", treatment.overlaps$Num.Genes, ")")

reduce.key <- apply(select(reduce.module.overlaps, contains("Overlap")), 1, `>=`, e2 = 10) %>% apply(2, any)
reduce.overlaps <- reduce.module.overlaps[reduce.key,]
heatmap.range <- c(0, 2)
gen.text.heatmap(treatment.logpvals, text.matrix.overlaps, color.values, heatmap.labels, "", "treatment_genelists", heatmap.range, colorRampPalette(c("red", "green", "white"))(n = 300))

#reduce.overlaps <- select(treatment.overlaps, Gene.List, Num.Genes, matches("\\.yellow"), matches("\\.green$")) %>% filter(Num.Genes >= 50)
#reduce.pvals <- select(reduce.overlaps, contains("P.value")) 
#reduce.pvals.adj <- apply(reduce.pvals, 2, p.adjust, method = "fdr", n = nrow(reduce.pvals) * ncol(reduce.pvals)) %>% as.matrix
#reduce.logpvals <- -log10(reduce.pvals.adj)
#reduce.count <- select(reduce.overlaps, contains("Overlap")) %>% as.matrix
#text.matrix.reduce <- str_c(reduce.count, '\n(', signif(reduce.pvals.adj, 2), ')')
#dim(text.matrix.reduce) <- dim(reduce.count)
#reduce.labels <- str_c(reduce.overlaps$Gene.List, " (", reduce.overlaps$Num.Genes, ")")
#gen.text.heatmap(reduce.logpvals, text.matrix.reduce, c("yellow", "green"), reduce.labels, "", "treatment_genelists_reduce", heatmap.range, colorRampPalette(c("red", "yellow", "white"))(n = 300))

gen.text.heatmap <- function(cor.dataset, text.matrix, x.names, y.names, maintitle, filename, zlim.plot = c(-1,1), color.scheme = greenWhiteRed(50))
{
    width.dynamic <- 3 + ncol(text.matrix)
    height.dynamic <- 3 + nrow(text.matrix)
    CairoPDF(filename, width = width.dynamic, height = 10)
    par(mar = c(8, 16, 3, 3))
    labeledHeatmap(Matrix = cor.dataset, xLabels = x.names, yLabels = y.names, ySymbols = y.names, yColorLabels = TRUE, colors = color.scheme, textMatrix = text.matrix, setStdMargins = F, cex.text = 0.5, zlim = zlim.plot, main = maintitle)
    dev.off()
}



green.only <- filter(modules.out, Module == "green")
pink.only <- filter(modules.out, Module == "pink")

#Final plots
green.gobiol.file <- "./enrichr/green/green_GO_Biological_Process_2015.xlsx"
green.gobiol <- read.xlsx(green.gobiol.file) 
green.gobiol$Num.Genes <- map(green.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
green.gobiol %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(green.gobiol, green.only$Symbol, file_path_sans_ext(green.gobiol.file))
green.gobiol$Database <- "GO Biological Process"
green.gobiol.final <- slice(green.gobiol, c(1, 25, 35))

green.molec <- read.xlsx("./enrichr/green/green_GO_Molecular_Function_2015.xlsx") %>% slice(14)
green.molec$Database <- "GO Molecular Function"
green.molec$Num.Genes <- map(green.molec$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)

green.reactome <- read.xlsx("./enrichr/green/green_Reactome_2016.xlsx") %>% slice(2)
green.reactome$Database <- "Reactome"
green.reactome$Num.Genes <- map(green.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)

green.enrichr <- rbind(green.gobiol.final, green.molec, green.reactome)
gen.enrichrplot(green.enrichr, "green.enrichr")

pink.gobiol.file <- "./enrichr/pink/pink_GO_Biological_Process_2015.xlsx"
pink.gobiol <- read.xlsx(pink.gobiol.file) 
pink.gobiol$Num.Genes <- map(pink.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
pink.gobiol %<>% filter(Num.Genes > 4) %>% filter(P.value < 0.01)
get.kappa.cluster(pink.gobiol, pink.only$Symbol, file_path_sans_ext(pink.gobiol.file))
pink.gobiol$Database <- "GO Biological Process"
pink.gobiol.final <- slice(pink.gobiol, c(4, 53, 41, 6))

pink.reactome <- read.xlsx("./enrichr/pink/pink_Reactome_2016.xlsx") %>% slice(1)
pink.reactome$Database <- "Reactome"
pink.reactome$Num.Genes <- map(pink.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)

pink.enrichr <- rbind(pink.gobiol.final, pink.reactome)

gen.enrichrplot(pink.gobiol.final, "pink.enrichr")

#PPI
biogrid.ppi <- read_tsv("~/Downloads/BIOGRID-ORGANISM-3.4.136.tab2/BIOGRID-ORGANISM-Homo_sapiens-3.4.136.tab2.txt") %>% data.frame
biogrid.ppi.reduce <- select(biogrid.ppi, contains("Official"))

inweb.ppi <- read_tsv("../../Dementia Project/mapt/InWeb3_HC_NonRed.txt", col_names = FALSE)
colnames(inweb.ppi) <- c("Interactor.A", "Interactor.B")
iw.hugo <- read_tsv("../../Dementia Project/mapt/IWtoHugo.txt", col_names = FALSE)
colnames(iw.hugo) <- c("IW.ID", "HUGO")

green.gobiol.resp <- str_split(green.gobiol.final[2,]$Genes, ",")[[1]]
pink.gobiol.toll <- str_split(pink.gobiol.final[4,]$Genes, ",")[[1]]
pink.gobiol.apop <- str_split(pink.gobiol.final[2,]$Genes, ",")[[1]]

ggr.ppi <- get.ppi(green.gobiol.resp)
pgt.ppi <- get.ppi(pink.gobiol.toll)
pga.ppi <- get.ppi(pink.gobiol.apop)

plot.ppi(adjacency.expr, green.gobiol.resp, ggr.ppi, "green_igraph", TRUE)
plot.ppi(adjacency.expr, pink.gobiol.apop, pga.ppi, "pga_igraph", TRUE)
#plot.ppi(adjacency.expr, pink.gobiol.toll, pgt.unique.final, "pgt_igraph", FALSE)

