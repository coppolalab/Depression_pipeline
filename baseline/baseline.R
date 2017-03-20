#Specific packages
library(R.utils) #for capitalize
library(tools) #for file_path_sans_ext
library(limma) #for DE 
library(lumiHumanIDMapping) #for converting Illumina probe IDs to nuIDs
library(lumiHumanAll.db) #for getting symbols for nuIDs
library(annotate) #for getting symbols for nuIDs
library(lumi) #main preprocessing package
library(WGCNA) #for fastcor and connecitivity stats
library(biomaRt) #for extra annotation of genes WARNING: BROKEN
library(irr) #for kappa clustering
library(sva) #for ComBat

#Plotting
library(ggplot2) #main plotting package
library(extrafont) #access extra fonts (redundant w/ fontconfig?)
library(Cairo) #Cairo graphics device
library(heatmap.plus) #main heatmap
library(gplots) #for heatmap.2
library(RColorBrewer) #for extra colors (unneeded?)

#Input/output manipulation
library(readr) #for reading csv, tsv
library(openxlsx) #best pacakage for reading writing xlsx
library(stringr) #much better string processing
library(broom) #for cleaning up output of lm

#Language enhancements
library(reshape2) #for reshaping data (should use tidyr)
library(plyr) #for reshaping data (should use dplyr)
library(magrittr) #for %>% operator
library(purrr) #makes R into Haskell ;)
library(dplyr) #for manipulating data
library(functional) #makes R into Haskell ;)
library(vadr) #makes R into Haskell ;)

source('../../FRDA project/common_functions.R')
source('../../code/GO/enrichr.R')

match.exact <- mkchain(map_chr(paste %<<<% "^" %<<% c("$", sep = "")), paste(collapse = "|"))

#Boxplot function
gen.boxplot <- function(filename, lumi.object, maintext, ylabtext)
{
    expr.df <- exprs(lumi.object) %>% t %>% data.frame
    dataset.addvars <- mutate(expr.df, Sample.Status = sampleNames(lumi.object), Batch = lumi.object$Batch)
    dataset.m <- melt(dataset.addvars, id = c("Sample.Status", "Batch"))
    p <- ggplot(dataset.m, aes(x = Sample.Status, y = value, fill = factor(Batch))) + geom_boxplot() + theme_bw()
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.0, size = 5))     
    p <- p + ggtitle(maintext) + ylab(ylabtext) + xlab("Sample") + theme(axis.text.x = element_text(size = 3))
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    ggsave(filename = filename, plot = p, family = "Oxygen", width = 20 , height = 8)
}

#Heatmap bar function
gen.heatmapbars <- function(batch.colorscheme, diagnosis.colorscheme, targetset)
{
    batch.heatmap <- data.frame("Batch" = seq(1:length(batch.colorscheme)), "Batch.Color" = batch.colorscheme)
    diagnosis.heatmap <- data.frame("Status" = levels(targetset$Status), "Diagnosis.Color" = diagnosis.colors)
    colorscheme <- data.frame("Batch" = targetset$Batch, "Status" = targetset$Status) %>% join(batch.heatmap) %>% join(diagnosis.heatmap)
    colorscheme <- as.matrix(subset(colorscheme, select = c("Batch.Color", "Diagnosis.Color")))
    return(colorscheme)
}

#Heatmap function
gen.heatmap <- function(filename, lumi.object, maintitle)
{
    intensities1.cor <- corFast(exprs(lumi.object))
    CairoPDF(filename, width = 10, height = 10)
    heatmap.plus(intensities1.cor, col = heat.colors(40), ColSideColors = cbind(lumi.object$Batch.Color, lumi.object$Diagnosis.Color), scale = "none", cexCol = 0.17, cexRow = 0.17, main = maintitle)
    dev.off()
}

gen.histogram <- function(filename, lumi.object)
{
    expr.df <- exprs(lumi.object) %>% t %>% data.frame
    dataset.addvars <- mutate(expr.df, Sample.Status = sampleNames(lumi.object), Group = lumi.object$Group)
    dataset.m <- melt(dataset.addvars, id = c("Sample.Status", "Group"))

    p <- ggplot(dataset.m, aes(value, group = Sample.Status, col = factor(Group))) + geom_density() + theme_bw()
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + ggtitle("Histogram of VST Expression") + ylab("Density") + xlab("VST Expression") 
    CairoPDF(filename, height = 5, width = 9)
    print(p)
    dev.off()
}

gen.mdsplot <- function(filename, lumi.object, variable)
{
    mds.vst <- exprs(lumi.object) %>% t %>% dist %>% cmdscale(eig = TRUE) 
    mds.vst.plot <- data.frame(Sample.Name = rownames(mds.vst$points), Group = lumi.object$Group, Batch = factor(lumi.object$Batch), Component.1 = mds.vst$points[,1], Component.2 = mds.vst$points[,2])

    p <- ggplot(mds.vst.plot, aes_string(x = "Component.1", y = "Component.2", col = variable)) + geom_point() 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + xlab("Component 1") + ylab("Component 2") + ggtitle("MDS of Group")
    CairoPDF(file = filename, height = 6, width = 7)
    print(p)
    dev.off()
}

gen.connectivityplot <- function(filename, dataset, maintitle) {
    norm.adj <- (0.5 + 0.5 * bicor(exprs(dataset)))
    colnames(norm.adj) <- dataset$External.ID
    rownames(norm.adj) <- dataset$External.ID
    net.summary <- fundamentalNetworkConcepts(norm.adj)
    net.connectivity <- net.summary$Connectivity
    connectivity.zscore <- (net.connectivity - mean(net.connectivity)) / sqrt(var(net.connectivity))

    connectivity.plot <- data.frame(Sample.Name = names(connectivity.zscore), Z.score = connectivity.zscore, Sample.Num = 1:length(connectivity.zscore))
    p <- ggplot(connectivity.plot, aes(x = Sample.Num, y = Z.score, label = Sample.Name) )
    p <- p + geom_text(size = 4, colour = "red")
    p <- p + geom_hline(aes(yintercept = -2)) + geom_hline(yintercept = -3) 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle(maintitle)
    CairoPDF(filename, width = 10, height = 10)
    print(p)
    dev.off()
    return(connectivity.zscore)
}

gen.decide <- function(test, fit.object, write.results) {
    results <- decideTests(fit.object, adjust.method = test[1], p = as.numeric(test[2]))
    if(write.results == TRUE)
    {
        write.fit(file = paste("./fit_", test[1], ".tsv", sep = ""), fit.object, adjust = test[1], results = results)
    }
    num.genes <- map_lgl(results, any, na.rm = T) %>% which %>% length  #Make this better
    mysum <- summary(factor(results))[-2] #Eliminate the row for no change in expression
    mysum[1] <- -(mysum[1])
    mysum.return <- data.frame(Test = paste(test[1], " p<", test[2], sep = ""), Num = paste(num.genes, "Genes", sep = " "), Direction = c("negative", "positive"), mysum) %>% data.frame
    return(mysum.return)
}

gen.ratios <- function(lumi.object, time.point) {
    CONT <- exprs(lumi.object[,lumi.object$Time.Point == time.point & lumi.object$Group == "CONT"])
    CONT.means <- rowMeans(CONT)
    ECT <- exprs(lumi.object[,lumi.object$Time.Point == time.point & lumi.object$Group == "ECT"])
    Diff <- ECT - CONT.means

    all.coefficients <- data.frame("Symbol" = rownames(Diff), Diff)
    all.samples <- data.frame("Symbol" = featureNames(lumi.object), exprs(lumi.object))
    colnames(all.samples)[2:length(all.samples)] %<>% str_c("expr", sep = ".") 
    ratio.exp <- merge(all.coefficients, all.samples)
    return(ratio.exp)
}

gen.workbook <- function(dataset, filename) {
    pval.cols <- colnames(dataset) %>% str_detect("p.value") %>% which
    coef.cols <- colnames(dataset) %>% str_detect("Coef") %>% which
    #colnames(dataset)[coef.cols] %<>% str_replace("Coef", "") %>% str_replace_all("_", " ")
    #dataset$Definition %<>% str_replace_all("Homo sapiens ", "") %>% str_replace_all("PREDICTED: ", "")

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = FALSE)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(dataset), rule = "<0.05", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = coef.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1, widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    setColWidths(wb, 1, cols = 3:8, widths = 10)
    setColWidths(wb, 1, cols = 9:14, widths = 15)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")

    saveWorkbook(wb, filename, overwrite = TRUE) 
}

#Create genelists
gen.tables <- function(dataset, lumi.object, ratio.exp, suffix) {
    treat.de <- data.frame("Symbol" = rownames(dataset), dataset)
    #colnames(treat.de)[str_detect(colnames(treat.de), "Genes")] %<>% str_replace("Genes\\.", "") %>% tolower %>% capitalize

    fitsel.ratio.all <- merge(treat.de, ratio.exp)

    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    bm.table <- getBM(attributes = c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', values = as.character(fitsel.ratio.all$Symbol), mart = ensembl)
    bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
    colnames(bm.table) <- c("Symbol", "Definition")
    fitsel.ratio.all <- join(fitsel.ratio.all, bm.table)

    fitsel.return.all <- select(fitsel.ratio.all, Symbol, Definition, dplyr::contains("Coef"), dplyr::contains("p.value"), F, F.p.value, dplyr::contains("ECT"), matches("^t"), A, dplyr::contains("DE")) %>% arrange(desc(F))

    anovalist <- map_lgl(fitsel.return.all$ECT_vs_CONT, any, na.rm = T) %>% which
    fitsel.return <- fitsel.return.all[anovalist,]

    gen.workbook(fitsel.return, paste("./significant_geneList_", suffix, ".xlsx", sep = ""))

    write.csv(fitsel.return.all, paste("./complete_genelist_", suffix, ".csv", sep = ""), row.names = FALSE)
    return(fitsel.return)
}

enrichr.submit <- function(dataset, enrichr.terms, subdir)
{
    dir.create(file.path("./enrichr", subdir), showWarnings = TRUE, recursive = TRUE)
    enrichr.data <- map(enrichr.terms, get.enrichrdata, dataset, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]

    names(enrichr.data) <- enrichr.names

    map(names(enrichr.data), enrichr.wkbk, enrichr.data, subdir)
    return(enrichr.data)
}

enrichr.wkbk <- function(database, full.df, subdir)
{
    dataset <- full.df[[database]]

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    
    dir.create(file.path("./enrichr", subdir), recursive = TRUE)
    filename = paste(file.path("./enrichr", subdir, database), ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

gen.enrichrplot <- function(enrichr.df, enrichr.expr, filename, plot.height = 5, plot.width = 8)
{
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map_int(Compose(unlist, length))
    enrichr.df$Log.pvalue <- -(log10(enrichr.df$P.value))
    enrichr.updown <- map(enrichr.df$Genes, get.updown, enrichr.expr) %>% reduce(rbind)
    colnames(enrichr.updown) <- c("Up", "Down")
    enrichr.df <- cbind(enrichr.df, enrichr.updown)
    enrichr.df$Log.Up <- enrichr.df$Log.pvalue * enrichr.df$Up / enrichr.df$Gene.Count
    enrichr.df$Log.Down <- enrichr.df$Log.pvalue * enrichr.df$Down / enrichr.df$Gene.Count
    enrichr.df$Term %<>% str_replace_all("\\ \\(.*$", "") %>% str_replace_all("\\_Homo.*$", "") %>% tolower #Remove any thing after the left parenthesis and convert to all lower case
    enrichr.df$Format.Name <- paste(enrichr.df$Database, ": ", enrichr.df$Term, " (", enrichr.df$Gene.Count, ")", sep = "")
    enrichr.df %<>% arrange(Log.pvalue)
    enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)
    enrichr.df.plot <- select(enrichr.df, Format.Name, Log.Up, Log.Down) %>% melt(id.vars = "Format.Name") 

    p <- ggplot(enrichr.df.plot, aes(Format.Name, value, fill = variable)) + geom_bar(stat = "identity") + geom_text(label = c(as.character(enrichr.df$Format.Name), rep("", nrow(enrichr.df))), hjust = "left", aes(y = 0.1)) + scale_fill_discrete(name = "Direction", labels = c("Up", "Down"))
    p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste('-', Log[10], ' P-value')))
    CairoPDF(filename, height = plot.height, width = plot.width)
    print(p)
    dev.off()
}

get.updown <- function(filter.vector, enrichr.expr)
{
    split.vector <- str_split(filter.vector, ",")[[1]]
    grep.vector <- match.exact(split.vector)
    enrichr.filter <- filter(enrichr.expr, grepl(grep.vector, Symbol))
    enrichr.vector <- c("Up" = length(which(sign(enrichr.filter$logFC) == 1)), "Down" = length(which(sign(enrichr.filter$logFC) == -1)))
    return(enrichr.vector)
}

#Read in metadata
sample.key <- read.xlsx("../ECT_GE_Data_042116.xlsx", sheet = 1) #Read first sheet with sample key
sample.key$slideID <- str_c(sample.key$general.array, sample.key$genexstripe.controling.stripe, sep = "_")
pheno.data <- read.xlsx("../ECT_GE_Data_042116.xlsx", sheet = 2) #Read second sheet with phenotype data
#colnames(pheno.data)[12] <- "Treatment"
sample.key.ect <- filter(sample.key, !grepl("KET", External.ID)) %>% filter(!(grepl("Excluded", Subject.ID)))
sample.key.ect$External.ID %<>% str_replace_all("\\-", "_")
sample.key.ect$Subject.ID <- str_sub(sample.key.ect$External.ID, 1, 9)
sample.key.ect[grepl("tube2", sample.key.ect$External.ID),]$Subject.ID %<>% paste("tube2", sep = "_")
pheno.data$Subject.ID %<>% str_replace(" 2nd tube", "_tube2")

pheno.data.join <- join(sample.key.ect, pheno.data)
pheno.tube2 <- filter(pheno.data.join, grepl("tube2", Subject.ID)) %>% getElement("Subject.ID") %>% str_sub(1, 9) 
pheno.data.filter <- filter(pheno.data.join, !is.na(SAMPLE)) %>% filter(!is.element(External.ID, pheno.tube2))
rownames(pheno.data.filter) <- pheno.data.filter$slideID

tsv.files <- list.files("../raw_data/", pattern = "*.txt", full.names = TRUE) %>% map(read_tsv) #Read the tab separated files
csv.files <- list.files("../raw_data/", pattern = "*.csv", full.names = TRUE) %>% map(read_csv) #Read the comma separated files
intensities.mat <- reduce(c(tsv.files, csv.files), merge) #Merge tables because they have different numbers of rows
write.table(intensities.mat, "../raw_data/all_batches.tsv", sep = "\t", row.names = FALSE) #Write merged table back to disk

lumi.raw <- lumiR("../raw_data/all_batches.tsv", lib.mapping = "lumiHumanIDMapping", QC = FALSE) #read in all sample probe profiles
slide.order <- match(pheno.data.filter$slideID, sampleNames(lumi.raw))
lumi.ect <- lumi.raw[,slide.order]
#pheno.data.filter$sampleNames <- pheno.data.filter$External.ID
pData(lumi.ect) <- pheno.data.filter
lumi.ectonly <- lumi.ect[,is.element(lumi.ect$Group, c("ECT", "CONT"))]
lumi.ectonly <- lumi.ectonly[,lumi.ectonly$RIN >= 6]
sampleNames(lumi.ectonly) <- lumi.ectonly$External.ID

lumi.vst <- lumiT(lumi.ectonly)
saveRDS(lumi.vst, "./save/lumi.vst.rda")

lumi.baseline2 <- lumi.vst[,lumi.vst$Time.Point == "TP2"]
gen.boxplot("baseline2_vst.jpg", lumi.baseline2, "VST Normalized Signal Intensity", "Intensity")
gen.histogram("baseline2_histogram", lumi.baseline2)
gen.mdsplot("baseline2_mdsplot_group", lumi.baseline2, "Group")
gen.mdsplot("baseline2_mdsplot_batch", lumi.baseline2, "Batch")

lumi.baseline2.norm <- lumiN(lumi.baseline2, "rsn")
lumi.baseline2.cutoff <- detectionCall(lumi.baseline2.norm)
lumi.baseline2.expr <- lumi.baseline2.norm[which(lumi.baseline2.cutoff > 0),]
symbols.lumi <- getSYMBOL(rownames(lumi.baseline2.expr), "lumiHumanAll.db") %>% is.na
lumi.baseline2.annot <- lumi.baseline2.expr[!symbols.lumi,]

model.combat <- model.matrix(~ factor(Group) + factor(Ethnicity) + factor(Gender) + Age + RIN, data = pData(lumi.baseline2.annot)) %>% data.frame
baseline2.combat.expr <- ComBat(dat = exprs(lumi.baseline2.annot), batch = factor(lumi.baseline2.annot$Batch), mod = model.combat) 
baseline2.combat <- lumi.baseline2.annot 
exprs(baseline2.combat) <- baseline2.combat.expr

gen.mdsplot("combat_mdsplot_group", baseline2.combat, "Group")
gen.mdsplot("combat_mdsplot_batch", baseline2.combat, "Batch")

tree.baseline2 <- exprs(baseline2.combat) %>% t %>% dist %>% hclust(method = "average")
CairoPDF("clustering_baseline2", width = 13, height = 10)
plot(tree.baseline2, main = "Hierarchical Clustering Sammples")
dev.off()

baseline2.connectivity <- gen.connectivityplot("baseline2_connectivity", baseline2.combat, "")
baseline2.outlier <- names(baseline2.connectivity[baseline2.connectivity < -3])
baseline2.names <- names(baseline2.connectivity)
baseline2.reps <- baseline2.names[grepl("rep", baseline2.names)]
baseline2.orig <- str_replace(baseline2.reps, "_rep.*$", "")
baseline2.reps.df <- data.frame(Orig = abs(baseline2.connectivity[match(baseline2.orig, names(baseline2.connectivity))]), Replicate = abs(baseline2.connectivity[match(baseline2.reps, names(baseline2.connectivity))]))
baseline2.whichmax <- apply(baseline2.reps.df, 1, which.max)
baseline2.reps.name <- data.frame(Orig = names(baseline2.connectivity[match(baseline2.orig, names(baseline2.connectivity))]), Replicate = names(baseline2.connectivity[match(baseline2.reps, names(baseline2.connectivity))]))
baseline2.rmreps.key <- baseline2.reps.name[cbind(seq_along(baseline2.whichmax), baseline2.whichmax)]

baseline2.remove <- c(baseline2.outlier, baseline2.rmreps.key)
baseline2.rmreps <- baseline2.combat[,!is.element(baseline2.combat$External.ID, baseline2.remove)]

gen.boxplot("baseline2_vst_rmreps.jpg", baseline2.rmreps, "VST Normalized Signal Intensity", "Intensity")
gen.histogram("baseline2_histogram_rmreps", baseline2.rmreps)
gen.mdsplot("baseline2_mdsplot_group_rmreps", baseline2.rmreps, "Group")
gen.mdsplot("baseline2_mdsplot_batch_rmreps", baseline2.rmreps, "Batch")

baseline2.rmreps$Group %<>% factor
baseline2.rmreps$Batch %<>% factor

lm.age <- lm(Age ~ Group, data = pData(baseline2.rmreps)) %>% anova %>% tidy
lm.rin <- lm(RIN ~ Group, data = pData(baseline2.rmreps)) %>% anova %>% tidy

lm.gender <- lm(Gender ~ Group, data = pData(baseline2.rmreps)) %>% anova %>% tidy
lm.ethnicity <- lm(Ethnicity ~ Gender, data = pData(baseline2.rmreps)) %>% anova %>% tidy
lm.batch <- lm(as.numeric(Batch) ~ Group, data = pData(baseline2.rmreps)) %>% anova %>% tidy

p <- ggplot(pData(baseline2.rmreps), aes(x = Group, y = Age)) + geom_boxplot()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank()) + ggtitle(paste("Age (p < ", round(lm.age$p.value[1], 3), ")", sep = ""))
CairoPDF("baseline2_age_boxplot", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(baseline2.rmreps), aes(x = Group, y = RIN)) + geom_boxplot()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank()) + ggtitle(paste("RIN (p < ", round(lm.age$p.value[1], 3), ")", sep = ""))
CairoPDF("baseline2_rin_boxplot", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(baseline2.rmreps), aes(x = Group, fill = factor(Gender))) + geom_bar()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.title = element_blank()) 
p <- p + ggtitle(paste("Gender (p < ", round(lm.gender$p.value[1], 3), ")", sep = "")) 
CairoPDF("baseline2_gender_barplot", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(baseline2.rmreps), aes(x = Group, fill = factor(Ethnicity))) + geom_bar()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.title = element_blank()) 
p <- p + ggtitle(paste("Ethnicity (p < ", round(lm.ethnicity$p.value[1], 3), ")", sep = "")) 
CairoPDF("baseline2_ethnicity_barplot", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(baseline2.rmreps), aes(x = Group, fill = factor(Batch))) + geom_bar()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.title = element_blank()) 
p <- p + ggtitle(paste("Batch (p < ", round(lm.batch$p.value[1], 3), ")", sep = "")) 
CairoPDF("baseline2_batch_barplot", height = 6, width = 6)
plot(p)
dev.off()

baseline2.symbols <- featureNames(baseline2.rmreps) %>% getSYMBOL('lumiHumanAll.db') %>% factor
baseline2.collapse.expr <- collapseRows(exprs(baseline2.rmreps), rowGroup = baseline2.symbols, rownames(exprs(baseline2.rmreps)))
baseline2.collapse <- ExpressionSet(assayData = baseline2.collapse.expr$datETcollapsed, phenoData = phenoData(baseline2.rmreps))
saveRDS.gz(baseline2.collapse, "./save/baseline2.collapse.rda")

model.baseline2 <- model.matrix(~ Group + Ethnicity + Gender + Age + RIN, data = pData(baseline2.collapse)) %>% data.frame
baseline2.fit <- lmFit(exprs(baseline2.collapse), model.baseline2) %>% eBayes
top.baseline2 <- topTable(baseline2.fit, coef = 2, n = Inf)
top.baseline2$Symbol <- rownames(top.baseline2)
write.xlsx(top.baseline2, "baseline2_genetable.xlsx")

lumi.baseline1 <- lumi.vst[,lumi.vst$Time.Point == "TP1"]
gen.boxplot("baseline1_vst.jpg", lumi.baseline1, "VST Normalized Signal Intensity", "Intensity")
gen.histogram("baseline1_histogram", lumi.baseline1)
gen.mdsplot("baseline1_mdsplot_group", lumi.baseline1, "Group")
gen.mdsplot("baseline1_mdsplot_batch", lumi.baseline1, "Batch")

lumi.baseline1.norm <- lumiN(lumi.baseline1, "rsn")
lumi.baseline1.cutoff <- detectionCall(lumi.baseline1.norm)
lumi.baseline1.expr <- lumi.baseline1.norm[which(lumi.baseline1.cutoff > 0),]
symbols.baseline1 <- getSYMBOL(rownames(lumi.baseline1.expr), "lumiHumanAll.db") %>% is.na
lumi.baseline1.annot <- lumi.baseline1.expr[!symbols.baseline1,]

model.combat1 <- model.matrix(~ factor(Group) + factor(Ethnicity) + factor(Gender) + Age + RIN, data = pData(lumi.baseline1.annot)) %>% data.frame
baseline1.combat.expr <- ComBat(dat = exprs(lumi.baseline1.annot), batch = factor(lumi.baseline1.annot$Batch), mod = model.combat1) 
baseline1.combat <- lumi.baseline1.annot 
exprs(baseline1.combat) <- baseline1.combat.expr

gen.mdsplot("baseline1_combat_mdsplot_group", baseline1.combat, "Group")
gen.mdsplot("baseline1_combat_mdsplot_batch", baseline1.combat, "Batch")

tree.baseline1 <- exprs(baseline1.combat) %>% t %>% dist %>% hclust(method = "average")
CairoPDF("clustering_baseline1", width = 13, height = 10)
plot(tree.baseline1, main = "Hierarchical Clustering Sammples")
dev.off()

baseline1.connectivity <- gen.connectivityplot("baseline1_connectivity", baseline1.combat, "")
baseline1.outlier <- names(baseline1.connectivity[baseline1.connectivity < -3])

baseline1.rmout <- baseline1.combat[,!is.element(baseline1.combat$External.ID, baseline1.outlier)]

gen.boxplot("baseline1_vst_rmout.jpg", baseline1.rmout, "VST Normalized Signal Intensity", "Intensity")
gen.histogram("baseline1_histogram_rmout", baseline1.rmout)
gen.mdsplot("baseline1_mdsplot_group_rmout", baseline1.rmout, "Group")
gen.mdsplot("baseline1_mdsplot_batch_rmout", baseline1.rmout, "Batch")

baseline1.rmout$Group %<>% factor
baseline1.rmout$Batch %<>% factor

lm.age.baseline1 <- lm(Age ~ Group, data = pData(baseline1.rmout)) %>% anova %>% tidy
lm.rin.baseline1 <- lm(RIN ~ Group, data = pData(baseline1.rmout)) %>% anova %>% tidy

lm.gender.baseline1 <- lm(Gender ~ Group, data = pData(baseline1.rmout)) %>% anova %>% tidy
lm.ethnicity.baseline1 <- lm(Ethnicity ~ Gender, data = pData(baseline1.rmout)) %>% anova %>% tidy
lm.batch.baseline1 <- lm(as.numeric(Batch) ~ Group, data = pData(baseline1.rmout)) %>% anova %>% tidy

p <- ggplot(pData(baseline1.rmout), aes(x = Group, y = Age)) + geom_boxplot()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank()) + ggtitle(paste("Age (p < ", round(lm.age$p.value[1], 3), ")", sep = ""))
CairoPDF("baseline1_age_boxplot", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(baseline1.rmout), aes(x = Group, y = RIN)) + geom_boxplot()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank()) + ggtitle(paste("RIN (p < ", round(lm.age$p.value[1], 3), ")", sep = ""))
CairoPDF("baseline1_rin_boxplot", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(baseline1.rmout), aes(x = Group, fill = factor(Gender))) + geom_bar()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.title = element_blank()) 
p <- p + ggtitle(paste("Gender (p < ", round(lm.gender$p.value[1], 3), ")", sep = "")) 
CairoPDF("baseline1_gender_barplot", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(baseline1.rmout), aes(x = Group, fill = factor(Ethnicity))) + geom_bar()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.title = element_blank()) 
p <- p + ggtitle(paste("Ethnicity (p < ", round(lm.ethnicity$p.value[1], 3), ")", sep = "")) 
CairoPDF("baseline1_ethnicity_barplot", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(baseline1.rmout), aes(x = Group, fill = factor(Batch))) + geom_bar()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.title = element_blank()) 
p <- p + ggtitle(paste("Batch (p < ", round(lm.batch$p.value[1], 3), ")", sep = "")) 
CairoPDF("baseline1_batch_barplot", height = 6, width = 6)
plot(p)
dev.off()

baseline1.symbols <- featureNames(baseline1.rmout) %>% getSYMBOL('lumiHumanAll.db') %>% factor
baseline1.collapse.expr <- collapseRows(exprs(baseline1.rmout), rowGroup = baseline1.symbols, rownames(exprs(baseline1.rmout)))
baseline1.collapse <- ExpressionSet(assayData = baseline1.collapse.expr$datETcollapsed, phenoData = phenoData(baseline1.rmout))
saveRDS.gz(baseline1.collapse, "./save/baseline1.collapse.rda")

model.baseline1 <- model.matrix(~ 0 + Group + Ethnicity + Gender + Age + RIN, data = pData(baseline1.collapse)) %>% data.frame
baseline1.fit <- lmFit(exprs(baseline1.collapse), model.baseline1) 
baseline1.contrasts <- makeContrasts(ECT_vs_CONT = GroupECT - GroupCONT, levels = model.baseline1)
baseline1.anova <- contrasts.fit(baseline1.fit, baseline1.contrasts) %>% eBayes

top.baseline1 <- topTable(baseline1.fit, coef = 2, n = Inf)
top.baseline1$Symbol <- rownames(top.baseline1)
write.xlsx(top.baseline1, "baseline1_genetable.xlsx")

decide <- list(c("fdr", 0.05), c("fdr", 0.1), c("none", 0.001), c("none", 0.005), c("none", 0.01)) #Specify significance cutoffs
#decide.baseline1 <- ldply(decide, gen.decide, baseline1.fit, FALSE) 
baseline1.decide.plot <- map_df(decide, gen.decide, baseline1.anova, FALSE) #Compute significance cutoffs
#gen.decideplot("./threshold_selection_baseline1", baseline1.decide.plot, height.plot = 10, width.plot = 7) #Plot different significance cutoffs

decide.final <- gen.decide(c("fdr", 0.1), baseline1.anova, TRUE) #%>% melt(id.vars = c("Test", "Num", "Direction"))
baseline1.ratio.exp <- gen.ratios(baseline1.collapse, "TP1")
baseline1.de.object <- read_tsv("./fit_fdr.tsv") #Read in unadjusted fit object
rownames(baseline1.de.object) <- featureNames(baseline1.collapse)
baseline1.fit.selection <- gen.tables(baseline1.de.object, ect.baseline1, baseline1.ratio.exp, "pLess005_baseline1") #create differential expression table for unadjusted fit
saveRDS.gz(baseline1.fit.selection, file = "./save/baseline1.fit.selection.rda")

gen.anova(baseline1.fit.selection, "none_baseline1", FALSE)

decide.baseline2 <- ldply(decide, gen.decide, baseline2.fit, FALSE) 
colnames(decide.baseline2) <- c("negative", "positive", "Test")
decide.baseline2$Time.Point <- "TP2"
decide.baseline2$negative %<>% as.numeric
decide.baseline2$positive %<>% as.numeric
decide.baseline2$negative <- -(decide.baseline2$negative)

decide.combined <- rbind(decide.baseline1, decide.baseline2)
decide.plot <- melt(decide.combined, id.vars = c("Test", "Time.Point"))
colnames(decide.plot)[3:4] <- c("Direction", "Num.Genes")

p <- ggplot()
p <- p + geom_bar(data = filter(decide.plot, Direction == "positive"),  aes(x = Time.Point, y = Num.Genes), stat = "identity", colour = "black", fill = "red", position = "dodge")   
p <- p + geom_text(data = filter(decide.plot, Direction == "positive"), stat = "identity", size = 4, aes(x = Time.Point, y = Num.Genes, ymax = max(Num.Genes) + 110, ymin = min(Num.Genes) - 110, label = Num.Genes), hjust = -0.3, position = position_dodge(width = 1))
p <- p + geom_bar(data = filter(decide.plot, Direction == "negative"),  aes(x = Time.Point, y = Num.Genes), stat = "identity", colour = "black", fill = "green", position = "dodge")   
p <- p + geom_text(data = filter(decide.plot, Direction == "negative"), stat = "identity", size = 4, aes(x = Time.Point, y = Num.Genes, ymax = max(Num.Genes) + 110, ymin = min(Num.Genes) - 110, label = Num.Genes), hjust = 1.3, position = position_dodge(width = 1))
p <- p + facet_grid(Test ~ .) 
p <- p + theme_bw() + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p <- p + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0))# + ylab("Differentially Expressed Genes")
CairoPDF("threshold_selection", width = 6, height = 7)
print(p)
dev.off()

top.fdr.baseline1 <- slice(top.baseline1, 1:500)
top.fdr.baseline2 <- slice(top.baseline2, 1:500)

enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "WikiPathways_2016", "Reactome_2016", "BioCarta_2016", "PPI_Hub_Proteins", "Humancyc_2016", "NCI-Nature_2016", "Panther_2016") 
baseline1.enrichr <- enrichr.submit(top.fdr.baseline1, enrichr.terms, "baseline1")
baseline2.enrichr <- enrichr.submit(top.fdr.baseline2, enrichr.terms, "baseline2")

baseline1.gobiol.file <- "./enrichr/baseline1/GO_Biological_Process_2015.xlsx"
baseline1.gobiol <- read.xlsx(baseline1.gobiol.file)
baseline1.gobiol$Database <- "GO Biological Process"
baseline1.gobiol$Num.Genes <- map(baseline1.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
baseline1.gobiol.filter <- filter(baseline1.gobiol, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(baseline1.gobiol.filter, top.fdr.baseline1$Symbol, file_path_sans_ext(baseline1.gobiol.file))
baseline1.gobiol.final <- slice(baseline1.gobiol, c(77, 56, 38, 1, 7, 43, 12, 11, 41))

baseline1.gomolec.file <- "./enrichr/baseline1/GO_Molecular_Function_2015.xlsx"
baseline1.gomolec <- read.xlsx(baseline1.gomolec.file)
baseline1.gomolec$Database <- "GO Molecular Function"
baseline1.gomolec$Num.Genes <- map(baseline1.gomolec$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
baseline1.gomolec.filter <- filter(baseline1.gomolec, Num.Genes > 4) %>% filter(P.value < 0.05)
baseline1.gomolec.final <- slice(baseline1.gomolec, 1)

baseline1.reactome.file <- "./enrichr/baseline1/Reactome_2016.xlsx"
baseline1.reactome <- read.xlsx(baseline1.reactome.file)
baseline1.reactome$Database <- "Reactome"
baseline1.reactome$Num.Genes <- map(baseline1.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
baseline1.reactome.filter <- filter(baseline1.reactome, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(baseline1.reactome.filter, top.fdr.baseline1$Symbol, file_path_sans_ext(baseline1.reactome.file))
baseline1.reactome.final <- slice(baseline1.reactome, c(37, 9, 1))

baseline1.final <- rbind(baseline1.gobiol.final, baseline1.gomolec.final, baseline1.reactome.final)
gen.enrichrplot(baseline1.final, top.fdr.baseline1, "baseline1.enrichr")

baseline2.gobiol.file <- "./enrichr/baseline2/GO_Biological_Process_2015.xlsx"
baseline2.gobiol <- read.xlsx(baseline2.gobiol.file)
baseline2.gobiol$Database <- "GO Biological Process"
baseline2.gobiol$Num.Genes <- map(baseline2.gobiol$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
baseline2.gobiol.filter <- filter(baseline2.gobiol, Num.Genes > 4) %>% filter(P.value < 0.01)
get.kappa.cluster(baseline2.gobiol.filter, top.fdr.baseline2$Symbol, file_path_sans_ext(baseline2.gobiol.file))
baseline2.gobiol.final <- slice(baseline2.gobiol, c(8, 23, 26, 2, 1))

baseline2.gomolec.file <- "./enrichr/baseline2/GO_Molecular_Function_2015.xlsx"
baseline2.gomolec <- read.xlsx(baseline2.gomolec.file)
baseline2.gomolec$Database <- "GO Molecular Function"
baseline2.gomolec$Num.Genes <- map(baseline2.gomolec$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
baseline2.gomolec.filter <- filter(baseline2.gomolec, Num.Genes > 4) %>% filter(P.value < 0.05)
baseline2.gomolec.final <- slice(baseline2.gomolec, c(1, 2, 7))

baseline2.reactome.file <- "./enrichr/baseline2/Reactome_2016.xlsx"
baseline2.reactome <- read.xlsx(baseline2.reactome.file)
baseline2.reactome$Database <- "Reactome"
baseline2.reactome$Num.Genes <- map(baseline2.reactome$Genes, str_split, ",") %>% map(getElement, 1) %>% map_int(length)
baseline2.reactome.filter <- filter(baseline2.reactome, Num.Genes > 4) %>% filter(P.value < 0.05)
get.kappa.cluster(baseline2.reactome.filter, top.fdr.baseline2$Symbol, file_path_sans_ext(baseline2.reactome.file))
baseline2.reactome.final <- slice(baseline2.reactome, c(1, 3, 5, 8))

baseline2.final <- rbind(baseline2.gobiol.final, baseline2.gomolec.final, baseline2.reactome.final)
gen.enrichrplot(baseline2.final, top.fdr.baseline2, "baseline2.enrichr")

#Overlap analysis
list.names <- list.files("../MDD_ReferenceGeneSets/GENE_LISTS/", full.names = TRUE)
gene.lists <- map(list.names, read_delim, delim = "\n")
saveRDS.gz(gene.lists, "./gene.lists.rda")

get.intersect <- function(gene.list, top.object) {
    gene.vector <- gene.list[[1]]
    shared.symbols <- intersect(top.object$Symbol, gene.vector)
    top.filter <- top.object$Symbol[1:length(shared.symbols)]    
    overlap <- intersect(shared.symbols, top.filter)
    overlap.pval <- phyper(length(overlap), length(shared.symbols), (nrow(top.object) - length(shared.symbols)), length(shared.symbols), lower.tail = FALSE)
    return(c(length(overlap), overlap.pval))
}

baseline1.overlaps <- map(gene.lists, get.intersect, top.baseline1) %>% reduce(rbind) %>% data.frame
colnames(baseline1.overlaps) <- c("Overlap.Baseline1", "P.value.Baseline1")

baseline2.overlaps <- map(gene.lists, get.intersect, top.baseline2) %>% reduce(rbind) %>% data.frame
colnames(baseline2.overlaps) <- c("Overlap.Baseline2", "P.value.Baseline2")

list.overlaps <- data.frame(Gene.List = map_chr(gene.lists, names), Num.Genes = map_int(gene.lists, nrow), baseline1.overlaps, baseline2.overlaps)
list.overlaps %<>% arrange(P.value.Baseline1, P.value.Baseline2)
saveRDS.gz(list.overlaps, "./save/list.overlaps.rda")

gen.text.heatmap <- function(cor.dataset, text.matrix, x.names, y.names, maintitle, filename, zlim.plot = c(-1,1), color.scheme = greenWhiteRed(50))
{
    width.dynamic <- 3 + ncol(text.matrix)
    height.dynamic <- 3 + nrow(text.matrix)
    CairoPDF(filename, width = 6, height = 10)
    par(mar = c(8, 16, 3, 3))
    labeledHeatmap(Matrix = cor.dataset, xLabels = x.names, yLabels = y.names, ySymbols = y.names, yColorLabels = TRUE, colors = color.scheme, textMatrix = text.matrix, setStdMargins = F, cex.text = 0.5, zlim = zlim.plot, main = maintitle)
    dev.off()
}

list.pvals <- select(list.overlaps, dplyr::contains("P.value")) %>% as.matrix
list.logpvals <- -log10(list.pvals)
list.count <- select(list.overlaps, dplyr::contains("Overlap")) %>% as.matrix
text.matrix.overlaps <- str_c(list.count, '\n(', signif(list.pvals, 2), ')')
dim(text.matrix.overlaps) <- dim(list.count)
heatmap.labels <- str_c(list.overlaps$Gene.List, " (", list.overlaps$Num.Genes, ")")

#reduce.key <- apply(select(reduce.module.overlaps, contains("Overlap")), 1, `>=`, e2 = 10) %>% apply(2, any)
#reduce.overlaps <- reduce.module.overlaps[reduce.key,]
heatmap.range <- c(0, 2)
gen.text.heatmap(list.logpvals, text.matrix.overlaps, c("Baseline1", "Baseline2"), heatmap.labels, "", "baseline_genelists", heatmap.range, colorRampPalette(c("red", "yellow", "white"))(n = 300))

timeresponse.overlaps <- readRDS.gz("../response_treatment/save/timeresponse.overlap.rda")
treatment.overlaps <- readRDS.gz("../treatment/save/treatment.overlaps.rda")

combined.overlaps <- join(list.overlaps, treatment.overlaps) %>% join(timeresponse.overlaps)
write.xlsx(combined.overlaps, "combined.overlaps.xlsx")
