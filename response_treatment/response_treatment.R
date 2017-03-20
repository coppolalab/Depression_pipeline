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
library(dplyr) #for manipulating data
library(magrittr) #for %>% operator
library(purrr) #makes R into Haskell ;)
library(functional) #makes R into Haskell ;)
library(vadr) #makes R into Haskell ;)
library(tidyr)

source('../../FRDA project/common_functions.R')
source('../../code/GO/enrichr.R')

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
    dataset.addvars <- mutate(expr.df, Sample.Status = sampleNames(lumi.object), Time.Point = lumi.object$Time.Point)
    dataset.m <- melt(dataset.addvars, id = c("Sample.Status", "Time.Point"))

    p <- ggplot(dataset.m, aes(value, group = Sample.Status, col = factor(Time.Point))) + geom_density() + theme_bw()
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + ggtitle("Histogram of VST Expression") + ylab("Density") + xlab("VST Expression") 
    CairoPDF(filename, height = 5, width = 9)
    print(p)
    dev.off()
}

gen.mdsplot <- function(filename, lumi.object, variable)
{
    mds.vst <- exprs(lumi.object) %>% t %>% dist %>% cmdscale(eig = TRUE) 
    mds.vst.plot <- data.frame(Sample.Name = rownames(mds.vst$points), Time.Point = lumi.object$Time.Point, Batch = factor(lumi.object$Batch), Component.1 = mds.vst$points[,1], Component.2 = mds.vst$points[,2])

    p <- ggplot(mds.vst.plot, aes_string(x = "Component.1", y = "Component.2", col = variable)) + geom_point() 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + xlab("Component 1") + ylab("Component 2") + ggtitle("MDS of Group")
    CairoPDF(file = filename, height = 6, width = 7)
    print(p)
    dev.off()
}

gen.connectivityplot <- function(filename, dataset, maintitle)
{
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

gen.decide <- function(test, fit.object, write.results)
{
    results <- decideTests(fit.object, adjust.method = test[1], p = as.numeric(test[2])) #Run test at specified cutoff
    if(write.results == TRUE)
    {
        write.fit(file = paste("./fit_", test[1], ".tsv", sep = ""), fit.object, adjust = test[1], results = results)
    }
    num.genes <- apply(results, 1, any, na.rm = T) %>% which %>% length  #Make this better
    mysum <- summary(results)[-2,] #Eliminate the row for no change in expression
    mysum[1,] <- -(mysum[1,])
    mysum.return <- data.frame(Test = paste(test[1], " p<", test[2], sep = ""), Num = paste(num.genes, "Genes", sep = " "), Direction = c("negative", "positive"), as.data.frame.matrix(mysum)) %>% data.frame
    return(mysum.return)
}

gen.decideplot <- function(filename, decide.plot, width.plot = 6, height.plot = 7)
{
    decide.plot$Comparison <- str_replace_all(decide.plot$Comparison, "_", " ")
    decide.plot$Test %<>% factor(levels = levels(factor(decide.plot$Test)))
    p <- ggplot()
    p <- p + geom_bar(data = subset(decide.plot, Direction == "positive"),  aes(x = Comparison, y = Count), stat = "identity", colour = "black", fill = "red", position = "dodge")   
    p <- p + geom_text(data = subset(decide.plot, Direction == "positive"), stat = "identity", size = 4, aes(x = Comparison, y = Count, ymax = max(Count) + 110, label = Count), hjust = -0.3, position = position_dodge(width = 1))
    p <- p + geom_bar(data = subset(decide.plot, Direction == "negative"),  aes(x = Comparison, y = Count), stat = "identity", colour = "black", fill = "green", position = "dodge") 
    p <- p + geom_text(data = subset(decide.plot, Direction == "negative"), stat = "identity", size = 4, aes(x = Comparison, y = Count, ymax = min(Count) - 110, label = abs(Count)), hjust = 1.3, position = position_dodge(width = 1))
    if (length(unique(decide.plot$Test)) > 1)
    {
        p <- p + facet_grid(Test + Num ~ .) 
        #p <- p + ggtitle("Threshold Selection")
    }
    #else
    #{
        #p <- p + ggtitle(paste(decide.plot$Test, "\n", decide.plot$Num))
    #}
    p <- p + theme_bw() + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    p <- p + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0))# + ylab("Differentially Expressed Genes")
    CairoPDF(filename, width = width.plot, height = height.plot)
    print(p)
    dev.off()
}

gen.ratios <- function(lumi.object)
{
    TP1.No <- exprs(lumi.object[,lumi.object$Time.Response == "TP1.No"])
    TP1.No.means <- rowMeans(TP1.No)
    TP2.No <- exprs(lumi.object[,lumi.object$Time.Response == "TP2.No"])
    TP2.No.means <- rowMeans(TP2.No)
    TP3.No <- exprs(lumi.object[,lumi.object$Time.Response == "TP3.No"])
    TP3.No.means <- rowMeans(TP3.No)

    TP1.Yes <- exprs(lumi.object[,lumi.object$Time.Response == "TP1.Yes"])
    TP2.Yes <- exprs(lumi.object[,lumi.object$Time.Response == "TP2.Yes"])
    TP3.Yes <- exprs(lumi.object[,lumi.object$Time.Response == "TP3.Yes"])

    TP1 <- TP1.Yes - TP1.No.means
    TP2 <- TP2.Yes - TP2.No.means
    TP3 <- TP3.Yes - TP3.No.means

    all.coefficients <- data.frame("Symbol" = rownames(TP1), TP1, TP2, TP3)
    all.samples <- data.frame("Symbol" = featureNames(lumi.object), exprs(lumi.object))
    colnames(all.samples)[2:length(all.samples)] %<>% str_c("expr", sep = ".") 
    ratio.exp <- merge(all.coefficients, all.samples)
    return(ratio.exp)
}

gen.workbook <- function(dataset, filename)
{
    pval.cols <- colnames(dataset) %>% str_detect("p.value.") %>% which
    coef.cols <- colnames(dataset) %>% str_detect("Coef.") %>% which
    colnames(dataset)[coef.cols] %<>% str_replace("Coef.", "") %>% str_replace_all("_", " ")
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
gen.tables <- function(dataset, lumi.object, ratio.exp, suffix)
{
    treat.de <- data.frame("Symbol" = rownames(dataset), dataset)
    #colnames(treat.de)[str_detect(colnames(treat.de), "Genes")] %<>% str_replace("Genes\\.", "") %>% tolower %>% capitalize

    fitsel.ratio.all <- merge(treat.de, ratio.exp)

    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    bm.table <- getBM(attributes = c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', values = as.character(fitsel.ratio.all$Symbol), mart = ensembl)
    bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
    colnames(bm.table) <- c("Symbol", "Definition")
    fitsel.ratio.all <- join(fitsel.ratio.all, bm.table)

    fitsel.return.all <- select(fitsel.ratio.all, Symbol, Definition, dplyr::contains("Coef."), dplyr::contains("p.value."), F, F.p.value, dplyr::contains("Res."), matches("^t."), A, dplyr::contains("DE")) %>% arrange(desc(F))

    anovalist <- apply(select(fitsel.return.all, dplyr::contains("Res.")), 1, any, na.rm = T) %>% which
    fitsel.return <- fitsel.return.all[anovalist,]
    coef.cols <- colnames(fitsel.return) %>% str_detect("Coef.") %>% which

    gen.workbook(fitsel.return, paste("./significant_geneList_", suffix, ".xlsx", sep = ""))

    write.csv(fitsel.return.all, paste("./complete_genelist_", suffix, ".csv", sep = ""), row.names = FALSE)
    return(fitsel.return.all)
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

#Make anova objects
gen.anova <- function(dataset, suffix)
{
    plot.TP1 <- filter(dataset, Res.Yes_vs_No_TP1 != 0) %>% select(matches("_01$"))
    plot.TP2 <- filter(dataset, Res.Yes_vs_No_TP2 != 0) %>% select(matches("_02$"))
    plot.TP3 <- filter(dataset, Res.Yes_vs_No_TP3 != 0) %>% select(matches("_03$"))

    gen.anova.heatmap(paste("./anova_heatmap_TP1", suffix, sep = "_"), plot.TP1, "Yes vs. No")
    gen.anova.heatmap(paste("./anova_heatmap_TP2", suffix, sep = "_"), plot.TP2, "Yes vs. No")
    gen.anova.heatmap(paste("./anova_heatmap_TP3", suffix, sep = "_"), plot.TP3, "Yes vs. No")

    #return(list(anova.pco, anova.pca))
}

#Generate anova heatmaps
gen.anova.heatmap <- function(filename, dataset, maintitle)
{ 
    CairoPDF(filename, width = 10, height = 10)
    heatmap.object <- heatmap.2(as.matrix(dataset), col = rev(redgreen(48)), breaks=c(seq(-1.2, 1.2, 0.05)), trace = "none", cexCol = 0.3, labRow = "", keysize = 0.9)
    dev.off()
    return(heatmap.object)
}

ect.collapse <- readRDS.gz("../treatment/save/ect.collapse.rda")

ect.collapse$Time.Response <- str_c(ect.collapse$Time.Point, ect.collapse$responder_to_treatment, sep = ".")
ect.timeresponse <- ect.collapse[,ect.collapse$Time.Point != "TP4"]
model.timeresponse <- model.matrix(~ 0 + Time.Response + Ethnicity + Gender + Age + RIN, data = pData(ect.timeresponse)) %>% data.frame

timeresponse.fit <- lmFit(exprs(ect.timeresponse), design = model.timeresponse) 

timeresponse.contrasts <- makeContrasts(Yes_vs_No_TP1 = Time.ResponseTP1.Yes - Time.ResponseTP1.No, Yes_vs_No_TP2 = Time.ResponseTP2.Yes - Time.ResponseTP2.No, Yes_vs_No_TP3 = Time.ResponseTP3.Yes - Time.ResponseTP3.No, levels = model.timeresponse)
timeresponse.anova <- contrasts.fit(timeresponse.fit, timeresponse.contrasts) %>% eBayes

decide <- list(c("fdr", 0.05), c("fdr", 0.1), c("none", 0.001), c("none", 0.005), c("none", 0.01)) #Specify significance cutoffs
timeresponse.decide.plot <- map_df(decide, gen.decide, timeresponse.anova, FALSE) %>% gather(Comparison, Count, -Test, -Num, -Direction) #Compute significance cutoffs
gen.decideplot("./threshold_selection_timeresponse", timeresponse.decide.plot, height.plot = 10, width.plot = 7) #Plot different significance cutoffs

decide.final <- gen.decide(c("none", 0.005), timeresponse.anova, TRUE) %>% melt(id.vars = c("Test", "Num", "Direction"))
timeresponse.ratio.exp <- gen.ratios(ect.timeresponse)
timeresponse.de.object <- read_tsv("./fit_none.tsv") #Read in unadjusted fit object
rownames(timeresponse.de.object) <- featureNames(ect.timeresponse)
timeresponse.fit.selection <- gen.tables(timeresponse.de.object, ect.timeresponse, timeresponse.ratio.exp, "pLess005_timeresponse") #create differential expression table for unadjusted fit
saveRDS.gz(timeresponse.fit.selection, file = "./save/timeresponse.fit.selection.rda")

gen.anova(timeresponse.fit.selection, "none_timeresponse")

#Overlap
gene.lists <- readRDS.gz("../baseline/gene.lists.rda")

top.iterate <- function(coef.num, fit.object, suffix) {
    top.object <- topTable(fit.object, coef = coef.num, n = Inf)
    top.object$Symbol <- rownames(top.object)
    top.overlaps <- map(gene.lists, get.intersect, top.object) %>% reduce(rbind) %>% data.frame
    suffix.num <- str_c(suffix, coef.num)
    colnames(top.overlaps) <- c(str_c("Overlap.", suffix.num), str_c("P.value.", suffix.num))
    return(top.overlaps)
}

get.intersect <- function(gene.list, top.object) {
    gene.vector <- gene.list[[1]]
    shared.symbols <- intersect(top.object$Symbol, gene.vector)
    top.filter <- top.object$Symbol[1:length(shared.symbols)]    
    overlap <- intersect(shared.symbols, top.filter)
    overlap.pval <- phyper(length(overlap), length(shared.symbols), (nrow(top.object) - length(shared.symbols)), length(shared.symbols), lower.tail = FALSE)
    return(c(length(overlap), overlap.pval))
}

gen.text.heatmap <- function(cor.dataset, text.matrix, x.names, y.names, maintitle, filename, zlim.plot = c(-1,1), color.scheme = greenWhiteRed(50))
{
    width.dynamic <- 3 + ncol(text.matrix)
    height.dynamic <- 3 + nrow(text.matrix)
    CairoPDF(filename, width = width.dynamic, height = 10)
    par(mar = c(8, 16, 3, 3))
    labeledHeatmap(Matrix = cor.dataset, xLabels = x.names, yLabels = y.names, ySymbols = y.names, yColorLabels = TRUE, colors = color.scheme, textMatrix = text.matrix, setStdMargins = F, cex.text = 0.5, zlim = zlim.plot, main = maintitle)
    dev.off()
}

timeresponse.overlaps.list <- map(1:ncol(timeresponse.contrasts), top.iterate, timeresponse.anova, "timeresponse") %>% reduce(cbind)
timeresponse.overlaps <- data.frame(Gene.List = map_chr(gene.lists, names), Num.Genes = map_int(gene.lists, nrow), timeresponse.overlaps.list)
colnames(timeresponse.overlaps) %<>% str_replace("timeresponse", "Yes_vs_No_TP")
saveRDS.gz(timeresponse.overlaps, "./save/timeresponse.overlap.rda")

timeresponse.pvals <- select(timeresponse.overlaps, dplyr::contains("P.value")) %>% as.matrix
timeresponse.logpvals <- -log10(timeresponse.pvals)
timeresponse.count <- select(timeresponse.overlaps, dplyr::contains("Overlap")) %>% as.matrix
text.matrix.overlaps <- str_c(timeresponse.count, '\n(', signif(timeresponse.pvals, 2), ')')
dim(text.matrix.overlaps) <- dim(timeresponse.count)
heatmap.labels <- str_c(timeresponse.overlaps$Gene.List, " (", timeresponse.overlaps$Num.Genes, ")")

heatmap.range <- c(0, 2)
gen.text.heatmap(timeresponse.logpvals, text.matrix.overlaps, colnames(timeresponse.contrasts), heatmap.labels, "", "timeresponse_genelists", heatmap.range, colorRampPalette(c("red", "yellow", "white"))(n = 300))
