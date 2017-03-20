#Specific packages
library(limma) #for DE 
library(lumiHumanIDMapping) #for converting Illumina probe IDs to nuIDs
library(lumiHumanAll.db) #for getting symbols for nuIDs
library(annotate) #for getting symbols for nuIDs
library(lumi) #main preprocessing package
library(WGCNA) #for fastcor and connecitivity stats
library(biomaRt) #for extra annotation of genes 
library(sva) #for ComBat
library(siggenes)

#Plotting
library(ggplot2) #main plotting package
library(Cairo) #Cairo graphics device
library(heatmap.plus) #main heatmap
library(gplots) #for heatmap.2

#Input/output manipulation
library(readr) #for reading csv, tsv
library(openxlsx) #best pacakage for reading writing xlsx
library(stringr) #much better string processing
library(broom) #for cleaning up output of lm

#Language enhancements
library(dplyr) #for manipulating data
library(magrittr) #for %>% operator
library(purrr) #makes R into Haskell ;)
library(tidyr)

source('../../FRDA project/common_functions.R')
source('../../code/GO/enrichr.R')

#Boxplot function
BoxPlot <- function(filename, lumi.object, maintext, ylabtext) {
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
gen.heatmapbars <- function(batch.colorscheme, diagnosis.colorscheme, targetset) {
    batch.heatmap <- data.frame("Batch" = seq(1:length(batch.colorscheme)), "Batch.Color" = batch.colorscheme)
    diagnosis.heatmap <- data.frame("Status" = levels(targetset$Status), "Diagnosis.Color" = diagnosis.colors)
    colorscheme <- data.frame("Batch" = targetset$Batch, "Status" = targetset$Status) %>% join(batch.heatmap) %>% join(diagnosis.heatmap)
    colorscheme <- as.matrix(subset(colorscheme, select = c("Batch.Color", "Diagnosis.Color")))
    return(colorscheme)
}

#Heatmap function
gen.heatmap <- function(filename, lumi.object, maintitle) {
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

gen.decide <- function(test, fit.object, write.results, suffix = "")
{
    results <- decideTests(fit.object, adjust.method = test[1], p = as.numeric(test[2])) #Run test at specified cutoff
    if(write.results == TRUE)
    {
        write.fit(file = paste("./fit_", test[1], "_", suffix, ".tsv", sep = ""), fit.object, adjust = test[1], results = results)
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

gen.ratios <- function(lumi.object, TP4 = TRUE)
{
    TP1 <- exprs(lumi.object[,lumi.object$Time.Point == "TP1"])
    TP1.means <- rowMeans(TP1)
    TP2 <- exprs(lumi.object[,lumi.object$Time.Point == "TP2"])
    TP2.means <- rowMeans(TP2)
    TP3 <- exprs(lumi.object[,lumi.object$Time.Point == "TP3"])
    TP3.means <- rowMeans(TP3)

    TP2.TP1 <- TP2 - TP1.means
    TP3.TP1 <- TP3 - TP1.means
    TP3.TP2 <- TP3 - TP2.means

    colnames(TP3.TP2) %<>% str_c("_vs_TP2")

    if (TP4 == TRUE) {
        TP4 <- exprs(lumi.object[,lumi.object$Time.Point == "TP4"])

        TP4.TP1 <- TP4 - TP1.means
        TP4.TP2 <- TP4 - TP2.means
        TP4.TP3 <- TP4 - TP3.means

        colnames(TP4.TP2) %<>% str_c("_vs_TP2")
        colnames(TP4.TP3) %<>% str_c("_vs_TP3")

        all.coefficients <- data.frame("Symbol" = rownames(TP2.TP1), TP2.TP1, TP3.TP1, TP3.TP2, TP4.TP1, TP4.TP2, TP4.TP3)
    } else {
        all.coefficients <- data.frame("Symbol" = rownames(TP2.TP1), TP2.TP1, TP3.TP1, TP3.TP2)
    }

    all.samples <- data.frame("Symbol" = featureNames(lumi.object), exprs(lumi.object))
    colnames(all.samples)[2:length(all.samples)] %<>% str_c("expr", sep = ".") 
    ratio.exp <- merge(all.coefficients, all.samples)
    return(ratio.exp)
}

gen.workbook <- function(dataset, filename)
{
    pval.cols <- colnames(dataset) %>% str_detect("P.Value") %>% which
    adj.pval.cols <- colnames(dataset) %>% str_detect("adj.P.Val") %>% which
    coef.cols <- colnames(dataset) %>% str_detect("logFC") %>% which
    #dataset$Definition %<>% str_replace_all("Homo sapiens ", "") %>% str_replace_all("PREDICTED: ", "")

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(dataset), rule = "<0.005", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = adj.pval.cols, rows = 1:nrow(dataset), rule = "<0.05", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = coef.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1, widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    setColWidths(wb, 1, cols = 3:8, widths = 18)
    setColWidths(wb, 1, cols = 9:20, widths = 19)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    #modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")

    saveWorkbook(wb, filename, overwrite = TRUE) 
}

#gen.workbook <- function(dataset, filename)
#{
    #pval.cols <- colnames(dataset) %>% str_detect("p.value.") %>% which
    #coef.cols <- colnames(dataset) %>% str_detect("Coef.") %>% which
    #colnames(dataset)[coef.cols] %<>% str_replace("Coef.", "") %>% str_replace_all("_", " ")
    ##dataset$Definition %<>% str_replace_all("Homo sapiens ", "") %>% str_replace_all("PREDICTED: ", "")

    #wb <- createWorkbook()
    #addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    #writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = FALSE)
    #sig.pvalues <- createStyle(fontColour = "red")
    #conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(dataset), rule = "<0.05", style = sig.pvalues)
    #conditionalFormatting(wb, 1, cols = coef.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    #setColWidths(wb, 1, cols = 1, widths = "auto")
    #setColWidths(wb, 1, cols = 2, widths = 45)
    #setColWidths(wb, 1, cols = 3:8, widths = 10)
    #setColWidths(wb, 1, cols = 9:14, widths = 15)
    #pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    #freezePane(wb, 1, firstRow = TRUE)
    #showGridLines(wb, 1, showGridLines = TRUE)
    #modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")

    #saveWorkbook(wb, filename, overwrite = TRUE) 
#}

#Create genelists
gen.tables <- function(dataset, lumi.object, ratio.exp, suffix)
{
    fitsel.ratio.all <- merge(dataset, ratio.exp)

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
gen.anova <- function(dataset, suffix, TP4 = TRUE)
{
    plot.TP2.TP1 <- filter(dataset, Res.TP2_vs_TP1 != 0) %>% select(matches("_02$"))
    plot.TP3.TP1 <- filter(dataset, Res.TP3_vs_TP1 != 0) %>% select(matches("_03$"))
    plot.TP3.TP2 <- filter(dataset, Res.TP3_vs_TP2 != 0) %>% select(matches("_03_vs_TP2$"))


    gen.anova.heatmap(paste("./anova_heatmap_TP2_TP1", suffix, sep = "_"), plot.TP2.TP1, "TP2 vs. TP1")
    gen.anova.heatmap(paste("./anova_heatmap_TP3_TP1", suffix, sep = "_"), plot.TP3.TP1, "TP3 vs. TP1")
    gen.anova.heatmap(paste("./anova_heatmap_TP3_TP2", suffix, sep = "_"), plot.TP3.TP2, "TP3 vs. TP2")

    if (TP4 == TRUE) {
        plot.TP4.TP1 <- filter(dataset, Res.TP4_vs_TP1 != 0) %>% select(matches("_04$"))
        plot.TP4.TP2 <- filter(dataset, Res.TP4_vs_TP2 != 0) %>% select(matches("_04_vs_TP2$"))
        plot.TP4.TP3 <- filter(dataset, Res.TP4_vs_TP3 != 0) %>% select(matches("_04_vs_TP3$"))

        gen.anova.heatmap(paste("./anova_heatmap_TP4_TP1", suffix, sep = "_"), plot.TP4.TP1, "TP4 vs. TP1")
        gen.anova.heatmap(paste("./anova_heatmap_TP4_TP2", suffix, sep = "_"), plot.TP4.TP2, "TP4 vs. TP2")
        gen.anova.heatmap(paste("./anova_heatmap_TP4_TP3", suffix, sep = "_"), plot.TP4.TP3, "TP4 vs. TP3")
    }
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

lumi.vst <- readRDS.gz("../baseline/save/lumi.vst.rda")

lumi.ect <- lumi.vst[,lumi.vst$Group == "ECT"]
BoxPlot("ect_vst.jpg", lumi.ect, "VST Normalized Signal Intensity", "Intensity")
gen.histogram("ect_histogram", lumi.ect)
gen.mdsplot("ect_mdsplot_timepoint", lumi.ect, "Time.Point")
gen.mdsplot("ect_mdsplot_batch", lumi.ect, "Batch")

lumi.ect.norm <- lumiN(lumi.ect, "rsn")
lumi.ect.cutoff <- detectionCall(lumi.ect.norm)
lumi.ect.expr <- lumi.ect.norm[which(lumi.ect.cutoff > 0),]
symbols.lumi.expr <- getSYMBOL(rownames(lumi.ect.expr), "lumiHumanAll.db") %>% is.na
lumi.ect.annot <- lumi.ect.expr[!symbols.lumi.expr,]

model.combat <- model.matrix(~ factor(Time.Point) + factor(Ethnicity) + factor(Gender) + Age + RIN, data = pData(lumi.ect.annot)) %>% data.frame
ect.combat.expr <- ComBat(dat = exprs(lumi.ect.annot), batch = factor(lumi.ect.annot$Batch), mod = model.combat) 
ect.combat <- lumi.ect.annot 
exprs(ect.combat) <- ect.combat.expr

gen.mdsplot("combat_mdsplot_timepoint", ect.combat, "Time.Point")
gen.mdsplot("combat_mdsplot_batch", ect.combat, "Batch")

tree.ect <- exprs(ect.combat) %>% t %>% dist %>% hclust(method = "average")
CairoPDF("clustering_ect", width = 13, height = 10)
plot(tree.ect, main = "Hierarchical Clustering Sammples")
dev.off()

ect.connectivity <- gen.connectivityplot("ect_connectivity", ect.combat, "")
ect.outlier <- names(ect.connectivity[ect.connectivity < -2])
ect.names <- names(ect.connectivity)
ect.reps <- str_subset(ect.names, "rep")
ect.orig <- str_replace(ect.reps, "_rep.*$", "")
ect.reps.df <- data.frame(Orig = abs(ect.connectivity[match(ect.orig, names(ect.connectivity))]), Replicate = abs(ect.connectivity[match(ect.reps, names(ect.connectivity))]))
ect.whichmax <- apply(ect.reps.df, 1, which.max)
ect.reps.name <- data.frame(Orig = names(ect.connectivity[match(ect.orig, names(ect.connectivity))]), Replicate = names(ect.connectivity[match(ect.reps, names(ect.connectivity))]))
ect.rmreps.key <- ect.reps.name[cbind(seq_along(ect.whichmax), ect.whichmax)]

ect.remove <- c(ect.outlier, ect.rmreps.key) %>% unique
ect.rmreps <- ect.combat[,!is.element(ect.combat$External.ID, ect.remove)]

BoxPlot("ect_vst_rmreps.jpg", ect.rmreps, "VST Normalized Signal Intensity", "Intensity")
gen.histogram("ect_histogram_rmreps", ect.rmreps)
gen.mdsplot("ect_mdsplot_timepoint_rmreps", ect.rmreps, "Time.Point")
gen.mdsplot("ect_mdsplot_batch_rmreps", ect.rmreps, "Batch")

ect.rmreps$Time.Point %<>% factor
ect.rmreps$Batch %<>% factor

lm.age <- lm(Age ~ Time.Point, data = pData(ect.rmreps)) %>% anova %>% tidy
lm.rin <- lm(RIN ~ Time.Point, data = pData(ect.rmreps)) %>% anova %>% tidy

lm.gender <- lm(Gender ~ Time.Point, data = pData(ect.rmreps)) %>% anova %>% tidy
lm.ethnicity <- lm(Ethnicity ~ Gender, data = pData(ect.rmreps)) %>% anova %>% tidy
lm.batch <- lm(as.numeric(Batch) ~ Time.Point, data = pData(ect.rmreps)) %>% anova %>% tidy

p <- ggplot(pData(ect.rmreps), aes(x = Time.Point, y = Age)) + geom_boxplot()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank()) + ggtitle(paste("Age (p < ", round(lm.age$p.value[1], 3), ")", sep = ""))
CairoPDF("ect_age_boxplot", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(ect.rmreps), aes(x = Time.Point, y = RIN)) + geom_boxplot()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank()) + ggtitle(paste("RIN (p < ", round(lm.age$p.value[1], 3), ")", sep = ""))
CairoPDF("ect_rin_boxplot", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(ect.rmreps), aes(x = Time.Point, fill = factor(Gender))) + geom_bar()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.title = element_blank()) 
p <- p + ggtitle(paste("Gender (p < ", round(lm.gender$p.value[1], 3), ")", sep = "")) 
CairoPDF("ect_gender_barplot", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(ect.rmreps), aes(x = Time.Point, fill = factor(Ethnicity))) + geom_bar()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.title = element_blank()) 
p <- p + ggtitle(paste("Ethnicity (p < ", round(lm.ethnicity$p.value[1], 3), ")", sep = "")) 
CairoPDF("ect_ethnicity_barplot", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(ect.rmreps), aes(x = Time.Point, fill = factor(Batch))) + geom_bar()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.title = element_blank()) 
p <- p + ggtitle(paste("Batch (p < ", round(lm.batch$p.value[1], 3), ")", sep = "")) 
CairoPDF("ect_batch_barplot", height = 6, width = 6)
plot(p)
dev.off()

ect.symbols <- featureNames(ect.rmreps) %>% getSYMBOL('lumiHumanAll.db') %>% factor
ect.collapse.expr <- collapseRows(exprs(ect.rmreps), rowGroup = ect.symbols, rownames(exprs(ect.rmreps)))
ect.collapse <- ExpressionSet(assayData = ect.collapse.expr$datETcollapsed, phenoData = phenoData(ect.rmreps))
saveRDS.gz(ect.collapse, "./save/ect.collapse.rda")

ect.notp4 <- ect.collapse[, ect.collapse$Time.Point != "TP4"]
ect.notp4$Time.Point %<>% droplevels
model.ect <- model.matrix(~ Ethnicity + Gender + Age + RIN + Time.Point + responder_to_treatment + Time.Point * responder_to_treatment, data = pData(ect.collapse)) %>% data.frame
ect.fit <- lmFit(exprs(ect.collapse), design = model.ect) %>% eBayes

#duplicate.correlation <- duplicateCorrelation(object = ect.collapse, design = model.ect, block = factor(ect.collapse$SUBJECT))
#ect.contrasts <- makeContrasts(TP2_vs_TP1 = Time.PointTP2 - Time.PointTP1, TP3_vs_TP1 = Time.PointTP3 - Time.PointTP1, TP3_vs_TP2 = Time.PointTP3 - Time.PointTP2, TP4_vs_TP1 = Time.PointTP4 - Time.PointTP1, TP4_vs_TP2 = Time.PointTP4 - Time.PointTP2, TP4_vs_TP3 = Time.PointTP4 - Time.PointTP3, levels = model.ect)
#fit.anova <- contrasts.fit(ect.fit, ect.contrasts) %>% eBayes

top.TP2_vs_TP1 <- topTable(fit.anova, coef = 1, n = Inf)
colnames(top.TP2_vs_TP1) %<>% str_c("TP2_vs_TP1", sep = ".")
top.TP2_vs_TP1$Symbol <- rownames(top.TP2_vs_TP1)

top.TP3_vs_TP1 <- topTable(fit.anova, coef = 2, n = Inf)
colnames(top.TP3_vs_TP1) %<>% str_c("TP3_vs_TP1", sep = ".")
top.TP3_vs_TP1$Symbol <- rownames(top.TP3_vs_TP1)

top.TP3_vs_TP2 <- topTable(fit.anova, coef = 3, n = Inf)
colnames(top.TP3_vs_TP2) %<>% str_c("TP3_vs_TP2", sep = ".")
top.TP3_vs_TP2$Symbol <- rownames(top.TP3_vs_TP2)

top.TP4_vs_TP1 <- topTable(fit.anova, coef = 4, n = Inf)
colnames(top.TP4_vs_TP1) %<>% str_c("TP4_vs_TP1", sep = ".")
top.TP4_vs_TP1$Symbol <- rownames(top.TP4_vs_TP1)

top.TP4_vs_TP2 <- topTable(fit.anova, coef = 5, n = Inf)
colnames(top.TP4_vs_TP2) %<>% str_c("TP4_vs_TP2", sep = ".")
top.TP4_vs_TP2$Symbol <- rownames(top.TP4_vs_TP2)

top.TP4_vs_TP3 <- topTable(fit.anova, coef = 6, n = Inf)
colnames(top.TP4_vs_TP3) %<>% str_c("TP4_vs_TP3", sep = ".")
top.TP4_vs_TP3$Symbol <- rownames(top.TP4_vs_TP3)

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "uswest.ensembl.org")
bm.table <- getBM(attributes = c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', values = as.character(top.joined$Symbol), mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Symbol", "Definition")

top.joined <- left_join(top.TP2_vs_TP1, top.TP3_vs_TP1) %>% left_join(top.TP3_vs_TP2) %>% left_join(top.TP4_vs_TP1) %>% left_join(top.TP4_vs_TP2) %>% left_join(top.TP4_vs_TP3)
top.annot <- left_join(top.joined, bm.table)
top.ordered <- select(top.annot, Symbol, Definition, dplyr::contains("logFC"), dplyr::contains("P.Value"), dplyr::contains("adj.P.Val"), dplyr::contains("AveExpr"), dplyr::matches("t\\."), dplyr::matches("B\\.")) %>% arrange(P.Value.TP3_vs_TP2)
gen.workbook(top.ordered, "genetable_nodetscore.xlsx")

top.detscore <- read.xlsx("./genetable_detscore.xlsx")
gen.workbook(top.detscore, "genetable_detscore.xlsx")

grey.genelist <- read_tsv("./myGeneList_darkgrey.txt")
grey.nodetscore <- intersect(grey.genelist$Symbol, top.ordered$Symbol)
grey.detscore <- intersect(grey.genelist$Symbol, top.detscore$Symbol)
grey.missing <- grey.genelist$Symbol[!(grey.genelist$Symbol %in% top.detscore$Symbol)]
grey.genelist$Present <- as.character(grey.genelist$Symbol %in% top.detscore$Symbol)
write.xlsx(grey.genelist, "grey.genelist.xlsx")

toptable.test <- topTable(ect.fit, coef = 8, n = Inf)
siggene.ebam <- limma2ebam(ect.fit, coef = 11)
ebam2excel(siggene.ebam, 0.8, "siggene.ebam.TP3_vs_TP2.csv")

tp3tp2 <- ect.collapse[,grepl("TP2|TP3", ect.collapse$Time.Point)]
tp3tp2.a0 <- find.a0(tp3tp2, as.integer(tp3tp2$Time.Point), B = 1000, rand = 12345, delta = 0.8)
ebam.tp3tp2 <- ebam(tp3tp2.a0, delta = 0.8)

decide <- list(c("fdr", 0.05), c("fdr", 0.1), c("none", 0.001), c("none", 0.005), c("none", 0.01)) #Specify significance cutoffs
decide.plot <- map_df(decide, gen.decide, fit.anova, FALSE) %>% gather(Comparison, Count, -Test, -Num, -Direction) #Compute significance cutoffs
gen.decideplot("./threshold_selection", decide.plot, height.plot = 10, width.plot = 7) #Plot different significance cutoffs

ratio.exp <- gen.ratios(ect.collapse)
decide.final <- gen.decide(c("none", 0.005), fit.anova, TRUE) %>% melt(id.vars = c("Test", "Num", "Direction"))

de.object <- read_tsv("./fit_none.tsv")  #Read in unadjusted fit object 
de.object$Symbol <- featureNames(ect.collapse)
fit.selection <- gen.tables(de.object, ect.collapse, ratio.exp, "pLess005") #create differential expression table for unadjusted fit
saveRDS.gz(fit.selection, file = "./save/fit.selection.rda")

gen.anova(fit.selection, "none")

ect.responder <- ect.collapse[,ect.collapse$responder_to_treatment == "Yes"]
model.responder <- model.matrix(~ 0 + Time.Point + Ethnicity + Gender + Age + RIN, data = pData(ect.responder)) %>% data.frame

responder.correlation <- duplicateCorrelation(object = ect.responder, design = model.responder, block = factor(ect.responder$SUBJECT))
responder.fit <- lmFit(exprs(ect.responder), design = model.responder, block = factor(ect.responder$SUBJECT), correlation = duplicate.correlation$consensus.correlation) 

responder.contrasts <- makeContrasts(TP2_vs_TP1 = Time.PointTP2 - Time.PointTP1, TP3_vs_TP1 = Time.PointTP3 - Time.PointTP1, TP3_vs_TP2 = Time.PointTP3 - Time.PointTP2, TP4_vs_TP1 = Time.PointTP4 - Time.PointTP1, TP4_vs_TP2 = Time.PointTP4 - Time.PointTP2, TP4_vs_TP3 = Time.PointTP4 - Time.PointTP3, levels = model.responder)
responder.anova <- contrasts.fit(responder.fit, responder.contrasts) %>% eBayes

responder.decide.plot <- map_df(decide, gen.decide, responder.anova, FALSE) %>% gather(Comparison, Count, -Test, -Num, -Direction) #Compute significance cutoffs
gen.decideplot("./threshold_selection_responder", responder.decide.plot, height.plot = 10, width.plot = 7) #Plot different significance cutoffs

responder.ratio.exp <- gen.ratios(ect.responder)
decide.final <- gen.decide(c("none", 0.005), responder.anova, TRUE, "responder") %>% melt(id.vars = c("Test", "Num", "Direction"))

responder.de.object <- read_tsv("./fit_none_responder.tsv") #Read in unadjusted fit object
rownames(responder.de.object) <- featureNames(ect.responder)
responder.fit.selection <- gen.tables(responder.de.object, ect.responder, responder.ratio.exp, "pLess005_responder") #create differential expression table for unadjusted fit
saveRDS.gz(responder.fit.selection, file = "./save/responder.fit.selection.rda")

gen.anova(responder.fit.selection, "none_responder")

ect.nonresponder <- ect.collapse[,ect.collapse$responder_to_treatment == "No" & ect.collapse$Time.Point != "TP4"]
ect.nonresponder$Time.Point %<>% droplevels
model.nonresponder <- model.matrix(~ 0 + Time.Point + Ethnicity + Gender + Age + RIN, data = pData(ect.nonresponder)) %>% data.frame

nonresponder.correlation <- duplicateCorrelation(object = ect.nonresponder, design = model.nonresponder, block = factor(ect.nonresponder$SUBJECT))
nonresponder.fit <- lmFit(exprs(ect.nonresponder), design = model.nonresponder, block = factor(ect.nonresponder$SUBJECT), correlation = duplicate.correlation$consensus.correlation) 

nonresponder.contrasts <- makeContrasts(TP2_vs_TP1 = Time.PointTP2 - Time.PointTP1, TP3_vs_TP1 = Time.PointTP3 - Time.PointTP1, TP3_vs_TP2 = Time.PointTP3 - Time.PointTP2, levels = model.nonresponder)
nonresponder.anova <- contrasts.fit(nonresponder.fit, nonresponder.contrasts) %>% eBayes

nonresponder.decide.plot <- map_df(decide, gen.decide, nonresponder.anova, FALSE) %>% gather(Comparison, Count, -Test, -Num, -Direction) #Compute significance cutoffs
gen.decideplot("./threshold_selection_nonresponder", nonresponder.decide.plot, height.plot = 10, width.plot = 7) #Plot different significance cutoffs

nonresponder.ratio.exp <- gen.ratios(ect.nonresponder, FALSE)
decide.final <- gen.decide(c("none", 0.005), nonresponder.anova, TRUE, "nonresponder") %>% melt(id.vars = c("Test", "Num", "Direction"))

nonresponder.de.object <- read_tsv("./fit_none_nonresponder.tsv") #Read in unadjusted fit object
rownames(nonresponder.de.object) <- featureNames(ect.nonresponder)
nonresponder.fit.selection <- gen.tables(nonresponder.de.object, ect.nonresponder, nonresponder.ratio.exp, "pLess005_nonresponder") #create differential expression table for unadjusted fit
saveRDS.gz(nonresponder.fit.selection, file = "./save/nonresponder.fit.selection.rda")

gen.anova(nonresponder.fit.selection, "none_nonresponder", FALSE)

top.fdr.baseline1 <- slice(top.baseline1, 1:500)
top.fdr.baseline2 <- slice(top.baseline2, 1:500)

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

treatment.overlaps <- map(1:ncol(ect.contrasts), top.iterate, fit.anova, "treatment") %>% reduce(cbind)
ect.overlaps <- data.frame(Gene.List = map_chr(gene.lists, names), Num.Genes = map_int(gene.lists, nrow), treatment.overlaps)
ect.overlapnames <- str_c("Overlap", colnames(ect.contrasts), sep = ".")
ect.pvaluenames <- str_c("P.value", colnames(ect.contrasts), sep = ".")
ect.colnames <- colnames(ect.overlaps)
colnames(ect.overlaps)[grepl("Overlap", ect.colnames)] <- ect.overlapnames
colnames(ect.overlaps)[grepl("P.value", ect.colnames)] <- ect.pvaluenames

responder.overlaps.list <- map(1:ncol(responder.contrasts), top.iterate, responder.anova, "treatment") %>% reduce(cbind)
responder.overlaps <- data.frame(Gene.List = map_chr(gene.lists, names), Num.Genes = map_int(gene.lists, nrow), responder.overlaps.list)
responder.overlapnames <- str_c("Overlap.Responder", colnames(responder.contrasts), sep = ".")
responder.pvaluenames <- str_c("P.value.Responder", colnames(responder.contrasts), sep = ".")
responder.colnames <- colnames(responder.overlaps)
colnames(responder.overlaps)[grepl("Overlap", responder.colnames)] <- responder.overlapnames
colnames(responder.overlaps)[grepl("P.value", responder.colnames)] <- responder.pvaluenames

nonresponder.overlaps.list <- map(1:ncol(nonresponder.contrasts), top.iterate, nonresponder.anova, "treatment") %>% reduce(cbind)
nonresponder.overlaps <- data.frame(Gene.List = map_chr(gene.lists, names), Num.Genes = map_int(gene.lists, nrow), nonresponder.overlaps.list)
nonresponder.overlapnames <- str_c("Overlap.Nonresponder", colnames(nonresponder.contrasts), sep = ".")
nonresponder.pvaluenames <- str_c("P.value.Nonresponder", colnames(nonresponder.contrasts), sep = ".")
nonresponder.colnames <- colnames(nonresponder.overlaps)
colnames(nonresponder.overlaps)[grepl("Overlap", nonresponder.colnames)] <- nonresponder.overlapnames
colnames(nonresponder.overlaps)[grepl("P.value", nonresponder.colnames)] <- nonresponder.pvaluenames

gen.text.heatmap <- function(cor.dataset, text.matrix, x.names, y.names, maintitle, filename, zlim.plot = c(-1,1), color.scheme = greenWhiteRed(50))
{
    width.dynamic <- 3 + ncol(text.matrix)
    height.dynamic <- 3 + nrow(text.matrix)
    CairoPDF(filename, width = width.dynamic, height = 10)
    par(mar = c(8, 16, 3, 3))
    labeledHeatmap(Matrix = cor.dataset, xLabels = x.names, yLabels = y.names, ySymbols = y.names, yColorLabels = TRUE, colors = color.scheme, textMatrix = text.matrix, setStdMargins = F, cex.text = 0.5, zlim = zlim.plot, main = maintitle)
    dev.off()
}

ect.pvals <- select(ect.overlaps, dplyr::contains("P.value")) %>% as.matrix
ect.logpvals <- -log10(ect.pvals)
ect.count <- select(ect.overlaps, dplyr::contains("Overlap")) %>% as.matrix
text.matrix.overlaps <- str_c(ect.count, '\n(', signif(ect.pvals, 2), ')')
dim(text.matrix.overlaps) <- dim(ect.count)
heatmap.labels <- str_c(ect.overlaps$Gene.List, " (", ect.overlaps$Num.Genes, ")")

heatmap.range <- c(0, 2)
gen.text.heatmap(ect.logpvals, text.matrix.overlaps, colnames(ect.contrasts), heatmap.labels, "", "ect_genelists", heatmap.range, colorRampPalette(c("red", "yellow", "white"))(n = 300))

responder.pvals <- select(responder.overlaps, dplyr::contains("P.value")) %>% as.matrix
responder.logpvals <- -log10(responder.pvals)
responder.count <- select(responder.overlaps, dplyr::contains("Overlap")) %>% as.matrix
text.matrix.overlaps <- str_c(responder.count, '\n(', signif(responder.pvals, 2), ')')
dim(text.matrix.overlaps) <- dim(responder.count)
heatmap.labels <- str_c(responder.overlaps$Gene.List, " (", responder.overlaps$Num.Genes, ")")

gen.text.heatmap(responder.logpvals, text.matrix.overlaps, colnames(responder.contrasts), heatmap.labels, "", "responder_genelists", heatmap.range, colorRampPalette(c("red", "yellow", "white"))(n = 300))

nonresponder.pvals <- select(nonresponder.overlaps, dplyr::contains("P.value")) %>% as.matrix
nonresponder.logpvals <- -log10(nonresponder.pvals)
nonresponder.count <- select(nonresponder.overlaps, dplyr::contains("Overlap")) %>% as.matrix
text.matrix.overlaps <- str_c(nonresponder.count, '\n(', signif(nonresponder.pvals, 2), ')')
dim(text.matrix.overlaps) <- dim(nonresponder.count)
heatmap.labels <- str_c(nonresponder.overlaps$Gene.List, " (", nonresponder.overlaps$Num.Genes, ")")

gen.text.heatmap(nonresponder.logpvals, text.matrix.overlaps, colnames(nonresponder.contrasts), heatmap.labels, "", "nonresponder_genelists", heatmap.range, colorRampPalette(c("red", "yellow", "white"))(n = 300))

overlaps.joined <- join(ect.overlaps, responder.overlaps) %>% join(nonresponder.overlaps)
saveRDS.gz(overlaps.joined, "./save/treatment.overlaps.rda")
