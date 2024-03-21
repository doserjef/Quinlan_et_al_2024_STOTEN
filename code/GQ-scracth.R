
dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=sample, 
                              design=~site + sampleDate, tidy = F)

keep_genes <- rowSums(counts(dds)) > 10
dds <- dds[ keep_genes, ]
dds <- DESeq(dds)

resSite <- as.data.frame(results(dds, contrast = c("site", "Heather's House", "PSU Arboretum"),  lfcThreshold = 1, alpha =0.05))
resSite[which(resSite$log2FoldChange > 1 & resSite$pvalue < 0.05),] # elastin & uncharacterized

resDate1 <- as.data.frame(results(dds, contrast = c("sampleDate", "7/28/23", "7/26/23"),  lfcThreshold = 1, alpha =0.05))
resDate1[which(resDate1$log2FoldChange > 1 & resDate1$padj < 0.05),] #protein Cep78 homolog & uncharacterized

resDate2 <- as.data.frame(results(dds, contrast = c("sampleDate", "7/28/23", "7/29/23"),  lfcThreshold = 1, alpha =0.05))
resDate2[which(resDate2$log2FoldChange > 1 & resDate2$padj < 0.05),]

resDate3 <- as.data.frame(results(dds, contrast = c("sampleDate", "7/28/23", "6/21/23"),  lfcThreshold = 1, alpha =0.05))
resDate3[which(resDate3$log2FoldChange > 1 & resDate3$padj < 0.05),]

####
resSite <- as.data.frame(results(dds, contrast = c("site", "Heather's House", "PSU Arboretum")))
resSite[which(resSite$log2FoldChange > 1 & resSite$padj < 0.05),] # elastin & uncharacterized

resDate1 <- as.data.frame(results(dds, contrast = c("sampleDate", "7/28/23", "7/26/23")))
resDate1[which(resDate1$log2FoldChange > 1 & resDate1$padj < 0.05),] #protein Cep78 homolog & uncharacterized

resDate2 <- as.data.frame(results(dds, contrast = c("sampleDate", "7/28/23", "7/29/23")))
resDate2[which(resDate2$log2FoldChange > 1 & resDate2$padj < 0.05),]

resDate3 <- as.data.frame(results(dds, contrast = c("sampleDate", "7/28/23", "6/21/23")))
resDate3[which(resDate3$log2FoldChange > 1 & resDate3$padj < 0.05),]

####

dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=sample, 
                              design=~site + treat, tidy = F) # colony + tissue + round +  trt

keep_genes <- rowSums(counts(dds)) > 10
dds <- dds[ keep_genes, ]
dds <- DESeq(dds)

resSite <- as.data.frame(results(dds, contrast = c("site", "Heather's House", "PSU Arboretum"),  lfcThreshold = 1, alpha =0.05))
resSite[which(resSite$log2FoldChange > 1 & resSite$padj < 0.05),] # elastin & uncharacterized

resDate1 <- as.data.frame(results(dds, contrast = c("treat", "pre", "post"),  lfcThreshold = 1, alpha =0.05))
resDate1[which(resDate1$log2FoldChange > 1 & resDate1$padj < 0.05),] #protein Cep78 homolog & uncharacterized

####
resSite <- as.data.frame(results(dds, contrast = c("site", "Heather's House", "PSU Arboretum")))
resSite[which(resSite$log2FoldChange > 1 & resSite$padj < 0.05),] # elastin & uncharacterized

resDate1 <- as.data.frame(results(dds, contrast = c("treat", "pre", "post")))
resDate1[which(resDate1$log2FoldChange > 1 & resDate1$padj < 0.05),] #protein Cep78 homolog & uncharacterized


# Plot results ----

# PCA for significant DEG --------  
defaultSig <- default[which(default$tiss=="all" & default$padj < 0.05),"gene"]
keep_genes <- rowSums(counts(dds)) > 10
keepDefault <- ifelse(names(keep_genes) %in% defaultSig & keep_genes==TRUE, TRUE, FALSE)
ddsDefault <- dds[ keepDefault, ] # two genes do not contail "LOC" these are probably the weird trna ala_1, etc. 
