library(stringr)
argv <- commandArgs(T)

judge <- argv[1]
All_group <- argv[2]
DEG_group <- argv[3]
All_sample_count <- argv[4]
output <- argv[5]

count_file <- read.table(All_sample_count,header=T,row.names=1,sep="\t")
A <- unlist(str_split(DEG_group,pattern="_vs_"))[1]
B <- unlist(str_split(DEG_group,pattern="_vs_"))[2]
All_list <- unlist(str_split(All_group,pattern=";"))

A_list <- unlist(str_split(All_list[grep(paste(A,":",sep=""),All_list)],pattern=":"))[2]
A_list <- unlist(str_split(A_list,pattern=","))
B_list <- unlist(str_split(All_list[grep(paste(B,":",sep=""),All_list)],pattern=":"))[2]
B_list <- unlist(str_split(B_list,pattern=","))

##计算差异
if(judge=="edgeR"){

  library(edgeR)
  library(limma)
  
  count <- count_file[,c(c(B_list),c(A_list))]
  group<-factor(rep(c(B,A),c(length(B_list),length(A_list))))
  dge = DGEList(counts = count)
  keep <- rowSums(cpm(dge) > 1) >= 0.5 * length(group)
  dge <- dge[keep, , keep.lib.sizes = TRUE]
  dge <- calcNormFactors(dge)

  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  comparison <- paste(A , B, sep = "-")
  contrast.matrix <- makeContrasts(contrasts = comparison,levels = design)
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit, contrast = contrast.matrix)
  DEGAll <- lrt$table
  DEGAll$FDR <- p.adjust(DEGAll$PValue, method = "fdr")
  DEGAll$Geneid <- rownames(DEGAll)
  DEGAll <- dplyr::select(DEGAll, Geneid, everything())
  write.table(DEGAll,paste(output,"/",A,"_vs_",B,"_raw_edgeR.xls",sep=""),sep="\t",quote=F,row.names=F)
}


if(judge=="DEseq2"){
  library(DESeq2)
  
  count_file <- count_file[rowMeans(count_file)>=1,]
  counts <- count_file[,c(c(B_list),c(A_list))]
  group <- data.frame("group"=c(rep(0,length(B_list)),rep(1,length(A_list))),"sample"=c(c(B_list),c(A_list)))
  colData <- data.frame(row.names=group$sample,"group"=group$group)
  group <- factor(group$group)

  dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds1 <- nbinomWaldTest(dds)
  # dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)

  sample_normalize <- counts(dds1,normalized=T)
  sample_normalize <- data.frame("Gene_ID"=rownames(sample_normalize),sample_normalize)

  res <- results(dds1)
  res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  res1 <- data.frame("Gene_ID"=rownames(res1),res1)
  res1 <- merge(res1,sample_normalize,by="Gene_ID")
  res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
  
  write.table(res1,paste(output,"/",A,"_vs_",B,"_raw_DEseq2.xls",sep=""),sep="\t",quote=F,row.names = F)
}

