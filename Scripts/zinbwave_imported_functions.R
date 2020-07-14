
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%% The functions negBinTestDESeq2_zinbweights, computeExactWeights, and normDESeq2
#%%% were taken from the repository https://github.com/mcalgaro93/sc2meta showing how 
#%%% to use the ZINB-WAVE method for microbiome data
#%%% I am not using the full file here as it features many additional packages

### performs negative binomial two-sample test of *DESeq2* to detect Diff. Abund. accounting for zero inflation through zinb weights
negBinTestDESeq2_zinbweights <- function(physeq, design = as.formula("~ grp"), IndepFilter = NULL,
                                         normFacts = c("TMM", "RLE", "poscounts", "ratio", "CSS", "UQ", "none", "TSS"), returnDispEsts = FALSE, weights)
{
  register(SerialParam())
  # require(DESeq2)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ### add 1 to zero counts
  # if (any(otu_table(physeq) == 0))
  # {
  #  otu_table(physeq) <- otu_table(physeq) + 0.0001L
  # } else {}
  
  dds <- phyloseq_to_deseq2(physeq, design = design)
  
  normFacts <- paste("NF", normFacts, sep = ".")
  NFs = get_variable(physeq, normFacts)
  
  if(normFacts == "NF.TMM"){
    NFs = NFs *sample_sums(physeq)
  }
  #   if(normFacts %in% c("NF.TMM","NF.TSS")){
  #     NFs = NFs/exp(mean(log(NFs)))
  #   }
  
  sizeFactors(dds) <- NFs/exp(mean(log(NFs))) #This has no impact on the results but facilitates fitting
  
  ### ZINB-WaVE weights
  counts <- as(otu_table(physeq), "matrix")
  weights[which(weights<1e-6)] <- 1e-06
  assays(dds)[["weights"]] = weights
  ### Run DESeq
  
  ddsRes <- DESeq(object = dds, test = "LRT", reduced = ~1, parallel = FALSE)
  dispEsts <- dispersions(ddsRes)
  ddsRes <- results(ddsRes, alpha = 0.05)#, cooksCutoff = FALSE)
  
  ### Independent Filtering, should be before everything
  if(!is.null(IndepFilter))
  {
    toKeep <- ddsRes$baseMean >= IndepFilter & !is.na(ddsRes$pvalue)
    ddsResFilt <- ddsRes[toKeep, ]
    ddsResFilt$padj <- p.adjust(ddsResFilt$pvalue, method = "BH")
    ddsRes <- as(ddsResFilt, "data.frame")
    ddsRes[order(ddsRes$padj), ]
  } else {}
  
  #  ddsRes$id <- rownames(ddsRes)
  pValMat <- as.matrix(ddsRes[, c("pvalue", "padj")])
  colnames(pValMat) <- c("rawP", "adjP")
  statInfo <- cbind("logFC" = ddsRes$log2FoldChange, "LRT" = ddsRes$stat)
  list("pValMat" = pValMat,"dispEsts" = dispEsts, "statInfo" = statInfo)
}# END - function: negBinTestDESeq2 + ZINBWaVE


computeExactWeights <- function (model, x) 
{
  mu <- getMu(model)
  pi <- getPi(model)
  theta <- getTheta(model)
  theta <- matrix(rep(theta, each = ncol(x)), ncol = nrow(x))
  nb_part <- dnbinom(t(x), size = theta, mu = mu)
  zinb_part <- pi * ( t(x) == 0 ) + (1 - pi) *  nb_part
  zinbwg <- ( (1 - pi) * nb_part ) / zinb_part
  zinbwg <- t(zinbwg)
  zinbwg[x > 0] <- 1
  zinbwg[zinbwg < 1e-15] <- 1e-15
  zinbwg
}

### function that apply different normalisations and build *DESeqDataSet* object
### for DESeq2 analysis
normDESeq2 <- function(physeq, whichOTUs = NULL, method = c("poscounts","ratio"))
{
  # require(DESeq2)
  method <- match.arg(method)
  
  ### Coerce count data to vanilla matrix of integers and check if there are zeroes
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ## select which OTUs to analyse
  if (!missing(whichOTUs) || !is.null(whichOTUs))
  {
    physeq <- prune_taxa(taxa_names(physeq)[whichOTUs], physeq)
  } else {}# END - if: whichOTUs
  
  #   otu_table(physeq) <- otu_table(otuTab, taxa_are_rows = TRUE)
  
  ## Calculate size factors
  if (method == "poscounts")
  {
    obj <- phyloseq_to_deseq2(physeq,design = ~grp)
    normFacts <- sizeFactors(DESeq2::estimateSizeFactors(obj,type = "poscounts"))
  } else {
    otuTab <- as(otu_table(physeq), "matrix")
    if (any(otuTab == 0))
    {
      otuTab <- otuTab + 1L
    } else {}
    normFacts <- DESeq2::estimateSizeFactorsForMatrix(otuTab)
  }
  
  physeq@sam_data@.Data <- c(physeq@sam_data@.Data, list(normFacts))
  aux <- physeq@sam_data@names
  aux[length(aux)] <- paste("NF", method, sep = ".")
  physeq@sam_data@names <- aux
  physeq
}# END - function: normDESeq2

#%%% END OF IMPORTED FUNCTIONS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%