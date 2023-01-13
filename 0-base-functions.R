###################################
# Base functions for analysis
###################################

# [TEMPLATE FOR NEW FUNCTIONS]
# Documentation: 
# Usage: 
# Description: 
# Args/Options: 
# 
# Returns: 
# Output: ...
# [TEMPLATE FOR NEW FUNCTIONS]

phyloseqTransformation <- function(ps, 
                                   c = .4){
  dds <- phyloseq_to_deseq2(ps, design = ~ 1)
  dds <-  estimateSizeFactors(dds, type = "poscounts")
  dds <- estimateDispersions(dds, fitType = "local")
  abund <- assay(dds)
  abund_temp <- matrix(nrow = nrow(abund), ncol = ncol(abund))
  
  # Column co-factor
  
  for(i in 1:dim(abund)[1]){
    for(j in 1:dim(abund)[2]){
      abund_temp[i,j] <- abund[i,j]/sizeFactors(dds)[j]
    }
  }
  
  # Row co-factor
  ind <- which(dispersions(dds) > 2)
  ind2 <- which(dispersions(dds) <= 2)
  
  for(i in ind){
    for(j in 1:dim(abund_temp)[2]){
      if(abund_temp[i,j] > 0){
        abund_temp[i,j] <- sqrt(dispersions(dds)[i]- .5)*asinh(sqrt((abund_temp[i, j]+ c)/(dispersions(dds)[i]-2*c)))
      }else{
        abund_temp[i,j] <- sqrt(dispersions(dds)[i]- .5)*asinh(sqrt(abund_temp[i, j]))
      }
      
    }
  }
  
  if(length(ind2) > 0){
    for(i in ind2){
      for(j in 1:dim(abund_temp)[2]){
        if(abund_temp[i,j] > 0){
          abund_temp[i,j] <- log2(abund_temp[i, j] + .5*dispersions(dds)[i])
        }
        
      }
    }
  }
  
  rownames(abund_temp) <- rownames(abund)
  colnames(abund_temp) <- colnames(abund)
  
  ps <- phyloseq(otu_table(abund_temp, taxa_are_rows = TRUE),
                 sample_data(ps), 
                 tax_table(ps), 
                 phy_tree(ps))
  
  return(ps)
  
}

get_scores <- function(pcoa_out, smp_data, axes = 1:2){
  scores <- pcoa_out$vectors[, axes]
  colnames(scores) <- paste0("PC", axes)
  # Combine sample scores and sample data
  sample_scores <- data.frame(scores, stringsAsFactors = FALSE) %>%
    rownames_to_column("id_16s") %>%
    left_join(smp_data) 
  return(sample_scores)
}



get_evals <- function(pcoa_out) {
  evals <- pcoa_out$values[,1]
  var_exp <- 100 * evals/sum(evals)
  return(list("evals" = evals, "variance_exp" = var_exp))
}

plot_ordination <- function(
  sample_scores, var_exp, group_column = "Subject", 
  color_continuous = FALSE, color_column = NULL, colors = NULL, shape = 21, 
  ax_names = c("PC1", "PC2")){
  if(is.null(color_column)){
    color_column <- "Color_Column"
    sample_scores$Color_Column <- "Color"
  }
  if(is.null(colors)){
    nCols <-length(unique(sample_scores[[color_column]]))
    colors <- RColorBrewer::brewer.pal(9, "Set1")
    colors <- colorRampPalette(colors)(nCols)
  }
  plt <- ggplot(
    sample_scores, 
    aes_string(
      x = ax_names[1], y = ax_names[2], fill = color_column)) +
    # geom_polygon(aes_string(group = group_column)) +
    # geom_path(aes_string(group = group_column, color = color_column)) +
    geom_point(size = 3, color = "grey10", shape = shape) +
    xlab(sprintf("PC1 [%s%% variance]", round(var_exp[1], 2))) +
    ylab(sprintf("PC2 [%s%% variance]", round(var_exp[2], 2))) +
    coord_fixed() +
    theme(text = element_text(size = 20)) 
  if (!color_continuous){ # if continuous
    plt <- plt + scale_fill_manual(values = colors, na.value = "grey") +
      scale_color_manual(values = colors, na.value = "grey") 
  }
  return(plt)
}




colorize <- function(x, color) {
  if (knitr::is_latex_output()) {
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else if (knitr::is_html_output()) {
    sprintf("<span style='color: %s;'>%s</span>", color,
            x)
  } else x
}

#--------------------------------------------
# define function for getting a list of species from a specific location/sample type
#--------------------------------------------
# Usage: used for making venn diagrams of species overlap between sample locations
# Description: Takes in a data.frame with variables 'species_id' and 'location' and filters 
#             the data.frame to obtain a list of 'species_ids' contained within the specified 'location'. 
# Args/Options: 
#   df: data.frame with variables named 'species_id' and 'location'
#   loc: location from which to obtain species list (should be one of the values contained in the variable 'location')
#   species_coverage: coverage within a specific sample required in order to keep in the species list.
# Returns:
  # A list of the 'species_ids' contained within the specified 'location'.
  
pull_species <- function(df, loc, species_coverage = 0) {
  df %>%
    filter(coverage > species_coverage) %>%
    filter(location == loc) %>%
    pull(species_id)
}

#--------------------------------------------
# define function for running a topic model on a phyloseq object
#--------------------------------------------
# Usage: 
# Description: Takes in a phyloseq object and computes k topics for 
#
# Args/Options: 
# physeq: phyloseq object
# k: number of topics; default is 15
# method: algorithm for computing the latent dirichlet allocation; default is variational expectation-maximization (VEM) 
# 
# Returns:
# A 


physeq_topic_model <- function(physeq, k = 15, method = "VEM") {
  dtm <- as.simple_triplet_matrix(physeq@otu_table)
  VEM_model = LDA(dtm, k = k, method = method, control = list(alpha = 0.1))
  return(VEM_model)
}

###########################
### Flattens correlation matrix
###########################

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat, diag = F)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

###########################
### Functions to easily manipulate correlation output.
### Copied from https://www.khstats.com/blog/corr-plots/corr-plots/
###########################
cors <- function(df) {
  M <- Hmisc::rcorr(as.matrix(df))
  Mdf <- map(M, ~data.frame(.x))
  return(Mdf)
}

formatted_cors <- function(df){
  cors(df) %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, "measure2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    rename(p = P) %>%
    mutate(sig_p = ifelse(p < .05, T, F),
           p_if_sig = ifelse(sig_p, p, NA),
           r_if_sig = ifelse(sig_p, r, NA)) 
}


###########################
### Splits phyloseq object
###########################
phyloseq_sep_variable <- function(physeq, variable, drop_zeroes = T){
  
  # require(phyloseq)
  # require(plyr)
  
  ## Check the input
  if(is.null(phyloseq::sample_data(physeq, errorIfNULL = F))){
    stop("Sample data is missing in the phyloseq-object.\n")
  }
  
  ## Extract sample meta-data
  mtd <- as(object = phyloseq::sample_data(physeq), Class = "data.frame")
  
  if(!variable %in% colnames(mtd)){
    stop("Grouping variable is missing from the sample data of phyloseq-object.\n")
  }
  
  if(class(mtd[, variable]) %in% c("integer", "numeric") ){
    if( length( unique(mtd[, variable]) ) > 5){
      stop("Groupping variable is numeric and it has too many levels. Consider transforming it to factor.\n")
    } else {
      warning("Groupping variable is numeric and it was coerced to factor.\n")
      mtd[, variable] <- factor(mtd[, variable])
    }
  }
  
  if(length(table(mtd[, variable])) == 1){
    cat("Warning: there is only one group of samples in the resulting list.\n")
  }
  
  ## Add sample IDs to the meta-data
  smp <- data.frame(
    SID = phyloseq::sample_names(physeq),
    mtd,
    stringsAsFactors = F)
  
  ## Extract sample names by the specified variable
  svv <- plyr::dlply(.data = smp, .variables = variable, .fun = function(z){ z$SID })
  
  ## Extract samples by groupping variable
  res <- plyr::llply(.data = svv, .fun = function(z){ phyloseq::prune_samples(z, x = physeq) })
  
  ## Remove taxa with zero abundance
  if(drop_zeroes == TRUE){
    res <- plyr::llply(.data = res, .fun = function(x){ phyloseq::prune_taxa(phyloseq::taxa_sums(x) > 0, x) })
  }
  
  return(res)
}



#--------------------------------------------
# define function for summarizing the number of unique bsh genes per sample
#--------------------------------------------
# Usage: 
# Description: Takes in a data.frame object, groups by subject, calculates summary statistics, and generates a plot
#
# Args/Options: 
# df: a data.frame

# Returns:
# A list containing: 
# df.summary: the grouped, summarized data.frame
# plot: the ggplot object 
unique_genes_per_sample <- function(df){
  df_summary <- df %>%
    group_by(Subject, Set, Type, meta_samplename, location) %>%
    summarise(n = n(), 
              n_norm = n/mean(reads_meta_filtered)*median_read_depth) %>%
    # filter(Set %in% c("Stool", "1", "2", "3", "4", "5")) %>%
    mutate(loc = ifelse(Set == "Stool", "Stool", "Small Intestine"),
           Subject = factor(Subject, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)))
  
  plot <- ggplot(df_summary,
                 aes(x = Subject, y = n_norm, color = loc)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_point(position = position_jitterdodge(0.2), alpha = 0.4) + 
    # ggsignif::geom_signif(comparisons = list(c("Small intestine", "Stool"))) +
    scale_color_manual(values = CapAndStoolColors) +
    labs(y = "Unique bsh genes (normalized by sequencing depth)", x = "Subject", color = "Location")
  
  plot
  
  return(list(df.summary = df_summary, plot = plot))
}



# tip agglomeration based on heirarchial clustering of the phylogenetic distance
clust_tip_glom <- function (physeq, h = 0.2, hcfun = cluster::agnes, ...) {
  dd <- as.dist(ape::cophenetic.phylo(phy_tree(physeq)))
  psclust <- cutree(as.hclust(hcfun(dd, ...)), h = h)
  return(psclust)
}

merge_tips <- function(physeq, psclust) {
  taxtab_collapse <- physeq@tax_table %>%
    as.data.frame() %>%
    rownames_to_column("SeqName") %>%
    select(SeqName, Genus, Species) %>%
    mutate(Tip = paste0("Tip", psclust)) %>%
    group_by(Tip) %>%
    mutate( 
      Org = paste(Genus, Species),
      Org = ifelse(is.na(Genus) & is.na(Species), NA, Org)
    ) %>%
    dplyr::summarise(
      OrgIncluded = paste0(unique(Org[!is.na(Org)]), collapse = "/"),
      SeqIncluded = paste0(SeqName[!is.na(SeqName)], collapse = "/"),
      Freq = n()
    ) %>%
    mutate(
      OrgIncluded = ifelse(OrgIncluded == "", "Unknown", OrgIncluded)
    )
  physeq@tax_table <- tax_table(
    cbind(Tip = paste0("Tip", psclust), physeq@tax_table[, 1:7]))
  freqtab <- data.frame(table(psclust))
  freqtab$Tip <-  paste0("Tip", freqtab$pclust)
  cliques <- freqtab %>% filter(Freq > 1) %>% 
    .[["psclust"]] %>% as.character()
  for (i in cliques) {
    physeq <- merge_taxa(physeq, 
                         eqtaxa = names(psclust)[psclust == i])
  }
  taxtab <- physeq@tax_table %>%
    as.data.frame() %>%
    left_join(taxtab_collapse) %>%
    column_to_rownames("Tip")
  physeq@tax_table <- tax_table(as(taxtab, "matrix"))
  taxa_names(physeq) <- rownames(taxtab)
  return(physeq)
} 


# filter out taxa not present in at least 'thresh' distinct subjects.
filter_subject_prevalence <- function(ps, thresh = 2) {
  smdata <- data.frame(ps@sam_data)
  seqtab <- t(as(ps@otu_table, "matrix"))
  
  subjects_prevalence <- sapply(1:ntaxa(ps), function(i) {
    subj.present <- unique(smdata[seqtab[i, ] > 0, "Subject"])
    length(unique(subj.present))
  })
  ps <- prune_taxa(rownames(seqtab)[subjects_prevalence >= thresh], ps)
  return(ps)
}

# limma DA
limma_fit <- function(seqtab, smpdata, 
                      dsgn, sizefac, block, 
                      taxtable, alpha = 0.05){
  mm <- model.matrix(dsgn, smpdata)
  colnames(mm) <- make.names(colnames(mm))
  v <- voom_ihs(seqtab, mm, plot = FALSE, 
                lib.size = sizefac)
  corfit <- duplicateCorrelation(v, mm, block = block)
  v <- voom_ihs(seqtab, mm, plot = FALSE, block = block,
                correlation = corfit$consensus,
                lib.size = sizefac)
  fit <- lmFit(v, mm, block = block, correlation = corfit$consensus)
  fit.ebayes <- eBayes(fit)
  res <-topTable(fit.ebayes, adjust="BH", n = Inf, coef = 2,
                   p.value = alpha) %>%
    rownames_to_column("SeqName") %>%
        left_join(taxtable)
  return(res)
}

## Voom_ihs
library(limma)
voom_ihs <- function (counts, design = NULL, lib.size = NULL, normalize.method = "none", 
                      span = 0.5, plot = FALSE, save.plot = FALSE, ...) 
{
  ihs <- function(x) {
    transformed <- log(x + sqrt(x ^ 2 + 1))
    return(transformed)
  }
  
  hs <- function(x) {
    y <- 0.5*exp(-x)*(exp(2*x)-1)
    return(y)
  }
  
  out <- list()
  if (is(counts, "DGEList")){
    out$genes <- counts$genes
    out$targets <- counts$samples
    if (is.null(design) && diff(range(as.numeric(counts$sample$group))) > 
        0) 
      design <- model.matrix(~group, data = counts$samples)
    if (is.null(lib.size)) 
      lib.size <- with(counts$samples, lib.size * norm.factors)
    counts <- counts$counts
  }else{
    isExpressionSet <- suppressPackageStartupMessages(is(counts, 
                                                         "ExpressionSet"))
    if (isExpressionSet) {
      if (length(Biobase::fData(counts))) 
        out$genes <- Biobase::fData(counts)
      if (length(Biobase::pData(counts))) 
        out$targets <- Biobase::pData(counts)
      counts <- Biobase::exprs(counts)
    }
    else {
      counts <- as.matrix(counts)
    }
  }
  
  n <- nrow(counts)
  
  if (n < 2L){stop("Need at least two genes to fit a mean-variance trend")}
  
  if (is.null(design)) {
    design <- matrix(1, ncol(counts), 1)
    rownames(design) <- colnames(counts)
    colnames(design) <- "GrandMean"
  }
  
  if (is.null(lib.size)){lib.size <- colSums(counts)}
  
  y <- t(ihs(t(counts)/lib.size))
  
  y <- normalizeBetweenArrays(y, method = normalize.method)
  
  fit <- lmFit(y, design, ...)
  
  if (is.null(fit$Amean)){fit$Amean <- rowMeans(y, na.rm = TRUE)} 
  
  #   fitted mean for each row (taxa)    
  sx <- fit$Amean 
  #   fitted standard deviation
  sy <- sqrt(fit$sigma)
  
  allzero <- rowSums(counts) == 0
  
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  
  l <- lowess(sx, sy, f = span)
  
  if (plot){
    plot(sx, sy, xlab = "mean(ihs(count/sizeFactor))", ylab = "Sqrt(standard deviation )", 
         pch = 16, cex = 0.25)
    title("voom_ihs: Mean-variance trend")
    lines(l, col = "red")
  }
  
  #   linearly interpolate given the smooth line from the lowess
  f <- approxfun(l, rule = 2)
  
  if(fit$rank < ncol(design)){
    j <- fit$pivot[1:fit$rank]
    fitted.values <- fit$coef[, j, drop = FALSE] %*% t(fit$design[, 
                                                                  j, drop = FALSE])
  }else{
    fitted.values <- fit$coef %*% t(fit$design)
  }
  
  fitted.cpm <- hs(fitted.values)
  
  fitted.count <- t(t(fitted.cpm)*(lib.size))
  
  fitted.ihs <- ihs(fitted.count)
  
  w <- 1/f(fitted.ihs)^4
  
  dim(w) <- dim(fitted.ihs)
  
  out$E <- y
  
  out$weights <- w
  
  out$design <- design
  
  if(is.null(out$targets)){
    out$targets <- data.frame(lib.size = lib.size)
  }else{
    out$targets$lib.size <- lib.size
  }
  
  if (save.plot) {
    out$voom.xy <- list(x = sx, y = sy, xlab = "mean(ihs(count/sizeFactor))", 
                        ylab = "Sqrt( standard deviation )")
    out$voom.line <- l
  }
  
  new("EList", out)
}


#We can plot normalized and transformed counts to see the trends.
get_plot_data <- function(seqtab, smpdata, size_factors, limma_res){
  sf <- size_factors[colnames(seqtab)]
  cnts_trans <- asinh(sweep(seqtab, 2, sf, "/"))
  
  cnts_trans.long <- data.frame(cnts_trans) %>%
    rownames_to_column("OrgName") %>%
    pivot_longer(cols = starts_with("S"), names_to = "id_16s", values_to = "count") %>%
    left_join(smpdata) %>%
    rename(SeqName = OrgName)
  
  DF <- cnts_trans.long %>%
    filter(SeqName %in% limma_res$SeqName) %>%
    left_join(limma_res)
  
  return(DF)
}


plot_percent <- function(DF, sbj_markersize = 3, smp_markersize = 1.5){
  subjDF <- DF %>%
    dplyr::select(
      SeqName, OrgName, logFC, adj.P.Val, 
      Subject, location, count) %>%
    group_by(SeqName, OrgName, logFC, adj.P.Val, Subject, 
             location) %>% 
    dplyr::summarise(count = mean(count))
  
  logfc_lim <- max(abs(DF$logFC))
  ggplot(
    DF,
    aes(x = location, y = count)) +
    geom_point(
      color = "grey55",
      size = smp_markersize, alpha = 0.5
    ) +
    geom_point(
      data = subjDF,
      aes(color = logFC), 
      size = sbj_markersize, shape = 17
    ) +
    scale_color_gradientn(
      colors = colorRampPalette(brewer.pal(9, "RdBu")[-(4:6)])(100),
      limits = c(-logfc_lim, logfc_lim)
    ) + 
    ylab("transformed counts") +
    xlab("Percent weight loss at 12 months [%]") +
    facet_wrap(~OrgName, scales = "free")
}


#--------------------------------------------
# Functions for correlating ba concentrations with ASVs
#--------------------------------------------


plot_ba_corr <- function(bile_acid,df2plot){

df2plot$bile_acid<-df2plot[,bile_acid]
df2plot <- df2plot %>% filter(bile_acid > 1) %>% filter(Abundance > 0)
# split the data by ASV
B <- split(df2plot, df2plot$OTU)
    
# Calculate the correlation in all data.frames using lapply 
M <- lapply(B, function(x) cor.test(x$Abundance, log10(x$bile_acid), method='pearson'))
            
# Get ASVs that are significant
pvals<-lapply(M, function(x) x$p.value)
pvals<-as.data.frame(pvals)
pvals.corrected<-p.adjust(pvals,method="hochberg") %>% as.data.frame() %>% filter(. < 0.05)
sig<-rownames(pvals.corrected)

# Slope positive for conjugated
#r<-lapply(M, function(x) x$estimate)         
#filtered<-data.frame(r) %>% t() %>% as.data.frame() %>% filter(cor > 0)
#sig.r<-rownames(filtered)         
    
# Filter dataframe to only include significant ASVs              
df2plot.filt <- df2plot %>%
    filter(OTU %in% sig) #%>%
#    filter(OTU %in% sig.r)

if (dim(df2plot.filt)[1] > 0) {
    a<-ggplot(df2plot.filt, aes(x=log10(bile_acid), y=(Abundance))) + 
        geom_point()  + geom_smooth(method=lm, se=FALSE) +
        labs(title=bile_acid,y=expression('log'[2]*' ASV count'), x=paste0(expression('log'[10]*' concentration of '),bile_acid)) +
        stat_cor(method = "pearson") +
        facet_wrap(~Family+Genus+OTU)+
        theme(strip.text = element_text(size=12),
        axis.text = element_text(size = 12)) 
    return(a)
} else {return(NULL)}
              
}
   
ba_corr_dataframe <- function(bile_acid,df2plot){

df2plot$bile_acid<-df2plot[,bile_acid]
df2plot <- df2plot %>% filter(bile_acid > 1)  %>% filter(Abundance > 0)
# split the data by ASV
B <- split(df2plot, df2plot$OTU)
    
# Calculate the correlation in all data.frames using lapply 
M <- lapply(B, function(x) cor.test(x$Abundance, log10(x$bile_acid), method='pearson'))
            
# Get ASVs that are significant
pvals<-lapply(M, function(x) x$p.value)
pvals<-as.data.frame(pvals)
pvals.corrected<-p.adjust(pvals,method="hochberg") %>% as.data.frame() %>% filter(. < 0.05)
sig<-rownames(pvals.corrected)
sig
# Slope 
r<-lapply(M, function(x) x$estimate)       
filtered<-data.frame(r) %>% t() %>% as.data.frame() 
filtered$OTU<-rownames(filtered)
sig.r <- filtered %>% filter(OTU %in% sig)

# Filter dataframe to only include significant ASVs              
df2plot.filt <- df2plot %>%
    filter(OTU %in% sig) %>%
    select(OTU, Phylum, Family, Genus, Species) %>%
    unique() %>%
    mutate(BA = bile_acid) %>%
    left_join(.,sig.r, by='OTU')

if (dim(df2plot.filt)[1] > 0) {
    return(df2plot.filt)
} else {return(NULL)}
              
}
  

plot_ba_corr_merge <- function(bile_acid,df2plot){


tmp<-df2plot.final %>% select(bile_acid) %>% mutate(bile_acid=rowSums(.))
df2plot$bile_acid<-tmp$bile_acid


df2plot <- df2plot %>% filter(bile_acid > 1)  %>% filter(Abundance > 0)
    
# split the data by ASV
B <- split(df2plot, df2plot$OTU)
    
# Calculate the correlation in all data.frames using lapply 
M <- lapply(B, function(x) cor.test(x$Abundance, log10(x$bile_acid), method='pearson'))
            
# Get ASVs that are significant
pvals<-lapply(M, function(x) x$p.value)
pvals<-as.data.frame(pvals)
pvals.corrected<-p.adjust(pvals,method="hochberg") %>% as.data.frame() %>% filter(. < 0.05)
sig<-rownames(pvals.corrected)

# Slope positive for conjugated
#r<-lapply(M, function(x) x$estimate)         
#filtered<-data.frame(r) %>% t() %>% as.data.frame() %>% filter(cor > 0)
#sig.r<-rownames(filtered)         
    
# Filter dataframe to only include significant ASVs              
df2plot.filt <- df2plot %>%
    filter(OTU %in% sig) #%>%
#    filter(OTU %in% sig.r)

if (dim(df2plot.filt)[1] > 0) {
    a<-ggplot(df2plot.filt, aes(x=log10(bile_acid), y=(Abundance))) + 
        geom_point()  + geom_smooth(method=lm, se=FALSE) +
        labs(title=paste(bile_acid, collapse=' and '),y=expression('log'[2]*' ASV count'), x=paste0(expression('log'[10]*' concentration of '),bile_acid)) +
        stat_cor(method = "pearson") +
        facet_wrap(~Family+Genus+OTU)+
        theme(strip.text = element_text(size=12),
        axis.text = element_text(size = 12)) 
    return(a)
} else {return(NULL)}
              
}
   

ba_corr_dataframe_merge <- function(bile_acid,df2plot){

tmp<-df2plot.final %>% select(bile_acid) %>% mutate(bile_acid=rowSums(.)) 
df2plot$bile_acid<-tmp$bile_acid


df2plot <- df2plot %>% filter(bile_acid > 1, Abundance > 0)
# split the data by ASV
B <- split(df2plot, df2plot$OTU)
    
# Calculate the correlation in all data.frames using lapply 
M <- lapply(B, function(x) cor.test(x$Abundance, log10(x$bile_acid), method='pearson'))
            
# Get ASVs that are significant
pvals<-lapply(M, function(x) x$p.value)
pvals<-as.data.frame(pvals)
pvals.corrected<-p.adjust(pvals,method="hochberg") %>% as.data.frame() %>% filter(. < 0.05)
sig<-rownames(pvals.corrected)

# Slope 
r<-lapply(M, function(x) x$estimate)       
filtered<-data.frame(r) %>% t() %>% as.data.frame() 
filtered$OTU<-rownames(filtered)
sig.r <- filtered %>% filter(OTU %in% sig)
          
new_name<-paste(bile_acid, collapse=' and ')

# Filter dataframe to only include significant ASVs              
df2plot.filt <- df2plot %>%
    filter(OTU %in% sig) %>%
    select(OTU, Phylum, Family, Genus, Species) %>%
    unique() %>%
    mutate(BA = new_name) %>%
    left_join(.,sig.r, by='OTU')

if (dim(df2plot.filt)[1] > 0) {
    return(df2plot.filt)
} else {return(NULL)}
              
}




## custom function to rarefy to a subject's minimum sequencing depth, 
## with functionality to record repeated iterations (i)
## This function is used in 'ExtendedData_Figure_3.R'
bySubj_rarefy <- function(s, i, phyloseq_object) {
   s = "1"
  ps_subj <- phyloseq_object %>%
    prune_samples(phyloseq_object@sam_data$Subject == s, .)
  
  min_sample_depth = min(sample_sums(ps_subj))
  
  ps_evendepth <- rarefy_even_depth(ps_subj, 
                                    sample.size = min_sample_depth, 
                                    rngseed = FALSE, 
                                    replace = T, 
                                    trimOTUs = T, 
                                    verbose = F)
  alpha_div_sample <- estimate_richness(ps_evendepth, measures = "Shannon") %>%
    mutate(Subject = s,
           level = "sample", 
           read_depth = min_sample_depth) %>%
    rownames_to_column("ID")
  
  ps_sampleType <- merge_samples(ps_evendepth, "Type", fun = sum)
  alpha_div_type <- estimate_richness(ps_sampleType, measures = "Shannon") %>%
    mutate(Subject = s, 
           level = "type",
           read_depth = min_sample_depth) %>%
    rownames_to_column("ID")
  
  ps_sampleSet <- merge_samples(ps_evendepth, "Set", fun = sum)
  alpha_div_set <- estimate_richness(ps_sampleSet, measures = "Shannon") %>%
    mutate(Subject = s, 
           level = "set",
           read_depth = min_sample_depth) %>%
    rownames_to_column("ID")
  
  alpha_div_allSubj <- alpha_div_sample %>%
    bind_rows(alpha_div_type) %>%
    bind_rows(alpha_div_set) %>%
    mutate(i = i)
  
  return(alpha_div_allSubj)
}





