#######################################
# Capscan Small Intestine Capsule Study 
# 
# Evaluation of capsule opening in  
# different locations along the small
# intestine from proximal to distal
#######################################

# configure data directories
# source base functions
# load libraries
library(here)
library(dplyr)
library(ggplot2)
library(vegan)
library(phyloseq)
library(tidyverse)
library(wesanderson)
library(reshape2)
library(DESeq2)
library(limma)
library(stringr)
library(plyr)
library(rstatix)
library(ggbeeswarm) 
library(Hmisc)
library(ggpubr)
library(cowplot)
library(ggrepel)
library(factoextra)
library(growthcurver)

#--------------------------------------------
# define directories
#--------------------------------------------
raw_data_dir = paste0(here::here(), "/1-data/raw_data/")
clean_data_dir = paste0(here::here(), "/1-data/clean_data/")
fig_dir = paste0(here::here(),"/3-figures/")

#--------------------------------------------
# define clean data paths
#--------------------------------------------

clean_phyloseq = paste0(clean_data_dir, 'clean_mouse_cultures_phlyoseq.RDS')

#--------------------------------------------
# Functions
#--------------------------------------------

source("voom_ihs.R")

# Helper function for converting 0s to 0.001
filterfun1 = function(x){
  x[(x / sum(x)) < (0.001)] <- 0.001
  return(x)
}

# Helper function for flattening a correlation matrix from rcorr
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

filter_subject_prevalence <- function(ps, thresh = 2) {
  smdata <- data.frame(ps@sam_data)
  seqtab <- t(as(ps@otu_table, "matrix"))
  
  subjects_prevalence <- sapply(1:ntaxa(ps), function(i) {
    subj.present <- unique(smdata[seqtab[i, ] > 0, "mouse_subj"])
    length(unique(subj.present))
  })
  ps <- prune_taxa(rownames(seqtab)[subjects_prevalence >= thresh], ps)
  return(ps)
}
            
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


get_scores <- function(pcoa_out, smp_data, axes = 1:2){
  scores <- pcoa_out$vectors[, axes]
  colnames(scores) <- paste0("PC", axes)
  # Combine sample scores and sample data
  sample_scores <- data.frame(scores, stringsAsFactors = FALSE) %>%
    rownames_to_column("SeqName") %>%
    left_join(smp_data) 
  return(sample_scores)
}


get_evals <- function(pcoa_out) {
  evals <- pcoa_out$values[,1]
  var_exp <- 100 * evals/sum(evals)
  return(list("evals" = evals, "variance_exp" = var_exp))
}
