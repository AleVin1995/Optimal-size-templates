## Compute Jaccard similarity between total and subsampled reference genes
Jaccard_Similarity <- function(Corr_FC,
                               es,
                               nones,
                               BAGEL_essential,
                               BAGEL_nonEssential,
                               index=1){
  
  JS <- vector(mode = 'numeric', length = ncol(Corr_FC))
  
  for (i in index){
    ## compute essential genes using the total reference sets
    geneFC <- Corr_FC[,i]
    names(geneFC) <- rownames(Corr_FC)
    ttp <- CRISPRcleanR::ccr.PrRc_Curve(geneFC,BAGEL_essential,BAGEL_nonEssential,FDRth = 0.05,display = F)
    EG_tot <- names(geneFC[which(geneFC<ttp$sigthreshold)])
    EG_tot <- setdiff(EG_tot, union(BAGEL_essential,BAGEL_nonEssential))
    
    ## compute essential genes using the subsampled reference sets
    geneFC_sub <- geneFC[setdiff(names(geneFC), setdiff(union(BAGEL_essential,BAGEL_nonEssential),
                                                        union(es,nones)))]
    ttp_sub <- CRISPRcleanR::ccr.PrRc_Curve(geneFC_sub,es,nones,FDRth = 0.05,display = F)
    EG_sub <- names(geneFC_sub[which(geneFC_sub<ttp_sub$sigthreshold)])
    EG_sub <- setdiff(EG_sub, union(es,nones))
    
    ## Jaccard Similarity
    JS[i] <- length(intersect(EG_sub,EG_tot))/length(union(EG_sub,EG_tot))
  }
  
  JS[is.na(JS)] <- 0
  JS <- JS[index]
  JS <- round(mean(JS), digits = 4)
  
  return(JS)
}

## get top % genes
top_genes <- function(ranking,freq){
  tot <- names(sort(ranking[,paste0(freq,'%')],decreasing = T))
  top <- tot[1:round(length(tot)*freq/100)]
  return(top)
}

## create validation matrix for FCs
matrix_FCs <- function(Corr_FC,
                       final_es,
                       final_nones,
                       BAGEL_essential,
                       BAGEL_nonEssential){
  
  out <- matrix(nrow = 1, ncol = 18, dimnames = list('HT-29', paste0(seq(10,95,5),'%')))
  
  for (i in seq(10,95,5)){
    es <- top_genes(final_es,i)
    nones <- top_genes(final_nones,i)
    
    out[,paste0(i,'%')] <- Jaccard_Similarity(Corr_FC,es,nones,BAGEL_essential,BAGEL_nonEssential)
  }
  return(out*100)
}

## Plot validation results on independent libraries
validation_ind_lib <- function(final_es,
                         final_nones,
                         BAGEL_essential,
                         BAGEL_nonEssential){
  
  ## Brunello library
  Brunello <- read.table('data/Libraries/Brunello_corrected_logFCs.tsv',stringsAsFactors = F, 
                             header = T,row.names = 1) 
  Brunello <- matrix_FCs(Brunello,final_es,final_nones,BAGEL_essential,BAGEL_nonEssential)
  
  ## GeCKOv2 library
  GeCKOv2 <- read.table(paste0('data/Libraries/GeCKOv2_corrected_logFCs.tsv'),stringsAsFactors = F,
                      header = T,row.names = 1)
  GeCKOv2 <- matrix_FCs(GeCKOv2,final_es,final_nones,BAGEL_essential,BAGEL_nonEssential)
  
  ## MinLib library
  MinLib <- read.table(paste0('data/Libraries/MinLib_corrected_logFCs.tsv'),stringsAsFactors = F, 
                       header = T,row.names = 1)
  MinLib <- matrix_FCs(MinLib,final_es,final_nones,BAGEL_essential,BAGEL_nonEssential)
  
  ## Whitehead library
  Whitehead <- read.table(paste0('data/Libraries/Whitehead_corrected_logFCs.tsv'),stringsAsFactors = F,
                          header = T,row.names = 1)
  Whitehead <- matrix_FCs(Whitehead,final_es,final_nones,BAGEL_essential,BAGEL_nonEssential)
  
  
  ## plot results
  par(mfrow=c(2,2), mar=c(2,4,2,2), pty='s')
  
  boxplot(Brunello, ylim = c(0,100),outline = F,frame.plot=F,xaxt='n',ylab = 'Jaccard Similarity %')
  title('Brunello', line = 1)
  axis(1, at = seq(18), labels=F)
  mtext('a', adj = -0.2, side = 3, line = 2, cex = 1)
  boxplot(GeCKOv2, ylim = c(0,100),outline = F,frame.plot=F,yaxt='n',xaxt='n')
  title('GeCKOv2', line = 1)
  axis(1, at = seq(18), labels=F)
  axis(2, labels=F)
  boxplot(MinLib,ylim = c(0,100),outline = F,frame.plot=F,ylab = 'Jaccard Similarity %')
  title('MinLibCas9', line = 1)
  boxplot(Whitehead, ylim = c(0,100), outline = F, frame.plot=F, yaxt='n')
  title('Whitehead', line = 1)
  axis(2, labels=F)
}
  