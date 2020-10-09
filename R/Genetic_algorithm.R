## Compute Jaccard similarity between total and subsampled reference genes
Jaccard_Similarity <- function(Corr_FC,
                               es,
                               nones,
                               BAGEL_essential,
                               BAGEL_nonEssential,
                               index){
  
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

## Shuffle fixed % of subsampled genes with held-out reference genes
Shuffle <- function(es,
                    nones,
                    BAGEL_essential,
                    BAGEL_nonEssential,
                    swap){
  
  ## nº subsampled genes to be swapped
  swap_es <- ceiling(length(es)*swap)
  swap_nones <- ceiling(length(nones)*swap)
  
  ## essential/non-essential reference genes repertoire to pick from
  pick_es <- setdiff(BAGEL_essential,es)
  pick_nones <- setdiff(BAGEL_nonEssential,nones)
  
  ## in case the predefined nº genes to be swapped is greater than the remaining ones
  ## this happens for high % of subsampling
  if (swap_es > length(pick_es) || swap_nones > length(pick_nones)){
    swap_es <- length(pick_es)
    swap_nones <- length(pick_nones)
  }
  
  ## shuffle subsampled reference genes
  es <- sample(es)
  nones <- sample(nones)
  
  es[1:swap_es] <- sample(pick_es, size = swap_es)
  nones[1:swap_nones] <- sample(pick_nones, size = swap_nones)
  
  return(list(essential = es, nonessential = nones))
}

## optimize essential and non-essential reference genes at a fixed subsampling %
Training <- function(Corr_FC, # genome-wide CRISPR-Cas9 screening matrix (values provided as logFCs or BFs)
                     BAGEL_essential,
                     BAGEL_nonEssential,
                     freq, # subsampling %
                     steps=15, # stop training nº consecutive iterations without improvement
                     swap=25, # % of subsampled essential/non-essential genes to be swapped
                     split_data=80, # % of cell lines to be used in the training phase
                     seed){ # set seed
  
  set.seed(seed)
  
  freq <- round(freq/100,digits = 2)
  swap <- round(swap/100,digits = 2)
  split_data <- round(split_data/100,digits = 2)
  
  ## remove non-screened genes from the reference sets
  screened_genes <- rownames(Corr_FC)
  BAGEL_essential <- intersect(BAGEL_essential,screened_genes)
  BAGEL_nonEssential <- intersect(BAGEL_nonEssential,screened_genes)
  
  ## initial subsampled reference sets
  es <- sample(BAGEL_essential,size = round(length(BAGEL_essential)*freq))
  nones <- sample(BAGEL_nonEssential,size = round(length(BAGEL_nonEssential)*freq))
  
  ## select training and test set
  cells <- seq(ncol(Corr_FC))
  training_set <- sample(cells,size = round(length(cells)*split_data))
  test_set <- setdiff(cells,training_set)
  
  ## training phase
  JS_training <- Jaccard_Similarity(Corr_FC,es,nones,BAGEL_essential,BAGEL_nonEssential,training_set)
  iteration <- 1
  GO <- T
  
  cat(paste0('Iteration ',iteration,': ',tail(JS_training,n=1),'\n'))
  
  while (GO){
    newSub <- Shuffle(es,nones,BAGEL_essential,BAGEL_nonEssential,swap=swap)
    newJS <- Jaccard_Similarity(Corr_FC,newSub$essential,newSub$nonessential,
                                BAGEL_essential,BAGEL_nonEssential,training_set)
    
    if (newJS>tail(JS_training,n=1)){
      es <- newSub$essential
      nones <- newSub$nonessential
      JS_training <- c(JS_training,newJS)
    } else {
      JS_training <- c(JS_training,tail(JS_training,n=1))
    }
    
    iteration <- iteration+1
    
    if (iteration%%5==0){
      cat(paste0('Iteration ',iteration,': ',tail(JS_training,n=1),'\n'))
    }
    
    if (length(JS_training)>=steps){
      if (tail(JS_training, n=1)==JS_training[length(JS_training)-steps+1]){
        GO <- F
        cat(paste0('Iteration ',iteration,': ',tail(JS_training,n=1),'\n'))
      }
    }
  }
  
  ## test phase
  JS_test <- Jaccard_Similarity(Corr_FC,es,nones,BAGEL_essential,BAGEL_nonEssential,test_set)
  
  return(list(JS_training=tail(JS_training,n=1),JS_test=JS_test,training_set=sort(training_set),
              test_set=sort(test_set),essential=es,nonessential=nones))
}