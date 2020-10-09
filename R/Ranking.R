## normalize scores in range 0-100
Normalize <- function(x){
  formula <- (x-min(x))*100/(max(x)-min(x))
  return(round(formula, digits = 2))
}

## library-specific gene rankings
ranking <- function(BAGEL_essential, ## set of BAGEL essential genes
                    BAGEL_nonEssential, ## set of BAGEL non-essential genes
                    lib, ## training library (KY or AVANA)
                    param){ ## score (BF or FC)
  es_ranking <- matrix(nrow = length(BAGEL_essential),ncol = 18)
  rownames(es_ranking) <- BAGEL_essential
  colnames(es_ranking) <- paste0(seq(10,95,5),'%')
  
  nones_ranking <- matrix(nrow = length(BAGEL_nonEssential),ncol = 18)
  rownames(nones_ranking) <- BAGEL_nonEssential
  colnames(nones_ranking) <- paste0(seq(10,95,5),'%')
  
  for (i in seq(10,95,5)){
    load(paste0('data/Binary/',lib,'_',param,'_',i,'_es.RData'))
    load(paste0('data/Binary/',lib,'_',param,'_',i,'_nones.RData'))
    load(paste0('data/Training/',lib,'_',param,'_test.RData'))
    
    es_prod <- get(paste0(lib,'_es'))%*%get(paste0(lib,'_test'))[,paste0(i,'%')]
    es_prod <- Normalize(es_prod[,1])
    
    nones_prod <- get(paste0(lib,'_nones'))%*%get(paste0(lib,'_test'))[,paste0(i,'%')]
    nones_prod <- Normalize(nones_prod[,1])
    
    es_ranking[,paste0(i,'%')] <- es_prod
    nones_ranking[,paste0(i,'%')] <- nones_prod
  }
  
  return(list(es_ranking=es_ranking,nones_ranking=nones_ranking))
}