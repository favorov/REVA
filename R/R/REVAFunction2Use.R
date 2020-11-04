#' @import Rcpp
#' @importFrom GSReg GSReg.kendall.tau.distance
##  ' @import Hmisc
##  ' @inport utils
#
#

REVA.p.value <- function(x,mydist, sampleNum2EstVar = 100, permsN = 10){
  
  reva <- AfsariCorcpp(x,mydist)
  if( nrow(mydist)  < sampleNum2EstVar ){
    medianCount <- AfsariCorVarcpp(mydist)
    medianCount2 <- AfsariCorVarcpp2(mydist)
    
    alpha <- medianCount/(3*choose(nrow(mydist),3)*choose(nrow(mydist)-2,2))
    
  }else{
    medianCountAll <- vector("list",length = permsN)
    medianCount2All <- medianCountAll
    
    for( i in 1:permsN){
      mysamples <-  sample(1:nrow(mydist), sampleNum2EstVar)
      
      medianCountAll[[i]] <- AfsariCorVarcpp(mydist[mysamples,mysamples])
      medianCount2All[[i]] <- AfsariCorVarcpp2(mydist[mysamples,mysamples])
    }
    
    medianCount <- Reduce("+",medianCountAll)/permsN
    medianCount2 <- Reduce("+",medianCount2All)/permsN
    
    alpha <- medianCount/(3*choose(sampleNum2EstVar,3)*choose(sampleNum2EstVar-2,2))
  }
  
  beta <- (1-3*alpha)/6
  gamma <- (1+3*alpha)/12
  
  covest <- (alpha*2/15)+(1/10)*beta*4+  (7/60)* gamma * 4 - 1/9
  covtest2 <- sum((rbind(c(2,1,1),c(1,2,1),c(1,1,2))/12 )*
                    medianCount2/sum(medianCount2))-1/9
  varapp <- covest*9/nrow(mydist)
  varapp2 <- varapp + covtest2*18/nrow(mydist)^2
  varapp3 <- varapp2 + 4/3/nrow(mydist)^3
  pvalue <- pnorm(reva ,mean = 1/3, sd = sqrt(max(varapp3,1e-6)),lower.tail = F) 
  return(list(reva = reva,
              varapp = varapp,
              varapp2 = varapp2,
              varapp3 = varapp3,
              medianCount = medianCount,
              medianCount2 = medianCount2,
              pvalue = pvalue))
  
}


### REVA matx must contain samples in the columns 
### and x must contain the corresponding scalars in the vector
## dist function must return a symmetric matrix with ncol = ncol matx 
## representing the distance between all pair of samples

REVA.user.distance <- function(  x, matx, 
                                 distfunc = GSReg:::GSReg.kendall.tau.distance, ...){
  
  intername <- intersect(names(x),colnames(matx))
  
  mydist <- distfunc( matx[,intername],...)
  
  
  
  revaval <- REVAPvalue( x[intername] , mydist[intername,intername] )
  revaval$dist <- mydist
  
  return(revaval)
  
}

REVA.user.distance.Geneset <- function( x, matxAll, genesetList = NULL,...){
  
  if(is.null(genesetList)){
    geneList <- list(AllGenes = rownames(matxAll))
  }
  
  z <- sapply( genesetList , FUN = function(y) REVA.user.distance(x,matxAll[y,],...) )
  
  return(z);
  
}

