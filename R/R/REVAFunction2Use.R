#' @import Rcpp
# import description ends
NULL

#' @title REVA.test
#' @description The main function of the REVA package. Calculates REVA statistics for a set of scalar values for each point (x) and the distances between points (mydist) 
#' @param x vector of scalar values associated with the points
#' @param mydist the distances between the points
#' @param sampleNum2EstVar if the number of samples (e.g. \code{lengths(x)}) is larger or equal to this value, the variance will be estimated by subsampling to \code{sampleNum2EstVar}
#' @param permsN number of subsampling estimations to average
#' @return a list with \code{reva}, \code{var} and \code{pvalue} fields for the REVA statistics values, the variance estimation and the p-value of REVA test
#' @export
REVA.test <- function(x, mydist, sampleNum2EstVar = 100, permsN = 10){
  
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
              var = varapp3,
              pvalue = pvalue))
}


### REVA matx must contain samples in the columns 
### and x must contain the corresponding scalars in the vector
## dist function must return a symmetric matrix with ncol = ncol matx 
## representing the distance between all pair of samples
#' @importFrom bioDist tau.dist
#
#

#' @export
REVA.user.distance <- function(  x, matx, 
                                 distfunc = function(m) 
                                   as.matrix(bioDist:::tau.dist(t(m))), 
                                 ...){
  
  intername <- intersect(names(x),colnames(matx))
  
  mydist <- distfunc( matx[,intername],...)
  
  
  
  revaval <- REVA.p.value( x[intername] , mydist[intername,intername] )
  revaval$dist <- mydist
  
  return(revaval)
  
}

#' @export
REVA.user.distance.Geneset <- function( x, matxAll, genesetList = NULL,...){
  
  if(is.null(genesetList)){
    geneList <- list(AllGenes = rownames(matxAll))
  }
  
  z <- sapply( genesetList , FUN = function(y) REVA.user.distance(x,matxAll[y,],...) )
  
  return(z);
  
}

