
#' Calculate test statistic for ANCOM
#' 
#' @description
#' Calculates test statistics for differences in OTU abundances between treatment groups.
#' 
#' @param otu_data the OTU dataset.
#' @param n_otu the number of OTUs.
#' @param alpha the significance level at which the tests are to be performed.
#' @param multcorr type of correction for multiple comparisons, see Details.
#' @param Wexact logical, should Wilcoxon tests return exact p-values?
#' @param ncore if ncore>1, then \pkg{doParallel} will be loaded and used.
#' 
#' @details
#' \code{multcorr} can take on values of 1 (no correction), 2 (a less stringent)
#' correction, or 3 (a more stringent correction).
#' 
#' @note
#' This function is intended to be called by \code{\link{ANCOM}}, see the documentation of
#' that function for details on using the method.
#' 
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom exactRankTests wilcox.exact
#' @importFrom coin kruskal_test
#' @importMethodsFrom coin pvalue
#' @export
#' 
#' 

ancom.detect <- function(otu_data, n_otu, alpha, multcorr, ncore){
  
  ## Detect whether the data are dependent or not
  if( ncol(otu_data) == n_otu+1  ){
    Group     <- otu_data[, ncol(otu_data) ]
    ID        <- rep( 1 , nrow(otu_data) )
    repeated <- FALSE
    fformula  <- formula("lr ~ Group")
  } else if( ncol(otu_data) == n_otu+2  ){
    Group     <- otu_data[, ncol(otu_data)-1 ]
    ID        <- otu_data[, ncol(otu_data)   ]
    repeated <- TRUE
    fformula  <- formula("lr ~ Group | ID")
    
  } else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  ## Detect which test to use, 
  ## Dependent data: Friedman test
  ## Independent data: Wilcoxon Rank-Sum or the Kruskal-Wallis
  ## exactRankTests::wilcox.exact is faster than stats::kruskal.test and stats::wilcox.test
  if( repeated==FALSE ){
    if( length(unique(Group))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  } else{
    tfun <- stats::friedman.test
  }

  
  
  ## Parallelized way to get the logratio.mat
  ## Doubles the number of computations to make, so only run the parallel
  ## version if there are multiple cores. Method may also add some computational
  ## overhead, so if only 2 cores, the nested for-loop shoud have advantage
  ## over the parallel loop (though I have not tested that).
  ## For some reason this is taking much longer, do not run the parallel loop as of now.
  if( FALSE ){
    registerDoParallel( cores=ncore )
  
    aa <- bb <- NULL
    logratio.mat <- foreach( bb = 1:n_otu, .combine='rbind', .packages="foreach" ) %:% 
      foreach( aa = 1:n_otu , .combine='c',  .packages="foreach" ) %dopar% {
        if( aa==bb ){
          p_out <- NA
        } else{
          data.pair <- otu_data[,c(aa,bb)]
          lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
          lr_dat <- data.frame( lr=lr, Group=Group, ID=ID )
          p_out  <- tfun(formula=fformula, data = lr_dat)$p.value
        }
        p_out
    }
    rownames(logratio.mat) <- colnames(logratio.mat) <- NULL
  } else{
    logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
    for(ii in 1:(n_otu-1)){
      for(jj in (ii+1):n_otu){
        data.pair <- otu_data[,c(ii,jj)]
        lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
        lr_dat <- data.frame( lr=lr, Group=Group, ID=ID )
        
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }
    }  
    ind <- lower.tri(logratio.mat)
    logratio.mat[ind] <- t(logratio.mat)[ind]
  }
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<alpha))
    })
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<alpha))
    })
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<alpha))
    })
  }
  return(W)
}

############################################################
############################################################


#' Run the ANCOM method
#' 
#' @description
#' Runs ANCOM to test for differences in OTU abundances between treatment groups.
#' 
#' @param OTUdat the OTU dataset. See Details for formatting instructions.
#' @param sig the significance level (or FDR) at which the tests are to be performed.
#' @param multcorr type of correction for multiple comparisons, see Details.
#' @param tau a tuning parameter in the Stepwise testing method. See Details.
#' @param theta a tuning parameter in the Stepwise testing method. See Details.
#' @param repeated logical determining whether the data have repeated measures (e.g., longitudinal design).
#' @details
#' The ANCOM method was developed and tested with default values of the two tuning parameters 
#' (\code{tau=0.02} and \code{theta=0.1}). For consistency, users are recommended to leave 
#' these tuning parameters at their default values, unless they wish to explore the performance
#'  of ANCOM for different values of the tuning parameters.
#' 
#' Data should be formatted as follows: each row is a subject, and each column is an OTU.
#' The final column should contain the grouping variable.
#' 
#' To adjust for multiple testing, \code{multcorr} may take on the following values:
#' \itemize{
#' \item{ \code{1}: }{ A stringent correction}
#' \item{ \code{2}: }{ A less stringent correction}
#' \item{ \code{3}: }{ No correction (default)}
#' }
#' The more stringent correction is not available in the shiny application.
#' 
#' @note
#' The function \code{\link{plot_ancom}} will produce plots for objects produced by \code{ANCOM}.
#' 
#' @return
#' The function produces a list with the following elements:
#' \itemize{
#' \item{ \code{W}: }{ values of the test statistics.}
#' \item{ \code{detected}: }{ names of OTUs detected.}
#' \item{ \code{dframe}: }{ the input dataframe.}
#' }
#'  
#' @export
#' 
#' @examples
#' 
#' \dontrun{
#' ## Create and run a small example
#' 
#' nn <- 10
#' pp <- 20
#' sim_otu <- matrix( 0, nrow=nn, ncol=pp+1 )
#' sim_otu <- data.frame(sim_otu)
#' colnames(sim_otu) <- c( paste0("OTU_", letters[1:pp] ), "Group" )
#' sim_otu[,pp+1]    <- c( rep("Control",nn/2), rep("Treatment",nn/2)  )
#' idx_trt <- sim_otu$Group=="Treatment"
#' 
#' for( ii in 1:pp ){
#'   sim_otu[,ii] <- rpois( nn, 1 )
#' }
#' 
#' # Create some significance
#' sim_otu[idx_trt,3] <- rpois( nn/2, 8)
#' sim_otu[idx_trt,7] <- rpois( nn/2, 8)
#' sim_otu[idx_trt,9] <- rpois( nn/2, 8)
#' 
#' ancom.out <- ANCOM( OTUdat = sim_otu, sig = 0.20, multcorr = 2 )
#' ancom.out$W
#' ancom.out$detected
#' }
#' 

ANCOM  <-  function(OTUdat, sig=0.05, multcorr=3, tau=0.02, theta=0.1, repeated=FALSE ){
  
  #OTUdat <-  read.delim(filepath,header=TRUE)
  
  num_col <- ncol( OTUdat )
  
  if( repeated==FALSE ){
    colnames(OTUdat)[ num_col ] <- "Group"    # rename last column as "Group"
    num_OTU      <- ncol(OTUdat) - 1
    
    sub_drop <- data.frame( nm_drop= "N/A" )
    sub_keep <- data.frame( nm_keep= "All subjects" )
    colnames(sub_drop) <- "Subjects removed"
    colnames(sub_keep) <- "Subjects retained"
    n_summary   <- paste0( "No subjects entirely removed (not a repeated-measures design)" )
    
  } else{
    colnames(OTUdat)[ num_col-1 ] <- "Group"    # rename 2nd last column as "Group"
    colnames(OTUdat)[ num_col   ] <- "ID"    # rename last column as "ID"
    OTUdat$ID    <- factor( OTUdat$ID )
    num_OTU      <- ncol(OTUdat) - 2
    
    ## Drop subjects if missing at a given time point
    crossTab <- table( OTUdat$Group , OTUdat$ID  )==0
    id_drop  <- apply( crossTab, 2, FUN=function(x) any(x)  )
    nm_drop  <- names( which( id_drop ) )
    idx_drop <- OTUdat$ID %in% nm_drop
    OTUdat   <- OTUdat[ idx_drop==FALSE, ]
    
    if( nrow(OTUdat)==0 ){ stop("Too many missing values in data, all subjects dropped") }    
    OTUdat$ID <- droplevels( OTUdat$ID )    
    num_dropped <- sum(id_drop)
    num_retain  <- length(id_drop) - num_dropped
    
    sub_drop <- data.frame( nm_drop=paste(nm_drop, collapse=", " ) )
    sub_keep <- data.frame( nm_keep= paste(levels(OTUdat$ID), collapse=", " ) )
    colnames(sub_drop) <- "Subjects removed"
    colnames(sub_keep) <- "Subjects retained"
    n_summary   <- paste0( "Analysis used ", num_retain, " subjects (", num_dropped, " were removed due to incomplete data)")
  }
  
  OTUdat$Group <- factor( OTUdat$Group )
  OTUdat       <- data.frame( OTUdat[ which(is.na(OTUdat$Group)==FALSE),],row.names=NULL )
  
  W.detected   <- ancom.detect(OTUdat, num_OTU, sig, multcorr, ncore=1 )
  W_stat       <- W.detected
  
  
  ## Per Shyamal (June 4, 2015): 
  ##  If number of OTUs is < 10, then use 'arbitrary' method, reject Ho if 
  ##  W > p - 1, where p = number of OTUs. If number of OTUs > 10, then use
  ##  the stepwise method. Consequently, only one output will be produced,
  ##  instead of detected_arbitrary and detected_stepwise, produce "detected"
  ## Rephrase "arbitrary", since it's not arbitrary, more of an empirical method.
  
  ## Detected using arbitrary cutoff
  # Previous code:
  # detected_arbitrary <- colnames(OTUdat)[ which( W.detected > num_OTU*theta ) ]
  
  if( num_OTU < 10 ){
    detected <- colnames(OTUdat)[which(W.detected > num_OTU-1 )]    
  } else{
    ## Detected using a stepwise mode detection
    if( max(W.detected)/num_OTU >= theta ){
      c.start <- max(W.detected)/num_OTU
      cutoff  <- c.start-c(0.05,0.10,0.15,0.20,0.25)
      
      prop_cut <- rep(0,length(cutoff))
      for(cut in 1:length(cutoff)){
        prop_cut[cut] <- length(which(W.detected>=num_OTU*cutoff[cut]))/length(W.detected)
      } 
      
      del <- rep(0,length(cutoff)-1)
      for( ii in 1:(length(cutoff)-1) ){
        del[ii] <- abs(prop_cut[ii]-prop_cut[ii+1])
      }
      
      if(       del[1]< tau & del[2]<tau & del[3]<tau ){ nu=cutoff[1]
      }else if( del[1]>=tau & del[2]<tau & del[3]<tau ){ nu=cutoff[2]
      }else if( del[2]>=tau & del[3]<tau & del[4]<tau ){ nu=cutoff[3]                                
      }else{ nu=cutoff[4] }
      
      up_point <- min(W.detected[ which( W.detected >= nu*num_OTU ) ])
      
      W.detected[W.detected>=up_point] <- 99999
      W.detected[W.detected<up_point]  <- 0
      W.detected[W.detected==99999]    <- 1
      
      detected <- colnames(OTUdat)[which(W.detected==1)]
      
    } else{
      W.detected <- 0
      detected   <- "No significant OTUs detected"
    }
    
  }
  
  #results_list <- list( W         = W_stat,
  #                      Arbitrary = detected_arbitrary,
  #                      Stepwise  = detected_stepwise )
  #idx0 <- lapply( results_list , FUN=length)
  #results_list[idx0==0] <- "No significant OTUs detected"
  
  results <- list( W=W_stat, detected=detected, dframe=OTUdat, repeated=repeated,
                   n_summary=n_summary, sub_drop=sub_drop, sub_keep=sub_keep)
  class(results) <- "ancom"
  
  return(results)  
  
}
#########################################################################





