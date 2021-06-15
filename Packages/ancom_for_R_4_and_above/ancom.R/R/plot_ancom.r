#' Plot of data from objects of class 'ancom'
#' 
#' @description
#' Produces comparison boxplots of data for objects of class \code{ancom}.
#' 
#' @param object object of class \code{ancom}.
#' @param ncols the number of columns for \code{ggplot} to produce using \code{facet_wrap}.
#'        If \code{ncol=-1}, then the function will attempt to 
#' @param ... space for additional arguments (none corrently)
#' 
#' @details
#' \code{plot_ancom} uses \pkg{ggplot} to produce graphics.
#' 
#' 
#' @import ggplot2
#' 
#' @export
#' 


plot_ancom <- function( object, ncols=-1, ... ){
  
  if( !(class(object)=="ancom") ){
    stop("'object' is not of class ancom")
  }
  
  repeated <- object$repeated
  
  Group  <- OTU <- ID <- NULL
  dframe <- object$dframe
  
  if( repeated==FALSE ){
    colnames(dframe)[ ncol(dframe) ] <- "Group"  
  } else{
    colnames(dframe)[ ncol(dframe)-1 ] <- "Group"
  }
    
  # OTUs that were detected
  Sig_OTU <- object$detected
  
  if( Sig_OTU[1] == "No significant OTUs detected" ){
    ## Plot a message so that SOMETHING is produced
    plot( 1:5 , 1:5 , col="white", xaxt='n', yaxt='n', xlab="", ylab="", frame.plot=FALSE )
    text( 3, 3, labels="No significant OTUs detected" )
    
  } else{
    
    if( repeated==FALSE ){
      W_check <- data.frame( colnames(dframe)[-ncol(dframe)], object$W, row.names=NULL)
      colnames(W_check) <- c("OTU_ID","W")  
    } else{
      W_check <- data.frame( colnames(dframe)[-c(ncol(dframe)-1, ncol(dframe) )], object$W, row.names=NULL)
      colnames(W_check) <- c("OTU_ID","W")  
    }
        
    W_check        <- W_check[which(  W_check$OTU_ID %in% Sig_OTU),]
    W_check        <- W_check[order(-W_check$W),]
    nplot          <- nrow(W_check)
    W_check$OTU_ID <- as.character(W_check$OTU_ID)
    
    # Get the DATA to plot
    for( ii in 1:nplot ){
      dsub <- dframe[ , colnames(dframe) %in% c(W_check$OTU_ID[ii],"Group") ]
      colnames(dsub)=c("OTU","Group")
      OTU_name <- rep( W_check$OTU_ID[ii] , nrow(dsub) )
      pltDat00 <- data.frame( OTU_name , dsub )              
      if( ii==1 ){
        pltDat <- pltDat00
      } else{
        pltDat <- rbind(pltDat, pltDat00 )
      } 
    }
    
    pltDat$OTU      <- log( pltDat$OTU + 1 )
    pltDat$OTU_name <- factor( pltDat$OTU_name , Sig_OTU )
    
    if( ncols<1 ){
      ncols <- min(3, nplot)
    }
    
    gplot <- ggplot( pltDat , aes(x=factor(Group), y=OTU ) ) + 
      facet_wrap( ~ OTU_name , ncol=ncols, scales="free_y") + 
      geom_boxplot()
    
    gplot + labs(x = "Grouping Factor" , y="Log of Abundance"  ) + theme(
      panel.background  = element_rect( fill="white" , colour="black"),
      panel.grid        = element_blank(),
      strip.text        = element_text( size=rel(1.25)),
      strip.background  = element_rect( fill="grey90" , color="black"),
      axis.title        = element_text( size=rel(1.25) , color="black"),
      axis.text         = element_text( size=rel(1.05) , color="black"),
      # legend.position   = c(0.90,0.10),
      legend.position   = "none",
      legend.background = element_rect(colour = "black", fill="white"),
      legend.key        = element_rect(fill = "white"),
      legend.title      = element_text(size = 15 ),
      legend.text       = element_text(size = 12 )
    )
    
  }
    
}

