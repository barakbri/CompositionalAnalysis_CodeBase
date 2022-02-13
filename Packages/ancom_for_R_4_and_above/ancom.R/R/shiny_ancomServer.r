#'
#' @rdname shiny_ancom
#' 
#' @importFrom DT renderDataTable
#' @importFrom DT dataTableOutput
#' @import shiny
#' @import ggplot2
#' @importFrom stringr str_pad
#' @importFrom stringr str_trim
#' 
#' @export
#' 

shiny_ancomServer <- function(input, output) {
    
  ## Function in import the data
  file_input <- function(file, header=TRUE, n_rows=-1 ){
    filename  <- file[1][[1]]
    fileloca  <- file[4][[1]]
    file_splt <- strsplit( paste(filename), "." , fixed=TRUE)[[1]]
    extn      <- file_splt[ length(file_splt)]
    
    if( extn=="csv" ){
      file_in  <- read.table( file=fileloca , sep=",", header=header, nrows=n_rows, stringsAsFactors=FALSE)
    }
    if( extn=="txt" ){
      file_in  <- read.table( file=fileloca , sep="\t", header=header, nrows=n_rows, stringsAsFactors=FALSE )
    }
    #if( extn=="xlsx" ){
    #  if( n_rows>0 ){
    #    m_rows <- 1:n_rows
    #  } else{
    #   m_rows <- NULL
    #  }
    #  file_in  <- read.xlsx( xlsxFile=fileloca , colNames=header, rows=m_rows )
    #}
    return( list( file_in=file_in, extn=extn ) )
  }
  
  
  
  ancom_out <- reactive({
    
    compute1 <- input$compute1
    
    ## Put all the code to run ANCOM inside this
    if( compute1 > 0 ){
      isolate({
        
        ## Extract the arguments
        file1a  <- input$file1
        file1   <- file1a[4]
        
        file_out <- input$file_out

        datafmt <- input$datafmt
        rept    <- input$repeated
        
        #ncores  <- input$ncores
        ncores <- 1 # <-- This is until I can fix the parallel
        
        tau     <- 0.02 # Hard-coded at least for shiny purposes
        theta   <- 0.1
        
        if( input$adjust == TRUE ){
          mc_type <- 2
          alpha   <- input$fdr
        } else{
          mc_type <- 3
          alpha   <- input$alpha
        }
        
        #colnames(results)= paste0("Significant OTU/taxa at FDR = ", sig )
        if( mc_type == 3 ){
          listTitle <- paste0("Significant OTU/taxa at ", alpha," level of significance" )
        } else{
          listTitle <- paste0("Significant OTU/taxa at FDR = ", alpha )
        }
        
        # dlmter  <- input$dlmter  ## For input some other time, to allow more flexible file formats     
        
        ## Load the data        
        if( datafmt=="wide"){
          #dat01  <- read.csv( paste(file1), header=TRUE, stringsAsFactors=FALSE )
          the_data <- file_input( file1a ) 
          dat01    <- the_data$file_in
          extn     <- the_data$extn
          
          if( is.character(dat01[,1]) ){
            stop("Specified subjects on rows, but first column is character.")
          }
          
          if( rept==FALSE ){
            pp <- ncol(dat01)
            group_name <- colnames(dat01)[]
            Group <- factor( str_trim( paste(dat01[,pp]), side="both" ) )
            dat01 <- data.frame( dat01[,-pp], Group )
          } else{
            pp <- ncol(dat01)
            group_name <- colnames(dat01)[]
            Group <- factor( str_trim( paste(dat01[,pp-1]), side="both" ) )
            ID    <- factor( str_trim( paste(dat01[,pp  ]), side="both" ) )
            dat01 <- data.frame( dat01[,-c(pp-2, pp)], Group, ID )
          }

          
        }
        if( datafmt=="tall" ){
          #dat_chk <- read.csv( paste(file1), header=TRUE, nrows=5 , stringsAsFactors=FALSE )
          the_data <- file_input( file1a , header=TRUE, n_rows=5 ) 
          dat_chk  <- the_data$file_in
          extn     <- the_data$extn
          
          if( is.numeric(dat_chk[,1]) ){
            stop("Specified OTUs on rows, but first column is not OTU names (character).")
          }
          #dat01a <- read.csv( paste(file1), header=FALSE, stringsAsFactors=FALSE  )
          the_data <- file_input( file1a , header=FALSE ) 
          dat01a   <- as.matrix(the_data$file_in)
          
          rnam   <- paste(dat01a[,1])

          if( rept==FALSE ){
            
            rnam   <- c( rnam[-1] )
            Group  <- factor( str_trim( paste(dat01a[1,])[-1], side="both" ) )
            dat01b <- as.data.frame( dat01a[-1,] )
            dat01b <- as.data.frame( dat01b[,-1] )
            for( ii in 1:ncol(dat01b) ){
              dat01b[,ii] <- as.numeric( paste(dat01b[,ii]))
            }
            dat01c <- as.data.frame(t(dat01b))
            colnames(dat01c) <- rnam
            dat01 <- cbind(dat01c, Group)
            
          } else{
            
            rnam   <- c( rnam[-c(1,2)] )
            Group  <- factor( str_trim( paste(dat01a[1,])[-1], side="both" ) )
            ID     <- factor( str_trim( paste(dat01a[2,])[-1], side="both" ) )
            
            dat01b <- as.data.frame( dat01a[-c(1,2),] )
            dat01b <- as.data.frame( dat01b[,-1] )
            for( ii in 1:ncol(dat01b) ){
              dat01b[,ii] <- as.numeric( paste(dat01b[,ii]))
            }
            dat01c <- as.data.frame(t(dat01b))
            colnames(dat01c) <- rnam
            dat01 <- cbind(dat01c, Group, ID)
            
          }
          
          rownames(dat01) <- NULL
          
        }
                
        ## Run the model and return the results
        ## May be good to NAME the arguments if there become too many
        
        withProgress(message = 'Status:', value = 1, {
          setProgress(1/2, detail = paste("Computing"))
          
          ancom.results <- ANCOM( OTUdat=as.data.frame(dat01) , sig=alpha, 
                                  multcorr=mc_type, tau=tau, theta=theta, repeated=rept )
        })
        
        
        ## Get the dataset of the detected OTUs and splice the output filename
        idx_detect   <- which( colnames(dat01) %in% ancom.results$detected )
        dat_selected <- t(dat01[ , c(ncol(dat01), idx_detect) ])
        
        sig_OTO               <- data.frame( names=ancom.results$detected )
        colnames(sig_OTO)     <- listTitle
        
        ancom.results$df_out  <- as.data.frame(dat_selected)
        ancom.results$sig_OTO <- sig_OTO
        
        if( nchar(file_out) >0 ){
          file_in  <- file1a[1][[1]]
          list_out <- paste0( file_out, "/OTU_list_" , file_in )
          data_out <- paste0( file_out, "/Selected_OTU_data_" , file_in )
          
          #file_parts   <- strsplit( file_out, "/")[[1]]
          if( extn=="csv" ){
            # file_parts[length(file_parts)] <- "Raw_data_detected_OTU.csv" 
            # data_out <- paste( file_parts, collapse="/")
            write.table(   sig_OTO,      list_out, row.names=FALSE, sep=",", quote=FALSE )
            write.table(   dat_selected, data_out, col.names=FALSE, sep=",", quote=FALSE )
          }
          if( extn=="txt" ){
            #file_parts[length(file_parts)] <- "Raw_data_detected_OTU.txt" 
            #data_out <- paste( file_parts, collapse="/")
            write.table( sig_OTO,      list_out, row.names=FALSE, sep="\t", quote=FALSE )
            write.table( dat_selected, data_out, col.names=FALSE, sep="\t", quote=FALSE )
          }
          #if( extn=="xlsx" ){
          #  #file_parts[length(file_parts)] <- "Raw_data_detected_OTU.xlsx" 
          #  #data_out <- paste( file_parts, collapse="/")
          #  write.xlsx(  sig_OTO,      list_out, row.names=FALSE )
          #  write.xlsx(  dat_selected, data_out, row.names=TRUE, col.names=FALSE )
          #}
        }
        
        return( ancom.results )
        
      })
      
    }
    
    
  }) ## End of object creation
  
  
  ##############################################################################
  
  output$etc <-  renderPrint({
    
    out1 <- ancom_out()
    print( out1$file_out )
    
    
  })
  
  ##############################################################################
  
  ##
  ## Make the SUMMARY output
  ##
  
  output$summary <-  renderPrint({
    out <-ancom_out()
    if( length(out) > 1 ){
      sig_OTO <- out$sig_OTO
      print.data.frame( sig_OTO , row.names=FALSE )
    }
  })
  
  output$summary02 <-  renderPrint({
    out <-ancom_out()
    if( length(out) > 1 ){
      n_summary <- out$n_summary
      cat( n_summary )
    }
  })
  
  output$summary03 <-  renderPrint({
    out <-ancom_out()
    if( length(out) > 1 ){
      sub_keep <- out$sub_keep
      cat( paste(sub_keep[1,]) )    }
  })
  
  output$summary04 <-  renderPrint({
    out <-ancom_out()
    if( length(out) > 1 ){
      sub_drop <- out$sub_drop
      cat( paste(sub_drop[1,]) )    }
  })
  
  
  ##############################################################################
  
  ##
  ## Make the DATSET output
  ##
  
  output$data_sel <-  DT::renderDataTable({
    out <-ancom_out()
    if( length(out) > 1 ){
      data_out <- cbind( Names=rownames(out$df_out), out$df_out )
      data_out
    }
  })
  
  ##############################################################################
  
  ##
  ## Make the GRAPH output
  ##
  
  output$plot <-  renderPlot({
    out     <- ancom_out()
    go_plot <- input$compute2 + input$compute1
    
    if( input$fixPlot ){
      ncols <- input$ncols
    } else{
      ncols <- -1
    }
    
    if( go_plot>1 ){
      isolate({
        
        plot_ancom( out, ncols=ncols )
        
      })
    } else{
      plot( 1:5 , 1:5 , col="white", xaxt='n', yaxt='n', xlab="", ylab="", frame.plot=FALSE )
    }
    
  })
  ##############################################################################
  
  ##
  ## Put all the output TOGETHER (so that we can manipulate the plot size)
  ##
  
  output$theTabset <- renderUI({
    
    if( input$fixPlot ){
      pltWidth  <- input$pltWidth
      pltHeight <- input$pltHeight
    } else{
      pltWidth  <- 100
      pltHeight <- 400
    }
    
    summary1 <- verbatimTextOutput("summary") 
    summary2 <- verbatimTextOutput("summary02") 
    summary3 <- verbatimTextOutput("summary03") 
    summary4 <- verbatimTextOutput("summary04") 
    plot1    <- plotOutput("plot", width = paste0(pltWidth, "%"), height = pltHeight)
    dat_out1 <- DT::dataTableOutput( "data_sel")
    etc1     <- verbatimTextOutput("etc")
    
    tabsetPanel(
      tabPanel("Summary", h4("Summary of Analysis"), summary1, summary2, 
               h5("Subjects retained"), summary3,
               h5("Subjects removed"),  summary4 ),
      tabPanel("Boxplots of Abundances", plot1 ),
      tabPanel("Data (detected OTU)", dat_out1 )
      #tabPanel("Etc", etc1 )
    )
    
  })
  
    
}

