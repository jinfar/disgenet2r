#' Plots the content of a \code{DataGeNET.DGN} object.
#' This functions allows to create a variety of plots for \code{DataGeNET.DGN} objects.
#'
#' @name plot
#' @rdname plot-methods
#' @aliases DataGeNET.DGN-methods
#'
#' @param x Object of class \code{DataGeNET.DGN}
#' @param y NOT USED
#' @param layout Function to design the location of the different nodes. By
#' default \code{layout.fruchterman.reingold} from \code{igraph} is used.
#' @param class Type of the drawn chart. By default it is \code{"disease"} but
#' it also can be \code{"gene"}, \code{"DiseaseClass"},
#' \code{"DiseaseDisease"}, \code{"Heatmap"}, \code{"DiseaseClassHeatmap"},
#' \code{"ProteinClass"}, \code{"ProteinClassHeatmap"}, \code{"Barplot"},
#' \code{"Venn"}.
#' @param verbose By default \code{FALSE}. If set to \code{TRUE} information
#' on the drawing process will be shown.
#' @param ... Passed to inner functions for different plots.
#' @return A plot for \code{DataGeNET.DGN}.
#' @examples
#' \dontrun{
#' #Being x an DataGeNET.DGN
#' qr <- plot(x) # Get number of unique diseases
#' }
#' @export plot

setMethod(
  f = "plot",
  signature = "DataGeNET.DGN",
  definition = function( x, y, ... ) {
    plot_datagenet_dis( x, ... )
  }
)


plot_datagenet_dis <- function( object, layout = igraph::layout.fruchterman.reingold, class = "disease", topLevel = FALSE, prop = 1, limit = 100, verbose = FALSE, cutoff = 0, ... ) {

  if( !class %in% c( "disease", "gene", "DiseaseClass", "DiseaseDisease", "Comorbidity", "Heatmap", "DiseaseClassHeatmap", "ProteinClass", "ProteinClassHeatmap", "Barplot", "Venn", "gda" ) ) {
    stop( "Invalid content of argument 'class'." )
  }

  if( object@type == "variant" | object@type == "variant-disease") {
    if (class == "disease" ){
      plot_variant_disease( search = object@search, type = object@type, input = object@qresult, layout = layout, verbose = verbose, prop = prop )
    }
    if (class == "gene" ){
      plot_variant_gene( search = object@search, type = object@type, input = object@qresult, layout = layout, verbose = verbose, prop = prop )
    }
    if (class == "gda" ){
      plot_variant_gda( search = object@search, type = object@type, input = object@qresult, layout = layout, verbose = verbose, prop = prop )
    }
    if (class == "Heatmap" ){
      plot_variant_disease_heatmap( type= object@type, search= object@search, input = object@qresult, cutoff = cutoff, verbose = verbose, limit = limit )
    }
    if (class == "DiseaseClass" ){
      plot_variant_disclass( search = object@search, type = object@type, input = object@qresult, nchars = nchars, layout = layout, verbose = verbose, prop = prop )
    }
    if (class == "DiseaseClassHeatmap" ){
      plot_dis_heatmapDiseaseClass( search = object@search, type = object@type, input = object@qresult,nchars=nchars, verbose = verbose)
    }
  }

  else if( object@type == "disease") {
    if(object@search == "dis-dis" & class == "disease"){
      plot_dis_comorbidity( search= object@search, type = object@type, input = object@qresult, layout = layout, verbose = verbose ,  prop = prop, limit = limit)
    }
    else if(object@search == "dis-dis" &  class == "DiseaseClass" ) {
      plot_dis_Class( search= object@search, type = object@type, input = object@qresult, layout = layout, verbose = verbose, nchars = nchars,prop = prop )
    }
    else if(object@search == "dis-dis" &  class == "Comorbidity" ) {
      plot_dis_dis_associations( search= object@search, type = object@type, input = object@qresult, layout = layout, verbose = verbose ,  prop = prop, limit = limit,  database= object@database )
    }
    else if( class == "disease" ) {
      plot_dis_disease( search = object@search, type = object@type, input = object@qresult, layout = layout, verbose = verbose, prop = prop , limit = limit)
    }  else if( class == "DiseaseDisease" ) {
      plot_dis_comorbidity( search= object@search, type = object@type, input = object@qresult, layout = layout, verbose = verbose, prop = prop, limit = limit )
    }  else if( class == "Heatmap" ) {
      plot_dis_heatmapDisease( type= object@type, search= object@search, input = object@qresult, cutoff = cutoff, nchars = nchars, verbose = verbose, limit = limit )
    } else if( class == "ProteinClass" ) {
      plot_dis_Panther(  search= object@search, input = object@qresult, layout = layout, verbose = verbose )
    }else if( class == "ProteinClassHeatmap" ) {
      plot_dis_heatmapPanther(  type= object@type, search= object@search, input = object@qresult,nchars= nchars, verbose = verbose )
    } else if( class == "Barplot" ) {
      plot_dis_barplot(  type= object@type, search= object@search, input = object@qresult, verbose = verbose, limit = limit )
    }  else if( class == "Venn" ) {
      plot_dis_venn(  type= object@type, search= object@search, input = object@qresult, verbose = verbose, database= object@database )
    }
    else{
      stop( "Invalid content of argument 'class' for type = 'disease'." )
    }
  }

  else if( object@type  == "gene" ) {
    if( class == "disease" ) {
      plot_dis_disease( search = object@search, type = object@type, input = object@qresult, layout = layout, verbose = verbose, prop = prop , limit = limit)
    } else if( class == "DiseaseClass" ) {
      plot_dis_Class( search= object@search, type = object@type, input = object@qresult, layout = layout, nchars= nchars, verbose = verbose, prop = prop )
    } else if( class == "Heatmap" ) {
      plot_dis_heatmapDisease( type= object@type, search= object@search, input = object@qresult, cutoff = cutoff, nchars= nchars, verbose = verbose, limit = limit )
    }  else if( class == "DiseaseClassHeatmap" ) {
      plot_dis_heatmapDiseaseClass( type= object@type, search= object@search, input = object@qresult, nchars= nchars, verbose = verbose, limit = limit)
    }  else{
      stop( "Invalid content of argument 'class' for type = 'gene'." )
    }
  }
  else{
    stop( "Invalid content of argument 'type'." )
  }
}


plot_variant_disease<- function( input, type, layout, prop, search, verbose ) {

  input <- unique(input[c("c0.score", "c0.diseaseId", "c2.name", "c0.snpId")])
  l <- length(unique(input$c0.snpId))
  input$score <- as.numeric(as.character(input$c0.score))
  input <-  aggregate(score ~ c0.diseaseId + c2.name + c0.snpId, input , max)
  if(l > 1 ) {
    ttl <- "Variant-Disease Network"
  }else if ( l == 1 ){
    ttl <- paste0( "Diseases associated to ", input$c0.snpId[ 1 ] )
  }

  edges <- data.frame( input[ , "c2.name"], input[ , "c0.snpId"] )
  netw  <- igraph::graph.data.frame( edges, directed = FALSE )
  netw  <- igraph::simplify( netw )
  lay   <- layout( netw )
  diseases <- unique( input$c2.name )
  igraph::plot.igraph( netw,
                       vertex.frame.color  = "white",
                       layout              = lay,
                       vertex.color        = ifelse( igraph::V( netw )$name %in% diseases, "royalblue2", "goldenrod3" ),
                       vertex.label.dist   = 0,      #puts the name labels slightly off the dots
                       vertex.frame.color  = 'blue', #the color of the border of the dots
                       vertex.label.color  = 'black',#the color of the name labels
                       vertex.label.font   = 0,      #the font of the name labels
                       vertex.label        = igraph::V( netw )$names, #specifies the lables of the vertices. in this case the 'name' attribute is used
                       edge.color          = "darkgrey",
                       edge.width          = 1 + ( prop * input$score ),
                       edge.arrow.size     = 0.5,
                       vertex.size         = 15,
                       vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
                       main                = ttl
  )
}

plot_variant_gene<- function( input, type, layout, prop, search, verbose ) {

  input <- unique(input[c("c0.snpId", "c1.DSI", "c1.DPI" , "c1.symbol")])
  l <- length(unique(input$c0.snpId))
  #input$score <- as.numeric(as.character(input$score))
  #input[input$c0.score < 0.01,"c0.score" ] <- 0.09
  #input <-  aggregate(score ~ geneId + geneSymbol  + snp, input , max)
  tt <- as.data.frame( transDiseaseClass ( input ) )
  row.names(tt) <- NULL
  tt[,"V4"] <- gsub("null", NA, tt[ ,"V4"] )
  nd <- length(unique(subset(tt, is.na(V4))$c0.snpId))
  input <- subset(tt, !is.na(V4))

  if( l > 1 ) {
    ttl <- "Variant-Gene Network"
  }else if ( l == 1 ){
    ttl <- paste0( "Genes associated to ", input$c0.snpId[ 1 ] )
  }
  if (nd > 0){
    show( paste0("warning: ", nd , " variant(s) not shown in the plot" ))
  }
  colnames(input)<- c("c0.snpId", "c1.DSI", "c1.DPI", "V4")
  edges <- data.frame( input[ , "V4"], input[ , "c0.snpId"] )
  netw  <- igraph::graph.data.frame( edges, directed = FALSE )
  netw  <- igraph::simplify( netw )
  lay   <- layout( netw )
  genes <- unique( input$V4 )
  igraph::plot.igraph( netw,
                       vertex.frame.color  = "white",
                       layout              = lay,
                       vertex.color        = ifelse( igraph::V( netw )$name %in% genes, "#ff349a" , "goldenrod3"),
                       vertex.label.dist   = 0,      #puts the name labels slightly off the dots
                       vertex.frame.color  = 'blue', #the color of the border of the dots
                       vertex.label.color  = 'black',#the color of the name labels
                       vertex.label.font   = 0,      #the font of the name labels
                       vertex.label        = igraph::V( netw )$names, #specifies the lables of the vertices. in this case the 'name' attribute is used
                       edge.color          = "darkgrey",
                       edge.width          = 1 + ( prop  ),
                       edge.arrow.size     = 0.5,
                       vertex.size         = 15,
                       vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
                       main                = ttl
  )
}

plot_variant_gda<- function( input, type, layout, prop, search, verbose ) {

  input <- unique(input[c("c0.score", "c0.snpId",  "c2.name", "c1.symbol")])
  input <- data.frame(lapply(input, as.character), stringsAsFactors=FALSE)
  l <- length(unique(input$c0.snpId))

  #input[input$c0.score < 0.01,"c0.score" ] <- 0.09
  tt <- as.data.frame( transDiseaseClass ( input ) )
  row.names(tt) <- NULL
  tt[,"V4"] <- gsub("null", NA, tt[ ,"V4"] )
  nd <- length(unique(subset(tt, is.na(V4))$c0.snpId))
  #input <- subset(tt, !is.na(V4))
  colnames(tt) <- colnames(input)
  input <-tt
  input$c0.score <- as.numeric(as.character(input$c0.score))
  aa <-  aggregate(c0.score ~  c1.symbol +c0.snpId , input , max)
  bb <-  aggregate(c0.score ~ c2.name +c0.snpId , input , max)
  cc <- unique(input[c("c2.name", "c1.symbol",  "c0.score")])
  colnames(aa) <- c("A1", "A2", "score")
  colnames(bb) <- c("A1", "A2", "score")
  colnames(cc) <- c("A1", "A2", "score")

  tt <- rbind(aa,bb)
  tt <- rbind(tt,cc)

  if( l > 1 ) {
    ttl <- "Variant-Gene-Disease Network"
  }else if ( l == 1 ){
    ttl <- paste0( "Variant-Gene-Disease Network for ", input$snp[ 1 ] )
  }
  tt <- subset( tt, !is.na(A1) & !is.na(A2))
  edges <- data.frame( tt[ ,1], tt[ , 2] )
  netw  <- igraph::graph.data.frame( edges, directed = FALSE )
  netw  <- igraph::simplify( netw )
  lay   <- layout( netw )
  genes <- unique( input$c1.symbol )
  diseases <- unique( input$c2.name )
  variants <- unique( input$c0.snpId )

  igraph::plot.igraph( netw,
                       vertex.frame.color  = "white",
                       layout              = lay,
                       vertex.color        = ifelse( igraph::V( netw )$name %in% genes, "#ff349a" , ifelse(igraph::V( netw )$name %in% diseases, "royalblue2", "goldenrod3")),
                       vertex.label.dist   = 0,      #puts the name labels slightly off the dots
                       vertex.frame.color  = 'blue', #the color of the border of the dots
                       vertex.label.color  = 'black',#the color of the name labels
                       vertex.label.font   = 0,      #the font of the name labels
                       vertex.label        = igraph::V( netw )$names, #specifies the lables of the vertices. in this case the 'name' attribute is used
                       edge.color          = "darkgrey",
                       edge.width          = 1 + ( prop * tt$score ),
                       edge.arrow.size     = 0.5,
                       vertex.size         = 15,
                       vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
                       main                = ttl
  )
}


plot_variant_disclass <- function( input, type, layout, nchars, prop, search, verbose ) {
  results <- unique(input[c("c0.score", "c0.snpId", "c2.name", "c2.diseaseClassName")])
  l <- length(unique(results$c0.snpId))
  #input$score <- as.numeric(as.character(input$score))
  if(l > 1 ) {
    stop( "Chart not possible for list. Use instead the plot function with DiseaseClassHeatmap class argument" )
  }

  tt <- as.data.frame( transDiseaseClass ( results ) )
  rownames(tt) <- NULL
  tt[,4] <- gsub (" ; ", ";", tt[,4] )
  tt[,4] <- gsub ( "null", NA, tt[,4] )
  tt[,4] <- gsub ( "Wounds and Injuries", NA, tt[,4] )
  tt[,4] <- gsub ( "Animal Diseases", NA, tt[,4] )

  nd <- length(unique(subset(tt, is.na(V4))$c2.name))
  tt <- subset(tt, !is.na(V4))
  for( i in 1:nrow(tt)){
    if(nchar(as.character(tt$V4[i])) > nchars){
      rename <- substr(tt$V4[i], 1, nchars)
      tt$V4[i] <- paste0(rename, "...")
    }
  }


  results <- tt

  freq <- as.data.frame ( table ( results$V4 ) )
  freq$perc <- ( freq$Freq / sum( freq$Freq ) ) * 100

  rw <- matrix ( c( as.character( results[ 1, "c0.snpId" ] ), 10, 10), ncol=3 )
  colnames( rw ) <- colnames( freq )
  freq <- rbind( freq, rw )

  sizes <- as.numeric( freq[ , 3 ] )
  names( sizes ) <- freq[ , 1 ]

  #plot
  edges <- data.frame( results[, "c0.snpId" ], results[ , "V4" ] )
  netw  <- igraph::graph.data.frame( edges, directed = FALSE )
  netw  <- igraph::simplify( netw )
  lay   <- layout( netw )

  if( verbose ) {
    message( "The network contains ", igraph::vcount( netw ), " nodes and ", igraph::ecount( netw ), " edges." )
  }

  diseases <- unique( results$V4 )
  ttl <- paste0("Disease Classes associated to ", results$c0.snpId[1] )
  if (nd > 0){
    show( paste0("warning: ", nd , " disease(s) not shown in the plot" ))
  }
  igraph::plot.igraph( netw,
                       vertex.frame.color  = "white",
                       layout              = lay,
                       vertex.color        = ifelse( igraph::V( netw )$name %in% diseases, "royalblue2", "goldenrod3" ),
                       vertex.label.dist   = 0,      #puts the name labels slightly off the dots
                       vertex.frame.color  = 'blue', #the color of the border of the dots
                       vertex.label.color  = 'black',#the color of the name labels
                       vertex.label.font   = 0,      #the font of the name labels
                       vertex.label        = igraph::V( netw )$names, #specifies the lables of the vertices. in this case the 'name' attribute is used
                       edge.color          = "darkgrey",
                       edge.arrow.size     = 0.5,
                       edge.width          = 0.5,
                       vertex.size         = as.numeric(sizes[ igraph::V( netw )$name ])*prop,
                       vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
                       main                = ttl
  )
}


plot_dis_disease <- function( input, type, layout, prop, search, verbose ,  limit) {

  if( search == "dis-dis" ) {
    #     print(limit)
    #     print(dim(input)[1])
    #     if (dim(input)[1]> limit){
    #       input <- input[c(1:limit),]
    #     }
    #     edges <- data.frame( input[ , 2 ], input[ , 4 ] )
    #     netw  <- igraph::graph.data.frame( edges, directed = FALSE )
    #     netw  <- igraph::simplify( netw )
    #     lay   <- layout( netw )
    #
    #     genNum <- input[ , c( 4, 6 ) ]
    #
    #     rw <- matrix ( c ( as.character ( input [ 1, 2 ] ), 10 ), ncol = 2 )
    #     colnames( rw ) <- colnames ( genNum )
    #     genNum <- rbind ( genNum, rw )
    #
    #     sizes <- as.numeric(genNum[ , 2 ])
    #     names( sizes ) <- genNum[ , 1 ]
    #
    #     if( verbose ) {
    #       message( "The network contains ", igraph::vcount( netw ), " nodes and ", igraph::ecount( netw ), " edges." )
    #     }
    #     #diseases <- unique( input$c1.name )
    #     #ttl <- paste0( "Diseases associated to ", input$c1.name[ 1 ] )
    #     ttl <- "Disease-Disease Network"
    #     igraph::plot.igraph( netw,
    #                          vertex.frame.color = "white",
    #                          layout              = lay,
    #                          vertex.color        = ifelse ( igraph::V( netw )$name %in% diseases, "royalblue2", "#9bcdff" ),
    #                          vertex.label.dist   = 0,      #puts the name labels slightly off the dots
    #                          vertex.frame.color  = 'blue', #the color of the border of the dots
    #                          vertex.label.color  = 'black',#the color of the name labels
    #                          vertex.label.font   = 0,      #the font of the name labels
    #                          vertex.label        = igraph::V( netw )$names, #specifies the lables of the vertices. in this case the 'name' attribute is used
    #                          edge.color          = "darkgrey",
    #                          edge.width          = 1 + ( prop * input$c0.score ),
    #                          edge.arrow.size     = 0.5,
    #                          vertex.size         = as.numeric( sizes[ igraph::V( netw )$name ] )*prop,
    #                          vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
    #                          main                = ttl
    #     )
    #
    #

  }else {

    if(search=="unify"){
      edges <- data.frame( input[ , 12 ], input[ , 6 ] )
      diseases <- unique( input$diseaseIdentifier )
    }else{
      edges <- data.frame( input[ , 2 ], input[ , 6 ] )
      diseases <- unique( input$c1.name )
    }

    netw  <- igraph::graph.data.frame( edges, directed = FALSE )
    netw  <- igraph::simplify( netw )
    lay   <- layout( netw )

    if( verbose ) {
      message( "The network contains ", igraph::vcount( netw ), " nodes and ", igraph::ecount( netw ), " edges." )
    }


    if( search == "list" ) {
      ttl <- "Gene-Disease Network"
    }else if ( search == "single" ){
      if ( type == "gene"){
        ttl <- paste0( "Diseases associated to ", input$c2.name[ 1 ]  )
      }else if ( type == "disease"){
        ttl <- paste0( "Genes associated to ", input$c1.name[ 1 ] )
      }
    }
    else if( search == "unify" ) {
      ttl <- paste0( "Genes associated to ", input$diseaseIdentifier[ 1 ] )
    }

    igraph::plot.igraph( netw,
                         vertex.frame.color  = "white",
                         layout              = lay,
                         vertex.color        = ifelse( igraph::V( netw )$name %in% diseases, "royalblue2", "#ff349a" ),
                         vertex.label.dist   = 0,      #puts the name labels slightly off the dots
                         vertex.frame.color  = 'blue', #the color of the border of the dots
                         vertex.label.color  = 'black',#the color of the name labels
                         vertex.label.font   = 0,      #the font of the name labels
                         vertex.label        = igraph::V( netw )$names, #specifies the lables of the vertices. in this case the 'name' attribute is used
                         edge.color          = "darkgrey",
                         edge.width          = 1 + ( prop * input$c0.score ),
                         edge.arrow.size     = 0.5,
                         vertex.size         = 15,
                         vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
                         main                = ttl
    )
  }
}

plot_dis_Class <- function ( input, type, layout, search, nchars, verbose, prop ) {

  if( search == "list" ) {
    stop( "Chart not possible for list. Use instead the plot function with DiseaseClassHeatmap class argument" )
  }
  else if( search == "dis-dis" ) {
    #stop( "Chart not possible for disease-disease assoc query" )
    if (length(unique(input[, 1])) == 1){

      input[,"diseaseClass2"] <- gsub ( "null", NA, input[,"diseaseClass2"] )
      input[,"diseaseClass2"] <- gsub ( "Wounds and Injuries", NA, input[,"diseaseClass2"]  )
      input[,"diseaseClass2"]  <- gsub ( "Animal Diseases", NA, input[,"diseaseClass2"] )
      nd <- length(unique(subset(input, is.na(diseaseClass2))$Disease2) )
      input <- subset(input, !is.na(diseaseClass2) )
      # this means there is a central disease
      results <- input [, c (1, 4,5,6,8,9 )]

      tt <- as.data.frame( transDiseaseClass ( results ) )
      row.names(tt) <- NULL
      tt$V4 <- as.character(tt$V4)
      tt <- unique(tt[c("Disease2", "V4")])

      for( i in 1:nrow(tt)){
        if(nchar(as.character(tt$V4[i])) > nchars){
          rename <- substr(tt$V4[i], 1, nchars)
          tt$V4[i] <- paste0(rename, "...")
        }
      }
      if (nd > 0){
        show( paste0("warning: ", nd , " disease(s) not shown in the plot" ))
      }

      results <- merge(results[, c(2,3,5,6)], tt, by.x = "Disease2")
      freq <- as.data.frame ( table ( results$V4 ) )
      freq$perc <- ( freq$Freq / sum( freq$Freq ) ) * 100
      rw <- matrix ( c( as.character( input[ 1,2 ] ),10, 10), ncol=3 )
      colnames( rw ) <- colnames( freq )
      freq <- rbind( freq, rw )

      sizes <- as.numeric( freq[ , 3 ] )
      names( sizes ) <- freq[ , 1 ]

      results <- merge( input, results[c(1,5)], by="Disease2")

      #plot
      edges <- data.frame( results[ , "NameDisease1" ], results[ , "V4" ] )
      netw  <- igraph::graph.data.frame( edges, directed = FALSE )
      netw  <- igraph::simplify( netw )
      lay   <- layout( netw )

      if( verbose ) {
        message( "The network contains ", igraph::vcount( netw ), " nodes and ", igraph::ecount( netw ), " edges." )
      }

      diseases <- unique( results$V4 )
      ttl <- paste0("Disease Classes associated to ", input[1,2] )

      igraph::plot.igraph( netw,
                           vertex.frame.color  = "white",
                           layout              = lay,
                           vertex.color        = ifelse( igraph::V( netw )$name %in% diseases, "#9bcdff", "royalblue2"),
                           vertex.label.dist   = 0,      #puts the name labels slightly off the dots
                           vertex.frame.color  = 'blue', #the color of the border of the dots
                           vertex.label.color  = 'black',#the color of the name labels
                           vertex.label.font   = 0,      #the font of the name labels
                           vertex.label        = igraph::V( netw )$name, #specifies the lables of the vertices. in this case the 'name' attribute is used
                           edge.color          = "darkgrey",
                           edge.arrow.size     = 0.5,
                           edge.width          = 0.5,
                           vertex.size         = as.numeric(sizes[ igraph::V( netw )$name ])*prop,
                           vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
                           main                = ttl
      )
    }
    else{
      stop( "Chart not possible for list. Use instead the plot function with DiseaseClassHeatmap class argument" )
    }
  }

  else if( search == "single" & type == "gene"  ) {
    results <- input [, c ( 6, 1, 2, 3 )]
    results[,4] <- gsub( "null", NA, results[ ,4 ] )

    results <- as.data.frame( transDiseaseClass ( results ) )
    nd <- dim(subset(results, is.na(V4)))[1]
    if (nd > 0){
      show( paste0("warning: ", nd , " disease(s) not shown in the plot" ))
    }
    results <- subset(results, !is.na(V4))
    freq <- as.data.frame ( table ( results$V4 ) )
    freq$perc <- ( freq$Freq / sum( freq$Freq ) ) * 100

    rw <- matrix ( c( as.character( results[ 1, 1 ] ), 10, 10), ncol=3 )
    colnames( rw ) <- colnames( freq )
    freq <- rbind( freq, rw )

    sizes <- as.numeric( freq[ , 3 ] )
    names( sizes ) <- freq[ , 1 ]

    #plot
    edges <- data.frame( results[ , 1 ], results[ , 4 ] )
    netw  <- igraph::graph.data.frame( edges, directed = FALSE )
    netw  <- igraph::simplify( netw )
    lay   <- layout( netw )

    if( verbose ) {
      message( "The network contains ", igraph::vcount( netw ), " nodes and ", igraph::ecount( netw ), " edges." )
    }

    diseases <- unique( results$V4 )
    ttl <- paste0("Disease Classes associated to ", results$c2.name[1] )

    igraph::plot.igraph( netw,
                         vertex.frame.color  = "white",
                         layout              = lay,
                         vertex.color        = ifelse( igraph::V( netw )$name %in% diseases, "royalblue2", "#ff349a" ),
                         vertex.label.dist   = 0,      #puts the name labels slightly off the dots
                         vertex.frame.color  = 'blue', #the color of the border of the dots
                         vertex.label.color  = 'black',#the color of the name labels
                         vertex.label.font   = 0,      #the font of the name labels
                         vertex.label        = igraph::V( netw )$names, #specifies the lables of the vertices. in this case the 'name' attribute is used
                         edge.color          = "darkgrey",
                         edge.arrow.size     = 0.5,
                         edge.width          = 0.5,
                         vertex.size         = as.numeric(sizes[ igraph::V( netw )$name ])*prop,
                         vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
                         main                = ttl
    )
  }
}

plot_dis_heatmapDisease <- function( type, search, input, cutoff=0,nchars, limit, verbose ) {

  if( search == "single" | search == "unify" ) {
    stop( "For this type of chart, a multiple query created with 'disgenetDisease'is required. For single disease or unify disease is not possible to apply this graphic." )
  }

  if( search == "list" ) {

    if( type == "disease" ) {
      tt <- ( input[ , c( "c1.name", "c2.symbol", "c0.score" ) ] )
      tt$c1.name <- as.character(tt$c1.name)
      tt <- tt [ tt$c0.score>= cutoff, ]
      tt <- tt [ with ( tt , order(-c0.score)), ]
      if ( dim( tt )[ 1 ] > limit ){
        tt <- tt[ 1:limit ,]
        show( paste0("Dataframe of ", nrow(input), " rows has been reduced to ", limit, " rows." ))

      }

      for( i in 1:nrow(tt)){
        if(nchar(as.character(tt$c1.name[i])) > nchars){
          rename <- substr(tt$c1.name[i], 1, nchars)
          tt$c1.name[i] <- paste0(rename, "...")
        }
      }

      casted_tt <- reshape::cast( tt, c2.symbol ~ c1.name, value= "c0.score", length )
      l <- dim( casted_tt )[ 2 ]
      casted_tt[ is.na( casted_tt ) ] <- 0
      res <- hclust( dist( casted_tt[ 2:l ] ),  method = "ward.D")
      casted_tt$clus <- cutree( res, k = l+1 )
      ordered <- casted_tt[order(casted_tt$clus), "c2.symbol"]

      tt$c2.symbol  <- factor(tt$c2.symbol , levels= as.factor(ordered))
      ggplot2::qplot( x = c1.name, y = c2.symbol, data= tt , fill = c0.score, geom = "tile" ) +
        ggplot2::ylab( "Genes" ) + ggplot2::xlab( "Disease Name" ) +
        ggplot2::ggtitle( "Gene-Disease Heatmap") +
        ggplot2::scale_fill_gradient2(low = "#0000FF", high ="#2f30c9",guide = "colourbar", "Score") +
        ggplot2::theme ( axis.text = ggplot2::element_text ( size = 9),
                         axis.title=ggplot2::element_text(size=10 ),
                         panel.background = ggplot2::element_rect(fill="white", colour = "black"),
                         text = ggplot2::element_text(size = 11),
                         axis.text.x = ggplot2::element_text(angle=45, size = 11, hjust = 1)
        )


    } else  if( type == "gene" ) {
      tt <- ( input[ , c( "c0.score", "c2.symbol", "c1.name" ) ] )
      tt$c2.symbol <- as.character( tt$c2.symbol )
      tt <- tt [ tt$c0.score>= cutoff, ]


      if ( dim( tt )[ 1 ] > limit ){
        tt <- tt[ 1:limit ,]
        show( paste0("warning: dataframe of ", nrow(input), " rows has been reduced to ", limit, " rows." ))

      }

      for( i in 1:nrow(tt)){

        if(nchar(as.character(tt$c1.name[i])) > nchars){
          tt$c1.name <- as.character(tt$c1.name )
          rename <- substr(tt$c1.name[i], 1, nchars)
          tt$c1.name[i] <- paste0(rename, "...")
        }
      }

      casted_tt <- reshape::cast(tt, c1.name ~ c2.symbol, value= "c0.score" , length )
      l <- dim(casted_tt)[2]
      casted_tt[is.na( casted_tt)] <- 0
      res <- hclust(dist(casted_tt[2:l]),  method = "ward.D")
      casted_tt$clus <- cutree(res, k = l +1)
      ordered <- casted_tt[order(casted_tt$clus), "c1.name"]

      tt$c1.name  <- factor(tt$c1.name , levels= as.factor(ordered))
      ggplot2::qplot( x = c2.symbol, y = c1.name, data= tt , fill = c0.score, geom = "tile" ) +
        ggplot2::ylab( "Diseases" ) + ggplot2::xlab( "Genes" ) + ggplot2::ggtitle( "Gene-Disease Heatmap") +
        ggplot2::scale_fill_gradient2(low = "#0000FF", high ="#2f30c9",guide = "colourbar", "Score") +
        ggplot2::theme ( axis.text = ggplot2::element_text ( size = 9), axis.title=ggplot2::element_text(size=10 ), panel.background = ggplot2::element_rect(fill="white", colour = "black"),
                         text = ggplot2::element_text(size = 11),
                         axis.text.x = ggplot2::element_text(angle=45, size = 11, hjust = 1))
    }
  }
}

plot_dis_heatmapPanther <- function( type, search,nchars, input, verbose ) {

  if( type == "gene" ) {
    stop( "For this type of chart, a multiple query created with 'disgenetDisease' is required." )
  } else if( type == "disease" ) {

    if( search == "single" | search == "unify" ) {
      stop( "For this type of chart, a multiple query created with 'disgenetDisease' is required.
            For single disease or unify disease is not possible to apply this graphic." )
    }else if ( search == "list" ) {
      tt <- ( input[ , c( 2, 6, 10, 9 ) ] )
      tt <- as.data.frame( transDiseaseClass ( tt ) )
      row.names(tt) <- NULL
      tt$V4 <- as.character(tt$V4)

      tt[,4] <- gsub ( "null", NA, tt[,4] )
      tt$c1.name <- as.character(tt$c1.name)

      nd <- length(unique(as.character(subset(tt, is.na(V4))$c2.symbol)))
      if (nd > 0){
        show( paste0("warning: ", nd , " gene(s) not shown in the plot" ))
      }
      tt <- subset(tt, !is.na(V4))

      for( i in 1:nrow(tt)){
        if(nchar(as.character(tt$V4[i])) > nchars){
          rename <- substr(tt$V4[i], 1, nchars)
          tt$V4[i] <- paste0(rename, "...")
        }
      }
      for( i in 1:nrow(tt)){
        if(nchar(as.character(tt$c1.name[i])) > nchars){
          rename <- substr(tt$c1.name[i], 1, nchars)
          tt$c1.name[i] <- paste0(rename, "...")
        }
      }

      for ( i in 1:length(unique(tt$c1.name))){
        disease1 <- tt[tt$c1.name==unique(tt$c1.name)[i],]
        freq <- as.data.frame ( table ( as.character(disease1$V4 ) ))
        freq$perc <- (freq$Freq/sum(freq$Freq))*100
        freq$dis <- unique(tt$c1.name)[i]

        if ( i == 1){
          rr <- freq
        } else if ( i > 1){
          rr <- rbind(rr, freq)
        }
      }

      rr <- rr [ with ( rr , order(-perc)), ]
      casted_rr <- reshape::cast(rr, Var1 ~ dis, value= "perc"  )
      l <- dim(casted_rr)[2]
      Rowm  <- rowMeans(casted_rr[2:l], na.rm = T)
      casted_rr<- cbind(casted_rr, Rowm)
      ordered<-  casted_rr[order(casted_rr["Rowm"]), "Var1" ]
      rr$Var1  <- factor(rr$Var1  , levels= as.factor(ordered))
      m <- max(rr$perc)
      p <- ggplot2::ggplot(rr, ggplot2::aes(dis, Var1)) + ggplot2::geom_tile(ggplot2::aes(fill = perc), colour = "white") + ggplot2::scale_fill_gradient(name="%",limits = c(0,m), low = "lightgrey",   high = "darkolivegreen", na.value = "white")
      p <- p + ggplot2::theme_grey(base_size = 11) + ggplot2::labs(title = "Protein Class Heatmap", x = "Diseases",  y = "Protein class") + ggplot2::scale_x_discrete(expand = c(0, 0))
      p <- p + ggplot2::theme( plot.margin=unit(x=c(5,15,5,15),units="mm") , axis.line = ggplot2::element_line(size = 0.7, color = "black"), text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(angle=45, size = 11, hjust = 1), panel.background = ggplot2::element_blank())
      p
      #qplot( x = dis, y = Var1, data= rr , fill = perc, geom = "tile" ) + ylab( "Panther Class" ) + xlab( "Disease Name" ) + ggtitle( "Panther class heatmap") + scale_fill_gradient2(low = "#0000FF", high ="#ffa244",guide = "colourbar") + theme ( axis.text = element_text ( size = 9) )
    }

  }
}

plot_dis_heatmapDiseaseClass <- function( type, search, input,nchars, limit, verbose ) {

  if( type == "disease" ) {
    if( search == "dis-dis" ){
      l <- length(unique( input$Disease1 ))
      if (l == 1){
        stop( "For this type of chart, a multiple query created with 'disgenetDisDis' is required." )
      }
      else{
        results <- input [, c (1,2, 5 ,6 )]

        tt <- as.data.frame( transDiseaseClass ( results ) )
        row.names(tt) <- NULL
        tt[,4] <- gsub (" ; ", ";", tt[,4] )
        tt[,4] <- gsub ( "null", NA, tt[,4] )
        tt[,4] <- gsub ( "Wounds and Injuries", NA, tt[,4] )
        tt[,4] <- gsub ( "Animal Diseases", NA, tt[,4] )
        tt$V4 <- as.character(tt$V4)
        tt <- unique(tt[c("NameDisease1", "NameDisease2",  "V4")])
        nd <- dim(subset(tt, is.na(V4)))[1]
        if (nd > 0){
          show( paste0("warning: ", nd , " disease(s) not shown in the plot" ))
        }
        for( i in 1:nrow(tt)){
          if(nchar(as.character(tt$V4[i])) > nchars){
            rename <- substr(tt$V4[i], 1, nchars)
            tt$V4[i] <- paste0(rename, "...")
          }
        }
        ##
        for ( i in 1:length(unique(tt$NameDisease1))){
          disease1 <- tt[tt$NameDisease1==unique(tt$NameDisease1)[i],]
          freq <- as.data.frame ( table ( as.character(disease1$V4 ) ))
          freq$perc <- (freq$Freq/sum(freq$Freq))*100
          freq$dis <- unique(tt$NameDisease1)[i]

          if ( i == 1){
            rr <- freq
          } else if ( i > 1){
            rr <- rbind(rr, freq)
          }
        }

        rr <- rr [ with ( rr , order(-perc)), ]
        casted_rr <- reshape::cast(rr, Var1 ~ dis, value= "perc"  )
        l <- dim(casted_rr)[2]
        Rowm  <- rowSums(casted_rr[2:l], na.rm = T)
        casted_rr<- cbind(casted_rr, Rowm)
        ordered<-  casted_rr[order(casted_rr["Rowm"]), "Var1" ]
        rr$Var1  <- factor(rr$Var1  , levels= as.factor(ordered))
        m <- max(rr$perc)
        p <- ggplot2::ggplot(rr, ggplot2::aes(dis, Var1)) + ggplot2::geom_tile(ggplot2::aes(fill = perc), colour = "white") + ggplot2::scale_fill_gradient(limits = c(0,m), low = "white",   high = "royalblue2", na.value = "white", "Percentage")
        p <- p + ggplot2::theme_grey(base_size = 10) + ggplot2::labs(title = "Gene-Disease Class Heatmap", x = "Genes",  y = "MeSH Disease Class") + ggplot2::scale_x_discrete(expand = c(0, 0))
        p <- p + ggplot2::theme( plot.margin=unit(x=c(5,15,5,15),units="mm") , axis.line = ggplot2::element_line(size = 0.7, color = "black"), text = ggplot2::element_text(size = 9),axis.title=ggplot2::element_text(size=10 ), axis.text.x = ggplot2::element_text(angle=45, size = 10, hjust = 1), panel.background = ggplot2::element_blank())
        p
      }
    }

  } else if( type == "gene" ) {

    if( search == "single" ) {
      stop( "For this type of chart, a multiple query created with 'disgenetGeneList' is required." )
    }else if ( search == "list" ) {
      tt <- ( input[ , c( 6, 5, 2, 3 ) ] )

      tt[,4] <- gsub (" ; ", ";", tt[,4] )
      tt[,4] <- gsub ( "null", NA, tt[,4] )
      tt[,4] <- gsub ( "Wounds and Injuries", NA, tt[,4] )
      tt[,4] <- gsub ( "Animal Diseases", NA, tt[,4] )

      tt <- as.data.frame( transDiseaseClass ( tt ) )
      row.names(tt) <- NULL
      tt$V4 <- as.character(tt$V4)

      nd <- length(unique(as.character(subset(tt, is.na(V4))$c1.name)))
      if (nd > 0){
        show( paste0("warning: ", nd , " disease(s) not shown in the plot" ))
      }
      tt <- subset(tt, !is.na(V4))

      for( i in 1:nrow(tt)){
        if(nchar(as.character(tt$V4[i])) > nchars){
          rename <- substr(tt$V4[i], 1, nchars)
          tt$V4[i] <- paste0(rename, "...")
        }
      }

      for ( i in 1:length(unique(tt$c2.symbol))){
        disease1 <- tt[tt$c2.symbol==unique(tt$c2.symbol)[i],]
        freq <- as.data.frame ( table ( as.character(disease1$V4 ) ))
        freq$perc <- (freq$Freq/sum(freq$Freq))*100
        freq$dis <- unique(tt$c2.symbol)[i]

        if ( i == 1){
          rr <- freq
        } else if ( i > 1){
          rr <- rbind(rr, freq)
        }
      }

      rr <- rr [ with ( rr , order(-perc)), ]
      casted_rr <- reshape::cast(rr, Var1 ~ dis, value= "perc"  )
      l <- dim(casted_rr)[2]
      Rowm  <- rowSums(casted_rr[2:l], na.rm = T)
      casted_rr<- cbind(casted_rr, Rowm)
      ordered<-  casted_rr[order(casted_rr["Rowm"]), "Var1" ]
      rr$Var1  <- factor(rr$Var1  , levels= as.factor(ordered))
      m <- max(rr$perc)
      p <- ggplot2::ggplot(rr, ggplot2::aes(dis, Var1)) + ggplot2::geom_tile(ggplot2::aes(fill = perc), colour = "white") + ggplot2::scale_fill_gradient(limits = c(0,m), low = "white",   high = "mediumorchid4", na.value = "white", "Percentage")
      p <- p + ggplot2::theme_grey(base_size = 10) + ggplot2::labs(title = "Gene-Disease Class Heatmap", x = "Genes",  y = "MeSH Disease Class") + ggplot2::scale_x_discrete(expand = c(0, 0))
      p <- p + ggplot2::theme( plot.margin=unit(x=c(5,15,5,15),units="mm") , axis.line = ggplot2::element_line(size = 0.7, color = "black"), text = ggplot2::element_text(size = 9),axis.title=ggplot2::element_text(size=10 ), axis.text.x = ggplot2::element_text(angle=45, size = 10, hjust = 1), panel.background = ggplot2::element_blank())
      p
    }
  }
  else if( type == "variant-disease"  ) {
    stop( "For this type of chart, a multiple query created with 'variant2disease' is required." )
  }
  else if( type == "variant"  ) {

    if( search == "single" ) {
      stop( "For this type of chart, a multiple query created with 'variant2disease' is required." )
    }

    tt <- unique( input[ , c( "c0.snpId" ,"c0.diseaseId" , "c2.name", "c2.diseaseClassName") ] )
    tt[,4] <- gsub (" ; ", ";", tt[,4] )
    tt[,4] <- gsub ( "null", NA, tt[,4] )
    tt[,4] <- gsub ( "Wounds and Injuries", NA, tt[,4] )
    tt[,4] <- gsub ( "Animal Diseases", NA, tt[,4] )

    nd <- length(unique(as.character(subset(tt, is.na(c2.diseaseClassName))$c0.diseaseId)))
    if (nd > 0){
      show( paste0("warning: ", nd , " disease(s) not shown in the plot" ))
    }
    tt <- subset(tt, !is.na(c2.diseaseClassName))
    tt <- as.data.frame( transDiseaseClass ( tt  ))
    row.names(tt) <- NULL
    tt$V4 <- as.character(tt$V4)
    for( i in 1:nrow(tt)){
      if(nchar(as.character(tt$V4[i])) > nchars){
        rename <- substr(tt$V4[i], 1, nchars)
        tt$V4[i] <- paste0(rename, "...")
      }
    }

    for ( i in 1:length(unique(tt$c0.snpId))){
      disease1 <- tt[tt$c0.snpId==unique(tt$c0.snpId)[i],]
      freq <- as.data.frame ( table ( as.character(disease1$V4 ) ))
      freq$perc <- (freq$Freq/sum(freq$Freq))*100
      freq$snp <- unique(tt$c0.snpId)[i]

      if ( i == 1){
        rr <- freq
      } else if ( i > 1){
        rr <- rbind(rr, freq)
      }
    }

    rr <- rr [ with ( rr , order(-perc)), ]
    casted_rr <- reshape::cast(rr, Var1 ~ snp, value= "perc"  )
    l <- dim(casted_rr)[2]
    Rowm  <- rowSums(casted_rr[2:l], na.rm = T)
    casted_rr<- cbind(casted_rr, Rowm)
    ordered<-  casted_rr[order(casted_rr["Rowm"]), "Var1" ]
    rr$Var1  <- factor(rr$Var1  , levels= as.factor(ordered))
    m <- max(rr$perc)
    p <- ggplot2::ggplot(rr, ggplot2::aes(snp, Var1)) + ggplot2::geom_tile(ggplot2::aes(fill = perc), colour = "white") + ggplot2::scale_fill_gradient(limits = c(0,m), low = "white",   high = "darkolivegreen", na.value = "white", "Percentage")
    p <- p + ggplot2::theme_grey(base_size = 10) + ggplot2::labs(title = "Variant-Disease Class Heatmap", x = "Variants",  y = "MeSH Disease Class") + ggplot2::scale_x_discrete(expand = c(0, 0))
    p <- p + ggplot2::theme( plot.margin=unit(x=c(5,15,5,15),units="mm") , axis.line = ggplot2::element_line(size = 0.7, color = "black"), text = ggplot2::element_text(size = 9),axis.title=ggplot2::element_text(size=10 ), axis.text.x = ggplot2::element_text(angle=45, size = 10, hjust = 1), panel.background = ggplot2::element_blank())
    p
  }
}

plot_dis_Panther <- function ( input, layout, search, verbose ) {

  if( search == "single"){
    results <- input [, c ( 1, 2, 6, 9 )]
    results <- as.data.frame( transDiseaseClass ( results  ))
    results[,4] <- gsub ( "null", NA, results[,4] )


    row.names(results) <- NULL
    nd <- dim(subset(results, is.na(V4)))[1]
    if (nd > 0){
      show( paste0("warning: ", nd , " gene(s) not shown in the plot" ))
    }
    results <- subset(results, !is.na(V4))

    freq <- as.data.frame ( table ( as.character(results$V4 ) ))
    freq$perc <- (freq$Freq/sum(freq$Freq))*100

    rw <- matrix(c(as.character(results[1,2]), 10, 10),ncol=3)
    colnames(rw) <- colnames(freq)
    freq <- rbind(freq,rw)

    sizes <- as.numeric(freq[ , 3 ])
    names( sizes ) <- freq[ , 1 ]

    #plot
    edges <- data.frame( results[ , 2 ], results[ , 4 ] )
    netw  <- igraph::graph.data.frame( edges, directed = FALSE )
    netw  <- igraph::simplify( netw )
    lay   <- layout( netw )

    if( verbose ) {
      message( "The network contains ", igraph::vcount( netw ), " nodes and ", igraph::ecount( netw ), " edges." )
    }

    pantherClass <- unique( results$V4 )
    ttl <- paste0("Protein Classes associated to ", results$c1.name[ 1 ] )

    igraph::plot.igraph( netw,
                         vertex.frame.color  = "white",
                         layout              = lay,
                         vertex.color        = ifelse( igraph::V( netw )$name %in% pantherClass, "darkolivegreen", "royalblue2" ),
                         vertex.label.dist   = 0.2,      #puts the name labels slightly off the dots
                         vertex.frame.color  = 'blue', #the color of the border of the dots
                         vertex.label.color  = 'black',#the color of the name labels
                         vertex.label.font   = 0,      #the font of the name labels
                         vertex.label        = igraph::V( netw )$names, #specifies the lables of the vertices. in this case the 'name' attribute is used
                         edge.color          = "darkgrey",
                         edge.arrow.size     = 0.5,
                         edge.width          = 0.5,
                         vertex.size         = sizes[ igraph::V( netw )$name ],
                         vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
                         main                = ttl
    )
  }
  if( search == "unify"){
    results <- input [, c ( 1, 12, 6, 9 )]

    results[,4] <- gsub ( "null", NA, results[,4] )
    results[,4] <- gsub ( "Unclassified", NA, results[,4] )
    results <- as.data.frame( transDiseaseClass ( results ) )
    row.names(results) <- NULL
    nd <- dim(subset(results, is.na(V4)))[1]
    if (nd > 0){
      show( paste0("warning: ", nd , " gene(s) not shown in the plot" ))
    }
    results <- subset(results, !is.na(V4))


    freq <- as.data.frame ( table ( as.character(results$V4 ) ))
    freq$perc <- (freq$Freq/sum(freq$Freq))*100

    rw <- matrix(c(as.character(results[1,2]), 10, 10),ncol=3)
    colnames(rw) <- colnames(freq)
    freq <- rbind(freq,rw)

    sizes <- as.numeric(freq[ , 3 ])
    names( sizes ) <- freq[ , 1 ]

    #plot
    edges <- data.frame( results[ , 2 ], results[ , 4 ] )
    netw  <- igraph::graph.data.frame( edges, directed = FALSE )
    netw  <- igraph::simplify( netw )
    lay   <- layout( netw )

    if( verbose ) {
      message( "The network contains ", igraph::vcount( netw ), " nodes and ", igraph::ecount( netw ), " edges." )
    }

    pantherClass <- unique( results$V4 )
    ttl <- paste0("Protein Classes associated to ", results$diseaseIdentifier[ 1 ] )

    igraph::plot.igraph( netw,
                         vertex.frame.color  = "white",
                         layout              = lay,
                         vertex.color        = ifelse( igraph::V( netw )$name %in% pantherClass, "darkolivegreen", "royalblue2" ),
                         vertex.label.dist   = 0,      #puts the name labels slightly off the dots
                         vertex.frame.color  = 'blue', #the color of the border of the dots
                         vertex.label.color  = 'black',#the color of the name labels
                         vertex.label.font   = 0,      #the font of the name labels
                         vertex.label        = igraph::V( netw )$names, #specifies the lables of the vertices. in this case the 'name' attribute is used
                         edge.color          = "darkgrey",
                         edge.arrow.size     = 0.5,
                         edge.width          = 0.5,
                         vertex.size         = sizes[ igraph::V( netw )$name ],
                         vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
                         main                = ttl
    )
  }
}

plot_dis_barplot <- function ( input, type, search, verbose, limit ) {

  if( search != "dis-dis" ) {
    stop( "For this type of chart, a query created with 'disgenetDisDis' is required." )
  }
  #input <- input[c(1:20),]
  input$NameDisease1 <- as.character(input$NameDisease1)
  input$NameDisease2 <- as.character(input$NameDisease2)
  #   for( i in 1:nrow(input)){
  #     if(nchar(as.character(input$NameDisease1[i])) > 50){
  #       rename1 <- substr(input$NameDisease1[i], 1, 50)
  #       input$NameDisease1[i] <- paste0(rename1, "...")
  #     }
  #     if(nchar(as.character(input$NameDisease2[i])) > 50){
  #       rename2 <- substr(input$NameDisease2[i], 1, 50)
  #       input$NameDisease2[i] <- paste0(rename2, "...")
  #     }
  #   }

  input$commonGenes <- as.numeric(input$commonGenes)
  input$JaccardIndex <- as.numeric(input$JaccardIndex)
  input <- subset(input, commonGenes>0)
  input <- input [ with ( input , order(-commonGenes)), ]

  if (dim(input)[1]> limit){
    input <- input[c(1:limit),]
  }

  input$pairs <- paste(input$NameDisease1, input$NameDisease2, sep = " - ")
  orderedDiseases <- input[order(-as.numeric(input$commonGenes)), "pairs"  ]
  input$pairs <- factor(input$pairs , levels= as.factor(orderedDiseases))

  p <- ggplot2::ggplot (input, ggplot2::aes ( x = pairs, y = as.numeric(commonGenes) ), order = as.numeric(commonGenes) ) +
    ggplot2::geom_bar ( stat = "identity", fill = "grey" ) +
    ggplot2::labs ( title = paste0 ( "Disease-Disease Barplot" ) , x = "Disease Pairs", y = "# of shared genes")


  p <- p + ggplot2::theme_classic( ) + ggplot2::theme( plot.margin = unit ( x = c ( 5, 15, 5, 15 ), units = "mm" ),
                                                       axis.line = ggplot2::element_line ( size = 0.7, color = "black" ), text = ggplot2::element_text ( size = 14 ) ,
                                                       axis.text.x = ggplot2::element_text ( angle = 45, size = 9, hjust = 1 ))


  p <- p + ggplot2::geom_point(ggplot2::aes(x = pairs , y =as.numeric(JaccardIndex) *100))
  p <- p + ggplot2::geom_text(ggplot2::aes(x = pairs , y =as.numeric(JaccardIndex) *100 + 1), label = c(round(as.numeric(input$JaccardIndex), digits = 3 )), size=3)
  p

  }

plot_dis_venn <- function( search, type, input, verbose, database ) {

  if( search != "dis-dis" ) {
    stop("For this type of chart, objects comming from disgenetDisDis function is required." )
  }

  tt <- input
  diseases <- c(unique( as.character(unique(tt[, "Disease1"]))) )

  if( length(diseases) > 2 ) {
    stop( "For this type of chart, the diseases list should have 2 diseases" )
  }

  gen1 <-  disgenetDisease(diseases[1], database = database  )
  gen1 <- gen1@qresult$c2.geneId

  gen2 <-  disgenetDisease(diseases[2], database = database )
  gen2 <- gen2@qresult$c2.geneId



  grid.newpage()
  venn.plot <- draw.pairwise.venn(area1        = length(gen1),
                                  area2        = length(gen2),
                                  cross.area   = length(intersect(gen1, gen2)),
                                  scaled       = T,
                                  category     = c( tt[tt$Disease1 == diseases[1], "NameDisease1" ][1],  tt[tt$Disease1 == diseases[2], "NameDisease1" ][1]),
                                  scales       = TRUE,
                                  fill         = c("blue", "red"),
                                  alpha        = 0.3,
                                  lty          = "blank",
                                  cex          = 2,
                                  cat.cex      = 1.5,
                                  cat.pos      = c(0.5 , 1),
                                  cat.dist     = 0.04,
                                  cat.just     = list(c(0.5, 1), c(0.5, 1)),
                                  ext.pos      = 30,
                                  ext.dist     = -0.5,
                                  ext.length   = 0.85,
                                  ext.line.lwd = 2,
                                  ext.line.lty = "dashed")
  grid.draw(venn.plot)

}

plot_dis_pathway <- function( search, type, input, verbose , limit = 10) {

  if( type != "pathway" ) {
    stop( "For this type of chart, a multiple query created with 'disease2pathway' is required." )
  }

  input <- subset(input, pathwayName !="null")
  if (dim(input)[1] > limit){
    input <- input[1:limit, ]
  }

  input$pathwayName <- ifelse( nchar(input$pathwayName) > 30,  paste0(substr(input$pathwayName, 1, 25), "..."),
                    as.character(input$pathwayName) )
  rr <- as.data.frame(table(input$pathwayName))
  #rr <- rr[order(-as.numeric(rr$Freq)),   ]
  rr$Var1 <-factor(rr$Var1, levels=rr[order(-rr$Freq), "Var1"])
  p <- ggplot2::ggplot (rr, ggplot2::aes ( x = Var1, y = Freq ), order = as.numeric(Freq) ) +
    ggplot2::geom_bar ( stat = "identity", fill = "grey" ) +
    ggplot2::labs ( title = paste0 ( "Disease-Pathway Barplot" ) , x = "Pathway", y = "# of genes")


  p <- p + ggplot2::theme_classic( ) + ggplot2::theme( plot.margin = unit ( x = c ( 5, 15, 5, 15 ), units = "mm" ),
                                                       axis.line = ggplot2::element_line ( size = 0.7, color = "black" ), text = ggplot2::element_text ( size = 14 ) ,
                                                       axis.text.x = ggplot2::element_text ( angle = 45, size = 9, hjust = 1 ))


  p

}

plot_dis_comorbidity <- function( input, type, layout, search, verbose, prop, limit ) {

  if( search != "dis-dis" ) {
    stop("For this type of chart, objects comming from disgenetDisDis function is required." )
  } else  if( type != "disease" ) {
    stop("For this type of chart, objects comming from disgenetDisDis function is required." )
  }

  #   disNum <- length ( unique ( input$c1.name) )
  #   diseases <- as.character(unique ( input$c1.name) )
  #   pairs <- ( disNum *  disNum - disNum ) / 2
  #
  #   tb <- as.data.frame(setNames(replicate(5,numeric(0), simplify = F), c("dis1", "dis2", "gen1", "gen2", "intr")))
  #   cnt <- 1
  #
  #   for ( i in 1:disNum){
  #     dis1 <- input[ input$c1.name == unique( input$c1.name )[ i ] ,]
  #     gen1 <- unique( dis1$c2.name)
  #
  #     for ( j in 2:disNum){
  #       if ( i < j){
  #         dis2 <- input[ input$c1.name == unique( input$c1.name )[ j ] ,]
  #         gen2 <- unique( dis2$c2.name)
  #         tb[ cnt, ] <- c(diseases[ i ], diseases[ j ], length( gen1 ), length( gen2 ), length( intersect( gen1, gen2 ) ))
  #         cnt <- cnt + 1
  #       }
  #     }
  #   }
  input <- subset(input, commonGenes>0)
  a <- nrow(input)
  if (dim(input)[1]> limit){
    input <- input[c(1:limit),]
    show( paste0("warning: dataframe of ", a , " rows has been reduced to ", limit, " rows." ))
  }
  edges <- data.frame( input[ , "NameDisease1" ], input[ , "NameDisease2" ] )
  netw  <- graph.data.frame( edges, directed = FALSE )
  netw  <- simplify( netw )
  lay   <- layout(netw)

  if( verbose ) {
    message( "The network contains ", igraph::vcount( netw ), " nodes and ", igraph::ecount( netw ), " edges." )
  }
  if( search == "single" ) {
    stop ( "This network is only available for objects comming from disgenetDisease function where search is a list of diseases" )
  }else if ( search == "list" ){
    if ( type == "gene"){
      stop ( "This network is only available for objects comming from disgenetDisease function where search is a list of diseases" )
    }
  }
  dn1 <- unique(input$NameDisease1)
  if (length(dn1) == 1){
    ttl <- paste0( "Disease-Disease Network for ", dn1 )
  }
  else{
    ttl <- paste0( "Disease-Disease Network" )
  }

  ds1 <- input[,c( "NameDisease1", "NGenes1" ) ]
  colnames( ds1 ) <- c("dis", "gen")
  ds2 <- input[, c("NameDisease2", "NGenes2") ]
  colnames( ds2 ) <- c("dis", "gen")

  fnl <- rbind ( ds1, ds2)
  fnl <- fnl[!duplicated(fnl), ]

  sizes <- as.numeric( fnl[ , 2 ] )
  names( sizes ) <- fnl[ , 1 ]

  igraph::plot.igraph( netw,
                       vertex.frame.color  = "white",
                       layout              = lay,
                       vertex.color        = ifelse( igraph::V( netw )$name %in% ds1$dis, "royalblue2", "#b4daff" ),
                       vertex.label.dist   = 0,      #puts the name labels slightly off the dots
                       vertex.frame.color  = 'blue', #the color of the border of the dots
                       vertex.label.color  = 'black',#the color of the name labels
                       vertex.label.font   = 0,      #the font of the name labels
                       vertex.label        = igraph::V( netw )$names, #specifies the lables of the vertices. in this case the 'name' attribute is used
                       edge.color          = "black",
                       edge.width          = prop * (as.numeric( input$commonGenes ))^1/2,
                       edge.arrow.size     = 0.5,
                       vertex.size         = prop *30*log2(as.numeric(sizes[ igraph::V( netw )$name ])+1),
                       vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
                       main                = ttl
  )
}

plot_dis_dis_associations <- function( input, type, layout, search, verbose, prop, limit, database ) {

  if( search != "dis-dis" ) {
    stop("For this type of chart, objects comming from disgenetDisDis function is required." )
  } else  if( type != "disease" ) {
    stop("For this type of chart, objects comming from disgenetDisDis function is required." )
  }
  l <- unique(input$Disease1)
  if (length(l) == 1){
    stop("For this type of chart, objects generated with disgenetDisDis function, and a disease list, are required ." )

  }

  input <- subset(input, commonGenes>0)
  disNum <- length ( unique ( input$Disease1) )
  diseases <- as.character(unique (  input$Disease1) )

  res <- disgenetDisease(disease = diseases,
                         database = as.character(database))
  res <- res@qresult
  pairs <- ( disNum *  disNum - disNum ) / 2
  tb <- as.data.frame(setNames(replicate(7,numeric(0), simplify = F), c("disease1","NameDisease1", "disease2","NameDisease2", "NGenes1", "NGenes2", "intr")))
  cnt <- 1
  for ( i in 1:disNum){
    dis1 <- diseases[i]
    nd1 <- as.character (  res[res[,1]==diseases[i],"c1.name"][1]  )
    gen1 <- unique( subset(res, c1.cui == dis1 )$c2.geneId)
    for ( j in 2:disNum){
      if ( i < j){
        dis2 <- diseases[j]
        nd2 <- as.character (  res[res[,1]==diseases[j],"c1.name"][1]  )
        gen2 <- unique( subset(res, c1.cui == dis2 )$c2.geneId)
        tb[ cnt, ] <- c(diseases[ i ], nd1, diseases[ j ],nd2, length( gen1 ), length( gen2 ), length( intersect( gen1, gen2 ) ))
        cnt <- cnt + 1
      }
    }
  }
  tb <- tb[tb[,"intr"]>0 ,]
  finaldiseases <- unique( c(as.character(tb$disease1), as.character(tb$disease2)))
  #   if (dim(input)[1]> limit){
  #     input <- input[c(1:limit),]
  #     show( paste0("warning: dataframe of ", a , " rows has been reduced to ", limit, " rows." ))
  #   }
  edges <- data.frame( tb[ , "NameDisease1" ], tb[ , "NameDisease2" ] )
  netw  <- graph.data.frame( edges, directed = FALSE )
  netw  <- simplify( netw )
  lay   <- layout(netw)

  if( verbose ) {
    message( "The network contains ", igraph::vcount( netw ), " nodes and ", igraph::ecount( netw ), " edges." )
  }
  if( search == "single" ) {
    stop ( "This network is only available for objects comming from disgenetDisease function where search is a list of diseases" )
  }else if ( search == "list" ){
    if ( type == "gene"){
      stop ( "This network is only available for objects comming from disgenetDisease function where search is a list of diseases" )
    }
  }
  ttl <- paste0( "Disease-Disease Network" )

  #input<- subset(input,Disease1 %in% finaldiseases & Disease2 %in% finaldiseases  )
  ds1 <- input[,c( "NameDisease1", "NGenes1" ) ]
  colnames( ds1 ) <- c("dis", "gen")
  ds2 <- input[, c("NameDisease2", "NGenes2") ]
  colnames( ds2 ) <- c("dis", "gen")

  fnl <- rbind ( ds1, ds2)
  fnl <- fnl[!duplicated(fnl), ]

  sizes <- as.numeric( fnl[ , 2 ] )
  names( sizes ) <- fnl[ , 1 ]

  igraph::plot.igraph( netw,
                       vertex.frame.color  = "white",
                       layout              = lay,
                       vertex.color        =  "royalblue2" ,
                       vertex.label.dist   = 0,      #puts the name labels slightly off the dots
                       vertex.frame.color  = 'blue', #the color of the border of the dots
                       vertex.label.color  = 'black',#the color of the name labels
                       vertex.label.font   = 0,      #the font of the name labels
                       vertex.label        = igraph::V( netw )$names, #specifies the lables of the vertices. in this case the 'name' attribute is used
                       edge.color          = "black",
                       edge.width          = prop * (as.numeric( input$commonGenes ))^1/2,
                       edge.arrow.size     = 0.5,
                       # vertex.size         = 10,
                       vertex.size         = prop *30*log2(as.numeric(sizes[ igraph::V( netw )$name ])+1),
                       vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
                       main                = ttl
  )
}

plot_variant_disease_heatmap <- function( type, search, input, cutoff=cutoff, limit, verbose ) {

  if( search == "single" ) {
    stop( "For this type of chart, a multiple query created with 'variant2disease'is required. " )
  }

  if( type == "variant-disease" ) {
    tt <- unique( input[ , c( "diseaseId","diseaseName", "snp","score" ) ] )
    tt$snp <- as.character(tt$snp)
    tt$diseaseName <- as.character(tt$diseaseName)

    tt <- tt [ tt$score>= cutoff, ]
    tt <- tt [ with ( tt , order(-score)), ]
    if ( dim( tt )[ 1 ] > limit ){
      tt <- tt[ 1:limit ,]
      show( paste0("Dataframe of ", nrow(input), " rows has been reduced to ", limit, " rows." ))
    }

    #       for( i in 1:nrow(tt)){
    #         if(nchar(as.character(tt$diseaseName[i])) > 30){
    #           rename <- substr(tt$diseaseName[i], 1, 30)
    #           tt$c1.name[i] <- paste0(rename, "...")
    #         }
    #       }

    casted_tt <- reshape::cast( tt, snp ~diseaseName, value= "score", length )
    l <- dim( casted_tt )[ 2 ]
    casted_tt[ is.na( casted_tt ) ] <- 0
    res <- hclust( dist( casted_tt[ 2:l ] ),  method = "ward.D")
    casted_tt$clus <- cutree( res, k = 5 )
    ordered <- casted_tt[order(casted_tt$clus), "snp"]

    #tt$diseaseName  <- factor(tt$diseaseName , levels= as.factor(ordered))
    tt$snp  <- factor(tt$snp , levels= as.factor(ordered))

    p<- ggplot2::qplot( x = snp, y = diseaseName, data= tt , fill = score, geom = "tile" ) +
      ggplot2::ylab( "Diseases" ) + ggplot2::xlab( "Variants" ) +
      ggplot2::ggtitle( "Variant-Disease Heatmap") +
      ggplot2::scale_fill_gradient2(low = "#0000FF", high ="#2f30c9",guide = "colourbar", "Score") +
      ggplot2::theme ( axis.text = ggplot2::element_text ( size = 9),
                       axis.title=ggplot2::element_text(size=10 ),
                       panel.background = ggplot2::element_rect(fill="white", colour = "black"),
                       text = ggplot2::element_text(size = 11),
                       axis.text.x = ggplot2::element_text(angle=45, size = 11, hjust = 1)
      )
    print(p)

  } else  if( type == "variant" ) {

    tt <- unique( input[ , c( "diseaseId","diseaseName", "snp","score" ) ] )
    tt$snp <- as.character(tt$snp)
    tt$diseaseName <- as.character(tt$diseaseName)

    tt <- tt [ tt$score>= cutoff, ]
    if ( dim( tt )[ 1 ] > limit ){
      tt <- tt[ 1:limit ,]
      show( paste0("warning: dataframe of ", nrow(input), " rows has been reduced to ", limit, " rows." ))
    }

    #       for( i in 1:nrow(tt)){
    #
    #         if(nchar(as.character(tt$c1.name[i])) > 30){
    #           rename <- substr(tt$c1.name[i], 1, 30)
    #           tt$c1.name[i] <- paste0(rename, "...")
    #         }
    #       }

    casted_tt <- reshape::cast( tt, diseaseName~snp, value= "score", length )

    l <- dim(casted_tt)[2]
    casted_tt[is.na( casted_tt)] <- 0
    res <- hclust(dist(casted_tt[2:l]),  method = "ward.D")
    casted_tt$clus <- cutree(res, k = 5)
    ordered <- casted_tt[order(casted_tt$clus), "diseaseName"]

    tt$diseaseName  <- factor(tt$diseaseName , levels= as.factor(ordered))
    p <- ggplot2::ggplot( tt,  ggplot2::aes ( x = diseaseName, y = snp ) )
    p <- p + ggplot2::geom_tile(ggplot2::aes(fill = score), colour = "white")
    p <- p + ggplot2::labs(title = "Variant-Disease Heatmap", x = "Diseases",  y = "Variants")
    p <- p + ggplot2::scale_fill_gradient2(low = "#0000FF", high ="#2f30c9",guide = "colourbar", "Score")
    p <- p + ggplot2::theme ( axis.text = ggplot2::element_text ( size = 9),
                     axis.title=ggplot2::element_text(size=10 ),
                     panel.background = ggplot2::element_rect(fill="white", colour = "black"),
                     text = ggplot2::element_text(size = 11),
                     axis.text.x = ggplot2::element_text(angle=45, size = 11, hjust = 1) )

    print(p)
   }
}
