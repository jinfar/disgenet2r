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
#' @param class Type of the drawn chart. By default it is \code{"network"} but
#' it also can be \code{"Piechart"}, \code{"DiseaseClass"}, \code{"Venn"}
#' \code{"Heatmap"}, \code{"ProteinClass"}, \code{"Barplot"}, or \code{"Lollipop"}.
#' @param verbose By default \code{FALSE}. If set to \code{TRUE} information
#' on the drawing process will be shown.
#' @param ... Passed to inner functions for different plots.
#' @return A plot for \code{DataGeNET.DGN}.
#' @examples
#' \dontrun{
#' Being x an DataGeNET.DGN
#' qr <- plot(x) # Get number of unique diseases
#' }
#' @export plot

setMethod(
  f = "plot",
  signature = "DataGeNET.DGN",
  definition = function( x, y, ... ) {
    plot_disgenet( x, ... )
  }
)


plot_disgenet <- function( object, layout = "layout.fruchterman.reingold", class = "Network",
                          prop = 1, limit = 50, verbose = FALSE, cutoff = 0, count = 1,
                          showGenes =F, nboot=100 , nchars=30, interactive=F, pvalue=F, ncores=1) {

  if(dim(object@qresult)[1] == 0 ) {
    stop( "Empty object!" )
  }

  if( !class %in% c( "Network",  "Heatmap", "DiseaseClass",  "ProteinClass", "Barplot", "Venn", "Diseasome", "Piechart", "Enrichment", "Lollipop" ) ) {
    stop( "Invalid content of argument 'class'." )
  }

  if (class == "Network" ){
    if( object@type %in%  c("gene-disease", "disease-gene")) {
      plot_gda_network( search = object@search, type = object@type, input = object@qresult, layout = layout, verbose = verbose, prop = prop , limit = limit, interactive= interactive)
    }
    else if( object@type %in%  c("variant-disease", "disease-variant")) {
      if (showGenes == T){
        plot_vgda_network( search = object@search, type = object@type, input = object@qresult, layout = layout, verbose = verbose, prop = prop , interactive = interactive )
      }
      else {
        plot_vda_network( search = object@search, type = object@type, input = object@qresult, layout = layout, verbose = verbose, prop = prop, interactive = interactive  )
      }
    }
    else if( object@type %in%  c("disease-disease-gene", "disease-disease-variant")) {
        # if(object@search == "single"){
        #   plot_dis_comorbidity( search= object@search, type = object@type, input = object@qresult, layout = layout, verbose = verbose ,  prop = prop, limit = limit)
        # }
        # else if(object@search == "list") {
          plot_dis_comorbidity( search= object@search, type = object@type, input = object@qresult, layout = layout, verbose = verbose ,  prop = prop, limit = limit, interactive = interactive )
        #}
    }
    else{
      stop( "Invalid content of argument 'type' for class = 'Network'." )
    }
  }
  if (class == "Heatmap" ){
    if( object@type %in%  c("variant-disease", "disease-variant") & object@search == "list") {
      plot_variant_disease_heatmap( type= object@type, search= object@search, input = object@qresult,nchars = nchars, cutoff = cutoff, verbose = verbose, limit = limit )
    }
    else if( object@type %in%  c("gene-disease", "disease-gene")) {
      plot_gda_heatmap( type= object@type, search= object@search, input = object@qresult, cutoff = cutoff, nchars=nchars, verbose = verbose, limit = limit )    }
    else if( object@type %in%  c("disease-disease-gene", "disease-disease-variant")) {
      plot_dda_network( search = object@search, type = object@type, input = object@qresult, layout = layout, verbose = verbose, prop = prop , limit = limit)
    }
    else{
      stop( "Invalid content of argument 'type' for class = 'Heatmap'." )
    }
  }
  if (class == "DiseaseClass" ){
    if( object@type %in%  c("variant-disease")) {
      if (object@search == "single"){
        plot_disease_class( search = object@search, type = object@type, input = object@qresult, nchars = nchars, layout = layout, verbose = verbose, prop = prop )
      }
      else if (object@search == "list"){
        plot_gda_heatmap_DiseaseClass( type= object@type, search= object@search, input = object@qresult, nchars=nchars, verbose = verbose, limit = limit)
      }
    }
    else if( object@type %in%  c("gene-disease")) {
      if (object@search == "single"){
        plot_disease_class( search= object@search, type = object@type, input = object@qresult, layout = layout, nchars=nchars,  verbose = verbose, prop = prop )
      }
      else if (object@search == "list"){
        plot_gda_heatmap_DiseaseClass( search = object@search, type = object@type, input = object@qresult, nchars=nchars, verbose = verbose)
      }
    }
    else if( object@type %in%  c("disease-disease-gene", "disease-disease-variant")) {
      plot_disease_class( search = object@search, type = object@type, input = object@qresult, layout = layout, nchars= nchars, verbose = verbose, prop = prop)
    }
    else{
      stop( "Invalid content of argument 'type' for class = 'DiseaseClass'." )
    }
  }
  if (class == "ProteinClass" ){
    if( object@type %in%  c("variant-disease")) {
      stop( "Invalid object for Protein class!" )
    }
    else if( object@type %in%  c("disease-gene")) {
      if (object@search == "single"){
        plot_gda_Panther(  search= object@search, input = object@qresult, layout = layout, verbose = verbose )
      }
      else if (object@search == "list"){
        plot_gda_heatmapPanther(  type= object@type, search= object@search, input = object@qresult, nchars =nchars, verbose = verbose )
      }
    }
    else if( object@type %in%  c("disease-disease-gene", "disease-disease-variant")) {
      stop( "Invalid object for Protein class!" )
    }
    else{
      stop( "Invalid content of argument 'type' for class = 'ProteinClass'." )
    }
  }
  else if( class == "Barplot" ) {
    if( ! object@type %in%  c("disease-disease-gene", "disease-disease-variant")) {
      stop( "Invalid object for Barplot!" )
    }
    plot_barplot(  type = object@type, search= object@search, input = object@qresult, verbose = verbose, limit = limit )
  }

  else if( class == "Venn" ) {
    if( ! object@type %in%  c("disease-disease-gene", "disease-disease-variant")) {
      stop( "Invalid object for Barplot!" )
    }
    plot_venn(  type= object@type, search= object@search, input = object@qresult, verbose = verbose, database= object@database )
  }
  else if (class == "Diseasome" ){
    if( ! object@type %in%  c("disease-disease-gene", "disease-disease-variant")) {
      stop( "Invalid object for Barplot!" )
    }
    plot_diseasome( search= object@search, type = object@type, input = object@qresult, layout = layout, verbose = verbose ,  prop = prop, limit = limit,  database= object@database, pvalue = pvalue, nboot = nboot, interactive = interactive , cutoff = cutoff, ncores = ncores )

  }
  else if (class == "Enrichment"){
    if( ! object@type %in%  c("disease-enrichment")) {
      stop( "Invalid object for plot Enrichment!" )
    }
    plot_enrichment( input = object@qresult, count = count,  cutoff = cutoff ,nchars = nchars, limit = limit )
  }
}

# plot_variant_gene( search = object@search, type = object@type, input = object@qresult, layout = layout, verbose = verbose, prop = prop )

######
plot_gda_network <- function( input, type, layout, prop, search, verbose ,  limit, interactive) {
  edges <- data.frame( input[ , c("gene_symbol", "disease_name") ])
  diseases <- unique( input$disease_name )
  netw  <- igraph::graph.data.frame( edges, directed = FALSE )
  netw  <- igraph::simplify( netw )
  #lay   <- layout( netw )
  lay <- apply_layout(layout, netw)
  if( verbose ) {
    message( "The network contains ", igraph::vcount( netw ), " nodes and ", igraph::ecount( netw ), " edges." )
  }


  if( search == "list" ) {
    ttl <- "Gene-Disease Network"
  }else if ( search == "single" ){
    if ( type == "disease-gene"){
      ttl <- paste0( "Genes associated to ", input$disease_name[ 1 ]  )
    }else if ( type == "gene-disease"   ){
      ttl <- paste0( "Diseases associated to ", input$gene_symbol[ 1 ] )
    }
    else{
      ttl <- ""
      show( paste0("warning: ",type , " not correct" ))
    }
  }

  if (interactive == F){
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
                       edge.width          = 1 + ( prop * input$score ),
                       edge.arrow.size     = 0.5,
                       vertex.size         = 15,
                       vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
                       main                = ttl
  )
  }else if (interactive == TRUE){
    igraph::tkplot(netw,  vertex.frame.color  = "white",
                   layout              = lay,
                   vertex.color        = ifelse( igraph::V( netw )$name %in% diseases, "royalblue2", "#ff349a" ),
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
}

plot_vda_network<- function( input, type, layout, prop, search, verbose, interactive ) {

  input <- unique(input[c("score", "diseaseid", "disease_name", "variantid")])
  l <- length(unique(input$variantid))
  input$score <- as.numeric(as.character(input$score))
  input <-  aggregate(score ~ diseaseid + disease_name + variantid, input , max)
  if(l > 1 ) {
    ttl <- "Variant-Disease Network"
  }else if ( l == 1 ){
    ttl <- paste0( "Diseases associated to ", input$variantid[ 1 ] )
  }

  edges <- data.frame( input[ , "disease_name"], input[ , "variantid"] )
  netw  <- igraph::graph.data.frame( edges, directed = FALSE )
  netw  <- igraph::simplify( netw )
  #lay   <- layout( netw )
  lay <- apply_layout(layout, netw)

  diseases <- unique( input$disease_name )
  if (interactive == F){
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
  } else {
    igraph::tkplot( netw,
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
                    main = ttl)
  }
}

plot_vgda_network<- function( input, type, layout, prop, search, verbose, interactive  ) {

  results <- unique(input[c("score", "variantid",  "disease_name", "gene_symbol")])
  #input <- data.frame(lapply(input, as.character), stringsAsFactors=FALSE)
  l <- length(unique(input$variantid))

  #input[input$c0.score < 0.01,"c0.score" ] <- 0.09
  results <- tidyr::separate_rows(results, gene_symbol, sep = ",")
  nd <- length(unique(subset(results, gene_symbol!="")$variantid))
  results <- subset(results, gene_symbol!="")

  results$score <- as.numeric(as.character(results$score))
  aa <-  aggregate(score ~  gene_symbol +variantid , results , max)
  bb <-  aggregate(score ~ disease_name +variantid , results , max)
  cc <- unique(input[c("disease_name", "gene_symbol",  "score")])
  colnames(aa) <- c("A1", "A2", "score")
  colnames(bb) <- c("A1", "A2", "score")
  colnames(cc) <- c("A1", "A2", "score")

  results <- rbind(aa,bb)
  results <- rbind(results,cc)

  if( l > 1 ) {
    ttl <- "Variant-Gene-Disease Network"
  }else if ( l == 1 ){
    ttl <- paste0( "Variant-Gene-Disease Network for ", results$variantid[ 1 ] )
  }
  results <- subset( results, !is.na(A1) & !is.na(A2))
  edges <- data.frame( results[ ,1], results[ , 2] )
  netw  <- igraph::graph.data.frame( edges, directed = FALSE )
  netw  <- igraph::simplify( netw )
  #lay   <- layout( netw )
  lay <- apply_layout(layout, netw)
  genes <- unique( input$gene_symbol )
  diseases <- unique( input$disease_name )
  variants <- unique( input$variantid )
 if(interactive == F){
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
                       edge.width          = 1 + ( prop * results$score ),
                       edge.arrow.size     = 0.5,
                       vertex.size         = 15,
                       vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
                       main                = ttl
  )
 } else {
   igraph::tkplot( netw,
                        vertex.frame.color  = "white",
                        layout              = lay,
                        vertex.color        = ifelse( igraph::V( netw )$name %in% genes, "#ff349a" , ifelse(igraph::V( netw )$name %in% diseases, "royalblue2", "goldenrod3")),
                        vertex.label.dist   = 0,      #puts the name labels slightly off the dots
                        vertex.frame.color  = 'blue', #the color of the border of the dots
                        vertex.label.color  = 'black',#the color of the name labels
                        vertex.label.font   = 0,      #the font of the name labels
                        vertex.label        = igraph::V( netw )$names, #specifies the lables of the vertices. in this case the 'name' attribute is used
                        edge.color          = "darkgrey",
                        edge.width          = 1 + ( prop * results$score ),
                        edge.arrow.size     = 0.5,
                        vertex.size         = 15,
                        vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
                        main                = ttl
   )
 }
}



plot_dis_comorbidity <- function( input, type, layout, search, verbose, prop, limit, interactive ) {

  # if( type != "disease-disease-gene" | type !="disease-disease-variant" ) {
  #   stop("For this type of chart, objects created with disease2disease_by_gene function is required." )
  # }
  if( type == "disease-disease-gene"){
    input$common <- input$ngenes
    ds1 <- input[,c( "disease1_name", "ngenes1" ) ]
    colnames( ds1 ) <- c("dis", "gen")
    ds2 <- input[, c("disease2_name", "ngenes2") ]
    colnames( ds2 ) <- c("dis", "gen")
  }
  else if(type == "disease-disease-variant"){
    input$common <- input$nvariants
    ds1 <- input[,c( "disease1_name", "nvariants1" ) ]
    colnames( ds1 ) <- c("dis", "gen")
    ds2 <- input[, c("disease2_name", "nvariants2") ]
    colnames( ds2 ) <- c("dis", "gen")
  }
  else {
    stop("For this type of chart, objects created with **disease2disease_by_gene** or **disease2disease_by_variant** functions is required.")
  }


  a <- nrow(input)
  if (dim(input)[1]> limit){
    input <- input[c(1:limit),]
    show( paste0("warning: dataframe of ", a , " rows has been reduced to ", limit, " rows." ))
  }
  edges <- data.frame( input[ c("disease1_name", "disease2_name")] )
  netw  <- igraph::graph.data.frame( edges, directed = FALSE )
  netw  <- igraph::simplify( netw )
  #lay   <- layout(netw)
  lay <- apply_layout(layout = layout, netw = netw)

  if( verbose ) {
    message( "The network contains ", igraph::vcount( netw ), " nodes and ", igraph::ecount( netw ), " edges." )
  }

  dn1 <- unique(input$disease1_name)
  if (length(dn1) == 1){
    ttl <- paste0( "Disease-Disease Network" )
  }   else{
    ttl <- paste0( "Disease-Disease Network" )
  }



  fnl <- rbind ( ds1, ds2)
  fnl <- fnl[!duplicated(fnl), ]

  sizes <- as.numeric( fnl[ , 2 ] )
  names( sizes ) <- fnl[ , 1 ]

  if (interactive == F){
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
  ) } else {
    igraph::tkplot( netw,
                         vertex.frame.color  = "white",
                         layout              = lay,
                         vertex.color        = ifelse( igraph::V( netw )$name %in% ds1$dis, "royalblue2", "#b4daff" ),
                         vertex.label.dist   = 0,      #puts the name labels slightly off the dots
                         vertex.frame.color  = 'blue', #the color of the border of the dots
                         vertex.label.color  = 'black',#the color of the name labels
                         vertex.label.font   = 0,      #the font of the name labels
                         vertex.label        = igraph::V( netw )$names, #specifies the lables of the vertices. in this case the 'name' attribute is used
                         edge.color          = "black",
                         edge.width          = 3,
                         edge.arrow.size     = 0.5,
                         vertex.size         = prop *30*log2(as.numeric(sizes[ igraph::V( netw )$name ])+1),
                         vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
                         main                = ttl
    )
  }
}

#
# plot_dis_dis_associations <- function( input, type, layout, search, verbose, prop, limit, database ) {
#
#   # if( search != "dis-dis" ) {
#   #   stop("For this type of chart, objects comming from disgenetDisDis function is required." )
#   # } else  if( type != "disease" ) {
#   #   stop("For this type of chart, objects comming from disgenetDisDis function is required." )
#   # }
#   l <- unique(input$disease1_name)
#   if (length(l) == 1){
#     stop("For this type of chart, objects generated with disease2disease_by_gene, or disease2disease_by_variant
#          function, and search = list, are required ." )
#   }
#
#   if( type == "disease-disease-gene"){
#     input$common <- input$ngenes
#     disNum <- length ( unique ( input$diseaseid1) )
#     diseases <- as.character(unique (  input$diseaseid1) )
#
#     res <- disease2gene(disease = diseases,
#                         database = as.character(database))
#     res <- res@qresult
#     pairs <- ( disNum *  disNum - disNum ) / 2
#     tb <- as.data.frame(setNames(replicate(7,numeric(0), simplify = F), c("disease1","NameDisease1", "disease2","NameDisease2", "NGenes1", "NGenes2", "intr")))
#     cnt <- 1
#     for ( i in 1:disNum){
#       dis1 <- diseases[i]
#       nd1 <- as.character (  res[res[,"diseaseid"]==diseases[i],"disease_name"][1]  )
#       gen1 <- unique( subset(res, diseaseid == dis1 )$geneid)
#       for ( j in 2:disNum){
#         if ( i < j){
#           dis2 <- diseases[j]
#           nd2 <- as.character (  res[res[,"diseaseid"]==diseases[j],"disease_name"][1]  )
#           gen2 <- unique( subset(res, diseaseid == dis2 )$geneid)
#           tb[ cnt, ] <- c(diseases[ i ], nd1, diseases[ j ],nd2, length( gen1 ), length( gen2 ), length( intersect( gen1, gen2 ) ))
#           cnt <- cnt + 1
#         }
#       }
#     }
#   }
#   else if(type == "disease-disease-variant"){
#     input$common <- input$nvariants
#     disNum <- length ( unique ( input$diseaseid1) )
#     diseases <- as.character(unique (  input$diseaseid1) )
#
#     res <- disease2variant(disease = diseases,
#                            database = as.character(database))
#     res <- res@qresult
#     pairs <- ( disNum *  disNum - disNum ) / 2
#     tb <- as.data.frame(setNames(replicate(7,numeric(0), simplify = F), c("disease1","NameDisease1", "disease2","NameDisease2", "NGenes1", "NGenes2", "intr")))
#     cnt <- 1
#     for ( i in 1:disNum){
#       dis1 <- diseases[i]
#       nd1 <- as.character (  res[res[,"diseaseid"]==diseases[i],"disease_name"][1]  )
#       gen1 <- unique( subset(res, diseaseid == dis1 )$variantid)
#       for ( j in 2:disNum){
#         if ( i < j){
#           dis2 <- diseases[j]
#           nd2 <- as.character (  res[res[,"diseaseid"]==diseases[j],"disease_name"][1]  )
#           gen2 <- unique( subset(res, diseaseid == dis2 )$variantid)
#           tb[ cnt, ] <- c(diseases[ i ], nd1, diseases[ j ],nd2, length( gen1 ), length( gen2 ), length( intersect( gen1, gen2 ) ))
#           cnt <- cnt + 1
#         }
#       }
#     }
#   }
#
#
#   tb <- tb[tb[,"intr"]>0 ,]
#   finaldiseases <- unique( c(as.character(tb$disease1), as.character(tb$disease2)))
#   #   if (dim(input)[1]> limit){
#   #     input <- input[c(1:limit),]
#   #     show( paste0("warning: dataframe of ", a , " rows has been reduced to ", limit, " rows." ))
#   #   }
#   edges <- data.frame( tb[ , "NameDisease1" ], tb[ , "NameDisease2" ] )
#   netw  <- graph.data.frame( edges, directed = FALSE )
#   netw  <- simplify( netw )
#   lay   <- layout(netw)
#
#   if( verbose ) {
#     message( "The network contains ", igraph::vcount( netw ), " nodes and ", igraph::ecount( netw ), " edges." )
#   }
#
#   ttl <- paste0( "Disease-Disease Network" )
#
#   if( type == "disease-disease-gene"){
#     ds1 <- input[,c( "disease1_name", "disease1_ngenes" ) ]
#     colnames( ds1 ) <- c("dis", "gen")
#     ds2 <- input[, c("disease2_name", "disease2_ngenes") ]
#     colnames( ds2 ) <- c("dis", "gen")
#   }
#   else if ( type == "disease-disease-variant"){
#     ds1 <- input[,c( "disease1_name", "disease1_nvariants" ) ]
#     colnames( ds1 ) <- c("dis", "gen")
#     ds2 <- input[, c("disease2_name", "disease2_nvariants") ]
#     colnames( ds2 ) <- c("dis", "gen")
#   }
#
#   fnl <- rbind ( ds1, ds2)
#   fnl <- fnl[!duplicated(fnl), ]
#
#   sizes <- as.numeric( fnl[ , 2 ] )
#   names( sizes ) <- fnl[ , 1 ]
#
#   igraph::plot.igraph( netw,
#                        vertex.frame.color  = "white",
#                        layout              = lay,
#                        vertex.color        =  "royalblue2" ,
#                        vertex.label.dist   = 0,      #puts the name labels slightly off the dots
#                        vertex.frame.color  = 'blue', #the color of the border of the dots
#                        vertex.label.color  = 'black',#the color of the name labels
#                        vertex.label.font   = 0,      #the font of the name labels
#                        vertex.label        = igraph::V( netw )$names, #specifies the lables of the vertices. in this case the 'name' attribute is used
#                        edge.color          = "black",
#                        edge.width          = prop * (as.numeric( input$ngenes ))^1/2,
#                        edge.arrow.size     = 0.5,
#                        # vertex.size         = 10,
#                        vertex.size         = prop *30*log2(as.numeric(sizes[ igraph::V( netw )$name ])+1),
#                        vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
#                        main                = ttl
#   )
# }

##################################################################
## HEATMAP
##################################################################

plot_variant_disease_heatmap <- function( type, search, input,nchars, cutoff, limit, verbose ) {

  if( search == "single" ) {
    stop( "For this type of chart, a multiple query created with 'variant2disease'is required. " )
  }

  if( type == "variant-disease" ) {
    results <- unique( input[ , c( "diseaseid","disease_name", "variantid","score" ) ] )

    results <- results [ results$score>= cutoff, ]
    results <- results [ with ( results , order(-score)), ]
    if ( dim( results )[ 1 ] > limit ){
      results <- results[ 1:limit ,]
      show( paste0("Dataframe of ", nrow(input), " rows has been reduced to ", limit, " rows." ))
    }
    results$disease_name <- as.character(results$disease_name)
    for( i in 1:nrow(results)){

      if(nchar(as.character(results$disease_name[i])) > nchars){
        rename <- substr(results$disease_name[i], 1, nchars)
        results$disease_name[i] <- paste0(rename, "...")
      }
    }
    casted_tt <- reshape::cast( results, variantid ~disease_name, value= "score", length )
    l <- dim( casted_tt )[ 2 ]
    casted_tt[ is.na( casted_tt ) ] <- 0
    res <- hclust( dist( casted_tt[ 2:l ] ),  method = "ward.D")
    casted_tt$clus <- cutree( res, k = 5 )
    ordered <- casted_tt[order(casted_tt$clus), "variantid"]

    #results$diseaseName  <- factor(results$diseaseName , levels= as.factor(ordered))
    results$variantid  <- factor(results$variantid , levels= as.factor(ordered))

    p<- ggplot2::qplot( x = variantid, y = disease_name, data= results , fill = score, geom = "tile" ) +
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

  } else  if( type == "disease-variant" ) {

    results <- unique( input[ , c( "diseaseid","disease_name", "variantid","score" ) ] )

    results <- results [ results$score>= cutoff, ]
    if ( dim( results )[ 1 ] > limit ){
      results <- results[ 1:limit ,]
      show( paste0("warning: dataframe of ", nrow(input), " rows has been reduced to ", limit, " rows." ))
    }

    #       for( i in 1:nrow(results)){
    #
    #         if(nchar(as.character(results$c1.name[i])) > 30){
    #           rename <- substr(results$c1.name[i], 1, 30)
    #           results$c1.name[i] <- paste0(rename, "...")
    #         }
    #       }

    casted_tt <- reshape::cast( results, disease_name~variantid, value= "score", length )

    l <- dim(casted_tt)[2]
    casted_tt[is.na( casted_tt)] <- 0
    res <- hclust(dist(casted_tt[2:l]),  method = "ward.D")
    mm<-max(res$order)
    if (mm>4){
      casted_tt$clus <- cutree(res, k = 5)
    } else {
      casted_tt$clus <- cutree(res, k = mm)
    }

    ordered <- casted_tt[order(casted_tt$clus), "disease_name"]

    results$disease_name  <- factor(results$disease_name , levels= as.factor(ordered))
    p <- ggplot2::ggplot( results,  ggplot2::aes ( x = disease_name, y = variantid ) )
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
##################################################################
## BARPLOTS
##################################################################
plot_barplot <- function ( input, type, search, verbose, limit ) {

  #input <- input[c(1:20),]
  input$disease1_name <- as.character(input$disease1_name)
  input$disease2_name <- as.character(input$disease2_name)
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

  if (type =="disease-disease-gene"){
    input$common <- as.numeric(input$ngenes)
    input$JaccardIndex <- as.numeric(input$jaccard_genes)
    input <- subset(input, common>0)
    input <- input [ with ( input , order(-common)), ]
  } else if (type =="disease-disease-variant"){
    input$common <- as.numeric(input$nvariants)
    input$JaccardIndex <- as.numeric(input$jaccard_variants)
    input <- subset(input, common>0)
    input <- input [ with ( input , order(-common)), ]
  }


  if (dim(input)[1]> limit){
    input <- input[c(1:limit),]
  }

  input$pairs <- paste(input$disease1_name, input$disease2_name, sep = " - ")
  orderedDiseases <- input[order(-as.numeric(input$common)), "pairs"  ]
  input$pairs <- factor(input$pairs , levels= as.factor(orderedDiseases))

  p <- ggplot2::ggplot (input, ggplot2::aes ( x = pairs, y = as.numeric(common) ), order = as.numeric(commonGenes) ) +
    ggplot2::geom_bar ( stat = "identity", fill = "grey" ) +
    ggplot2::labs ( title = paste0 ( "Disease-Disease Barplot" ) , x = "Disease Pairs", y = "# of shared genes")


  p <- p + ggplot2::theme_classic( ) + ggplot2::theme( plot.margin = grid::unit( x = c ( 5, 15, 5, 15 ), units = "mm" ),
                                                       axis.line = ggplot2::element_line ( size = 0.7, color = "black" ), text = ggplot2::element_text ( size = 14 ) ,
                                                       axis.text.x = ggplot2::element_text ( angle = 45, size = 9, hjust = 1 ))


  p <- p + ggplot2::geom_point(ggplot2::aes(x = pairs , y =as.numeric(JaccardIndex) *100))
  p <- p + ggplot2::geom_text(ggplot2::aes(x = pairs , y =as.numeric(JaccardIndex) *100 + 1), label = c(round(as.numeric(input$JaccardIndex), digits = 3 )), size=3)
  p

}
##################################################################
## VENN
##################################################################
plot_venn <- function( search, type, input, verbose, database ) {

  # if( type != "" ) {
  #   stop("For this type of chart, objects comming from disgenetDisDis function is required." )
  # }

  results <- input
  diseases <- c(unique( as.character(unique(results[, "diseaseid1"]))) )

  if( length(diseases) > 5 | length(diseases) ==1 ) {
    stop( "For this type of chart, the diseases list should have 2-5 diseases" )
  }

  else if (length(diseases) == 2){
    gen1 <- disease2gene(disease = diseases[1] )
    name1 <-as.character(unique(gen1@qresult$disease_name))
    gen1 <- unique(gen1@qresult$gene_symbol)
    gen2 <- disease2gene(disease = diseases[2] )
    name2 <-as.character(unique(gen2@qresult$disease_name))
    gen2 <- unique(gen2@qresult$gene_symbol)

    grid::grid.newpage()
    venn.plot <- VennDiagram::draw.pairwise.venn(
      area1        = length(gen1),
      area2        = length(gen2),
      cross.area   = length(intersect(gen1, gen2)),
      scaled       = T,
      category     = c(name1, name2),
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
  }
  else if (length(diseases) == 3) {

    gen1 <- disease2gene(disease = diseases[1] )
    name1 <-as.character(unique(gen1@qresult$disease_name))
    gen1 <- unique(gen1@qresult$gene_symbol)
    gen2 <- disease2gene(disease = diseases[2] )
    name2 <-as.character(unique(gen2@qresult$disease_name))
    gen2 <- unique(gen2@qresult$gene_symbol)
    gen3 <- disease2gene(disease = diseases[3] )
    name3 <-unique(gen3@qresult$disease_name)
    gen3 <- unique(gen3@qresult$gene_symbol)

    grid::grid.newpage()
    venn.plot <- VennDiagram::draw.triple.venn(
      area1           = length(gen1),
      area2           = length(gen2),
      area3           = length(gen3),
      n12             = length(intersect(gen1, gen2)),
      n23             = length(intersect(gen2, gen3)),
      n13             = length(intersect(gen1, gen3)),
      n123            = length(intersect(intersect(gen1, gen2), gen3)),
      category        = c(name1, name2, name3),
      fill            = c('red', 'blue', 'green'),
      alpha = c(0.4,0.4,0.4),
      cat.col         = c('red', 'blue', 'green'),
      cat.pos = 1,
      #cex             = c(1/2,2/2,3/2,4/2,5/2,6/2,7/2),
      cex             = 2,
      cat.cex         = c(1.5,1.5,1.5),
      euler           = TRUE,
      euler.d  = TRUE,
      scaled          = TRUE
    )
  }
  else if (length(diseases) == 4) {

    gen1 <- disease2gene(disease = diseases[1] )
    name1 <-as.character(unique(gen1@qresult$disease_name))
    gen1 <- unique(gen1@qresult$gene_symbol)
    gen2 <- disease2gene(disease = diseases[2] )
    name2 <-as.character(unique(gen2@qresult$disease_name))
    gen2 <- unique(gen2@qresult$gene_symbol)
    gen3 <- disease2gene(disease = diseases[3] )
    name3 <-as.character(unique(gen3@qresult$disease_name))
    gen3 <- unique(gen3@qresult$gene_symbol)
    gen4 <- disease2gene(disease = diseases[4] )
    name4 <-as.character(unique(gen4@qresult$disease_name))
    gen4 <- unique(gen4@qresult$gene_symbol)
    grid::grid.newpage()
    venn.plot <- VennDiagram::draw.quad.venn(
      area1  = length(gen1),
      area2  = length(gen2),
      area3  = length(gen3),
      area4  = length(gen4),
      n12    = length(intersect(gen1, gen2)),
      n23    = length(intersect(gen2, gen3)),
      n13    = length(intersect(gen1, gen3)),
      n14 = length(intersect(gen1, gen4)),
      n24= length(intersect(gen2, gen4)),
      n34= length(intersect(gen3, gen4)),
      n123 = length(intersect(intersect(gen1, gen2), gen3)),
      n124 = length(intersect(intersect(gen1, gen2), gen4)),
      n134 = length(intersect(intersect(gen1, gen3), gen4)),
      n234 = length(intersect(intersect(gen2, gen3), gen4)),
      n1234 = length(intersect(intersect(intersect(gen1, gen2), gen3),gen4)),

      category        = c(name1, name2, name3, name4),

      fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
      cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
      cat.cex = 2,
      cat.pos = c(-20, 20, 60, 120),
      margin = 0.3,
      cex = 1.5,
      euler           = TRUE,
      scaled          = FALSE
    )
  }
  else if (length(diseases) == 5){

    gen1 <- disease2gene(disease = diseases[1] )
    name1 <-as.character(unique(gen1@qresult$disease_name))
    gen1 <- unique(gen1@qresult$gene_symbol)
    gen2 <- disease2gene(disease = diseases[2] )
    name2 <-as.character(unique(gen2@qresult$disease_name))
    gen2 <- unique(gen2@qresult$gene_symbol)
    gen3 <- disease2gene(disease = diseases[3] )
    name3 <-as.character(unique(gen3@qresult$disease_name))
    gen3 <- unique(gen3@qresult$gene_symbol)
    gen4 <- disease2gene(disease = diseases[4] )
    name4 <-as.character(unique(gen4@qresult$disease_name))
    gen4<- unique(gen4@qresult$gene_symbol)
    gen5 <- disease2gene(disease = diseases[5] )
    name5 <-as.character(unique(gen5@qresult$disease_name))
    gen5<- unique(gen5@qresult$gene_symbol)

    grid::grid.newpage()
    venn.plot <- VennDiagram::draw.quintuple.venn(
      area1  = length(gen1),
      area2  = length(gen2),
      area3  = length(gen3),
      area4  = length(gen4),
      area5 = length(gen5),
      n12 = length(intersect(gen1, gen2)),
      n13 = length(intersect(gen1, gen3)),
      n14 = length(intersect(gen1, gen4)),
      n15 = length(intersect(gen1, gen5)),
      n23 = length(intersect(gen2, gen3)),
      n24 = length(intersect(gen2, gen4)),
      n25 = length(intersect(gen2, gen5)),
      n34 = length(intersect(gen3, gen4)),
      n35 = length(intersect(gen3, gen5)),
      n45 = length(intersect(gen4, gen5)),
      n123 = length(intersect(intersect(gen1, gen2), gen3)),
      n124 = length(intersect(intersect(gen1, gen2), gen4)),
      n125 = length(intersect(intersect(gen1, gen2), gen5)),
      n134 = length(intersect(intersect(gen1, gen3), gen4)),
      n135 = length(intersect(intersect(gen1, gen3), gen5)),
      n145 = length(intersect(intersect(gen1, gen4), gen5)),
      n234 = length(intersect(intersect(gen2, gen3), gen4)),
      n235 = length(intersect(intersect(gen2, gen3), gen5)),
      n245 = length(intersect(intersect(gen2, gen4), gen5)),
      n345 = length(intersect(intersect(gen3, gen4), gen5)),
      n1234 = length(intersect(intersect(intersect(gen1, gen2), gen3),gen4)),
      n1235 = length(intersect(intersect(intersect(gen1, gen2), gen3),gen5)),
      n1245 = length(intersect(intersect(intersect(gen1, gen2), gen4),gen5)),
      n1345 = length(intersect(intersect(intersect(gen1, gen3), gen4),gen5)),
      n2345 = length(intersect(intersect(intersect(gen2, gen3), gen4),gen5)),
      n12345 = length(intersect(intersect(intersect(intersect(gen1, gen2), gen3),gen4), gen5)),
      category        = c(name1, name2, name3, name4, name5),
      fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
      cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
      cat.cex = 2,
      margin = 0.05,
      cat.pos =  1,
     # cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
      #        1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
     cex = 1.5,

      ind = TRUE
    )
  }
  else {
    stop("venn for more than 5 groups not possible! ")
  }
  grid::grid.draw(venn.plot)

}

##################################################################
##################################################################
# plot_variant_gene<- function( input, type, layout, prop, search, verbose ) {
#
#   results <- unique(input[c("variantid",  "gene_symbol")])
#   l <- length(unique(results$variantid))
#   #input <-  aggregate(score ~ geneId + geneSymbol  + snp, input , max)
#   results <- separate_rows(results, gene_symbol, sep = ",")
#   nd <- length(unique(subset(results, gene_symbol!="")$variantid))
#   results <- subset(results, gene_symbol!="")
#
#   if( l > 1 ) {
#     ttl <- "Variant-Gene Network"
#   }else if ( l == 1 ){
#     ttl <- paste0( "Genes associated to ", results$variantid[ 1 ] )
#   }
#   if (nd > 0){
#     show( paste0("warning: ", nd , " variant(s) not shown in the plot" ))
#   }
#
#   edges <- data.frame( results[c("variantid", "gene_symbol")] )
#   netw  <- igraph::graph.data.frame( edges, directed = FALSE )
#   netw  <- igraph::simplify( netw )
#   lay   <- layout( netw )
#   genes <- unique( results$gene_symbol )
#   igraph::plot.igraph( netw,
#                        vertex.frame.color  = "white",
#                        layout              = lay,
#                        vertex.color        = ifelse( igraph::V( netw )$name %in% genes, "#ff349a" , "goldenrod3"),
#                        vertex.label.dist   = 0,      #puts the name labels slightly off the dots
#                        vertex.frame.color  = 'blue', #the color of the border of the dots
#                        vertex.label.color  = 'black',#the color of the name labels
#                        vertex.label.font   = 0,      #the font of the name labels
#                        vertex.label        = igraph::V( netw )$names, #specifies the lables of the vertices. in this case the 'name' attribute is used
#                        edge.color          = "darkgrey",
#                        edge.width          = 1 + ( prop  ),
#                        edge.arrow.size     = 0.5,
#                        vertex.size         = 15,
#                        vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
#                        main                = ttl
#   )
# }


plot_disease_class <- function ( input, type, layout, search, nchars, verbose, prop ) {

  if( search == "list" ) {
    stop( "Chart not possible for list. Use instead the plot function with DiseaseClassHeatmap class argument" )
  }
  else if( type == "disease-disease-gene" | type == "disease-disease-variant" ) {

    if (length(unique(input[, "diseaseid1"])) == 1){
      disats <- disease2attribute(disease = unique(c(input$diseaseid1 , input$diseaseid2)))
      disats <- disats@qresult[c("diseaseid", "diseaseclassname")]
      disats$disease_class_name <- disats$diseaseclassname
      disats <- tidyr::separate_rows(disats, disease_class_name, sep = "; ")
      disats$disease_class_name <- gsub ( "Wounds and Injuries", NA, disats$disease_class_name )
      disats$disease_class_name  <- gsub ( "Animal Diseases", NA, disats$disease_class_name)
      disats$disease_class_name <-ifelse(disats$disease_class_name == "", NA, as.character(disats$disease_class_name))
      nd <- length(unique(subset(disats, is.na(disease_class_name))$diseaseid))
      input<- merge(input, disats,by.x = "diseaseid1", by.y = "diseaseid")
      input <- merge(input, disats,by.x = "diseaseid2", by.y = "diseaseid")

      if(dim(input[is.na(input$disease_class_name.x), ])[1]> 0){
        stop("index diseases has no disease class")
      }
      results <- subset(input, !is.na(disease_class_name.y) )
      # this means there is a central disease
      results <- unique(results[c("diseaseid2", "disease_class_name.y")])

      for( i in 1:nrow(results)){
        if(nchar(as.character(results$disease_class_name.y[i])) > nchars){
          rename <- substr(results$disease_class_name.y[i], 1, nchars)
          results$disease_class_name.y[i] <- paste0(rename, "...")
        }
      }
      if (nd > 0){
        show( paste0("warning: ", nd , " disease(s) not shown in the plot" ))
      }

      #results <- merge(results[, c(2,3,5,6)], results, by.x = "Disease2")
      freq <- as.data.frame ( table ( results$disease_class_name.y ) )
      freq$perc <- ( freq$Freq / sum( freq$Freq ) ) * 100
      rw <- matrix ( c( as.character( input[ 1,"disease1_name" ] ),10, 10), ncol=3 )
      colnames( rw ) <- colnames( freq )
      freq <- rbind( freq, rw )

      sizes <- as.numeric( freq[ , 3 ] )
      names( sizes ) <- freq[ , 1 ]

      results <- merge( input, results, by="diseaseid2")

      #plot
      edges <- data.frame( results[ , "disease1_name" ], results[ , "disease_class_name.y.y" ] )
      netw  <- igraph::graph.data.frame( edges, directed = FALSE )
      netw  <- igraph::simplify( netw )
      #lay   <- layout( netw )
      lay <- apply_layout(layout, netw)
      if( verbose ) {
        message( "The network contains ", igraph::vcount( netw ), " nodes and ", igraph::ecount( netw ), " edges." )
      }

      diseases <- unique( results$disease_class_name.y.y )
      ttl <- paste0("Disease Classes associated to ", input[1,"disease1_name"] )

      igraph::plot.igraph( netw,
                           vertex.frame.color  = "white",
                           layout              = lay,
                           vertex.color        = ifelse( igraph::V( netw )$name %in% diseases, "#9bcdff", "royalblue2"),
                           vertex.label.dist   = 0,      #puts the name labels slightly off the dots
                           vertex.frame.color  = 'blue', #the color of the border of the dots
                           vertex.label.color  = 'black',#the color of the name labels
                           vertex.label.font   = 0,      #the font of the name labels
                           vertex.label        = igraph::V( netw )$name, #specifies the labels of the vertices. in this case the 'name' attribute is used
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

  else if( search == "single" & type == "gene-disease" ) {
    results <- as.data.frame(input [, c ("gene_symbol", "diseaseid", "disease_name", "disease_class_name")])
    # results$gene_symbol <- as.character(results$gene_symbol)
    results$disease_class_name <- gsub ( "Wounds and Injuries", NA, results$disease_class_name )
    results$disease_class_name  <- gsub ( "Animal Diseases", NA, results$disease_class_name )
    results$disease_class_name <- gsub(";\\s+", ";", results$disease_class_name)
    results$disease_class_name <- gsub("^\\s+", "", results$disease_class_name)
    results <- as.data.frame (  tidyr::separate_rows(results, disease_class_name, sep = ";"))
    results$disease_class_name <-ifelse(results$disease_class_name == "", NA, as.character(results$disease_class_name))

    nd <- dim(subset(results, is.na(disease_class_name)))[1]
    if (nd > 0){
      show( paste0("warning: ", nd , " disease(s) not shown in the plot" ))
    }
    results <- subset(results, !is.na(disease_class_name))
    freq <- as.data.frame ( table ( results$disease_class_name ) )
    freq$perc <- ( freq$Freq / sum( freq$Freq ) ) * 100

    rw <- matrix ( c( as.character( results[ 1, 1 ] ), 10, 10), ncol=3 )
    colnames( rw ) <- colnames( freq )
    freq <- rbind( freq, rw )

    sizes <- as.numeric( freq[ , 3 ] )
    names( sizes ) <- freq[ , 1 ]

    #plot
    edges <- data.frame( results[ , c("gene_symbol" ,"disease_class_name") ] )
    netw  <- igraph::graph.data.frame( edges, directed = FALSE )
    netw  <- igraph::simplify( netw )
    #lay   <- layout( netw )
    lay <- apply_layout(layout, netw)

    if( verbose ) {
      message( "The network contains ", igraph::vcount( netw ), " nodes and ", igraph::ecount( netw ), " edges." )
    }

    diseases <- unique( results$disease_class_name )
    ttl <- paste0("Disease Classes associated to ", results$gene_symbol[1] )

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

  else if (search == "single" & type == "variant-disease"){
    results <- unique(input[c("score", "variantid", "disease_name", "disease_class_name")])
    l <- length(unique(results$variantid))

    if(l > 1 ) {
      stop( "Chart not possible for list. Use instead the plot function with DiseaseClassHeatmap class argument" )
    }

    results$disease_class_name <- gsub ( "Wounds and Injuries", NA, results$disease_class_name )
    results$disease_class_name  <- gsub ( "Animal Diseases", NA, results$disease_class_name )
    results <- data.frame( tidyr::separate_rows(results, disease_class_name, sep = "; "))
    results$disease_class_name <-ifelse(results$disease_class_name == "", NA, as.character(results$disease_class_name))

    nd <- dim(subset(results, is.na(disease_class_name)))[1]
    if (nd > 0){
      show( paste0("warning: ", nd , " disease(s) not shown in the plot" ))
    }
    results <- subset(results, !is.na(disease_class_name))

    for( i in 1:nrow(results)){
      if(nchar(as.character(results$disease_class_name[i])) > nchars){
        rename <- substr(results$disease_class_name[i], 1, nchars)
        results$disease_class_name[i] <- paste0(rename, "...")
      }
    }

    freq <- as.data.frame ( table ( results$disease_class_name ) )
    freq$perc <- ( freq$Freq / sum( freq$Freq ) ) * 100

    rw <- matrix ( c( as.character( results[ 1, "variantid" ] ), 10, 10), ncol=3 )
    colnames( rw ) <- colnames( freq )
    freq <- rbind( freq, rw )

    sizes <- as.numeric( freq[ , 3 ] )
    names( sizes ) <- freq[ , 1 ]

    #plot
    edges <- data.frame( results[, "variantid" ], results[ , "disease_class_name" ] )
    netw  <- igraph::graph.data.frame( edges, directed = FALSE )
    netw  <- igraph::simplify( netw )
    #lay   <- layout( netw )
    lay <- apply_layout(layout, netw)

    if( verbose ) {
      message( "The network contains ", igraph::vcount( netw ), " nodes and ", igraph::ecount( netw ), " edges." )
    }

    diseases <- unique( results$disease_class_name )
    ttl <- paste0("Disease Classes associated to ", results$variantid[1] )
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
  else {
    stop("disease class plot not for this type of object")
  }

}

plot_gda_heatmap <- function( type, search, input, nchars, cutoff=0, limit, verbose ) {

  if( search == "single"  ) {
    stop( "For this type of chart, a multiple query created with 'disease2gene', or 'gene2disease' is required. For single disease is not possible to plot this graphic." )
  }
  if( search == "list" ) {
    if( type == "disease-gene" ) {
      results <- ( input[ , c( "disease_name", "gene_symbol", "score" ) ] )
      results$disease_name <- as.character(results$disease_name)
      results <- results [ results$score>= cutoff, ]
      results <- results [ with ( results , order(-score)), ]
      if ( dim( results )[ 1 ] > limit ){
        results <- results[ 1:limit ,]
        show( paste0("Dataframe of ", nrow(input), " rows has been reduced to ", limit, " rows." ))

      }

      for( i in 1:nrow(results)){
        if(nchar(as.character(results$disease_name[i])) > nchars){
          rename <- substr(results$disease_name[i], 1, nchars)
          results$disease_name[i] <- paste0(rename, "...")
        }
      }

      casted_tt <- reshape::cast( results, gene_symbol ~ disease_name, value= "score", length )
      l <- dim( casted_tt )[ 2 ]
      casted_tt[ is.na( casted_tt ) ] <- 0
      res <- hclust( dist( casted_tt[ 2:l ] ),  method = "ward.D")
      casted_tt$clus <- cutree( res, k = l+1 )
      ordered <- casted_tt[order(casted_tt$clus), "gene_symbol"]

      results$gene_symbol  <- factor(results$gene_symbol , levels= as.factor(ordered))
      p<- ggplot2::qplot( x = disease_name, y = gene_symbol, data= results , fill = score, geom = "tile" ) +
        ggplot2::ylab( "Genes" ) + ggplot2::xlab( "Disease Name" ) +
        ggplot2::ggtitle( "Gene-Disease Heatmap") +
        # ggplot2::scale_fill_gradient2(low = "#0000FF", high ="#2f30c9",guide = "colourbar", "Score") +
        ggplot2::scale_fill_gradient(low = "#0000FF", high ="#2f30c9",  guide = "colourbar", "Score", breaks = c(0,1))  +
        ggplot2::theme ( axis.text = ggplot2::element_text ( size = 9),
                         axis.title=ggplot2::element_text(size=10 ),
                         panel.background = ggplot2::element_rect(fill="white", colour = "black"),
                         text = ggplot2::element_text(size = 11),
                         axis.text.x = ggplot2::element_text(angle=45, size = 11, hjust = 1)
        )
      print(p)

    }

    else  if( type == "gene-disease" ) {

      results <- ( input[ , c( "score", "gene_symbol", "disease_name" ) ] )
      results$gene_symbol <- as.character( results$gene_symbol )
      results <- results [ results$score>= cutoff, ]


      if ( dim( results )[ 1 ] > limit ){
        results <- results[order(-results$score), ]
        results <- results[ 1:limit ,]
        show( paste0("warning: dataframe of ", nrow(input), " rows has been reduced to ", limit, " rows." ))

      }

      for( i in 1:nrow(results)){

        if(nchar(as.character(results$disease_name[i])) > nchars){
          results$disease_name <- as.character(results$disease_name )
          rename <- substr(results$disease_name[i], 1, nchars)
          results$disease_name[i] <- paste0(rename, "...")
        }
      }

      casted_tt <- reshape::cast(results, disease_name ~ gene_symbol, value= "score" , length )
      l <- dim(casted_tt)[2]
      casted_tt[is.na( casted_tt)] <- 0
      res <- hclust(dist(casted_tt[2:l]),  method = "ward.D")

      casted_tt$clus <- cutree(res, k = 5)
      ordered <- casted_tt[order(casted_tt$clus), "disease_name"]

      results$disease_name  <- factor(results$disease_name , levels= as.factor(ordered))
      p<- ggplot2::qplot( x = gene_symbol, y = disease_name, data= results , fill = score, geom = "tile" ) +
        ggplot2::ylab( "Diseases" ) + ggplot2::xlab( "Genes" ) + ggplot2::ggtitle( "Gene-Disease Heatmap") +
        #ggplot2::scale_fill_gradient2(low = "#0000FF", high ="#2f30c9", breaks=c(0,1), guide = "colourbar", "Score") +
         ggplot2::scale_fill_gradient(low = "#0000FF", high ="#2f30c9",  guide = "colourbar", "Score", breaks = c(0,1))  +
        ggplot2::theme ( axis.text = ggplot2::element_text ( size = 9), axis.title=ggplot2::element_text(size=10 ), panel.background = ggplot2::element_rect(fill="white", colour = "black"),
                         text = ggplot2::element_text(size = 11),
                         axis.text.x = ggplot2::element_text(angle=45, size = 11, hjust = 1))
      #plot(p)
      print(p)
    }
  }
}

plot_gda_heatmapPanther <- function( type, search, input,nchars, verbose ) {
  results <- input [, c ("gene_symbol", "diseaseid", "disease_name", "protein_class_name")]
  results <- tidyr::separate_rows(results, protein_class_name, sep = "; ")
  results$protein_class_name <-ifelse(results$protein_class_name == "", NA, as.character(results$protein_class_name))
  nd <- length(unique(as.character(subset(results, is.na(protein_class_name))$gene_symbol)))
  if (nd > 0){
    show( paste0("warning: ", nd , " gene(s) not shown in the plot" ))
  }
  results <- subset(results, !is.na(protein_class_name))

  # for( i in 1:nrow(results)){
  #   if(nchar(as.character(results$disease_class_name)) > 25){
  #     rename <- substr(results$disease_class_name[i], 1, 25)
  #     results$disease_class_name[i] <- paste0(rename, "...")
  #   }
  # }
  for( i in 1:nrow(results)){
    if(nchar(as.character(results$disease_name[i])) > nchars){
      rename <- substr(results$disease_name[i], 1,nchars)
      results$disease_name[i] <- paste0(rename, "...")
    }
  }

  for ( i in 1:length(unique(results$disease_name))){
    disease1 <- results[results$disease_name==unique(results$disease_name)[i],]
    freq <- as.data.frame ( table ( as.character(disease1$protein_class_name ) ))
    freq$perc <- (freq$Freq/sum(freq$Freq))*100
    freq$dis <- unique(results$disease_name)[i]

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
  p <- p + ggplot2::theme( plot.margin=grid::unit(x=c(5,15,5,15),units="mm") , axis.line = ggplot2::element_line(size = 0.7, color = "black"), text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(angle=45, size = 11, hjust = 1), panel.background = ggplot2::element_blank())
  p<- p +   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  print(p)

}

plot_gda_heatmap_DiseaseClass <- function( type, search, input, limit,nchars, verbose ) {

  # if( type == "disease-disease-gene" ) {
  #
  #   results <- input [, c (1,2, 5 ,6 )]
  #
  #   results <- as.data.frame( transDiseaseClass ( results ) )
  #   row.names(results) <- NULL
  #   results[,4] <- gsub (" ; ", ";", results[,4] )
  #   results[,4] <- gsub ( "null", NA, results[,4] )
  #   results[,4] <- gsub ( "Wounds and Injuries", NA, results[,4] )
  #   results[,4] <- gsub ( "Animal Diseases", NA, results[,4] )
  #   results$V4 <- as.character(results$V4)
  #   results <- unique(results[c("NameDisease1", "NameDisease2",  "V4")])
  #   nd <- dim(subset(results, is.na(V4)))[1]
  #   if (nd > 0){
  #     show( paste0("warning: ", nd , " disease(s) not shown in the plot" ))
  #   }
  #   for( i in 1:nrow(results)){
  #     if(nchar(as.character(results$V4[i])) > 25){
  #       rename <- substr(results$V4[i], 1, 25)
  #       results$V4[i] <- paste0(rename, "...")
  #     }
  #   }
  #   ##
  #   for ( i in 1:length(unique(results$NameDisease1))){
  #     disease1 <- results[results$NameDisease1==unique(results$NameDisease1)[i],]
  #     freq <- as.data.frame ( table ( as.character(disease1$V4 ) ))
  #     freq$perc <- (freq$Freq/sum(freq$Freq))*100
  #     freq$dis <- unique(results$NameDisease1)[i]
  #
  #     if ( i == 1){
  #       rr <- freq
  #     } else if ( i > 1){
  #       rr <- rbind(rr, freq)
  #     }
  #   }
  #
  #   rr <- rr [ with ( rr , order(-perc)), ]
  #   casted_rr <- reshape::cast(rr, Var1 ~ dis, value= "perc"  )
  #   l <- dim(casted_rr)[2]
  #   Rowm  <- rowSums(casted_rr[2:l], na.rm = T)
  #   casted_rr<- cbind(casted_rr, Rowm)
  #   ordered<-  casted_rr[order(casted_rr["Rowm"]), "Var1" ]
  #   rr$Var1  <- factor(rr$Var1  , levels= as.factor(ordered))
  #   m <- max(rr$perc)
  #   p <- ggplot2::ggplot(rr, ggplot2::aes(dis, Var1)) + ggplot2::geom_tile(ggplot2::aes(fill = perc), colour = "white") + ggplot2::scale_fill_gradient(limits = c(0,m), low = "white",   high = "royalblue2", na.value = "white", "Percentage")
  #   p <- p + ggplot2::theme_grey(base_size = 10) + ggplot2::labs(title = "Gene-Disease Class Heatmap", x = "Genes",  y = "MeSH Disease Class") + ggplot2::scale_x_discrete(expand = c(0, 0))
  #   p <- p + ggplot2::theme( plot.margin=grid::unit(x=c(5,15,5,15),units="mm") , axis.line = ggplot2::element_line(size = 0.7, color = "black"), text = ggplot2::element_text(size = 9),axis.title=ggplot2::element_text(size=10 ), axis.text.x = ggplot2::element_text(angle=45, size = 10, hjust = 1), panel.background = ggplot2::element_blank())
  #   p +   ggplot2::theme(plot.title = element_text(hjust = 0.5))
  #
  # }
  # else

  if( type == "gene-disease" ) {

    results <- ( input[ , c("geneid", "gene_symbol", "diseaseid", "disease_name", "disease_class_name") ] )
    results$disease_class_name <- gsub ( "Wounds and Injuries", NA, results$disease_class_name )
    results$disease_class_name  <- gsub ( "Animal Diseases", NA, results$disease_class_name )
    results <- tidyr::separate_rows(results, disease_class_name, sep = "; ")
    results$disease_class_name <-ifelse(results$disease_class_name == "", NA, as.character(results$disease_class_name))


    nd <- length(unique(as.character(subset(results, is.na(disease_class_name))$disease_name)))
    if (nd > 0){
      show( paste0("warning: ", nd , " disease(s) not shown in the plot" ))
    }
    results <- subset(results, !is.na(disease_class_name))

    for( i in 1:nrow(results)){
      if(nchar(as.character(results$disease_class_name[i])) > nchars){
        rename <- substr(results$disease_class_name[i], 1, nchars)
        results$disease_class_name[i] <- paste0(rename, "...")
      }
    }

    for ( i in 1:length(unique(results$gene_symbol))){
      disease1 <- results[results$gene_symbol==unique(results$gene_symbol)[i],]
      freq <- as.data.frame ( table ( as.character(disease1$disease_class_name ) ))
      freq$perc <- (freq$Freq/sum(freq$Freq))*100
      freq$dis <- unique(results$gene_symbol)[i]

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
    p <- p + ggplot2::theme( plot.margin=grid::unit(x=c(5,15,5,15),units="mm") , axis.line = ggplot2::element_line(size = 0.7, color = "black"), text = ggplot2::element_text(size = 9),axis.title=ggplot2::element_text(size=10 ), axis.text.x = ggplot2::element_text(angle=45, size = 10, hjust = 1), panel.background = ggplot2::element_blank())
    p<-  p +   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    print(p)

  }

  else if( type == "variant-disease"  ) {

    if( search == "single" ) {
      stop( "For this type of chart, a multiple query created with 'variant2disease' is required." )
    }

    results <- unique( input[ , c( "variantid" ,"diseaseid" , "disease_name", "disease_class_name") ] )
    results$disease_class_name <- gsub ( "Wounds and Injuries", NA, results$disease_class_name )
    results$disease_class_name  <- gsub ( "Animal Diseases", NA, results$disease_class_name )
    results <- tidyr::separate_rows(results, disease_class_name, sep = "; ")
    results$disease_class_name <-ifelse(results$disease_class_name == "", NA, as.character(results$disease_class_name))


    nd <- length(unique(as.character(subset(results, is.na(disease_class_name))$disease_name)))
    if (nd > 0){
      show( paste0("warning: ", nd , " disease(s) not shown in the plot" ))
    }
    results <- subset(results, !is.na(disease_class_name))

    for( i in 1:nrow(results)){
      if(nchar(as.character(results$disease_class_name[i])) > nchars){
        rename <- substr(results$disease_class_name[i], 1, nchars)
        results$disease_class_name[i] <- paste0(rename, "...")
      }
    }


    for ( i in 1:length(unique(results$variantid))){
      disease1 <- results[results$variantid==unique(results$variantid)[i],]
      freq <- as.data.frame ( table ( as.character(disease1$disease_class_name ) ))
      freq$perc <- (freq$Freq/sum(freq$Freq))*100
      freq$snp <- unique(results$variantid)[i]

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
    p <- p + ggplot2::theme_grey(base_size = 10) + ggplot2::labs(title = "Variant-Disease Class Heatmap", x = "Variants",  y = "Disease Class") + ggplot2::scale_x_discrete(expand = c(0, 0))
    p <- p + ggplot2::theme( plot.margin=grid::unit(x=c(5,15,5,15),units="mm") , axis.line = ggplot2::element_line(size = 0.7, color = "black"), text = ggplot2::element_text(size = 9),axis.title=ggplot2::element_text(size=10 ), axis.text.x = ggplot2::element_text(angle=45, size = 10, hjust = 1), panel.background = ggplot2::element_blank())
    p<- p +   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    print(p)
  }
}

plot_gda_Panther <- function ( input, layout, search, verbose ) {

  if( search == "single"){
    results <- input [, c ("diseaseid", "disease_name", "geneid", "gene_symbol", "protein_class_name" )]

    results <- data.frame( tidyr::separate_rows(results, protein_class_name, sep = "; ") )
    results$protein_class_name <-ifelse(results$protein_class_name == "", NA, as.character(results$protein_class_name))


    nd <- length(unique(as.character(subset(results, is.na(protein_class_name))$gene_symbol)))
    if (nd > 0){
      show( paste0("warning: ", nd , " gene(s) not shown in the plot" ))
    }
    results <- subset(results, !is.na(protein_class_name))

    freq <- as.data.frame ( table ( as.character(results$protein_class_name ) ))
    freq$perc <- (freq$Freq/sum(freq$Freq))*100

    rw <- matrix(c(as.character(results[1,2]), 10, 10),ncol=3)
    colnames(rw) <- colnames(freq)
    freq <- rbind(freq,rw)

    sizes <- as.numeric(freq[ , 3 ])
    names( sizes ) <- freq[ , 1 ]

    #plot
    edges <- data.frame( results[c("disease_name", "protein_class_name") ] )
    netw  <- igraph::graph.data.frame( edges, directed = FALSE )
    netw  <- igraph::simplify( netw )
    #lay   <- layout( netw )
    lay <- apply_layout(layout, netw)

    if( verbose ) {
      message( "The network contains ", igraph::vcount( netw ), " nodes and ", igraph::ecount( netw ), " edges." )
    }

    pantherClass <- unique( results$protein_class_name )
    ttl <- paste0("Protein Classes associated to ", results$disease_name[ 1 ] )

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

}


plot_diseasome <- function( input, type, layout, search, verbose, prop, limit, cutoff , database, pvalue, nboot, ncores, interactive ) {

  indexdiseases <- unique(c(as.character(input$disease1_name)))
  diseases<- unique(c(as.character(input$diseaseid1), as.character(input$diseaseid2)))

  title <- ""

  if( type == "disease-disease-gene"  ) {
    title <- "DDA network via shared genes"

    input$pvalue_jaccard_genes <- as.numeric(as.character(input$pvalue))
    input  <- input[input$pvalue_jaccard_genes < cutoff,]
    res <- disease2gene(disease = diseases, database = database)

    res <- res@qresult
    pairs <- ( length(diseases) *   length(diseases) -  length(diseases) ) / 2
    tb <- as.data.frame(setNames(replicate(8,numeric(0), simplify = F), c("diseaseid1","disease1_name", "diseaseid2","disease2_name", "ngenes1", "ngenes2", "common", "jaccard_index")))
    cnt <- 1
    for ( i in 1:length(diseases)){
      dis1 <- diseases[i]
      nd1 <- as.character (  res[res[,"diseaseid"]==diseases[i],"disease_name"][1]  )
      gene1 <- unique( subset(res, diseaseid == dis1 )$geneid)
      for ( j in 2:length(diseases)){
        if ( i < j){
          dis2 <- diseases[j]
          nd2 <- as.character (  res[res[,"diseaseid"]==diseases[j],"disease_name"][1]  )
          gene2 <- unique( subset(res, diseaseid == dis2 )$geneid)
          ji <- round(length(intersect(gene1, gene2))/length(union(gene1, gene2)), 5)

          tb[ cnt, ] <- c(diseases[ i ], nd1, diseases[ j ],nd2, length( gene1 ), length( gene2 ), length( intersect( gene1, gene2 ) ), ji )
          cnt <- cnt + 1
        }
      }
    }
    tb <- tb[tb[,"common"]>0 ,]

    finaldiseases <- unique( c(as.character(tb$diseaseid1), as.character(tb$diseaseid2)))
    #   if (dim(input)[1]> limit){
    #     input <- input[c(1:limit),]
    #     show( paste0("warning: dataframe of ", a , " rows has been reduced to ", limit, " rows." ))
    #   }

    if(pvalue == T){
      if( verbose ) {
        message(paste0("performing ", nboot, " randomizations. Please wait..."))
      }
      tb <- pvalue_estimation(tb, database = "CURATED", ncores = ncores, nboot = nboot)
    }
    tb <- subset(tb, pvalue < 0.05)
    edges <- data.frame( tb[ , "disease1_name" ], tb[ , "disease2_name" ] )
    netw  <- igraph::graph.data.frame( edges, directed = FALSE )
    netw  <- igraph::simplify( netw )
    lay <- apply_layout(layout = layout, netw = netw)


    ds1 <- tb[,c( "disease1_name", "ngenes1" ) ]
    colnames( ds1 ) <- c("dis", "gen")
    ds2 <- tb[, c("disease2_name", "ngenes2") ]
    colnames( ds2 ) <- c("dis", "gen")

    fnl <- rbind ( ds1, ds2)
    fnl <- fnl[!duplicated(fnl), ]

    sizes <- as.numeric( fnl[ , 2 ] )
    names( sizes ) <- fnl[ , 1 ]
  }
   else if( type == "disease-disease-variant" ) {
    title <- "DDA network via shared variants"
    input$pvalue_jaccard_variants <- as.numeric(as.character(input$pvalue))

    res <- disease2variant(disease = diseases, database = database)

    res <- res@qresult
    pairs <- ( length(diseases) *   length(diseases) -  length(diseases) ) / 2
    tb <- as.data.frame(setNames(replicate(8,numeric(0), simplify = F), c("diseaseid1","disease1_name", "diseaseid2","disease2_name", "nvariants1", "nvariants2", "common", "jaccard_index")))
    cnt <- 1
    for ( i in 1:length(diseases)){
      dis1 <- diseases[i]
      nd1 <- as.character (  res[res[,"diseaseid"]==diseases[i],"disease_name"][1]  )
      gene1 <- unique( subset(res, diseaseid == dis1 )$variantid)
      for ( j in 2:length(diseases)){
        if ( i < j){
          dis2 <- diseases[j]
          nd2 <- as.character (  res[res[,"diseaseid"]==diseases[j],"disease_name"][1]  )
          gene2 <- unique( subset(res, diseaseid == dis2 )$variantid)
          ji <- round(length(intersect(gene1, gene2))/length(union(gene1, gene2)), 5)

          tb[ cnt, ] <- c(diseases[ i ], nd1, diseases[ j ],nd2, length( gene1 ), length( gene2 ), length( intersect( gene1, gene2 ) ), ji )
          cnt <- cnt + 1
        }
      }
    }
    tb <- tb[tb[,"common"]>0 ,]
    finaldiseases <- unique( c(as.character(tb$diseaseid1), as.character(tb$diseaseid2)))
    #   if (dim(input)[1]> limit){
    #     input <- input[c(1:limit),]
    #     show( paste0("warning: dataframe of ", a , " rows has been reduced to ", limit, " rows." ))
    #   }

    if(pvalue == T){
      if( verbose ) {
        message(paste0("performing ", nboot, " randomizations. Please wait..."))
      }
      tb <- pvalue_estimation(tb, database = "CURATED", ncores = ncores, nboot = nboot)
    }
    tb <- subset(tb, pvalue < 0.05)

    edges <- data.frame( tb[ , "disease1_name" ], tb[ , "disease2_name" ] )
    netw  <- igraph::graph.data.frame( edges, directed = FALSE )
    netw  <- igraph::simplify( netw )
    lay <- apply_layout(layout = layout, netw = netw)


    ds1 <- tb[,c( "disease1_name", "nvariants1" ) ]
    colnames( ds1 ) <- c("dis", "gen")
    ds2 <- tb[, c("disease2_name", "nvariants2") ]
    colnames( ds2 ) <- c("dis", "gen")

    fnl <- rbind ( ds1, ds2)
    fnl <- fnl[!duplicated(fnl), ]

    sizes <- as.numeric( fnl[ , 2 ] )
    names( sizes ) <- fnl[ , 1 ]

  }

  if( verbose ) {
    message( "The network contains ", igraph::vcount( netw ), " nodes and ", igraph::ecount( netw ), " edges." )
  }

  diseaseColor = "blue"
  comorColor = "lightblue"
  if(interactive == FALSE){
    igraph::plot.igraph( netw,
                         vertex.frame.color = "white",
                         layout              = lay,
                         vertex.color        = ifelse(igraph::V(netw)$name %in% indexdiseases, diseaseColor,comorColor),
                         vertex.label.dist   = 0,      #puts the name labels slightly off the dots
                         vertex.frame.color  = 'blue', #the color of the border of the dots
                         vertex.label.color  = 'black',#the color of the name labels
                         vertex.label.font   = 0,      #the font of the name labels
                         vertex.label        = igraph::V( netw )$name, #specifies the lables of the vertices. in this case the 'name' attribute is used
                         #edge.label          = round(as.numeric(as.character(tb[,"jaccard_index"])), 3),
                         edge.color          = "darkgrey",
                         edge.width          =  1 + as.numeric(as.character(tb[,"jaccard_index"])) *3,
                         edge.arrow.size     = 0.5,
                         vertex.size         = as.numeric( sizes[ igraph::V( netw )$name ] * prop),
                         vertex.size         = 0.5,
                         vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
                         main                = title
    )

  }else if (interactive == TRUE){
    igraph::tkplot(netw,
                           vertex.frame.color = "white",
                           layout              = lay,
                           vertex.color        = ifelse(igraph::V(netw)$name %in% indexdiseases, diseaseColor,comorColor),
                           vertex.label.dist   = 0,      #puts the name labels slightly off the dots
                           vertex.frame.color  = 'blue', #the color of the border of the dots
                           vertex.label.color  = 'black',#the color of the name labels
                           vertex.label.font   = 0,      #the font of the name labels
                           vertex.label        = igraph::V( netw )$name, #specifies the lables of the vertices. in this case the 'name' attribute is used
                           #edge.label          = round(as.numeric(as.character(tb[,"jaccard_index"])), 3),
                           edge.color          = "darkgrey",
                           edge.width          =  1 + as.numeric(as.character(tb[,"jaccard_index"])) *3,
                           edge.arrow.size     = 0.5,
                           vertex.size         = as.numeric( sizes[ igraph::V( netw )$name ] * prop),
                           vertex.size         = 0.5,
                           vertex.label.cex    = 0.8,    #specifies the size of the font of the labels
                           main                = title
    )
  }
}


plot_enrichment <- function( input, verbose, cutoff , count, nchars, limit ) {

  input<- subset(input, FDR < cutoff & Count > count  )
  if ( dim( input )[ 1 ] > limit ){
    input <- input[ 1:limit ,]
    show( paste0("warning: dataframe of ", nrow(input), " rows has been reduced to ", limit, " rows." ))
  }

  idx <- order(input$gg, decreasing = T)
  input$Description <- ifelse(nchar(input$Description) < nchars, as.character(input$Description), paste0(substr(input$Description, 0, nchars), "..."))
  input$Description<-factor(input$Description, levels=rev(unique(input$Description[idx])))
  title <- "DisGeNET enrichment"
  p <- ggplot2::ggplot(input, ggplot2::aes_string(x= input$gg, y="Description", size=input$Count, color=input$FDR)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_continuous(low="red", high="blue", name = "FDR", guide=ggplot2::guide_colorbar(reverse=TRUE)) +
    ggplot2::xlab("Ratio") + ggplot2::ylab(NULL) +
    ggplot2::ggtitle(title) + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(colour = "black",  size = 14, vjust = 1),
                   axis.text.y = ggplot2::element_text(colour = "black", size = 14, hjust = 1),
          axis.title = ggplot2::element_text(margin = ggplot2::margin(10, 5, 0, 0), color = "black", size = 12),
          axis.title.y = ggplot2::element_text(angle = 90)) + ggplot2::scale_size(range=c(3, 8))
   print(p)
}

