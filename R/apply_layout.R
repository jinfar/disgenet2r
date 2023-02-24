apply_layout <- function( layout, netw){
  if(layout == "layout.fruchterman.reingold"){
    lay   <- igraph::layout.fruchterman.reingold( netw )
  }else if(layout == "layout.auto"){
    lay   <- igraph::layout.auto( netw )
  }else if(layout == "layout.random"){
    lay   <- igraph::layout.random( netw )
  }else if(layout == "layout.circle"){
    lay   <- igraph::layout.circle( netw )
  }else if(layout == "layout.sphere"){
    lay   <- igraph::layout.sphere( netw )
  }else if(layout == "layout.fruchterman.reingold"){
    lay   <- igraph::layout.fruchterman.reingold( netw )
  }else if(layout == "layout.kamada.kawai"){
    lay   <- igraph::layout.kamada.kawai( netw )
  }else if(layout == "layout.spring"){
    lay   <- igraph::layout.spring( netw )
  }else if(layout == "layout.reingold.tilford"){
    lay   <- igraph::layout.reingold.tilford( netw )
  }else if(layout == "layout.fruchterman.reingold.grid"){
    lay   <- igraph::layout.fruchterman.reingold.grid( netw )
  }else if(layout == "layout.lgl"){
    lay   <- igraph::layout.lgl( netw )
  }else if(layout == "layout.svd"){
    lay   <- igraph::layout.svd( netw )
  }else if(layout == "layout.graphopt"){
    lay   <- igraph::layout.graphopt( netw )
  }else if(layout == "layout.norm"){
    lay   <- igraph::layout.norm( netw )
  }else{
    stop
    message("Check if your layout is included in the igraph layout options:
            layout.auto, layout.random, layout.circle, layout.sphere,
            layout.fruchterman.reingold, layout.kamada.kawai, layout.spring,
            layout.reingold.tilford, layout.fruchterman.reingold.grid,
            layout.lgl, layout.graphopt, layout.svd, layout.norm"
    )
  }
  return(lay)
}
