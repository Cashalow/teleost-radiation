

library(phylobase)


crop_clades <- function(tree, polyploids){
  node <- MRCA(phylo4(tree), polyploids)
  cropped_tree <- extract.clade(tree, ancestors(phylo4(tree), node, "parent"))
  poly <- extract.clade(tree, node)
  sister_node <- setdiff(children(phylo4(tree), ancestors(phylo4(tree), node, "parent")), node)
  if (sister_node > Ntip(tree)){di <- extract.clade(tree, sister_node)}
  else {di <- NULL}
  return(list(poly=poly, di=di, total=cropped_tree))
}

checking_monophyly_polyploids <- function(consensus, trees, node){
  polyploids <- names(descendants(phylo4(consensus), node, type = "tips"))
  anc <- names(descendants(phylo4(trees), MRCA(phylo4(trees), polyploids), type="tips"))
  if (!(identical(setdiff(polyploids, anc), character(0))) || !(identical(setdiff(anc, polyploids), character(0)))) {return(FALSE)}
  return(TRUE)
}


extract_nodes_from_trees <- function(consensus, trees, node){
  nodes <- matrix(ncol=3, nrow=length(trees), dimnames=list(1:length(trees), c("poly", "di", "ancestor")))
  polyploids <- names(descendants(phylo4(consensus), node, type = "tips"))
  for (i in 1:length(trees)){
    if (checking_monophyly_polyploids(consensus, trees[[i]], node)){
      poly <- MRCA(phylo4(trees[[i]]), polyploids)
      ancestor <- ancestor(phylo4(trees[[i]]), poly)
      di <-  setdiff(children(phylo4(trees[[i]]), ancestor), poly)
      nodes[i, "poly"] <- poly
      nodes[i, "ancestor"] <- ancestor
      nodes[i, "di"] <- di
    }
  }
  return(nodes)
}




plot_polyploid_diploid_clade <- function(cons, node, clade, mean_time, col=red){
  k <- extract.clade(cons, ancestor(phylo4(cons), node))
  cons_time <- branching.times(k)[1]
  k$edge.length <- k$edge.length*mean_time/cons_time
  l <- MRCA(phylo4(k), names(descendants(phylo4(cons), node, "tips")))
  edg=rep("black", Nedge(k))
  edg <- setNames(edg, k$edge[,2])
  edg[as.character(descendants(phylo4(k), l, "all"))] <- col
  k$node.label <- rep("", Nnode(k))
  k$node.label[l-Ntip(k)] <- capwords(gsub("_", "", clade))
  plot(k, edge.color=edg, show.tip.label=FALSE, show.node.label=TRUE)
  axisPhylo()
  return(branching.times(k)[1])
}


mean_time_of_ancestor <- function(trees, cons, node){
  polyploids <- names(descendants(phylo4(cons), node, "tips"))
  times <- c()
  for (i in trees){
    tree_node <- MRCA(phylo4(i), polyploids)
    times <- c(times, branching.times(i)[as.character(ancestor(phylo4(i), tree_node))])
  }
  print(times[!is.na(times)])
  return(mean(times[!is.na(times)]))
}



capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s,1,1)),
                           {s <- substring(s,2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}


#for (i in 1:length(fuchsia)){
# CairoSVG(file=as.character(i))
 #fuchsia[[i]]$node.label <- 1:(Ntip(fuchsia[[i]])-1)+Ntip(fuchsia[[i]])
 #plot(fuchsia[[i]], show.node.label = TRUE)
 #dev.off()
 #}
