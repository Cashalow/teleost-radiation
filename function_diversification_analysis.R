library(picante)


source("reading_trees_polyploid_nodes.R")
source("extract_polyploid_clades.R")
source("run_tp_until_not_significant.R")
source("diversification-functions-from-RPANDA/fit.ln.timevar_errap_log.R")
source("diversification-functions-from-RPANDA/getLikelihood.ln.timevar_errap_log.R")
source("diversification-functions-from-RPANDA/Phi_timevar_errap.R")
source("diversification-functions-from-RPANDA/Psi_timevar_errap.R")





read_results<- function(folder){
  res <- list()
  for (i in list.files(folder, full.names = FALSE)){
    load(paste(folder, i, sep="/"))
    name <- paste(strsplit(i, "_")[[1]][1], strsplit(i, "_")[[1]][2], sep="_")
    if (!is.list(res[[name]])) res[[name]] <- list()
    for (j in 1:length(results)){
      if (!is.null(results[[j]])) res[[name]][[j]] <- results[[j]]
    }
  }
  return(res)
}



linear <- function(x,y) {y[1]+y[2]*x}
constant <- function(x,y) {y[1]}
expo <- function(x,y) {y[1]*exp(y[2]*x)}


extract_mean_mu_lam <- function(results){
  poly_mu <- c()
  poly_lam <- c()
  di_mu <- c()
  di_lam <- c()
  for (i in 1:length(results)){
    res_poly <- return_mu_lam_tree(results[[i]]$res_poly[[2]])
    poly_lam <- c(poly_lam, res_poly[[1]])
    poly_mu <- c(poly_mu, res_poly[[2]])
    res_di <- return_mu_lam_tree(results[[i]]$res_di[[2]])
    di_lam <- c(di_lam, res_di[[1]])
    di_mu <- c(di_mu, res_di[[2]])
  }
  return(list(mpl=mean(poly_lam), mpm=mean(poly_mu), mdl=mean(di_lam), mdm=mean(di_mu)))
}


return_mu_lam_tree <- function(single_res){
  nb_shifts <- length(single_res)-1
  div <- single_res[[nb_shifts]][(nb_shifts+2):(nb_shifts*2+1)]
  turn <- single_res[[nb_shifts]][2:(nb_shifts+1)]
  return(list(div/(1-turn), div/(1-turn)-div))
}


loop_over_clades <- function(clades, func, starting_point, sampling){
  r1 <- func(clades$poly, list(starting_point$mpl, starting_point$mpm), sampling[1])
  if (is.null(clades$di)) {r2 <- NULL}
  else if (Ntip(clades$di) > 2) { r2 <- func(clades$di, list(starting_point$mdl, starting_point$mdm), sampling[2])}
  else {r2 <- NULL}
  return( list(res_poly=r1, res_di=r2))
}



run_diversification_tp <- function(tree, starting_point, sampling){
  step <- branching.times(tree)[1]/50
  print(sampling)
  return(run_tp_looking_for_rateshift_until_not_significant(getx(tree), 0, step, branching.times(tree)[1], 0.99, sampling))
}



run_diversification_rpanda <- function(clade, starting_point, sampling){
  funcs <- c("constant", "linear", "expo")
  starting_points_lam <- list(starting_point[[1]], c(starting_point[[1]], 0), c(starting_point[[1]], 0))
  starting_points_mu <- list(starting_point[[2]], c(starting_point[[2]], 0), c(starting_point[[2]], 0))
  print(starting_points_lam)
  print(starting_points_mu)
  models <- list()
  for (i in 1:length(funcs)){
    for (j in 1:length(funcs)){
      k <- paste(funcs[i], funcs[j], sep="_")
      print(k)
      models[[k]] <- fit.ln.timevar_errap_log(list(clade), branching.times(clade)[1], get(funcs[i]), get(funcs[j]), starting_points_lam[[i]], starting_points_mu[[j]])
    }
  }
  return(models)
}


loop_over_trees <- function(trees, poly_node, cons, name_clade, name_poly, func, starting_point=NULL, begin=0, end=500, dest="./"){
  print(func)
  results <- list()
  if (name_clade=="salmoniformes") sampling <- c(.28, .69)
  else sampling <- c(1, 1)
  for (i in begin:end){
    if (checking_monophyly_polyploids(cons, trees[[i]], poly_node)){
      clades <- crop_clades(trees[[i]], names(descendants(phylo4(cons), poly_node, "tips")))
      results[[i]] <- loop_over_clades(clades, get(func), starting_point, sampling)
      save(results, file=paste(dest, as.character(name_clade), "_", as.character(name_poly), "_", func, "_", begin, "_", end, "_results.R", sep=""))
    }
 }
}






