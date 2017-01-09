
library(plotrix)
library(Cairo)
library(RSvgDevice)

source("extract_polyploid_clades.R")
source("reading_trees_polyploid_nodes.R")
source("analyse_results_tp.R")


return_diversification_function <- function(model, name){
  print(name)
  model <- model[[name]]
  print(model$lamb_par)
  print(model$mu_par)
  if (name == "constant_constant") { return( function(x) model$lamb_par-model$mu_par+x*0)}
  else if (name == "constant_linear") { return( function(x) model$lamb_par-(model$mu_par[1]+x*model$mu_par[2]))}
  else if (name == "constant_expo") { return( function(x) model$lamb_par - (model$mu_par[1]*exp(x*model$mu_par[2])))}
  else if (name == "linear_constant") { return( function(x) model$lamb_par[1] + x*model$lamb_par[2] - model$mu_par)}
  else if (name == "linear_linear") { return( function(x) model$lamb_par[1]+x*model$lamb_par[2]-(model$mu_par[1]+model$mu_par[2]*x))}
  else if (name =="linear_expo") { return( function(x) model$lamb_par[1]+x*model$lamb_par[2]-(model$mu_par[1]*exp(x*model$mu_par[2])))}
  else if (name =="expo_constant"){ return( function(x) model$lamb_par[1]*exp(x*model$lamb_par[2])-model$mu_par)}
  else if (name =="expo_linear"){ return( function(x) model$lamb_par[1]*exp(x*model$lamb_par[2])-(model$mu_par[1]+x*model$mu_par[2]))}
  else if (name =="expo_expo"){return( function(x) model$lamb_par[1]*exp(x*model$lamb_par[2])-model$mu_par[1]*exp(x*model$mu_par[2]))}
}



select_best_model <- function(res, nb_param){
  index <- sapply(res, "[[", "aicc")
  prob <- exp((min(index)-index)/2)
  most_probable <- nb_param[names(prob[prob>0.05])]
  if (length(which(most_probable==min(most_probable)))==1) {return(names(which(most_probable==min(most_probable))))}
  else {
    return(names(which(prob==1)))
  }
}


plot_list_diversification_values <- function(res, nb_param, timemax, bound, cons, node, trees, clade, color, mean_time){
  par(mfrow = c(2,1), xaxs = "i")
  par(mar=c(1,7,0,0))
  timemax <- plot_polyploid_diploid_clade(cons, node, strsplit(clade, "_")[[1]][1], mean_time, color)
  par(mar=c(5,7,1,4))
  revaxis(x=c(0, timemax), y=c(0, 500), yrev=FALSE, xrev=TRUE, col="#ffffff00", xlab="Time before present", ylab = "Diversification rate\n\n\n", yside=2)
  ratemax <- c()
  ratemin <- c()
  func_poly <- list()
  func_di <- list()
  poly <- names(descendants(phylo4(cons), node, "tips"))
  for (i in 1:length(res)){
    if (checking_monophyly_polyploids(cons, trees[[i]], node)){
      lol <- crop_clades(trees[[i]], poly)
      if (!is.null(res[[i]]$res_poly)){
	print("poly")
        k <- select_best_model(res[[i]]$res_poly, nb_param)
        func_poly[[i]] <- return_diversification_function(res[[i]]$res_poly, k)
        max <- branching.times(lol$poly)[1]
        l <- rev(func_poly[[i]](seq(0, max, length=200)))
        lines(seq(-max, 0, length=200), l, type="l", col=color)
      }
      if (!is.null(res[[i]]$res_di)){
	print("di")
        k2 <- select_best_model(res[[i]]$res_di, nb_param)
	func_di[[i]] <- return_diversification_function(res[[i]]$res_di, k2)
        if (!is.null(lol$di))  max <- branching.times(lol$di)[1]
        else max <- timemax
        l <- rev(func_di[[i]](seq(0, max, length=200)))
        lines(seq(-max, 0, length=200), l, type="l", col="black")
      }
    }
  }
}


write_rpanda_results <- function(file, dest, boundaries, means){
  results <- read_results(file)
  print(names(results))
  num_param_models <- setNames(c(2, 3, 3, 3, 4, 4, 3, 4, 4), c("constant_constant", "constant_linear", "constant_expo", "linear_constant", "linear_linear", "linear_expo", "expo_constant", "expo_linear", "expo_expo"))
  for (i in names(results)[7]){
    print(i)
    devSVG(file=paste(dest, i, ".svg", sep=""), width=12, height=8)
    plot_list_diversification_values(results[[i]], num_param_models, branching.times(get(paste(strsplit(i, "_")[[1]][1], "_cons", sep="")))[1], boundaries[[i]], get(paste(strsplit(i, "_")[[1]][1], "_cons", sep="")), get(paste("polyploid_node_", strsplit(i, "_")[[1]][1], "_cons", sep=""))[as.integer(strsplit(i, "_")[[1]][2])], get(strsplit(i, "_")[[1]][1]), i, rgb(95, 128, 77, max=255), means[[i]])#rgb(95, 128, 77, max=255))
    dev.off()
  }
}


result_rpanda_folder <- "results_rpanda/"



write_rpanda_results(result_rpanda_folder, plot_folder, boundaries, means)


