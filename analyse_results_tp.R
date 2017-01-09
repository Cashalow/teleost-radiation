library(plotrix)
library(phylobase)
library(RSvgDevice)
library(ape)


source("extract_polyploid_clades.R")
source("reading_stored_results.R")
source("reading_trees_polyploid_nodes.R")



plot_single_clade_results <- function(clade_res){
  timemax <- max(clade_res[[3]])
  print(timemax)
  nb_shifts <- length(clade_res[[2]])-1
  estimates <- clade_res[[2]][[nb_shifts]][(nb_shifts+2):(nb_shifts*2+1)]
  time <- tail(clade_res[[2]][[nb_shifts]], nb_shifts-1)
  time <- sort(c(time, time, 0, timemax))
  div <- c()
  for (j in estimates) div <- c(div, j, j)
  return(list(time=-time, div=div))
}


plot_list_tp_diversification_values <- function(results, trees, cons, node, clade, color, mean_time){
  par(mfrow = c(2,1), xaxs = "i")
  par(mar=c(1,7,0,0))
  timemax <- plot_polyploid_diploid_clade(cons, node, clade, mean_time, color)
  line_poly <- list()
  line_di <- list()
  ratemax <- c()
  ratemin <- c()
  par(mar=c(5,7,1,4))
  for (i in 1:length(results)) {
    line_poly[[i]] <- plot_single_clade_results(results[[i]]$res_poly)
    line_di[[i]] <- plot_single_clade_results(results[[i]]$res_di)
    ratemax <- max(line_poly[[i]]$div, line_di[[i]]$div, ratemax)
    ratemin <- min(line_poly[[i]]$div, line_di[[i]]$div, ratemin)
  }
  revaxis(c(0, timemax ), c(ratemin, ratemax), col = "#ffffff00", xlab = "Time before present", ylab = "Diversification rate\n\n\n", frame.plot = FALSE, yrev=FALSE, xrev=TRUE, yside=2)
  for (i in 1:length(results)){
    lines(line_poly[[i]]$time, line_poly[[i]]$div, col=color)
    lines(line_di[[i]]$time, line_di[[i]]$div, col="black")
  }
  return(c(ratemax, ratemin))
}

write_tp_results <- function(file, dest, polyploid_color, means){
  all_tp_results <- read_results(file)
  boundaries <- list()
  for (i in names( all_tp_results)){
    print(i)
    devSVG(file=paste(dest, i, "_tp.svg", sep=""), width=12, height=8)
    boundaries[[i]] <- plot_list_tp_diversification_values( all_tp_results[[i]], get(strsplit(i, "_")[[1]][1]), get(paste(strsplit(i, "_")[[1]][1], "_cons", sep="")), get(paste("polyploid_node_", strsplit(i, "_")[[1]][1], "_cons", sep=""))[as.integer(strsplit(i, "_")[[1]][2])], i, polyploid_color, means[[i]])
    dev.off()
  }
  return(boundaries)
}

result_tp_folder <- "results_tp/" #where the results for the tp analysis files are stored
plot_folder <- "plots/" #where the plots will be stored
color = rgb(95, 128, 77, max=255) #color for polyploids


boundaries <- write_tp_results(result_tp_folder, plot_folder, color, means)




