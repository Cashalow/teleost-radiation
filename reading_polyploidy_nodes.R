library(ape)


data_folder <- "./"

means <- list()



cyprininae <- read.tree(paste(data_folder, "Cyprininae.mb_trees.shuffled.pruned", sep=""))
cyprininae_cons <- read.tree(paste(data_folder, "Cyprininae.mb_trees.pruned.con.newick", sep=""))
polyploid_node_cyprininae_cons <- c(429, 507, 527, 580, 588, 362, 479)


for (i in 1:length(polyploid_node_cyprininae_cons)){
  print(i)
  means[[paste("cyprininae", as.character(i), sep="_")]] <- mean_time_of_ancestor(cyprininae, cyprininae_cons, polyploid_node_cyprininae_cons[i])
}




botiidae_cons <- read.tree(paste(data_folder, "Botiidae.mb_trees.pruned.con.newick", sep=""))
botiidae <- read.tree(paste(data_folder, "Botiidae.mb_trees.shuffled.pruned", sep=""))
polyploid_node_botiidae_cons <- 37

means[["botiidae_1"]] <- mean_time_of_ancestor(botiidae, botiidae_cons, polyploid_node_botiidae_cons)



acipenseridae <- read.tree(paste(data_folder, "Acipenseridae.mb_trees.shuffled.pruned", sep=""))
acipenseridae_cons <- read.tree(paste(data_folder, "Acipenseridae_MULTILOCI.mb_trees.pruned.con.newick", sep=""))
polyploid_node_acipenseridae_cons <- c(40, 32)


for (i in 1:length(polyploid_node_acipenseridae_cons)){
  means[[paste("acipenseridae", as.character(i), sep="_")]] <- mean_time_of_ancestor(acipenseridae, acipenseridae_cons, i)
}




salmoniformes <- read.tree(paste(data_folder, "Salmoniformes.mb_trees.shuffled.pruned", sep=""))
salmoniformes_cons <- read.nexus(paste(data_folder, "Actinopterygians/Clades/Salmoniformes.mb_trees.pruned.con", sep=""))
polyploid_node_salmoniformes_cons <- 79


means[["salmoniformes_1"]] <- mean_time_of_ancestor(salmoniformes, salmoniformes_cons, polyploid_node_salmoniformes_cons)

