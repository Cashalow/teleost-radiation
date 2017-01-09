
#reads the stored results of the diversification analysis, only R files should be in the folder


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
