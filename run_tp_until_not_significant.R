library(TreePar)


run_tp_looking_for_mee_until_not_significant <- function(br.ti, start, step, mx.ti, pvalue, sampling=1, pdiv=FALSE){
  output <- bd.shifts.optim(br.ti, c(sampling, 1), step, start, mx.ti, ME=TRUE, all=FALSE, posdiv=pdiv)
  i <- length(output[[2]])
  while ((pchisq(2*(output[[2]][[i-1]][1]-output[[2]][[i]][1]), df=2)) > pvalue)
    {
      surv_rates <- output[[2]][[i]][4:(4+i-2)]
      i <- i+1
      output <- bd.shifts.optim(br.ti, c(sampling, surv_rates, 1), step, start, mx.ti, ME=TRUE, all=FALSE, miniall=output[[2]], posdiv=pdiv)
    }
  return(output)
}



run_tp_looking_for_rateshift_until_not_significant <- function(br.ti, start, step, mx.ti, pvalue, sampling, pdiv=FALSE){
  output <- bd.shifts.optim(br.ti, c(sampling, 1), step, start, mx.ti, posdiv=pdiv)
  i <- length(output[[2]])
  while ((pchisq(2*(output[[2]][[i-1]][1]-output[[2]][[i]][1]), df=3)) > pvalue)
    {
      output <- bd.shifts.optim(br.ti, c(sampling, rep(1,i)), step, start, mx.ti, miniall=output[[2]], posdiv=pdiv)
      i <- i+1
    }
  return(output)
}
