library(ape)
library(picante)

linear <- function(x,y) {y[1]+y[2]*x}
constant <- function(x,y) {y[1]}
expo <- function(x,y) {y[1]*exp(y[2]*x)}




source("fit.ln.timevar_errap_log.R")
source("getLikelihood.ln.timevar_errap_log.R")
source("Phi_timevar_errap.R")
source("Psi_timevar_errap.R")

cetacean <- read.tree("../../Data/Mammals/cetacean_phylogeny.txt")

fit.ln.timevar_errap_log(list(cetacean), branching.times(cetacean)[1], expo, constant, c(.1,.1), 0, cst.mu = T)

poaceae <- read.nexus("../../Data/Angiosperms/poaceae.nx")

poa_div_exp <- fit.ln.timevar_errap_log(list(poaceae), branching.times(poaceae)[1], expo, constant, c(.1,.1), 0, cst.mu = T)
poa_div_linear<- fit.ln.timevar_errap_log(list(poaceae), branching.times(poaceae)[1], linear, constant, c(.1,.1), 0, cst.mu = T)

poa_div <- fit.ln.timevar_errap_log(list(poaceae), branching.times(poaceae)[1], expo, constant, c(.1,.1), 0, cst.mu = T)

cetacean$node.label <- c(88:173)
plot(cetacean, show.node.label=TRUE)

delphinidae <- extract.clade(cetacean, 140)
phocoenidae <- extract.clade(cetacean, 135)
ziphiidae <- extract.clade(cetacean, 109)
balaenopteridae <- extract.clade(cetacean, 95)

other_cetacean <- drop.tip(cetacean, c(delphinidae$tip.label, phocoenidae$tip.label, ziphiidae$tip.label, balaenopteridae$tip.label))

best_fit_other <- fit.ln.timevar_errap_log(list(other_cetacean), branching.times(other_cetacean)[1], constant, expo, .1, c(0.1,0.1), cst.mu = F, fix.mu = F, cst.lamb = T)

fit.ln.timevar_errap_log(list(other_cetacean), branching.times(other_cetacean)[1], constant, linear, .1, c(0, .5), cst.mu = F, fix.mu = F, cst.lamb = T)

curve(best_fit_other$mu_par[1]*exp(x*best_fit_other$mu_par[2]), 0, 35)
abline(b=0, a=best_fit_other$lamb_par[1])

fit.ln.timevar_errap_log(list(balaenopteridae), branching.times(balaenopteridae)[1], expo, constant, c(.1,.1), 0, cst.mu = T, fix.mu = T, cst.lamb = F)

best_fit_bala <- fit.ln.timevar_errap_log(list(balaenopteridae), branching.times(balaenopteridae)[1], expo, constant, c(.1,.1), 0, cst.mu = T, fix.mu = T, cst.lamb = F)

curve(best_fit_bala$lamb_par[1]*exp(x*best_fit_bala$lamb_par[2]), 0, 15)
