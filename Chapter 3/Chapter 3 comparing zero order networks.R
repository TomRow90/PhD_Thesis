#######################################################################
#
#        Section 3.5.4: Simulations and toy examples for comparing zero
#                        order networks
#
#######################################################################

#load required packages
library(simDAG)
library(qgraph)
library(bootnet)
library(dagitty)
library(welfareNet)

# set wd
setwd("/Users/tomrowland/Documents/PhD/Thesis write up/Chapter 3")


#---------------- required functions ---------------------#

# function to simulate gamma distributed response
node_gamma <- function(data, parents, betas, intercept, phi) {
  
  data <- as.data.frame(data)
  
  s <- nrow(data)
  
  if(length(parents) == 0) {
    
    log.mu <- rep(as.numeric(intercept), s)
    
  } else {
    
    X <- as.matrix(data[, parents, drop = FALSE]) 
    log.mu <- as.numeric(intercept + X %*% betas)  
    
  }
  
  
  mu <- exp(log.mu) 
  
  shape <- 1 / phi
  
  rate <- 1 / (phi*mu)
  
  out <- rgamma(s, shape=shape, rate=rate)
  
  return(out)
  
}

# function to simulate beta distributed response
node_beta <- function(data, parents, betas, intercept, phi) {
  
  data <- as.data.frame(data)
  
  s <- nrow(data)
  
  vars <- length(parents)
  
  if(length(parents) == 0) {
    
    logit.mu <- rep(as.numeric(intercept), s)
    
  } else {
    
    X <- as.matrix(data[, parents, drop = FALSE]) 
    logit.mu <- as.numeric(intercept + X %*% betas)  
    
  }
  
  
  phi <- phi     # this is shape1+shape2
  
  mu <- 1 / (1 + exp(-logit.mu)) # this is shape1/(shape1+shape2)
  
  shape1 <- mu * phi
  
  shape2 <- phi - shape1
  
  out <- rbeta(s, shape1=shape1, shape2=shape2)
  
  return(out)
  
}

#------------------------------------------------#
#--- checking error control in permutation test--#
#--- comparison of distance correlation networks-#
#--- for section 3.5.4 of Chapter 3             -#
#------------------------------------------------#

# first generate parameters for random networks
random_parameters <- expand.grid(
  nodes = c(10,15,20),
  density = c(0.1,0.2,0.3)
)


random_networks <- list()

# generate reproducible GGMs 
set.seed(123)
for(i in 1:nrow(random_parameters)) {
  
  random_networks[[i]] <- genGGM(Nvar = random_parameters[i,1],
                                 p = random_parameters[i, 2],
                                 parRange = c(0.5,1),
                                 constant = 1.5,
                                 propPositive = 0.5,
                                 graph = "random")
  
  diag(random_networks[[i]]) <- 1
}



#marginal probability distributions to transform MVN into
node_margins_10 <- c("beta","norm","exp","pois","gamma",
                     "nbinom", "beta","norm","exp","pois")
node_margins_15 <- c("beta","norm","exp","pois","gamma",
                     "nbinom", "beta","norm","exp","pois",
                     "gamma","nbinom","beta","norm","exp")
node_margins_20 <- c("beta","norm","exp","pois","gamma",
                     "nbinom", "beta","norm","exp","pois",
                     "gamma","nbinom","beta","norm","exp",
                     "pois","gamma","nbinom","beta","norm")

#parameters of marginal distributions
param_margins_10 <- list(
  list(shape1 = 6, shape2 = 2),     
  list(mean = 0, sd = 1),      
  list(rate=0.1),      
  list(lambda=20),
  list(shape=4, rate=1),     
  list(size=1, prob=0.1),
  list(shape1 = 2, shape2 = 5),
  list(mean = 0, sd = 1),
  list(rate = 0.3),
  list(lambda = 10))


param_margins_15 <- list(
  list(shape1 = 6, shape2 = 2),     
  list(mean = 0, sd = 1),      
  list(rate=0.1),      
  list(lambda=20),
  list(shape=4, rate=1),     
  list(size=1, prob=0.1),
  list(shape1 = 2, shape2 = 5),
  list(mean = 0, sd = 1),
  list(rate = 0.3),
  list(lambda = 10),
  list(shape = 2, rate = 0.5),
  list(size = 1, prob = 0.2),
  list(shape1 = 1, shape2 = 7),
  list(mean = 0, sd = 1),
  list(rate = 0.2))

param_margins_20 <- list(
  list(shape1 = 6, shape2 = 2),     
  list(mean = 0, sd = 1),      
  list(rate=0.1),      
  list(lambda=20),
  list(shape=4, rate=1),     
  list(size=1, prob=0.1),
  list(shape1 = 2, shape2 = 5),
  list(mean = 0, sd = 1),
  list(rate = 0.3),
  list(lambda = 10),
  list(shape = 2, rate = 0.5),
  list(size = 1, prob = 0.2),
  list(shape1 = 1, shape2 = 7),
  list(mean = 0, sd = 1),
  list(rate = 0.2),
  list(lambda = 20),
  list(shape = 1, rate = 2),
  list(size = 1, prob = 0.3),
  list(shape1 = 5, shape2 = 5),
  list(mean = 0, sd = 1))

# this function generates 2 mixed datasets from a true network and
# performs the permutation teston it
sim.function <- function(true.net, sample.size, margins, param.margins, perms) {
  
  #transform GGM to zero order cors
  cors <- corpcor::pcor2cor(true.net)
  
  #in vector form for copula package
  cop_params <- cors[lower.tri(cors)]
  
  # Create the Gaussian (normal) copula.
  cop <- copula::normalCopula(param = cop_params, dim = ncol(true.net), dispstr = "un")
  
  # Simulate from the copula (uniform marginals)
  u1 <- copula::rCopula(sample.size, cop)
  u2 <- copula::rCopula(sample.size, cop)
  
  # Match margin names to quantile function names
  qfuns <- paste0("q", margins)
  
  # Generate data 1
  data1 <- as.data.frame(mapply(function(qfun, col, params) {
    do.call(qfun, c(list(p = col), params))
  }, qfuns, as.data.frame(u1), param.margins))
  
  colnames(data1) <- paste0("V", seq_len(ncol(data1)))
  
  
  # Generate data 2
  data2 <- as.data.frame(mapply(function(qfun, col, params) {
    do.call(qfun, c(list(p = col), params))
  }, qfuns, as.data.frame(u2), param.margins))
  
  colnames(data2) <- paste0("V", seq_len(ncol(data2)))
  
  #estimate dcor networks
  net1 <- welfareNet::dcor.network(data1, select = "saturated")
  
  net2 <- welfareNet::dcor.network(data2, select="saturated")
  
  #run permutation test
  comparison <- welfareNet::network.comparison.analysis(data1=data1, data2=data2, 
                                                        network1 = net1, network2 = net2, 
                                                        permutations=perms)
  
  #returns P values for network metrics - these are what we will check are uniformly distributed
  return(list(frob_norm = comparison$dissimilarity.results[1, "P value"],
              density = comparison$global.results[1, "P value"],
              strength = comparison$strength.results[1, "P value"],
              edge_weight = comparison$edge.differences[1, "P value"]))
  
}      

# number of samples
reps=1000
sample.size = c(10,40,80)

# use parSim to implement monte carlo simulation across all sample sizes and
# network structures, calling them sim.function above within it
permutation_test_sim <- parSim(
  sample.size = sample.size,
  reps = reps,
  nCores = 10,
  export = c("sim.function","random_networks",
             "node_margins_10","param_margins_10",
             "node_margins_15","param_margins_15",
             "node_margins_20","param_margins_20"),
  
  expression = {
    
    if(sample.size == 10) perms = 700
    if(sample.size == 40) perms = 325
    if(sample.size == 80) perms = 263
    
    results <- list(
      list(node_number = 10, density = 0.1,perms = perms,
           out = sim.function(random_networks[[1]], sample.size,
                              node_margins_10, param_margins_10, perms = perms)),
      
      list(node_number = 15, density = 0.1,perms = perms,
           out = sim.function(random_networks[[2]], sample.size,
                              node_margins_15, param_margins_15, perms = perms)),
      
      list(node_number = 20, density = 0.1,perms = perms,
           out = sim.function(random_networks[[3]], sample.size,
                              node_margins_20, param_margins_20, perms = perms)),
      
      list(node_number = 10, density = 0.2,perms = perms,
           out = sim.function(random_networks[[4]], sample.size,
                              node_margins_10, param_margins_10, perms = perms)),
      
      list(node_number = 15, density = 0.2,perms = perms,
           out = sim.function(random_networks[[5]], sample.size,
                              node_margins_15, param_margins_15, perms = perms)),
      
      list(node_number = 20, density = 0.2,perms = perms,
           out = sim.function(random_networks[[6]], sample.size,
                              node_margins_20, param_margins_20, perms = perms)),
      
      list(node_number = 10, density = 0.3,perms = perms,
           out = sim.function(random_networks[[7]], sample.size,
                              node_margins_10, param_margins_10, perms = perms)),
      
      list(node_number = 15, density = 0.3, perms = perms,
           out = sim.function(random_networks[[8]], sample.size,
                              node_margins_15, param_margins_15, perms = perms)),
      
      list(node_number = 20, density = 0.3,perms = perms,
           out = sim.function(random_networks[[9]], sample.size,
                              node_margins_20, param_margins_20, perms = perms))
    )
    
    do.call(rbind, lapply(results, function(x) {
      data.frame(
        node_number = x$node_number,
        density = x$density,
        frob_norm = x$out$frob_norm,
        weighted_density = x$out$density,
        strength = x$out$strength,
        edge_weight = x$out$edge_weight,
        perms = perms
      )
    }))
  }
)

# save given how long it takes to run this (~10h or so)
#saveRDS(object = random_network_sim, file = "random_network_sim.rds")
permutation_test_sim <- readRDS("permutation_test_sim.rds") # file is accessible in Github thesis repository under Chapter 3 to download so set wd accordingly

# summaries for table
aggregate(frob_norm ~ sample.size + density + node_number, data = permutation_test_sim, FUN = function(x) mean(x < 0.05))
aggregate(weighted_density ~ sample.size + density + node_number, data = permutation_test_sim, FUN = function(x) mean(x < 0.05))
aggregate(strength ~ sample.size + density + node_number, data = permutation_test_sim, FUN = function(x) mean(x < 0.05))
aggregate(edge_weight ~ sample.size + density + node_number, data = permutation_test_sim, FUN = function(x) mean(x < 0.05))

# function to generate histogram plots of P value distributions
nca_histograms_plot <- function(data, x) {
  
  p <- ggplot(data, aes(x = {{ x }})) +
    geom_histogram(
      binwidth = 0.1,
      boundary = 0
    ) +
    facet_wrap(~ sample.size + density + node_number) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.1)
    ) +
    xlab(label = "") + 
    theme_classic() +
    theme(
      strip.text = element_text(size = 5),
      text = element_text(size = 5)
    )
  
  return(p)
  
}


# plots not shown in chapter main text
frob_norm_p_hist <- nca_histograms_plot(data = permutation_test_sim,
                                        x = frob_norm)

frob_norm_p_hist

density_p_hist <- nca_histograms_plot(data = permutation_test_sim,
                                      x = weighted_density)

density_p_hist

strength_p_hist <- nca_histograms_plot(data = permutation_test_sim,
                                       x = strength)

strength_p_hist

edge_p_hist <- nca_histograms_plot(data = permutation_test_sim,
                                   x = edge_weight)

edge_p_hist

#-------------------------------------------------#
#.       Fuzzy example for section 3.5.4
#-------------------------------------------------#


# structural causal model for network 1 from mixed data
network1 <- empty_dag() +
  node("A", type="rbinom", size=1, prob=0.5) +
  node("B", type="rbeta", shape1=2, shape2=6) +
  node("C", type="negative_binomial", parents="A", betas=0.6, intercept=1.5, theta=1) +
  node("D", type="gaussian", parents="A", betas=20, intercept=100, error=15) +
  node("E", type="poisson", parents="B", betas=1.5, intercept=1.5) +
  node("F", type="poisson", parents="C", betas=0.02, intercept=2) +
  node("G", type="gaussian", parents=c("F"), betas=c(0.4), intercept=25, error=6) +
  node("H", type="gaussian", parents=c("D"), betas=c(0.25), intercept=3, error=1) +
  node("I", type=node_beta, parents="E", betas=0.25, intercept=-1, phi=25) +
  node("J", type=node_beta, parents="H", betas=-0.025, intercept=3, phi=40) +
  node("K", type="gaussian", parents=c("B"), betas=c(15), intercept=20, error=3)

# structural causal model for network 2 from mixed data
network2 <- empty_dag() +
  node("A", type="rbinom", size=1, prob=0.5) +
  node("B", type=node_beta, parents=c("A"), betas=0.25, intercept=-1, phi=25) +
  node("C", type="negative_binomial", parents="A", betas=1.2, intercept=1.5, theta=1) +
  node("D", type="gaussian", parents="A", betas=20, intercept=100, error=25) +
  node("E", type="poisson", parents=c("B","C"), betas=c(2.5,0.01), intercept=1.5) +
  node("F", type="poisson", parents=c("A","C"), betas=c(0.5,0.02), intercept=2) +
  node("G", type="gaussian", parents=c("F","C"), betas=c(0.4,0.4), intercept=25, error=6) +
  node("H", type="gaussian", parents=c("D","E"), betas=c(0.25,2), intercept=3, error=1) +
  node("I", type=node_beta, parents="E", betas=0.25, intercept=-1, phi=25) +
  node("J", type=node_beta, parents="H", betas=-0.025, intercept=3, phi=40) +
  node("K", type="gaussian", parents=c("B","I"), betas=c(25,25), intercept=20, error=3)

# convert to matrix 
dag1 <- dag2matrix(network1)
dag2 <- dag2matrix(network2)

node.names <- c("A","B","C","D","E","F","G","H","I","J","K")

# get implied zero order dependence structure
zero.order.network1 <- get.zero.order.structure(dag1, node.names = node.names)
zero.order.network2 <- get.zero.order.structure(dag2, node.names = node.names)

#true density difference
network1.density <- NetworkToolbox::conn(zero.order.network1)$total
network2.density <- NetworkToolbox::conn(zero.order.network2)$total
network2.density - network1.density

# take a n = 50 sample from each SCM
n <- 50
set.seed(123)
data1 <- as.data.frame(sim_from_dag(dag = network1, n_sim = n))
data2 <- as.data.frame(sim_from_dag(dag = network2, n_sim = n))

# check data doesn't look pathological/unrealistic
par(mfrow = c(1,1))
hist(data1$A)
hist(data1$B)
hist(data1$C)
hist(data1$D)
hist(data1$E)
hist(data1$F)
hist(data1$G)
hist(data1$H)
hist(data1$I)
hist(data1$J)
hist(data1$K)


hist(data2$A)
hist(data2$B)
hist(data2$C)
hist(data2$D)
hist(data2$E)
hist(data2$F)
hist(data2$G)
hist(data2$H)
hist(data2$I)
hist(data2$J)
hist(data2$K)


# estimate fuzzy zero order dcor networks on the sample
est.net1 <- zero.order.network(data = data1, 
                               method="dcor", 
                               select = "sig", 
                               threshold = 0.05, 
                               dcor.permutations = 10000, 
                               adjust = "none")

est.net2 <- zero.order.network(data = data2, 
                               method="dcor", 
                               select = "sig", 
                               threshold = 0.05, 
                               dcor.permutations = 10000, 
                               adjust = "none")

# convert to fuzzy networks
fuzzy.net1 <- estimate.fuzzy.network(est.net1$p.values, prior = 0.5)
fuzzy.net2 <- estimate.fuzzy.network(est.net2$p.values, prior = 0.5)

# compute density distributions of each
density1 <- fuzzy.density.distribution(fuzzy.net1$fuzzy.network)
density2 <- fuzzy.density.distribution(fuzzy.net2$fuzzy.network)

# JSD of the density distributions
philentropy::JSD(rbind(density1$density.pmf$probability, density2$density.pmf$probability))

# compute the delta distribution
diff.dist <- fuzzy.difference.distribution(prob.x = density2$density.pmf$probability, 
                                           prob.y = density1$density.pmf$probability)


# plot layout
layout.mat <- matrix(c(1,2,
                       3,4,
                       5,6,
                       7,7), nrow = 4, ncol = 2, byrow = TRUE)


# plot
pdf("Figure 3.3.pdf", width=8.5, height=14)

layout(layout.mat)
qgraph(dag1, layout = 'circle', theme = 'colorblind', colFactor = 2, asize = 10, label.cex = 1.5)
title(main = "A) Causal DAG for network 1", cex.main = 1.5, adj = 0, font.main = 2, line = 2)

qgraph(dag2, layout = 'circle', theme = 'colorblind', colFactor = 2, asize = 10, label.cex = 2)
title(main = "B) Causal DAG for network 2", cex.main = 1.5, adj = 0, font.main = 2, line = 2)

qgraph(zero.order.network1, layout = 'circle', theme = 'colorblind', label.cex = 2)
title(main = "C) Dag implied zero-order network 1", cex.main = 1.5, adj = 0, font.main = 2, line = 2)

qgraph(zero.order.network2, layout = 'circle', theme = 'colorblind',label.cex = 2)
title(main = "D) Dag implied zero-order network 2", cex.main = 1.5, adj = 0, font.main = 2, line = 2)

qgraph(fuzzy.net1$fuzzy.network, layout = 'circle', theme = 'colorblind', colFactor = 2, edge.labels = TRUE, label.cex = 2)
title(main = "E) Estimated fuzzy network 1", cex.main = 1.5, adj = 0, font.main = 2, line = 2)

qgraph(fuzzy.net2$fuzzy.network, layout = 'circle', theme = 'colorblind', colFactor = 2, edge.labels = TRUE, label.cex = 2)
title(main = "F) Estimated fuzzy network 2", cex.main = 1.5, adj = 0, font.main = 2, line = 2)

plot(diff.dist$delta, diff.dist$prob, ylim = c(0,0.15), pch = 19, cex = 0.5, xlab = "Density difference", ylab = "Probability", cex.axis = 1.5, cex.lab = 1.5,)
title(main = "G) Density delta distribution", cex.main = 1.5, adj = 0, font.main = 2, line = 2)
lines(diff.dist$delta, diff.dist$prob)
abline(v = 0, lty = 2)

dev.off()

# for fuzzy differential network we need to the p values for the edge weight differences
# increase permutations for more accurate P values
nca <- compare.zero.order.networks(data1 = data1,
                                   data2 = data2, 
                                   network1 = est.net1, 
                                   network2 = est.net2, 
                                   permutations = 10000, 
                                   adjust = "none")

# estimate fuzzy differential
fuzzy.differential <- estimate.fuzzy.network(nca$differential.p.values, prior = 0.25)

# get degree distributions
differential.degrees <- fuzzy.degree.distributions(fuzzy.differential$fuzzy.network)

differential.degrees

# colours for plotting
cols <- RColorBrewer::brewer.pal(name = "Set3", n = 12)


pdf("Figure 3.4.pdf", width=11, height=5)

par(mfrow = c(1,2))
qgraph(fuzzy.differential$fuzzy.network, 
       layout = 'circle', 
       theme = 'colorblind', 
       colFactor = 2, 
       edge.labels = TRUE,
       label.cex = 2.5,
       vsize=4)
title(main = "A) Fuzzy differential network", font.main = 2, adj = 0, cex.main=1.25, line = 2)

plot(0:10, differential.degrees$degree.pmf[1,-1], type = "l", lwd = 2,
     col = cols[1],
     xlab = "Node Degree", 
     ylab = "Probability",
     xlim = c(0,10),
     ylim = c(0,0.35))
title(main = "B) Differential node degrees", font.main = 2, adj = 0, cex.main = 1.25, line = 2)
lines(0:10, differential.degrees$degree.pmf[2,-1], lwd = 2, col = cols[3])
lines(0:10, differential.degrees$degree.pmf[3,-1], lwd = 2, col = cols[4])
lines(0:10, differential.degrees$degree.pmf[4,-1], lwd = 2, col = cols[5])
lines(0:10, differential.degrees$degree.pmf[5,-1], lwd = 2, col = cols[6])
lines(0:10, differential.degrees$degree.pmf[6,-1], lwd = 2, col = cols[7])
lines(0:10, differential.degrees$degree.pmf[7,-1], lwd = 2, col = cols[8])
lines(0:10, differential.degrees$degree.pmf[8,-1], lwd = 2, col = cols[9])
lines(0:10, differential.degrees$degree.pmf[9,-1], lwd = 2, col = cols[10])
lines(0:10, differential.degrees$degree.pmf[10,-1], lwd = 2, col = cols[11])
lines(0:10, differential.degrees$degree.pmf[11,-1], lwd = 2, col = cols[12])

legend("topright", 
       legend = c("Node A",
                  "Node B",
                  "Node C",
                  "Node D",
                  "Node E",
                  "Node F",
                  "Node G",
                  "Node H",
                  "Node I",
                  "Node J",
                  "Node K"), 
       col = cols[-2],
       lty = 1,
       lwd = 4,
       cex = 0.8)

dev.off()

# P(degree >=6) 
sum(differential.degrees$degree.pmf[8,7:11]) # node H
sum(differential.degrees$degree.pmf[11,7:11]) # node K
sum(differential.degrees$degree.pmf[1,7:11]) # node A


#########
# END.  #
#########





