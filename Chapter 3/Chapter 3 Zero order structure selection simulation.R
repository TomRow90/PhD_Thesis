#######################################################################
#
#        Section 3.5.3: Structure selection examples from Chapter 3
#
#######################################################################

#load required packages
library(simDAG)
library(qgraph)
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


#------------------------------------------------------------------------------#


# structural causal model for mixed data for our toy model example
# uses simDAG and the created functions for gamma and beta above
linear.mixed.dag <- empty_dag() +
  node("A", type="rbinom", size=1, prob=0.5) +
  node("B", type="rbeta", shape1=2, shape2=6) +
  node("C", type="negative_binomial", parents="A", betas=1.2, intercept=1.5, theta=1) +
  node("D", type="gaussian", parents="A", betas=20, intercept=100, error=15) +
  node("E", type="poisson", parents="B", betas=2.5, intercept=1.5) +
  node("F", type="poisson", parents="C", betas=0.02, intercept=2) +
  node("G", type="gaussian", parents=c("F","C"), betas=c(0.4,0.4), intercept=25, error=6) +
  node("H", type="gaussian", parents=c("D","E"), betas=c(2,4), intercept=3, error=1) +
  node("I", type=node_beta, parents="E", betas=0.25, intercept=-1, phi=25) +
  node("J", type=node_beta, parents="H", betas=-0.025, intercept=3, phi=40) +
  node("K", type="gaussian", parents=c("B","I"), betas=c(25,25), intercept=20, error=3)


nodes <- c("A","B","C","D","E","F","G","H","I","J","K")

# causal dag
true.dag <- dag2matrix(linear.mixed.dag)

# causal skeleton
true.skeleton <- dag2matrix(linear.mixed.dag)
true.skeleton[lower.tri(true.skeleton)] <- t(true.skeleton)[lower.tri(true.skeleton)]

# PMRF
true.pmrf <- gRbase::moralize(dag2matrix(linear.mixed.dag))

# marginal/zero-order structure
true.marginal <- get.zero.order.structure(dag = true.dag, node.names = nodes)

# not figures in chapter 3 but view here
par(mfrow=c(2,2))
g <- qgraph::qgraph(true.dag, layout='spring', esize=2, asize=5, title="A) Causal DAG")
qgraph::qgraph(true.skeleton, layout=g$layout, esize=2, title="B) Causal skeleton")
qgraph::qgraph(true.pmrf, layout=g$layout, esize=2, title="C) Pairwise Markov Random Field")
qgraph::qgraph(true.marginal, esize=2, curve = c(rep(0,4),0.7,rep(0,50)), title="D) Marginal dependence structure")


# simulate a data set to produce Figure 3.2, and set seed so reproducible
set.seed(123)
data <- as.data.frame(sim_from_dag(linear.mixed.dag, n = 50))

# estimate dcor network
dcor.net <- zero.order.network(data = data,
                               method = "dcor", 
                               select = "sig", 
                               dcor.permutations = 100000,
                               threshold = 0.05, 
                               adjust = "none", 
                               npn = FALSE)

# function to apply different P value thresholds given p value matrix from dcor.net
threshold.function <- function(p, threshold) {
  
  adj <- matrix(0, nrow = ncol(p), ncol = ncol(p))
  idx <- which(p < threshold)
  adj[idx] <- 1
  diag(adj) <- 0
  
  return(adj)
}

# apply different P value thresholds
high.threshold <- threshold.function(p = dcor.net$p.values, threshold = 0.01)
dimnames(high.threshold) <- list(nodes,nodes)
mod.threshold <- threshold.function(p = dcor.net$p.values, threshold = 0.05)
dimnames(mod.threshold) <- list(nodes,nodes)
low.threshold <- threshold.function(p = dcor.net$p.values, threshold = 0.10)
dimnames(low.threshold) <- list(nodes,nodes)

# estimate fuzzy network using welfareNet
dcor.fuzzy <- estimate.fuzzy.network(dcor.net$p.values, prior = 0.5)

# get estimates of degree of node H for the example
fuzzy.degree <- fuzzy.degree.distributions(dcor.fuzzy$fuzzy.network)
sum(fuzzy.degree$degree.pmf[8,8:11]) #P(degree >=6)
sum(dcor.fuzzy$fuzzy.network[8,-8]) # expected degree for node H

degree(high.threshold)["H"]
degree(mod.threshold)["H"]
degree(low.threshold)["H"]


# layout for plotting
layout.mat <- matrix(c(0,1,2,0,
                       3,4,5,6,
                       7,8,9,10), nrow = 3, ncol = 4, byrow = TRUE)

# get layout based on this graph
g <- qgraph::qgraph(true.marginal, esize=2, curve = c(rep(0,4),0.7,rep(0,50)), label.cex = 2)


pdf("Figure 3.2.pdf", width=18, height=10)

layout(layout.mat)

qgraph(true.dag, layout=g$layout, esize=2, asize=10, label.cex = 2)
title(main = "A) True causal DAG", cex.main = 1.5, adj = 0, font.main = 2, line = 3) 

g <- qgraph::qgraph(true.marginal, esize=2, curve = c(rep(0,4),0.7,rep(0,50)), label.cex = 2)
title(main = "B) True zero-order dependence structure", cex.main = 1.5, adj = 0, font.main = 2, line = 3) 

qgraph(high.threshold, layout = g$layout, theme = "colorblind", maximum = 1, label.cex = 2)
title(main = "C) High threshold (P < 0.01)", cex.main = 1.5, adj = 0, font.main = 2, line = 3) 

qgraph(mod.threshold, layout = g$layout, theme = "colorblind", maximum = 1, label.cex = 2)
title(main = "D) Moderate threshold (P < 0.05)", cex.main = 1.5, adj = 0, font.main = 2, line = 3) 

qgraph(low.threshold, layout = g$layout, theme = "colorblind", maximum = 1, label.cex = 2)
title(main = "E) Low threshold (P < 0.10)", cex.main = 1.5, adj = 0, font.main = 2, line = 3) 

qgraph(dcor.fuzzy$fuzzy.network, layout=g$layout, theme = "colorblind", maximum = 1, colFactor = 2, label.cex = 2)
title(main = "F) Fuzzy network", cex.main = 1.5, adj = 0, font.main = 2, line = 3) 

plot.new()
text(x = 0.5, y = 1, labels = "Node H degree = 3", cex = 1.5)

plot.new()
text(x = 0.5, y = 1, labels = "Node H degree = 7", cex = 1.5)

plot.new()
text(x = 0.5, y = 1, labels = "Node H degree = 8", cex = 1.5)

plot(0:10, fuzzy.degree$degree.pmf[8,-1], ylab = "Probability", xlab = "Degree of node H", cex.lab = 1.5, cex.axis = 1.2)
lines(0:10, fuzzy.degree$degree.pmf[8,-1])
abline(v = sum(dcor.fuzzy$fuzzy.network[8,-8]), lty = 2)
title(main = "G) Fuzzy degree distribution node H", cex.main = 1.5, adj = 0, font.main = 2, line = 3) 

dev.off()


############
#.  END    #
############

