# uncomment line below if you need to install devtools package 
# install.packages(“remotes”)
# then install welfareNet from github
remotes::install_github("TomRow90/welfareNet")

# load required packages
# qgraph, igraph, dagitty should install as dependencies in welfareNet
# if not they can be installed from CRAN
library(welfareNet)
library(qgraph)
library(igraph)
library(dagitty)
library(tidyr)
library(ggplot2)

# sample size
n <- 30

# set seed for reproducibility
set.seed(123)

# structural causal model for mixed data for our toy model example
A <- rnorm(n, mean = 0, sd = 1)
B <- rnorm(n, mean = 0, sd = 1)
C <- 0.8*A + rnorm(n, mean = 0, sd = 1)
D <- 0.8*C + rnorm(n, mean = 0, sd = 1)
E <- 0.8*D + 0.8*C + rnorm(n, mean = 0, sd = 1)
F <- 0.8*A + 0.8*B + 0.8*C + rnorm(n, mean = 0, sd = 1)

# combined into dataframe
data <- as.data.frame(cbind(A,B,C,D,E,F))

# this corresponds to the following DAG
dag <- matrix(c(0,0,1,0,0,1,
                0,0,0,0,0,1,
                0,0,0,1,1,1,
                0,0,0,0,1,0,
                0,0,0,0,0,0,
                0,0,0,0,0,0), nrow = 6, ncol = 6, byrow = TRUE)
node.names <- c("A","B","C","D","E","F")
dimnames(dag) <- list(node.names,node.names)


# get the implied zero order dependence structure using a welfareNet function
zero.order.structure <- get.zero.order.structure(dag = dag, 
                                                 node.names = node.names)

# view dag and implied zero order structure
par(mfrow = c(1,2))
qgraph(dag, layout = 'circle', asize =5)
qgraph(zero.order.structure, layout = 'circle')

# estimate distance correlation network using welfareNet
est.net <- zero.order.network(data = data, 
                              method = "dcor", 
                              select = "sig", 
                              dcor.permutations = 10000, 
                              threshold = 0.05, 
                              adjust = "none")

# Returned networks from zero.order.network function
par(mfrow = c(1,3))
qgraph(est.net$structure, 
       layout='circle', 
       title = "A) Structure")

qgraph(est.net$selected.net, 
       theme = 'colorblind', 
       layout = 'circle', 
       edge.labels = TRUE,
       edge.label.cex = 2, 
       title= "B) Weighted selected network")

qgraph(est.net$saturated.net, 
       theme = 'colorblind', 
       layout = 'circle', 
       edge.labels = TRUE,
       edge.label.cex = 2, 
       title = "C) Weighted saturated network")

# unweighted density
unweighted.density(est.net$structure)

# For degree, we use welfareNet::degree because there is a conflict with igraph
welfareNet::degree(est.net$structure)


# Bayesian fuzzy network with 0.5 prior on all edges
est.fuzzy <- estimate.fuzzy.network(est.net$p.values, prior = 0.5)

# visualise fuzzy
par(mfrow = c(1,1))
qgraph(est.fuzzy$fuzzy.network, 
       theme = 'colorblind',
       layout = 'circle',
       edge.labels = TRUE,
       edge.label.cex = 2)

# get distribution for fuzzy estimate of density
fuzzy.density <- fuzzy.density.distribution(est.fuzzy$fuzzy.network)

# get summary (note we supply the fuzzy network not the fuzzy.density object)
density.summaries <- fuzzy.summarise.density(est.fuzzy$fuzzy.network)
density.summaries

# get prior implied distribution of density and summary
prior.fuzzy.density <- fuzzy.density.distribution(est.fuzzy$priors)
prior.density.summaries <- fuzzy.summarise.density(est.fuzzy$priors)

# number of possible edges
nodes <- 6
edges <- nodes * (nodes-1) / 2

# example plotting using base R graphics
# first add posterior
plot(fuzzy.density$density.pmf$density, 
     fuzzy.density$density.pmf$probability, 
     pch = 19, 
     xlab = "Density", 
     ylab = "Probability")
lines(fuzzy.density$density.pmf$density, 
      fuzzy.density$density.pmf$probability, 
      col= "black",
      lty = 1)
abline(v = density.summaries$expected.value, lty = 1) # expected value
# now add prior to plot
points(prior.fuzzy.density$density.pmf$density, 
       prior.fuzzy.density$density.pmf$probability,
       pch = 21)
lines(prior.fuzzy.density$density.pmf$density, 
      prior.fuzzy.density$density.pmf$probability,
      col= "black",lty = 2)
abline(v = prior.density.summaries$expected.value, lty = 2) # expected value

# P(density >= 9)
# idx will be the value of interest + 1 because density can be 0
idx <- which(fuzzy.density$density.pmf$density == 9) 
sum(fuzzy.density$density.pmf$probability[idx:nrow(fuzzy.density$density.pmf)])

# get degree distributions for fuzzy network
fuzzy.degree <- fuzzy.degree.distributions(network = est.fuzzy$fuzzy.network)
fuzzy.degree

# convert to long format for easy plotting with ggplot
degree.pmf.data <- fuzzy.degree$degree.pmf %>%
                   pivot_longer(cols = !node,
                               names_to = "degree",
                               values_to = "probability")
degree.pmf.data$node <- as.factor(degree.pmf.data$node)

#plot
ggplot(degree.pmf.data, aes(x = degree, y = probability, group = node, col= node)) +
  geom_point() +
  geom_line() +
  theme_minimal()

# P(degree 4 or greater)
idx <- which(colnames(fuzzy.degree$degree.pmf) == "4")
sum(fuzzy.degree$degree.pmf[1, idx:ncol(fuzzy.degree$degree.pmf)]) # node A: 0.95
sum(fuzzy.degree$degree.pmf[2, idx:ncol(fuzzy.degree$degree.pmf)]) # node B: 0.26
sum(fuzzy.degree$degree.pmf[3, idx:ncol(fuzzy.degree$degree.pmf)]) # node C: 0.996
sum(fuzzy.degree$degree.pmf[4, idx:ncol(fuzzy.degree$degree.pmf)]) # node D: 0.93
sum(fuzzy.degree$degree.pmf[5, idx:ncol(fuzzy.degree$degree.pmf)]) # node E: 0.98
sum(fuzzy.degree$degree.pmf[6, idx:ncol(fuzzy.degree$degree.pmf)]) # node F: 0.98


# ------------- network comparisons ------------------#



# structural causal model for network 2
A2 <- rnorm(n, mean = 0, sd = 1)
B2 <- rnorm(n, mean = 0, sd = 1)
C2 <- 0.8*A2 + rnorm(n, mean = 0, sd = 1)
D2 <- rnorm(n, mean = 0, sd = 1)
E2 <- 0.8*B2 + rnorm(n, mean = 0, sd = 1)
F2 <- 0.8*D2 + rnorm(n, mean = 0, sd = 1)

# combined into dataframe
data2 <- as.data.frame(cbind(A2,B2,C2,D2,E2,F2))
colnames(data2) <- c("A","B","C","D","E","F")

# this corresponds to the following DAG
dag2 <- matrix(c(0,0,1,0,0,0,
                 0,0,0,0,1,0,
                 0,0,0,0,0,0,
                 0,0,0,0,0,1,
                 0,0,0,0,0,0,
                 0,0,0,0,0,0), nrow = 6, ncol = 6, byrow = TRUE)
dimnames(dag2) <- list(node.names,node.names)

# get the implied zero order dependence structure using a welfareNet function
zero.order.structure2 <- get.zero.order.structure(dag = dag2, node.names = node.names)

# view dag and zero order structure
par(mfrow = c(1,2))
qgraph(dag2, layout = 'circle', asize = 5)
qgraph(zero.order.structure2, layout = 'circle')

# estimate distance correlation network from data2
est.net2 <- zero.order.network(data = data2, 
                              method = "dcor", 
                              select = "sig", 
                              dcor.permutations = 10000, 
                              threshold = 0.05, 
                              adjust = "none")



# perform network comparison analysis permutation test
nca <- compare.zero.order.networks(data1 = data,
                                   data2 = data2,
                                   network1 = est.net,
                                   network2 = est.net2,
                                   permutations = 10000,
                                   adjust = "none")

# view results of permutation test at different levels of network description
nca$dissimilarity.results
nca$global.results
nca$strength.results
nca$edge.differences

# construct differential network showing significant edge weight differences
frequentist.differential <- differential.network(net = nca,
                                                  select = "sig",
                                                  threshold = 0.05,
                                                  adjust = FALSE)
# plot differential network
par(mfrow = c(1,1))
qgraph(frequentist.differential$selected.net, 
       theme = 'colorblind', 
       layout = 'circle', 
       edge.labels =TRUE, 
       edge.label.cex =2)


# ---------- fuzzy network comparisons ------------- #

# fuzzy network for network 2
est.fuzzy2 <- estimate.fuzzy.network(est.net2$p.values, prior = 0.5)

# density distribution for network 2
fuzzy.density2 <- fuzzy.density.distribution(est.fuzzy2$fuzzy.network)

# delta distributions
density.delta.dist <- fuzzy.difference.distribution(prob.x = fuzzy.density$density.pmf$probability,
                                                    prob.y = fuzzy.density2$density.pmf$probability)

prior.density.delta.dist <- fuzzy.difference.distribution(prob.x = prior.fuzzy.density$density.pmf$probability,
                                                          prob.y = prior.fuzzy.density$density.pmf$probability)


# plot density delta distribution
ggplot(density.delta.dist, aes(x = delta, y = prob)) +
  geom_point() +
  geom_line() +
  geom_point(data = prior.density.delta.dist, shape = 21) +
  geom_line(data = prior.density.delta.dist, linetype = "dotted") +
  theme_minimal()

# Probabilities
density.probs <- fuzzy.delta.probabilities(delta = density.delta.dist)
density.probs

# vector of priors to use for prior sensitivity
priors <- seq(0.1,0.5,0.05)

# perform density delta prior sensitivity
delta.density.prior.sensitivity <- fuzzy.density.prior.sensitivity(
  p.value.matrix1 = est.net$p.values,
  p.value.matrix2 = est.net2$p.values,
  priors = priors,
  epsilon = 0)

# view 
delta.density.prior.sensitivity

# plot delta density prior sensitivity 
ggplot(data = delta.density.prior.sensitivity, aes(x = prior, y = `P(delta > epsilon)`)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  theme_minimal()

# get degree distributions for fuzzy network 2
fuzzy.degree2 <- fuzzy.degree.distributions(network = est.fuzzy2$fuzzy.network)
fuzzy.degree2

# get distributions for node D
node.D.degree1 <- as.numeric(fuzzy.degree$degree.pmf[4,-1])
node.D.degree2 <- as.numeric(fuzzy.degree2$degree.pmf[4,-1])

# compute delta distributions of difference
node.D.degree.delta <- fuzzy.difference.distribution(prob.x = node.D.degree1,
                                                     prob.y = node.D.degree2)

node.D.degree.delta

# plot density delta distribution
ggplot(node.D.degree.delta, aes(x = delta, y = prob)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(limits = c(-5,5), breaks = seq(-5,5,1)) +
  theme_minimal()

# delta node D degree probabilities
node.D.degree.probs <- fuzzy.delta.probabilities(delta = node.D.degree.delta)
node.D.degree.probs
