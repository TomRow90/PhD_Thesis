################################################################################
#            R script for re-producing simulations in Chapter 1                #
################################################################################

#load required packages
library('corpcor')
library('qgraph')
library('bootnet')
library('GGMnonreg')
library('IsingSampler')
library('NetworkComparisonTest')

#------------------------------------------------------------------------------#
#                             SECTION 1.3: FIGURE 1.1
#------------------------------------------------------------------------------#


# network 1
unweighted_undirected <- matrix(c(0,1,1,0,
                                  1,0,0,1,
                                  1,0,0,1,
                                  0,1,1,0), nrow=4, ncol = 4, byrow = TRUE)
dimnames(unweighted_undirected) <- list(c("A","B","C","D"),
                                        c("A","B","C","D"))


g1 <- qgraph(unweighted_undirected, layout="circle", esize = 2)


# network 2
weighted_undirected <- matrix(c(0,0.25,0.15,0,
                                0.25,0,0,0.3,
                                0.15,0,0,0.1,
                                0,0.3,0.1,0), nrow=4, ncol = 4, byrow = TRUE)
dimnames(weighted_undirected) <- list(c("A","B","C","D"),
                                      c("A","B","C","D"))


g2 <- qgraph(weighted_undirected, layout="circle", theme = 'colorblind')


# network 3
unweighted_directed <- matrix(c(0,1,1,0,
                                0,0,0,1,
                                0,0,0,1,
                                0,0,0,0), nrow=4, ncol = 4, byrow = TRUE)
dimnames(unweighted_directed) <- list(c("A","B","C","D"),
                                      c("A","B","C","D"))


g3 <- qgraph(unweighted_directed, layout="circle", esize = 2, asize = 9)


# network 4
weighted_directed <- matrix(c(0,0.25,0.15,0,
                              0,0,0,0.3,
                              0,0,0,0.1,
                              0,0,0,0), nrow=4, ncol = 4, byrow = TRUE)
dimnames(weighted_directed) <- list(c("A","B","C","D"),
                                    c("A","B","C","D"))


g4 <- qgraph(weighted_directed, layout="circle", theme = "colorblind", asize = 9)


pdf("Figure 1.1.pdf", width=10, height=8)

par(mfrow=c(2,2))

qgraph(g1, title = "A) Undirected unweighted")
qgraph(g2, title = "B) Undirected weighted")
qgraph(g3, title = "C) Directed unweighted")
qgraph(g4, title = "D) Directed weighted")

dev.off()


#--------------------------------------------------------------#
#          SECTION 1.4: PLOT OF NETWORK THEORY MODEL
#--------------------------------------------------------------#

# network adjacency matrix
example_network <- matrix(c(0,1,0,0,0,0,0,0,
                            1,0,1,1,0,0,0,0,
                            0,1,0,1,0,0,1,0,
                            0,1,1,0,1,1,0,1,
                            0,0,0,1,0,1,0,0,
                            0,0,0,1,1,0,0,0,
                            0,0,1,0,0,0,0,0,
                            0,0,0,1,0,0,0,0), nrow=8, ncol = 8, byrow = TRUE)
dimnames(example_network) <- list(c("I1","I2","I3","I4","I5","I6","E1","E2"),
                                  c("I1","I2","I3","I4","I5","I6", "E1","E2"))

# plotting positions
node_positions <- matrix(c(0,0.9,
                           0.2,0.6,
                          -0.2,0.4,
                           0,0.1,
                           0.2,-0.1,
                          -0.1,-0.3,
                          -0.1,-0.8,
                           0.1,-0.8), nrow = 8, ncol = 2, byrow=TRUE)

# node shapes
shapes <- c(rep("circle",6), "square","square")

par(mfrow = c(1,1))

pdf("Figure 1.2.pdf", width=5, height=5)

qgraph(example_network,
       layout = node_positions,
       shape = shapes,
       edge.color = "black",
       labels = colnames(example_network),
       label.scale = FALSE)

# get plotting limits
usr <- par("usr")

# shade lower region
rect(usr[1], usr[3], usr[2], -0.6,
     col = adjustcolor("grey", alpha.f = 0.3),
     border = NA)

# shade upper region
rect(-1.3, -0.6, 1.3, 1.3,
     col = adjustcolor("grey90", alpha.f = 0.3),
     border = NA)

abline(h = -0.6, lty = 2)
text(x = -0.9, y = 1, labels = "INTERNAL FIELD", cex = 0.7, font = 2)
text(x = -0.9, y = -1.2, labels = "EXTERNAL FIELD", cex = 0.7, font = 2)

dev.off()

#------------------------------------------------------------------------------#
#-------------- SECTION 1.5: STATISTICAL ANIMAL WELFARE NETWORKS --------------#
#------------------------------------------------------------------------------#

#set seed to make below reproducible
set.seed(2)

#Generate data for the 3 node network 
NOV <- rnorm(500)
COR <- NOV + rnorm(500, 0, 2)
COG <- -NOV + rnorm(500, 0, 1)

#create data matrix
mat <- cbind(NOV,COR,COG)

#correlation matrix
corMat <- cor(mat)
corMat

#partial correlation matrix
parMat <- cor2pcor(corMat)
dimnames(parMat) <- list(c("NOV", "COR", "COG"),c("NOV", "COR", "COG"))
parMat

#Figure 2A & 2B

tiff("Figure 2.tiff", width=7.48, height=3.15, units='in', res=1200)

par(mfrow=c(1,2), oma=c(0,0,1,0))

qgraph(corMat, layout=matrix(c(2,2.5,1.25,1,2.75,1), nrow=3, ncol=2, byrow=T), node.width=1.1, node.height=1.1, 
       edge.labels=T, labels=c("NOV", "COR", "COG"), label.cex=1, edge.label.cex=1.2, theme="colorblind")
title(main="A: Correlation Structure", line=3.6, adj=0, cex.main=0.8)

qgraph(parMat, layout=matrix(c(2,2.5,1.25,1,2.75,1), nrow=3, ncol=2, byrow=T), node.width=1.1, node.height=1.1, 
       edge.labels=T, labels=c("NOV", "COR", "COG"), label.cex=1, edge.label.cex=1.2, theme="colorblind")
title(main="B: Partial Correlation Structure", line=3.6, adj=0, cex.main=0.8)

dev.off()

#------------------------------------------------------------------------------#
#--- SECTION 1.6: NETWORK MODELS AS A TOOL TO IDENTIFY INTERVENTION TARGETS ---#
#------------------------------------------------------------------------------#

#set sample size for simulated data
n2 <- 5000

#define network structure simulated data should re-produce: edge weights derived from Cramer et al (2016) "major depression as complex dynamic system"
#data obtained from: https://figshare.com/projects/Major_depression_as_a_complex_dynamic_system_accepted_for_publication_in_PLoS_ONE_/17360
network <- matrix(data=c(0.000,	  2.1407,	0.7232,	0.2041,	1.1296,	0.5217,	1.0530,	0.9409,	0.7484,	0.6849,	1.0979,	1.8733,	1.0211,	2.0693,
                         2.1407,	0.0000,	0.1766,	0.2811,	0.5763,	0.2392,	0.4273,	0.5311,	0.4459,	0.6564,	0.5070,	0.6826,	0.8178,	0.4986,
                         0.7232,	0.1766,	0.0000,-0.6082,	2.984,	0.0000,	0.2045,	0.0000,	0.0000,	0.0000,	0.0000,	0.3772,	0.1063,	0.2520,
                         0.2041,	0.2811,-0.6082,	0.0000,-0.5389,	3.1650,	0.2672,	0.2041,	0.0000,	0.0000,	0.4112,	0.5226,	0.0000,	0.0000,
                         1.1296,	0.5763,	2.9840,-0.5389,	0.0000,	0.0000,	0.7033,	0.4020,	0.4724,	0.2219,	0.2284,	0.1203,	0.4177,	0.1198,
                         0.5217,	0.2392,	0.0000,	3.1650,	0.0000,	0.0000,	0.0000,	0.5475,	0.4890,	0.2914,	0.4546,	0.0000,	0.0000,	0.1620,
                         1.0530,	0.4273,	0.2045,	0.2672,	0.7033,	0.0000,	0.0000,-0.5009,	1.2951,	0.0000,	0.8279,	0.0000,	0.2585,	0.4514,
                         0.9409,	0.5311,	0.0000,	0.2041,	0.4020,	0.5475,-0.5009,	0.0000,	0.0000,	0.4048,	1.4768,	0.2708,	0.0597,	0.2151,
                         0.7484,	0.4459,	0.0000,	0.0000,	0.4724,	0.4890,	1.2951,	0.0000,	0.0000,	0.0000,	0.3751,	0.3893,	0.9414,	0.1939,
                         0.6849,	0.6564,	0.0000,	0.0000,	0.2219,	0.2914,	0.0000,	0.4048,	0.0000,	0.0000,	1.5718,	0.3491,	0.7233,	0.1407,
                         1.0979,	0.5070,	0.0000,	0.4112,	0.2284,	0.4546,	0.8279,	1.4768,	0.3751,	1.5718,	0.0000,	0.2362,	0.4935,	0.0000,
                         1.8733,	0.6826,	0.3772,	0.5226,	0.1203,	0.0000,	0.0000,	0.2708,	0.3893,	0.3491,	0.2362,	0.0000,	0.6660,	1.4769,
                         1.0211,	0.8178,	0.1063,	0.0000,	0.4177,	0.0000,	0.2585,	0.0597,	0.9414,	0.7233,	0.4935,	0.6660,	0.0000,	0.2156,
                         2.0693,	0.4986,	0.2520,	0.0000,	0.1198,	0.1620,	0.4514,	0.2151,	0.1939,	0.1407,	0.0000,	1.4769,	0.2156,	0.0000),
                  nrow=14, ncol=14, byrow=T, dimnames=list(c("NOV","uAGG","BARK","NOISE","FRU","LAME","ARB","COP","ACT","BCS","APP","uFEAR","FEC","aFEAR"),c("NOV","uAGG","BARK","NOISE","FRU","LAME","ARB","COP","ACT","BCS","APP","uFEAR","FEC","aFEAR")))

#check graph is symmetrical with no errors in input
isSymmetric(network)

#threshold parameters also derived from Cramer et al (2016)
thresholds <- c(-2.3129,-3.1946,-4.3092,-3.8332,-3.9153,-3.9012,-3.0246,-4.448,-3.1753,-4.3372,-2.8269,-4.4272,-4.0421,-5.8303)

#temperature parameter
beta <- 1

#simulate data to reproduce network/graph
set.seed(123)
data <- IsingSampler(n=n2, graph=network, thresholds=thresholds, beta=beta, nIter=100, responses=c(0L, 1L), method="MH")
colnames(data) <- c("NOV","uAGG","BARK","NOISE","FRU","LAME","ARB","COP","ACT","BCS","APP","uFEAR","FEC","aFEAR")

#fit unconstrained network model (Ising model) to simulated data
fit <- estimateNetwork(data=data, default="IsingSampler", method='uni')
dimnames(fit$graph) <- list(c("NOV","uAGG","BARK","NOISE","FRU","LAME","ARB","COP","ACT","BCS","APP","uFEAR","FEC","aFEAR"),c("NOV","uAGG","BARK","NOISE","FRU","LAME","ARB","COP","ACT","BCS","APP","uFEAR","FEC","aFEAR"))

#not in main article but can view centrality plot
centralityPlot(fit$graph)
#store strength centrality in vector to be used later for altering node size based on centrality
centrality <- centralityTable(fit$graph)[29:42,5]

#these lines run for loop to give each node a size as a function of its strength centrality 
nodeSize <- rep(0, ncol(data))

for (i in 1:ncol(data)) {
  if(centrality[i] <= 0) {
    nodeSize[i] <- (8*exp(-ncol(data)/80) + 2*centrality[i]) 
  }
  else {
    nodeSize[i] <- (8*exp(-ncol(data)/80) + 2*centrality[i])
  }
}

#colour matrix to use as input to qgraph  
colMatrix2.1 <- matrix("gray", nrow=nrow(network), ncol=ncol(network))

#find which indices are >= 0.5, and which are <=-0.5 for edge colouring 
posIndices <- which(fit$graph >= 0.5, arr.ind=T)
negIndices <- which(fit$graph <= -0.5, arr.ind=T)

#colour edges >= 0.5 blue, and <= -0.5 red, leave the rest gray
colMatrix2.1[posIndices[,1:2]] <- "blue"
colMatrix2.1[negIndices[,1:2]] <- "red"

#plot network
tiff("Figure 1.4.tiff", width=7.48, height=4.73, units='in', res=1200)

par(mfrow=c(1,1), mar=c(0,0,0,0))

qgraph(fit$graph, 
       layout='spring',
       vsize=nodeSize, 
       theme="colorblind", 
       edge.color=colMatrix2.1)

dev.off()

#---------------centrality simulation using Cramer et al (2016) dynamical model------------------#

#set number of time steps
time <- 10000

#create matrix for data to be stored in
nodes <- matrix (0, nrow=14, ncol=time)
#set all nodes to be inactive (0) at start of simulation
nodes[,1] <- rep(0,14)

#start simulation loop 
for(t in 2:time) {
  if(t <= 5000) {
    for(n in 1:nrow(nodes)) {
      
      A <- sum(1.6*(fit$graph[n,] %*% nodes[,t-1])) #this is the total activitation function i.e. total amount of activation symptom n receives at time t
      b <- abs(fit$intercepts[n]) #this is the threshold of symptom n 
      p <- 1/(1+exp(b - A)) #this is the probability function for symptom n being active i.e. taking a value of 1
      nodes[n,t] <- rbinom(n=1, size=1, prob=p) #this uses p to assign a value of 0 (inactive) or 1 (active) to the nodes matrix
    }
  }
  else { # halfway through time steps, intervene on the system
    for(n in 1:nrow(nodes)) {
      
      A <- sum(1.6*(fit$graph[n,] %*% nodes[,t-1])) 
      b <- abs(fit$intercepts[n]) 
      p <- 1/(1+exp(b - A)) 
      nodes[n,t] <- rbinom(n=1, size=1, prob=p)
      nodes[1,t] <- 0 #this is the simulated intervention on highly central NOV node i.e. set fear of novelty to always be 0
    }
  }
}

# now repeat simulation but this time constraining low strength centrality node ACT to 0 

nodes2 <- matrix (0, nrow=14, ncol=time)
nodes2[,1] <- rep(0,14)

for(t in 2:time) {
  if(t <= 5000) {
    for(n in 1:nrow(nodes2)) {
      
      A <- sum(1.6*(fit$graph[n,] %*% nodes2[,t-1])) 
      b <- abs(fit$intercepts[n]) 
      p <- 1/(1+exp(b - A)) 
      nodes2[n,t] <- rbinom(n=1, size=1, prob=p) 
    }
  }
  else {
    for(n in 1:nrow(nodes2)) {
      
      A <- sum(1.6*(fit$graph[n,] %*% nodes2[,t-1])) 
      b <- abs(fit$intercepts[n]) 
      p <- 1/(1+exp(b - A)) 
      nodes2[n,t] <- rbinom(n=1, size=1, prob=p) 
      nodes2[9,t] <- 0 #intervene on low centrality node ACT by setting to always be 0
    }
  }
}

#calculate total activity of system at each time i.e. number of nodes active at each time point
activity <- colSums(nodes)
activity2 <- colSums(nodes2)

#these plots are not in the article. They show the number of nodes active at each time point and how this changes after the system is intervened on 
par(mfrow=c(2,1), mar=c(5.1, 4.1, 4.1, 2.1))
plot(activity[1:10000], type="l", col="blue", ylim=c(0,14), ylab="No. of active nodes", xlab="Time", main="A: High centrality intervention")
abline(v=5000, lty=2, col="green")
plot(activity2[1:10000], type="l", col="red", ylim=c(0,14), ylab="No. of active nodes", xlab="Time", main="B: Low centrality intervention")
abline(v=5000, lty=2, col="green")

#calculate the median number of nodes active in each simulation after intervention
median(colSums(nodes[,5001:10000]))
median(colSums(nodes2[,5001:10000]))

#The following produces networks shown in Figure 1.5
#create function to calculate probability node is active i.e. number of time steps active / total number of time steps 
fun <- function (x) sum(x)/5000

#use above function to calcualte probabilities
hcBaselineProb <- apply(nodes[,1:5000], 1, fun)
hcInterProb <- apply(nodes[ ,5001:10000], 1, fun)
lcBaselineProb <- apply(nodes2[,1:5000], 1, fun) #this will be almost identical to hcBaselineProb so isn't used in following plots
lcInterProb <- apply(nodes2[ ,5001:10000], 1, fun)

#node shape vectors for input to qgraph plot
shapes1 <- c("square", "circle","circle","circle","circle","circle","circle","circle","circle","circle","circle","circle","circle","circle")
shapes2 <- c("circle", "circle","circle","circle","circle","circle","circle","circle","square","circle","circle","circle","circle","circle")

par(mfrow=c(1,1))
#networks used in Figure 4, note: layout in publication created in powerpoint
qgraph(fit$graph, layout='spring', vsize=nodeSize, edge.color=colMatrix2.1, label.cex=1.1, pie=hcBaselineProb, pieBorder=1, pieColor="gold",color="white", theme="colorblind")
qgraph(fit$graph, layout='spring', shape=shapes1, vsize=nodeSize, edge.color=colMatrix2.1, label.cex=1.1, pie=hcInterProb, pieBorder=1, pieColor="gold",color="white",theme="colorblind")
qgraph(fit$graph, layout='spring', shape=shapes2, vsize=nodeSize, edge.color=colMatrix2.1, label.cex=1.1, pie=lcInterProb, pieBorder=1, pieColor="gold",color="white",theme="colorblind")

#This plot is not shown in article but is useful to see. It plots the networks above each phase in the simulation shown as no. of nodes active at each time step
layout(matrix(c(2,3,1,1,4,5,6,6), nrow=4, byrow=T))
plot(activity[1:10000], type="l", col="green", ylim=c(0,14), ylab="No. of active nodes", xlab="Time")
qgraph(fit$graph, layout='spring', theme='colorblind', vsize=nodeSize, label.cex=1, pie=hcBaselineProb, pieBorder=1, pieColor="gold",color="white")
qgraph(fit$graph, layout='spring', theme='colorblind', shape=shapes1, vsize=nodeSize, label.cex=1, pie=hcInterProb, pieBorder=1, pieColor="gold",color="white")
qgraph(fit$graph, layout='spring', theme='colorblind', vsize=nodeSize, label.cex=1, pie=lcBaselineProb, pieBorder=1, pieColor="gold",color="white")
qgraph(fit$graph, layout='spring', theme='colorblind', shape=shapes2, vsize=nodeSize, label.cex=1.2, pie=lcInterProb, pieBorder=1, pieColor="gold",color="white")
plot(activity2[1:10000], type="l", col="red", ylim=c(0,14), ylab="No. of active nodes", xlab="Time")

#------------------------------------------------------------------------------#
#------- SECTION 1.7: COMPARING NETWORKS TO AID WELFARE INFERENCE--------------#
#------------------------------------------------------------------------------#


#--------- perform network analysis for group 1 ---------#

# set sample size to realistic value
n=100

#create matrix of partial correlations that we want to roughly reproduce with simulated data
network1 <- matrix(data = c(1.00, 0.00,  0.00, 0.00, 0.00, 
                            0.00, 1.00,  0.30, 0.30, 0.00,  
                            0.00, 0.30,  1.00, 0.30, 0.00,  
                            0.00, 0.30,  0.30, 1.00, 0.00, 
                            0.00, 0.00,  0.00, 0.00, 1.00) ,nrow=5,ncol=5)
dimnames(network1) <- list(c("AGG","BOLD","ACT","EXP","SOC"),
                             c("AGG","BOLD","ACT","EXP","SOC"))

#check matrix is symmetrical & no human error entering matrix above
isSymmetric(network1)

#following simulates data from multivariate normal distribution - taken from ggmGenerator source code in bootnet package
graph <- network1
intercepts <- rep(0, ncol(graph))

#standardize:
if (!all(diag(graph) == 0 | diag(graph) == 1)){
  graph <- cov2cor(graph)
}

#Remove diag:
diag(graph) <- 0

#True sigma (covariance matrix):
if (any(eigen(diag(ncol(graph)) - graph)$values < 0)){
  stop("Precision matrix is not positive semi-definite")
}

Sigma <- cov2cor(solve(diag(ncol(graph)) - graph))

# Generate data:
set.seed(1)
data1 <- mvtnorm::rmvnorm(n, sigma = Sigma)
colnames(data1) <- c("AGG","BOLD","ACT","EXP","SOC")

#estimate non-regularised partial correlation network with model selection based on fisher Z transformations
fit1 <- ggm_inference(Y=as.matrix(data1), 
                      alpha=0.05, 
                      control_precision = TRUE,
                      boot = FALSE, 
                      method = "pearson", 
                      progress = TRUE)
dimnames(fit1$wadj) <- list(c("AGG","BOLD","ACT","EXP","SOC"),
                            c("AGG","BOLD","ACT","EXP","SOC"))


# colour matrix to use as input to qgraph  
colMatrix1 <- matrix("gray", nrow=nrow(network1), ncol=ncol(network1))

# find significant node indices, as well as positive and negative partial cors 
sigInd1 <- which(fit1$wadj < 0 | fit1$wadj > 0, arr.ind=T)
posInd1 <- which(fit1$pcors >= 0, arr.ind=T)
negInd1 <- which(fit1$pcors < 0, arr.ind=T)

# run for loop to overwrite non-significant colours with different colours for sig. edges
for (i in 1:nrow(sigInd1)) {
  
  if (fit1$wadj[sigInd1[i,1], sigInd1[i,2]] > 0) {
    colMatrix1[sigInd1[i,1], sigInd1[i,2]] <- "blue"
  }
  else {
    if(fit1$wadj[sigInd1[i,1],sigInd1[i,2]] < 0) {
      colMatrix1[sigInd1[i,1], sigInd1[i,2]] <- "red"
    }
    else {
      
    }
  }
  
}

#plot
qgraph(fit1$wadj, layout='circle', edge.labels=T, edge.color=colMatrix1)

#------------now repeat for group 2 network ---------------------#

# create matrix of 'ideal' partial correlations for group 2 in our conceptual example
# for this, we make the partial correlations much lower or 0 compared to group 1
# create matrix of partial correlations that we want to roughly reproduce with simulated data
network2 <- matrix(data = c(1.00, 0.50,  0.00, 0.00,-0.42, 
                            0.50, 1.00,  0.42, 0.42, 0.00,  
                            0.00, 0.42,  1.00, 0.42, 0.00,  
                            0.00, 0.42,  0.42, 1.00, 0.00, 
                           -0.42, 0.00,  0.00, 0.00, 1.00) ,nrow=5,ncol=5)
dimnames(network2) <- list(c("AGG","BOLD","ACT","EXP","SOC"),
                           c("AGG","BOLD","ACT","EXP","SOC"))

#check matrix is symmetrical & no human error entering matrix above
isSymmetric(network2)

#following simulates data from multivariate normal distribution - taken from ggmGenerator source code in bootnet package
graph <- network2
intercepts <- rep(0, ncol(graph))

#standardize:
if (!all(diag(graph) == 0 | diag(graph) == 1)){
  graph <- cov2cor(graph)
}

#Remove diag:
diag(graph) <- 0

#True sigma (covariance matrix):
if (any(eigen(diag(ncol(graph)) - graph)$values < 0)){
  stop("Precision matrix is not positive semi-definite")
}

Sigma <- cov2cor(solve(diag(ncol(graph)) - graph))

# Generate data:
set.seed(1)
data2 <- mvtnorm::rmvnorm(n, sigma = Sigma)
colnames(data2) <- c("AGG","BOLD","ACT","EXP","SOC")

#estimate non-regularised partial correlation network with model selection based on fisher Z transformations
fit2 <- ggm_inference(Y=as.matrix(data2), 
                      alpha=0.05, 
                      control_precision = TRUE,
                      boot = FALSE, 
                      method = "pearson", 
                      progress = TRUE)
dimnames(fit2$wadj) <- list(c("AGG","BOLD","ACT","EXP","SOC"),
                            c("AGG","BOLD","ACT","EXP","SOC"))

#----the following sets edges to be coloured if significant according to fit1.1, and sets other edge weights identified in fit1.2 to gray---#

#colour matrix to use as input to qgraph  
colMatrix2 <- matrix("gray", nrow=nrow(network2), ncol=ncol(network2))

#find significant node indices, as well as positive and negative partial cors 
sigInd2 <- which(fit2$wadj < 0 | fit2$wadj > 0, arr.ind=T)
posInd2 <- which(fit2$pcors >= 0, arr.ind=T)
negInd2 <- which(fit2$pcors < 0, arr.ind=T)

#run for loop to overwrite non-significant colours with different colours for sig. edges
for (i in 1:nrow(sigInd2)) {
  
  if (fit2$wadj[sigInd2[i,1], sigInd2[i,2]] > 0) {
    colMatrix2[sigInd2[i,1], sigInd2[i,2]] <- "blue"
  }
  else {
    if(fit2$wadj[sigInd2[i,1],sigInd2[i,2]] < 0) {
      colMatrix2[sigInd2[i,1], sigInd2[i,2]] <- "red"
    }
    else {
      
    }
  }
  
}

#plot
qgraph(fit2$wadj, layout='circle', edge.labels=T, edge.color=colMatrix2)

# colour nodes with diff centralities grey
node.colors <- c(rep("grey",2), rep("white",3))

#following produces Figure 6

pdf("Figure 1.6.pdf", width=7.48, height=3.15)

par(mfrow=c(1,2), oma=c(0,0,1,0))

qgraph(fit1$pcors, 
       layout='circle', 
       labels = c("AGG","BOLD","ACT","EXP","SOC"), 
       edge.labels=T, 
       edge.color=colMatrix1, 
       color = node.colors,
       maximum=0.5) #both networks need to be scaled to same maximum to be comparable so set maximum to 0.5
title(main="A) Enriched housing group", line=3.6, adj=0, cex.main=0.8)

qgraph(fit2$pcors, 
       layout='circle', 
       labels = c("AGG","BOLD","ACT","EXP","SOC"), 
       edge.labels=T, 
       edge.color=colMatrix2, 
       color = node.colors,
       maximum=0.5)
title(main="B) Unenriched housing group", line=3.6, adj=0, cex.main=0.8)

dev.off()

#-------------------------network comparisons-------------------------------#

# compare groups network connectivity (sum of the absolute partial correlations) using unconstrained model
# connectivity of each group, divided by 2 due to symmetric matrix
sum(abs(fit1$pcors[lower.tri(fit1$pcors)])); sum(abs(fit2$pcors[lower.tri(fit2$pcors)]))
#difference in connectivity
sum(abs(fit1$pcors[lower.tri(fit1$pcors)])) - sum(abs(fit2$pcors[lower.tri(fit2$pcors)]))

#following performs the permutation testing

#define function to estimate network at each iteration. Here we use unconstrained partial correlations network
fun <- function (d) { 
  obj <- ggm_inference(Y=as.matrix(d), 
                       alpha=0.05, 
                       control_precision = FALSE,
                       boot = FALSE, 
                       method = "pearson", 
                       progress = TRUE)
  return(obj$pcors)
}

#run NetworkComparisonTest
nct <- NCT(data1 = data1, 
           data2 = data2,
           it=10000, 
           test.edges=TRUE, 
           test.centrality = TRUE, 
           p.adjust.methods="fdr", 
           estimator=fun)

#diff in connectivity permuted data, real data, & pvalue
nct$glstrinv.sep
nct$glstrinv.real 
nct$glstrinv.pval

#strength centrality - corrected using fdr
nct$diffcen.real
nct$diffcen.pval

#individual edge differences - corrected using fdr
nct$edges.tested
nct$einv.real
nct$einv.pvals

#not shown in article but histograms of permuted statistics with real data point estimate on
par(mfrow=c(1,1))
plot(nct, what="strength")

#plot network of significant edge differences 
#first create matrix of differences
diff.matrix <- matrix(0, nrow=5, ncol=5)
diff.matrix[1,2] <- abs(fit1$pcors[1,2] - fit2$pcors[1,2])
diff.matrix[1,5] <- abs(fit1$pcors[1,5] - fit2$pcors[1,5]) 

diff.matrix[lower.tri(diff.matrix)] <- t(diff.matrix)[lower.tri(diff.matrix)]

diff.matrix

#network in Figure 1.7 showing significantly different edge weights between the 2 groups
pdf("Figure 1.7.pdf", width=3.54, height=3.15)

qgraph(diff.matrix, 
       layout='circle', 
       labels=c("AGG","BOLD","ACT","EXP","SOC"), 
       edge.labels=T, 
       theme='colorblind',
       maximum = 1)

dev.off()

#------------------------------------------------------------------------------#
#SECTION 1.8.NETWORK MODELS AS A TOOL TO IDENTIFY SPECIFIC DIRECT ASSOCIATIONS # 
#               BETWEEN STRESSORS AND COMPONENTS OF WELFARE                    #
#------------------------------------------------------------------------------#

#sample size
n4 <- 600

#partial cor network to reproduce with simulated data
network4 <- matrix(data=c(0.00,	0.40,	0.40,	0.00,	0.00,	0.00,	0.00,	0.00,
                          0.40,	0.00,	0.00,	0.20,	0.30,	0.00,	0.00,	0.00,
                          0.40,	0.00,	0.00, 0.00,	0.30,	0.20,	0.20,	0.00,	
                          0.00,	0.20, 0.00,	0.00, 0.20,	0.20,	0.00,	0.00,	
                          0.00,	0.30,	0.30, 0.20,	0.00,	0.00,	0.00,	0.00,	
                          0.00,	0.00,	0.20,	0.20,	0.00,	0.00,	0.00,	0.15,	
                          0.00,	0.00,	0.20,	0.00,	0.00,	0.00,	0.00, 0.00,	
                          0.00,	0.00,	0.00,	0.00,	0.00,	0.15, 0.00,	0.00),
                   nrow=8, ncol=8, byrow=T, dimnames=list(c("AVER","SLE","AGG","NOV","PLA","APP","EXP","STO"),c("AVER","SLE","AGG","NOV","PLA","APP","EXP","STO")))

#check matrix is symmetrical & no human error entering matrix above
isSymmetric(network4)

#generates a function that simulates ordinal data (from bootnet package)
ggmFunction<- ggmGenerator(ordinal=T, nLevels=5)
#simulate ordinal survery data using function just generated
set.seed(123)
data4 <- ggmFunction(n=n4, input=network4)
colnames(data4) <- c("AVER","SLE","AGG","NOV","PLA","APP","EXP","STO")

#estimate constrained network with fisher z method of thresholding
fit4.1 <- ggm_inference(Y=data4, alpha=0.05, control_precision=TRUE, boot=FALSE)
dimnames(fit4.1$wadj) <- list(c("AVER","SLE","AGG","NOV","PLA","APP","EXP","STO"),c("AVER","SLE","AGG","NOV","PLA","APP","EXP","STO"))

#estimate unconstrained network without any significance thresholding
fit4.2 <- estimateNetwork(data=data4, default="pcor")

#the following sets edges to be coloured if significant according to fit 1, and sets other edge weights identified in fit 2 to gray
#colour matrix to use as input to qgraph  
colMatrix4 <- matrix("gray", nrow=nrow(network4), ncol=ncol(network4))

#find significant node indices, as well as positive and negative partial cors 
sigInd4 <- which(fit4.1$wadj < 0 | fit4.1$wadj > 0, arr.ind=T)
sigInd4 
posInd4 <- which(fit4.1$graph >= 0, arr.ind=T)
negInd4 <- which(fit4.1$graph < 0, arr.ind=T)

#run for loop to overwrite non-significant colours with different colours for sig. edges
for (i in 1:nrow(sigInd4)) {
  
  if (fit4.1$wadj[sigInd4[i,1], sigInd4[i,2]] > 0) {
    colMatrix4[sigInd4[i,1], sigInd4[i,2]] <- "blue"
  }
  else {
    if(fit4.1$wadj[sigInd4[i,1],sigInd4[i,2]] < 0) {
      colMatrix4[sigInd4[i,1], sigInd4[i,2]] <- "red"
    }
    else {
      
    }
  }
  
}

nodeColours4 <- c("grey","white","white","white","white","white","white","white")
nodeShapes4 <- c("square","circle","circle","circle","circle","circle","circle","circle")

tiff("Figure 1.8.tiff", width=5.51, height=3.15, units='in', res=1200)

par(mfrow=c(1,1), oma=c(0,0,0,0))

qgraph(fit4.2$graph, 
       layout='spring', 
       edge.color=colMatrix4, 
       label.cex=1.2,edge.labels=T, 
       color=nodeColours4, 
       shape=nodeShapes4)

dev.off()


##################################################################################
#                                End of code                                     #
##################################################################################