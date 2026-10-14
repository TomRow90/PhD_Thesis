################################################################################
#                                                                              #
#                             Chapter 2:                                       #                             
#       Using Network Analysis to Detect Direct Associations Between Pain      #
#                          and Behaviour in Dogs                               #
#                                                                              #
#                                                                              #
#           2. Syntax file for accuracy, robustness & sensitivity              #
#                                                                              #  
################################################################################

library(mgm)

# set wd
setwd("/Users/tomrowland/Documents/PhD/Thesis write up/Chapter 2")

# load data
load("Chapter 2 Analysis.RData")

#----------------------------------------#
#    Accuracy Analysis
#----------------------------------------#

# visualise mean & sd of saturated rank based Gaussian Copula for parameter accuracy estimate

pdf("Figure A3.pdf", width=7.5, height=8.75)

plot(summary(mod), size= 1) + 
  theme(text = element_text(size=4)) + 
  coord_flip() + 
  scale_y_continuous(limits = c(-0.3,0.5), breaks = seq(-0.25,0.6,0.05)) +
  theme(axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10)) + 
  labs(x="Edge", y="Edge Weight / Partial Correlation")

dev.off()

range(summary(mod)$dat_results$Post.sd)

#-------------------------------------------#
#        Robustness & Sensitivity analysis
#-------------------------------------------#

#------effect of prior and BF evidence thresholds---------#


# we will write a function that takes data and a range of priors to perform this analysis
priorFunction <- function(data, priors) {
  
  # loop to estimate models & select graphs for range of supplied priors
  n <- length(priors)
  
  for(i in 1:n) {
    
    pri <- priors[i]
    
    mod <- explore(data, type="mixed", prior_sd=pri, iter=50000, progress=FALSE, impute=TRUE)
    graph <- select(mod, BF_cut=3, alternative="exhaustive")
    
    
    # store results
    if(i == 1) {
      g <- graph$post_prob
    } else {
      g <- rbind(g, graph$post_prob)
    }
    
  }
  
  # to make plotting easier outside function, we will have this function return 3 data frames for the 3 hypotheses tested
  
  # start with sensitivity to prior for positive edges 
  # first create data frame with just that info in
  positiveRobustness <- matrix(0,nrow=190,ncol=4)
  colnames(positiveRobustness) <- c("edge","median","lower","upper")
  positiveRobustness <- as.data.frame(positiveRobustness)
  
  # fill in 1st column with edge list
  positiveRobustness[,1] <- g[1:190,1]
  
  # fill in data frame with data specific for positive edges
  for(i in 1:length(positiveRobustness[,1])) {
    
    ind <- which(g$edge == positiveRobustness[i,1])
    positiveRobustness[i,2] <- median(g[ind,3])
    positiveRobustness[i,3] <- range(g[ind,3])[1]
    positiveRobustness[i,4] <- range(g[ind,3])[2]
    
  }
  
  
  
  # repeat as above but this time for negative edges
  negativeRobustness <- matrix(0,nrow=190,ncol=4)
  colnames(negativeRobustness) <- c("edge","median","lower","upper")
  negativeRobustness <- as.data.frame(negativeRobustness)
  
  negativeRobustness[,1] <- g[1:190,1]
  
  for(i in 1:length(negativeRobustness[,1])) {
    
    ind <- which(g$edge == negativeRobustness[i,1])
    negativeRobustness[i,2] <- median(g[ind,4])
    negativeRobustness[i,3] <- range(g[ind,4])[1]
    negativeRobustness[i,4] <- range(g[ind,4])[2]
    
  }
  
  
  # repeat again but this time for conditional independence edges
  ciRobustness <- matrix(0,nrow=190,ncol=4)
  colnames(ciRobustness) <- c("edge","median","lower","upper")
  ciRobustness <- as.data.frame(ciRobustness)
  
  ciRobustness[,1] <- g[1:190,1]
  
  for(i in 1:length(ciRobustness[,1])) {
    
    ind <- which(g$edge == ciRobustness[i,1])
    ciRobustness[i,2] <- median(g[ind,2])
    ciRobustness[i,3] <- range(g[ind,2])[1]
    ciRobustness[i,4] <- range(g[ind,2])[2]
    
  }
  
  
  # return objects which will be used for plotting outside function
  return <- list(edgeProbs = g, 
                 positiveRobustness = positiveRobustness,
                 negativeRobustness = negativeRobustness, 
                 condIndRobust = ciRobustness)
  
}

# range of priors we will test sensitivity/robustness to
priors <- seq(0.1,0.6,0.05)

# run prior sensitivity function
priorSim <- priorFunction(data=netData, priors = priors)


# robustness plot for hypothesis edges > 0
positivePlot <- ggplot(priorSim$positiveRobustness, 
                       mapping = aes(x = reorder(edge,median),
                                     y = median, 
                                     ymin =lower, 
                                     ymax =upper)) +
                      geom_pointrange(size=0.25, fatten=1) + 
                      coord_flip() + 
                      labs(x="", y="Posterior Probability") +
                      geom_hline(yintercept=0.75, linetype="dotted") + 
                      geom_hline(yintercept=0.91, linetype="dotted") +
                      theme(axis.text.y = element_text(size=4))

pdf("pos edge prior sense.pdf", width=8, height=11)

positivePlot 

dev.off()


# robustness plot for hypothesis edges < 0
negativePlot <- ggplot(priorSim$negativeRobustness, 
                       mapping = aes(x = reorder(edge,median),
                                     y = median, 
                                     ymin =lower, 
                                     ymax =upper)) +
  geom_pointrange(size=0.25, fatten=1) + 
  coord_flip() + 
  labs(x="", y="Posterior Probability") +
  geom_hline(yintercept=0.75, linetype="dotted") + 
  geom_hline(yintercept=0.91, linetype="dotted") +
  theme(axis.text.y = element_text(size=4))

pdf("negative edge prior sens.pdf", width=8, height=11)

  negativePlot 

dev.off()


# robustness plot for hypothesis of conditional independence
ciPlot <- ggplot(priorSim$condIndRobust, 
                 mapping = aes(x = reorder(edge,median),
                               y = median, 
                               ymin =lower, 
                               ymax =upper)) +
  geom_pointrange(size=0.25, fatten=1) + 
  coord_flip() + 
  labs(x="", y="Posterior Probability") +
  geom_hline(yintercept=0.75, linetype="dotted") + 
  geom_hline(yintercept=0.91, linetype="dotted") +
  theme(axis.text.y = element_text(size=4))

pdf("ci prior sens.pdf", width=8, height=11)

  ciPlot 
  
dev.off()


# take a look at pain edges in detail

ind <- which(priorSim$edgeProbs$edge == "age--pain")
age_pain <- priorSim$edgeProbs$prob_greater[ind]

ind <- which(priorSim$edgeProbs$edge == "pain--Energy")
pain_Energy <- priorSim$edgeProbs$prob_less[ind]


pdf("Figure A1.pdf", width=7.5, height=4)

par(mfrow=c(1,2))

plot(priors,age_pain, 
     xlim=c(0.1,0.6), 
     ylim=c(0,1), 
     pch = 19,
     cex = 0.5,
     xlab="Prior SD", 
     ylab="posterior probability that > 0", 
     main="Pain-Age")
lines(priors, age_pain)
abline(h=0.75, lty=2)
abline(h=0.91, lty=2)

plot(priors,pain_Energy, 
     xlim=c(0.1,0.6), 
     ylim=c(0,1),
     pch = 19,
     cex = 0.5,
     xlab="Prior SD", 
     ylab="posterior probability that < 0", main="Pain-Energy")
lines(priors, pain_Energy)
abline(h=0.75, lty=2)
abline(h=0.91, lty=2)

dev.off()


#---------------robustness to UVA procedure------------------#


# robustness to choice of item removal for redundant items
colnames(panas)

# data set with removals present instead 
itemRedData <- as.matrix(panas[,c("age","pain","Fright","NervousStartle",
                                  "PersistPlay","Interest","ExciteWalk","Phobia",
                                  "Escape","Energy","FrightNoises","NotRelaxed",
                                  "LittleEncouragement","DifficultyAdapting",
                                  "AfraidAppliances","PersistNaughty","NotCalmCrowds",
                                  "UnsettledRoutine")])

# estimate rank-based Gaussian Copula network model
mod2 <- explore(itemRedData, type="mixed", prior_sd=0.25, iter=50000, progress=TRUE, impute=TRUE)

# select graph using exhaustive hypothesis testing with evidence set at BF=3
graph2 <- select(mod2, BF_cut=3, alternative="exhaustive")

# use graph object to create conditional dependence adjacency (with + & - pcors in) 
posPcors2 <- graph2$pcor_mat * graph2$pos_mat
negPcors2 <- graph2$pcor_mat * graph2$neg_mat
pcors2 <- posPcors2 + negPcors2
dimnames(pcors2) <- list(colnames(itemRedData), colnames(itemRedData))

# view conditional dependence adjacency 
round(pcors2,2)

# correlated adjacencies
cor(as.vector(pcors[lower.tri(pcors)]),as.vector(pcors2[lower.tri(pcors2)]))

# plot network partial correlations in scatter plot
plot(as.vector(pcors[lower.tri(pcors)]),as.vector(pcors2[lower.tri(pcors2)]), xlab="primary network partial correlations",
     ylab="redundancy robustness network partial correlations")


# Plotting networks for supporting info

# network visualisation objects for qgraph input
names2 <- c("age","pain","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16")
nodeNames2 <- c("","","Fright","EasilyStartled/NervousStartle","PersistPlay","Interest","Excite/ExciteWalk",
                "Phobia","Escape","Energy","FrightNoises","NotRelaxed","NotLazy/LittleEncouragement",
                "DifficultyAdapting","AfraidAppliances","PersistNaughty","NotCalmUnfamiliar/NotCalmCrowds",
                "UnsettledRoutine")

# get graph layout
g <- qgraph(pcors, 
       layout='spring', 
       esize=5, 
       edge.width=2, 
       label.cex=1.5, 
       shape=shape, 
       color=colour, 
       labels=names, 
       theme='colorblind', 
       maximum=0.5,
       vsize=4, 
       GLratio=2.5, 
       layoutOffset=c(0,-0.1), 
       groups=groups, 
       nodeNames=nodeNames, 
       legend=FALSE, 
       title="A) Original network")

layout <- g$layout
layout[2,2] <- -0.4

# plot
pdf("Figure A10.pdf", width=8.75, height=4)

layout(matrix(c(1,2),1,2, byrow=TRUE), widths = c(3.25, 5.5))
par(mar=c(5,2,2,1)+0.1)

qgraph(pcors, 
       layout=g$layout, 
       esize=5, 
       edge.width=2, 
       label.cex=1.5, 
       shape=shape, 
       color=colour, 
       labels=names, 
       theme='colorblind', 
       maximum=0.5,
       vsize=5, 
       GLratio=2.5, 
       layoutOffset=c(0,-0.1), 
       groups=groups, 
       nodeNames=nodeNames, 
       legend=FALSE, 
       title="A) Primary network")

qgraph(pcors2, 
       layout=g$layout, 
       esize=5, 
       edge.width=2, 
       label.cex=1.5, 
       shape=shape, 
       color=colour, 
       labels=names2, 
       theme='colorblind', 
       maximum=0.5,
       vsize=3.75, 
       GLratio=2, 
       layoutOffset=c(0,-0.1), 
       groups=groups, 
       nodeNames=nodeNames2, 
       legend=TRUE,
       legend.cex = 0.3,
       title="B) Redundancy network")

dev.off()

# perform prior & BF robustness on this redundancy robustness network
# before doing this I need to remove some objects from global environment
# otherwise I run into memory storage problems
rm(priorSim) # remove the original prior sensitivity object
rm(mod) # dont use primary model again so remove
rm(graph) # dont use primary graph selected again so remove

priorSim2 <- priorFunction(data=itemRedData, priors=priors)

# robustness plot for hypothesis edges > 0
positivePlot2 <- ggplot(priorSim2$positiveRobustness, 
                       mapping = aes(x = reorder(edge,median),
                                     y = median, 
                                     ymin =lower, 
                                     ymax =upper)) +
  geom_pointrange(size=0.25, fatten=1) + 
  coord_flip() + 
  labs(x="", y="Posterior Probability") +
  geom_hline(yintercept=0.75, linetype="dotted") + 
  geom_hline(yintercept=0.91, linetype="dotted") +
  theme(axis.text.y = element_text(size=4))

pdf("redundancy positives plot.pdf", width=8, height=11)

  positivePlot2
  
dev.off()


# robustness plot for hypothesis edges < 0
negativePlot2 <- ggplot(priorSim2$negativeRobustness, 
                       mapping = aes(x = reorder(edge,median),
                                     y = median, 
                                     ymin =lower, 
                                     ymax =upper)) +
  geom_pointrange(size=0.25, fatten=1) + 
  coord_flip() + 
  labs(x="", y="Posterior Probability") +
  geom_hline(yintercept=0.75, linetype="dotted") + 
  geom_hline(yintercept=0.91, linetype="dotted") +
  theme(axis.text.y = element_text(size=4))

pdf("reundancy negative plot.pdf", width=8, height=11)

  negativePlot2

dev.off()
  
# robustness plot for hypothesis edges = 0
ciPlot2 <- ggplot(priorSim2$condIndRobust, 
                       mapping = aes(x = reorder(edge,median),
                                     y = median, 
                                     ymin =lower, 
                                     ymax =upper)) +
  geom_pointrange(size=0.25, fatten=1) + 
  coord_flip() + 
  labs(x="", y="Posterior Probability") +
  geom_hline(yintercept=0.75, linetype="dotted") + 
  geom_hline(yintercept=0.91, linetype="dotted") +
  theme(axis.text.y = element_text(size=4))

pdf("reundancy cond ind plot.pdf", width=8, height=11)

  ciPlot2

dev.off()
  
# take a look at the Pain - Node edges in detail, are the additional edges robust or not?
ind <- which(priorSim2$edgeProbs$edge == "age--pain")
age_pain <- priorSim2$edgeProbs$prob_greater[ind]

ind <- which(priorSim2$edgeProbs$edge == "pain--Energy")
pain_Energy <- priorSim2$edgeProbs$prob_less[ind]

ind <- which(priorSim2$edgeProbs$edge == "pain--FrightNoises")
pain_frightNoises <- priorSim2$edgeProbs$prob_less[ind]

ind <- which(priorSim2$edgeProbs$edge == "pain--DifficultyAdapting")
pain_difficultyAdapting <- priorSim2$edgeProbs$prob_less[ind]

# visualise these items robustness in the redundancy robustness network
par(mfrow=c(2,2))

plot(priors,age_pain, xlim=c(0.1,0.6), ylim=c(0,1), ylab="posterior probability", main="Pain-Age")
lines(priors, age_pain)
abline(h=0.75, lty=2)
abline(h=0.91, lty=2)
abline(h=0.50, lty=2)

plot(priors,pain_Energy, xlim=c(0.1,0.6), ylim=c(0,1), ylab="posterior probability", main="Pain-Energy")
lines(priors, pain_Energy)
abline(h=0.75, lty=2)
abline(h=0.91, lty=2)
abline(h=0.50, lty=2)

plot(priors,pain_frightNoises, xlim=c(0.1,0.6), ylim=c(0,1), ylab="posterior probability", main="Pain-FrightNoises")
lines(priors, pain_frightNoises)
abline(h=0.75, lty=2)
abline(h=0.91, lty=2)
abline(h=0.50, lty=2)

plot(priors,pain_difficultyAdapting, xlim=c(0.1,0.6), ylim=c(0,1), ylab="posterior probability", main="Pain-DifficultyAdapting")
lines(priors, pain_difficultyAdapting)
abline(h=0.75, lty=2)
abline(h=0.91, lty=2)
abline(h=0.50, lty=2)


#------alternative modelling using mgm-----------

# convert back to data frame & remove missing values (MGM cannot handle NAs)
mgmData <- as.data.frame(netData)
mgmData <- na.omit(mgmData)

# argument for mgm function
type <-  c("g","c",rep("g",16))
level <- c(1,3,rep(1,16))

# estimate mgm
mgmMod <- mgm(as.matrix(mgmData), type=type, level=level, lambdaSel = "EBIC", lambdaGam = 0)

# create unweighted adjancecy matric
edgeInd <- which(mgmMod$pairwise$wadj > 0)
adj <- mgmMod$pairwise$wadj
adj[edgeInd] <- 1

# some plotting arguments
curve <- matrix(0, nrow=18,ncol=18)
curve[1,8] <- 0.7
curve[4,12] <- 0.9

# plot for supporting information
pdf("Figure A3.pdf", width=8.75, height=4)

par(mfrow=c(1,1))

qgraph(adj, 
       layout='spring', 
       weighted=FALSE, 
       label.cex=1, 
       shape=shape, 
       color=colour, 
       labels=names, 
       theme='colorblind', 
       curve=curve, 
       vsize=4, 
       GLratio=1.8, 
       groups=groups, 
       nodeNames=nodeNames, 
       legend=TRUE, 
       legend.cex=0.35, 
       title="Unweighted MGM")

dev.off()

##############################
#           END              #
##############################
