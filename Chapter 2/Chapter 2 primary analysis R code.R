################################################################################
# Chapter 2: Using Network Analysis to Detect Direct Associations Between Pain #
#                         and Behaviour in Dogs                                #
#                                                                              #
#         1. Syntax file for data preprocessing, descriptive statistics,       #
#            Unique Variable Analysis, & primary network analysis              #
#                                                                              #  
################################################################################

#----------------------------------------------
# loading required packages & data-------------
#----------------------------------------------
library(BGGM)
library(qgraph)
library(ggplot2)
library(ggpubr)
library(wTO)
library(bayesplot)
library(e1071)

setwd("/Users/tomrowland/Documents/PhD/Thesis write up/Chapter 2")

# load data
panas <- read.table("panasFinal.txt", header=T)

#----------------------------------------------
# data processing 
#----------------------------------------------

# remove id col
panas <- panas[,-1]

# rename variables/columns
colnames(panas) <- c("age","pain","Fright","ExciteWalk","EasilyStartled","PersistPlay","Interest","NervousStartle","Excite","Phobia","Escape","NotCalmCrowds",
                     "Energy","FrightNoises","NotRelaxed","NotLazy","DifficultyAdapting","AfraidAppliances","LittleEncouragement","PersistNaughty",
                     "NotCalmUnfamiliar","Boisterous","UnsettledRoutine")

# sample size for each pain category
length(which(panas[,2]==1)) #725
length(which(panas[,2]==2)) #96
length(which(panas[,2]==3)) #121

# missing data
for(i in 1:ncol(panas)) {
  
  v <- colnames(panas)[i]
  missing <- ((length(which(is.na(panas[,i])))) / nrow(panas)) * 100
  print(c("variable:", v, "missing data = ", round(missing,2)))
  
}

#---------------------------------------------------------------
# Identify redundancies with weighted topological overlap (wTO)
#---------------------------------------------------------------

# redundancy function
redundancyDetection <- function(data, DoNotPlot = FALSE) {
  
  d <- as.matrix(data)
  
  #step 1: estimate spCGGM
  spCGGM <- explore(d, type="mixed", iter = 50000, prior_sd=0.25, impute=TRUE, progress=TRUE)
  
  #step 2: select conditional dependence graph
  graph <- select(spCGGM, BF_cut=3, alternative="exhaustive")
  
  #use graph object to create conditional dependence adjacency (with + & - pcors in) 
  posPcors <- graph$pcor_mat * graph$pos_mat
  negPcors <- graph$pcor_mat * graph$neg_mat
  pcors <- posPcors + negPcors
  
  dimnames(pcors) <- list(colnames(d), colnames(d))
  
  #step 3: estimate wTo
  wTo <- wTO(pcors, sign="sign")
  
  #visualise wTo
  print(round(wTo, 4))
  shape <- c(rep("circle",ncol(d)))
  names <- seq(1,ncol(d),1)
  nodeNames <- c(colnames(d))
  title <- "Redundancy Chains"
  
  #visualise graph thresholded at 0.2 to show redundancies detected
  network <- qgraph(wTo, 
                    threshold=0.25, 
                    vsize=3, 
                    edge.labels=TRUE, 
                    edge.label.cex=0.6, 
                    cut=0.25,
                    label.cex=1.8,
                    maximum=0.3, 
                    GLratio=1.2, 
                    legend.cex=0.6, 
                    title=title, 
                    labels=names, 
                    nodeNames=nodeNames, 
                    legend=TRUE,
                    DoNotPlot = DoNotPlot)
  
  #return in case needed
  return(list(wTo=wTo, redundancyChain=network, graph=graph, qgraph_obj = network))
  
}

# iteration 1
# Run it 4 times 'chains' and check get same result
par(mfrow = c(2,2))

for(i in 1:4) {
  
  it1 <- redundancyDetection(as.matrix(panas[,-c(1,2)]), DoNotPlot = TRUE)
  
  qgraph(it1$wTo, 
         threshold=0.25, 
         vsize=3, 
         edge.labels=TRUE, 
         edge.label.cex=1.5, 
         cut=0.25,
         label.cex=0.8,
         maximum=0.3, 
         legend.cex=0.6, 
         labels=colnames(panas[,-c(1,2)]), 
         label.scale = FALSE)  
  
}


# redundancies detected:
# EasilyStartled - NervousStartle
# NotCalmCrowds - NotCalmUnfamiliar
# Little encouragement - Not Lazy
# based on more general item content, retain Easilystartled, NotCalmUnfamiliar, and Not Lazy

# create new dataset with redundant items removed
rm.idx <- which(colnames(panas) %in% c("NervousStartle","NotCalmCrowds","LittleEncouragement"))
netData <- as.matrix(panas[,-rm.idx])
colnames(netData) <- colnames(panas[,-rm.idx])

# remove it1 object otherwise encounter problems with memory
rm(it1)


# iteration 2 of UVA to check for remaining redundancies 
# again, run it 4 times 'chains' and check get same result
par(mfrow = c(2,2))

for(i in 1:4) {
  
  it2 <- redundancyDetection(as.matrix(netData[,-c(1,2)]), DoNotPlot = TRUE)
  qgraph(it2$wTo, 
         threshold=0.25, 
         vsize=3, 
         edge.labels=TRUE, 
         edge.label.cex=1.5, 
         cut=0.25,
         label.cex=0.8,
         maximum=0.3, 
         legend.cex=0.6, 
         labels=colnames(netData[,-c(1,2)]), 
         label.scale = FALSE)   
}

# redundancies detected:
# Excite - ExciteWalk
# based on more general item content, retain Excite

# remove it2 object otherwise encounter problems with memory
rm(it2)

# create new dataset with redundant items removed
rm.idx <- which(colnames(netData) == "ExciteWalk")
netData <- as.matrix(netData[,-rm.idx])

# iteration 3: repeat UVA to check for remaining redundancies 
# Run it 4 times 'chains' and check get same result
par(mfrow = c(2,2))

for(i in 1:4) {
  
  it3 <- redundancyDetection(as.matrix(netData[,-c(1,2)]), DoNotPlot = TRUE)
  qgraph(it3$wTo, 
         threshold=0.25, 
         vsize=3, 
         edge.labels=TRUE, 
         edge.label.cex=1.5, 
         cut=0.25,
         label.cex=0.8,
         maximum=0.3, 
         legend.cex=0.6, 
         labels=colnames(netData[,-c(1,2)]), 
         label.scale = FALSE)    
}

# redundancies detected:
# Boisterous - Excite
# based on more general item content, retain Excite

# data set with redundancies removed
rm.idx <- which(colnames(netData) == "Boisterous")
netData <- as.matrix(netData[,-rm.idx])

# remove it3 object other encounter problems with memory
rm(it3)


# iteration 4: repeat UVA to check for remaining redundancies 
# Run it 4 times 'chains' and check get same result
par(mfrow = c(2,2))

for(i in 1:4) {
  
  it4 <- redundancyDetection(as.matrix(netData[,-c(1,2)]), DoNotPlot = TRUE)
  qgraph(it4$wTo, 
         threshold=0.25, 
         vsize=3, 
         edge.labels=TRUE, 
         edge.label.cex=1.5, 
         cut=0.25,
         label.cex=0.8,
         maximum=0.3, 
         legend.cex=0.6, 
         labels=colnames(netData[,-c(1,2)]), 
         label.scale = FALSE)    
}

# while two redundancies suggested, theoretically I am happy for them to
# remain in the graph because they can be argued to represent different constructs

# remove it4 object other encounter problems with memory
rm(it4)

#-------------------------------
# descriptives
#--------------------------------

for(i in 1:ncol(netData)) {
  hist(netData[,i], main=colnames(netData)[i])
}

for(i in 1:ncol(netData)) {
  
  na <- round((length(which(is.na(netData[,i]))) / nrow(netData)) * 100, 2)
  print(c("Variable",colnames(netData)[i],"missing data =",na, "%"))
}

#-------------------------------------
# Primary network analysis with BGGM
#-------------------------------------

# estimate rank-based Gaussian Copula network model
mod <- explore(netData, type="mixed", prior_sd=0.25, iter=50000, progress=TRUE, impute=TRUE)
# select graph using exhaustive hypothesis testing with evidence set at BF=3
graph <- select(mod, BF_cut=3, alternative="exhaustive")

# get posterior probs to correlation against chain 2
mod1.null <- graph$post_prob$prob_zero
mod1.pos <- graph$post_prob$prob_greater
mod1.neg <- graph$post_prob$prob_less

# estimate again to check posterior probabilities are stable
mod2 <- explore(netData, type="mixed", prior_sd=0.25, iter=50000, progress=TRUE, impute=TRUE)
# select graph using exhaustive hypothesis testing with evidence set at BF=3
graph2 <- select(mod2, BF_cut=3, alternative="exhaustive")

# get posterior probs to correlation against chain 1
mod2.null <- graph2$post_prob$prob_zero
mod2.pos <- graph2$post_prob$prob_greater
mod2.neg <- graph2$post_prob$prob_less

# plot probs against each other from each chain
# null probabilities
par(mfrow = c(1,1))
plot(mod1.null, mod2.null,
     ylim = c(0,1),
     xlim = c(0,1),
     xlab = "chain 1 null probabilities",
     ylab = "chain 2 null probabilities")
abline(a = 0, b = 1, col = "red", lwd = 2)

cor(mod1.null, mod2.null)

# positive probabilities
plot(mod1.pos, mod2.pos,
     ylim = c(0,1),
     xlim = c(0,1),
     xlab = "chain 1 positive probabilities",
     ylab = "chain 2 positive probabilities")
abline(a = 0, b = 1, col = "red", lwd = 2)

cor(mod1.pos, mod2.pos)

# neg probabilities
plot(mod1.neg, mod2.neg,
     ylim = c(0,1),
     xlim = c(0,1),
     xlab = "chain 1 negative probabilities",
     ylab = "chain 2 negative probabilities")
abline(a = 0, b = 1, col = "red", lwd = 2)

cor(mod1.neg, mod2.neg)

# they all look pretty good - very high, use mod1 here on out
# remove mod2 in case of memory problems
rm(mod2)

#get posterior mean and sd
summary(mod)

#-----------------------------------------
#     POSTERIOR PREDICTIVE CHECKS
# before interpreting model, perform PPC
#-----------------------------------------

# generate samples from posterior predictive distribution
posteriorPredictions <- posterior_predict(mod, iter=5000, progress=TRUE)

# due to missing data, create empirical set with listwise deletion, and delete corresponding indices from PPD
ind <- which(is.na(netData), arr.ind=TRUE)
ind <- ind[,1]
ind <- unique(ind)
netData2 <- netData[-ind, ]
yrep <- posteriorPredictions[-ind, , ]  

point.size <- 0.1
line.width <- 0.5
iter <- 5000

# create ppc_bar plots for each variable for qualitative PPC
agePPC <- ppc_bars(y = netData2[,1], t(yrep[,1,1:iter]), size=point.size, linewidth=line.width) + ggtitle("Age") + theme(plot.title = element_text(hjust = 0.5, size=8))
painPPC <- ppc_bars(y = netData2[,2], t(yrep[,2,1:iter]), size=point.size, linewidth=line.width) + scale_x_discrete(limits=c("1","2","3")) + theme(axis.title.y = element_blank()) + ggtitle("Pain") + theme(plot.title = element_text(hjust = 0.5, size=8)) 
frightPPC <- ppc_bars(y = netData2[,3], t(yrep[,3,1:iter]), size=point.size, linewidth=line.width) + theme(axis.title.y = element_blank()) + ggtitle("Fright") + theme(plot.title = element_text(hjust = 0.5, size=8))
easilyStartledPPC <- ppc_bars(y = netData2[,4], t(yrep[,4,1:iter]), size=point.size, linewidth=line.width) + ggtitle("Easily Startled") + theme(plot.title = element_text(hjust = 0.5, size=8))
persistPlayPPC <- ppc_bars(y = netData2[,5], t(yrep[,5,1:iter]), size=point.size, linewidth=line.width) + theme(axis.title.y = element_blank()) + ggtitle("Persist Play") + theme(plot.title = element_text(hjust = 0.5, size=8))
interestPPC <- ppc_bars(y = netData2[,6], t(yrep[,6,1:iter]), size=point.size, linewidth=line.width) + theme(axis.title.y = element_blank()) + ggtitle("Interest") + theme(plot.title = element_text(hjust = 0.5, size=8))
excitePPC <- ppc_bars(y = netData2[,7], t(yrep[,7,1:iter]), size=point.size, linewidth=line.width) + theme(axis.title.y = element_blank()) + ggtitle("Excite") + theme(plot.title = element_text(hjust = 0.5, size=8))
phobiaPPC <- ppc_bars(y = netData2[,8], t(yrep[,8,1:iter]), size=point.size, linewidth=line.width) + ggtitle("Phobia") + theme(plot.title = element_text(hjust = 0.5, size=8))
escapePPC <- ppc_bars(y = netData2[,9], t(yrep[,9,1:iter]), size=point.size, linewidth=line.width) + theme(axis.title.y = element_blank()) + ggtitle("Escape") + theme(plot.title = element_text(hjust = 0.5, size=8))
energyPPC <- ppc_bars(y = netData2[,10], t(yrep[,10,1:iter]), size=point.size, linewidth=line.width) + theme(axis.title.y = element_blank()) + ggtitle("Energy") + theme(plot.title = element_text(hjust = 0.5, size=8))
frightNoisesPPC <- ppc_bars(y = netData2[,11], t(yrep[,11,1:iter]), size=point.size, linewidth=line.width) + theme(axis.title.y = element_blank()) + ggtitle("Fright Noises") + theme(plot.title = element_text(hjust = 0.5, size=8))
notRelaxedPPC <- ppc_bars(y = netData2[,12], t(yrep[,12,1:iter]), size=point.size, linewidth=line.width) + ggtitle("Not Relaxed") + theme(plot.title = element_text(hjust = 0.5, size=8))
notLazyPPC <- ppc_bars(y = netData2[,13], t(yrep[,13,1:iter]), size=point.size, linewidth=line.width) + theme(axis.title.y = element_blank()) + ggtitle("Not Lazy") + theme(plot.title = element_text(hjust = 0.5, size=8))
difficultyAdaptingPPC <- ppc_bars(y = netData2[,14], t(yrep[,14,1:iter]), size=point.size, linewidth=line.width) + theme(axis.title.y = element_blank()) + ggtitle("Difficulty Adapting") + theme(plot.title = element_text(hjust = 0.5, size=8))
afraidAppliancesPPC <- ppc_bars(y = netData2[,15], t(yrep[,15,1:iter]), size=point.size, linewidth=line.width) + theme(axis.title.y = element_blank()) + ggtitle("Afraid Appliances") + theme(plot.title = element_text(hjust = 0.5, size=8))
persistNaughtyPPC <- ppc_bars(y = netData2[,16], t(yrep[,16,1:iter]), size=point.size, linewidth=line.width)  + ggtitle("Persist Naughty") + theme(plot.title = element_text(hjust = 0.5, size=8))
notCalmUnfamiliarPPC <- ppc_bars(y = netData2[,17], t(yrep[,17,1:iter]), size=point.size, linewidth=line.width) + theme(axis.title.y = element_blank()) + ggtitle("Not Calm Unfamiliar") + theme(plot.title = element_text(hjust = 0.5, size=8))
unsettledRoutinePPC <- ppc_bars(y = netData2[,18], t(yrep[,18,1:iter]), size=point.size, linewidth=line.width) + theme(axis.title.y = element_blank()) + ggtitle("Unsettled Routine") + theme(plot.title = element_text(hjust = 0.5, size=8))


#save ppc_bar plots for supplementary matrial
pdf("Figure A4.pdf", width=7.5, height=8.75)

ggarrange(agePPC, painPPC, frightPPC, easilyStartledPPC, persistPlayPPC, interestPPC, 
          excitePPC, phobiaPPC, escapePPC, energyPPC, frightNoisesPPC,
          notRelaxedPPC, notLazyPPC, difficultyAdaptingPPC, afraidAppliancesPPC, 
          persistNaughtyPPC, notCalmUnfamiliarPPC, unsettledRoutinePPC, 
          ncol=4, nrow=5, align="h", common.legend = TRUE)

dev.off()

# note: qualitative inspection of PPD suggest good reproduction of empirical data

# PPC to examine how well model captures skew of original data
ppcSkewList <- list()

for(i in 1:ncol(netData)) {
  
  ppcSkewList[[i]] <- ppc_stat(y = netData2[,i], yrep=t(yrep[,i,1:iter]), stat=skewness)
  ppcSkewList[[i]]  <- ppcSkewList[[i]] + ggtitle(colnames(netData2)[i]) + theme(plot.title = element_text(hjust = 0.5, size=8))
}

for(i in 1:ncol(netData)) {
  
  print(ppcSkewList[[i]])
}


# save ppc_bar plots for appendix A
pdf("Figure A5.pdf", width=7.5, height=8.75)

ggarrange(ppcSkewList[[1]],ppcSkewList[[2]],ppcSkewList[[3]],ppcSkewList[[4]],ppcSkewList[[5]],
          ppcSkewList[[6]],ppcSkewList[[7]],ppcSkewList[[8]],ppcSkewList[[9]],ppcSkewList[[10]],
          ppcSkewList[[11]],ppcSkewList[[12]],ppcSkewList[[13]],ppcSkewList[[14]],ppcSkewList[[15]],
          ppcSkewList[[16]],ppcSkewList[[17]],ppcSkewList[[18]], ncol=4, nrow=5, align="h", common.legend = TRUE)
dev.off()


# note: PPC checks qualitatively adequate - move on to inference on copula model

#---------------------------------------------
# Graph selections
#---------------------------------------------

# use graph object to create conditional dependence adjacency (with + & - pcors in) 
posPcors <- graph$pcor_mat * graph$pos_mat
negPcors <- graph$pcor_mat * graph$neg_mat
pcors <- posPcors + negPcors
dimnames(pcors) <- list(colnames(netData), colnames(netData))

# view conditional dependence adjacency 
round(pcors,2)

# view conditional independence adjacency
graph$null_mat

# create ambiguous edge set adjacency matrix
fullMatrix <- matrix(1, nrow=18, ncol=18)
ambAdj <- fullMatrix - (graph$pos_mat + graph$neg_mat + graph$null_mat)

# view ambiguous adjacency
ambAdj

#----------Marginal Correlation Network Structure for Indirect effects---------------------#

zeroOrderNet <- pcor_to_cor(mod, iter=50000)
zeroOrderNet$R_mean

# select adjacency matrix of marginal correlation structure based on 95% cred. int. 
cors <- zeroOrderNet$R[,,1:50000]

adj <- ifelse(apply(cors, 1:2, quantile, 0.025) < 0 &
                apply(cors, 1:2, quantile, 0.975) > 0, 0, 1)

# find indices of which dependencies were selected
ind <- which(adj == 1)
zeroOrderWadj <- zeroOrderNet$R_mean
dimnames(zeroOrderWadj) <- list(colnames(netData), colnames(netData))

# replace entries whose CI crosses 0 in the adjacency matrix with 0 to obtain final graph
zeroOrderWadj[-ind] <- 0

# find 1sd of marginal cor posterior distributions
sds <- apply(cors, 1:2, sd)
dimnames(sds) <- list(colnames(netData), colnames(netData))
round(zeroOrderWadj,2)
round(sds,2)

#------------------------------------------
# Plotting networks for publication
#------------------------------------------

# network visualisation objects for qgraph input
shape <- c("square","square",rep("circle",16))
colour <- c("grey","green4","orange2")
names <- c("age","pain","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16")
nodeNames <- c("","",colnames(netData[,-c(1,2)]))
groups <- list("External system variables" = c(1,2), "Positive Activation"=c(5,6,7,9,10,13,16), "Negative Activation"=c(3,4,8,11,12,14,15,17,18))


pdf("Figure 2.1.pdf", width=5.5, height=7)

  layout(matrix(c(1,4,2,4,3,4),3,2, byrow=TRUE))

  par(mar=c(5,2,2,1)+0.1)

  qgraph(zeroOrderWadj, 
         layout='circle', 
         esize=5, 
         edge.width=2,
         label.cex=1.5, 
         shape=shape, 
         color=colour, 
         labels=names, 
         theme='colorblind', 
         maximum=0.5,
         curve=curve, 
         vsize=9, 
         GLratio=2.5, 
         layoutOffset=c(0,-0.1), 
         groups=groups, 
         legend=FALSE, 
         nodeNames=nodeNames, 
         title="A) Zero-order Dependencies")
  
  qgraph(pcors, 
         layout='circle', 
         esize=5, 
         edge.width=2, 
         label.cex=1.5, 
         shape=shape, 
         color=colour, 
         labels=names, 
         theme='colorblind', 
         maximum=0.5,
         vsize=9, 
         GLratio=2.5, 
         layoutOffset=c(0,-0.1), 
         groups=groups, 
         nodeNames=nodeNames, 
         legend=FALSE, 
         title="B) Conditional Dependencies")
  
  qgraph(graph$null_mat, 
         layout='circle', 
         label.cex=1.5, 
         shape=shape, 
         color=colour, 
         labels=names,
         vsize=9, 
         GLratio=2.5, 
         layoutOffset=c(0,-0.1), 
         groups=groups, 
         nodeNames=nodeNames, 
         legend=FALSE, 
         title="C) Conditional Independencies")
  
  
  plot.new()
  
  legend(.1,0.75, legend=c("3: PersistPlay","4: Interest","5: Excite","7: Escape","8: Energy","11: NotLazy","14: PersistNaughty"), title=as.expression(bquote(bold("Positive Activation Items"))), col="green4", title.adj=0, pch=16, cex=1.2, bty="n")
  
  legend(.13,0.45, legend=c("1: Fright","2: EasilyStartle","6: Phobia","9:FrightNoises","10:NotRelaxed","12: DifficultyAdapting","13: AfraidAppliances","15: NotCalmUnfamiliar","16: UnsettledRoutine"), title.adj=0, title=as.expression(bquote(bold("Negative Activation Items"))), col="orange2", pch=16, cex=1.2, bty="n")

dev.off()

# maximum values referred to in figure legend
pos.idx <- which(zeroOrderWadj[lower.tri(zeroOrderWadj)] > 0)
max((zeroOrderWadj[lower.tri(zeroOrderWadj)])[pos.idx])

neg.idx <- which(zeroOrderWadj[lower.tri(zeroOrderWadj)] < 0)
min((zeroOrderWadj[lower.tri(zeroOrderWadj)])[neg.idx])

pos.idx <- which(pcors[lower.tri(pcors)] > 0)
max((pcors[lower.tri(pcors)])[pos.idx])

neg.idx <- which(pcors[lower.tri(pcors)] < 0)
min((pcors[lower.tri(pcors)])[neg.idx])

#------------------------------------------------
#     Conditional probability estimation
# from the posterior predictive distribution
#------------------------------------------------

# let's do this 4x and check estimates and bounds are stable
# with just 1000 iter. More iter and we run into a problem with simulated
# datasets having 0 entries for <6m. So let's set a seed and use 1000 iter
# but check the estimates and bounds are stable across 4 separate runs

seeds <- c(11,12,133,1274)
conditional.table.list <- list()

for(i in 1:4) {
  
  set.seed(seeds[i])
  posteriorPredictions <- posterior_predict(mod, iter=1000, progress=TRUE)

  #Create conditional probability table for P(Pain | ......)
  conditional.table.list[[i]] <- as.data.frame(matrix(0, nrow=11, ncol=4))
  colnames(conditional.table.list[[i]]) <- c("conditioning","probability","lowerBound","upperBound")

  # age is 6 levels, energy is 5 levels, we want to condition in each level of these
  conditional.table.list[[i]][,1] <- c(1:6, 1:5)


  #fill in remaining columns with mean and credible bounds
  for(j in 1:11) {
  
    tryCatch(if(j <= 6) {
      p <- predicted_probability(posteriorPredictions, Y = netData, outcome="pain", age=j)
      conditional.table.list[[i]][j,2] <- median(p$collect[,3])
      conditional.table.list[[i]][j,3] <- quantile(p$collect[,3], prob=0.025)
      conditional.table.list[[i]][j,4] <- quantile(p$collect[,3], prob=0.975)
    } else {
      if(j > 6 & j <=11) {
        p <- predicted_probability(posteriorPredictions, Y = netData, outcome="pain", Energy=j-6)
        conditional.table.list[[i]][j,2] <- median(p$collect[,3])
        conditional.table.list[[i]][j,3] <- quantile(p$collect[,3], prob=0.025)
        conditional.table.list[[i]][j,4] <- quantile(p$collect[,3], prob=0.975)
      } 
    }, error=function(e) {NA})
  }

}

conditional.table.list

#---------Dot plot of point est. & CIs---------

agePlots <- list()

for(i in 1:4) {
  
  agePlots[[i]] <- ggplot(conditional.table.list[[i]][1:6,], aes(x = conditioning, y = probability)) +
                geom_pointrange(aes(y = probability, ymin = lowerBound, ymax = upperBound), 
                                size=0.3, 
                                fatten=2) +
                ylim(c(0,1)) + 
                scale_x_discrete(limits=c("<6m","6m-1y","1y-2y","2y-6y","6y-10y",">10y")) +
                ggtitle("A)") +
                labs(x="Item response", y="P(pain | Age)") + 
                theme_classic(base_size = 12) +
                theme(axis.text.x = element_text(size=10, angle=90), 
                      axis.text.y = element_text(size=10)) 
}

# view plots to check estimates and intervals stable
ggarrange(agePlots[[1]], agePlots[[2]],
          agePlots[[3]], agePlots[[4]], nrow = 2, ncol = 2)


# same again for P(pain | energy)
energyPlots <- list()

for(i in 1:4) {
  
  energyPlots[[i]] <- ggplot(conditional.table.list[[i]][7:11,], aes(x = conditioning, y = probability)) +
    geom_pointrange(aes(y = probability, ymin = lowerBound, ymax = upperBound), 
                    size=0.3, 
                    fatten=2) +
    ylim(c(0,1)) + 
    scale_x_discrete(limits=c("1","2","3","4","5")) +
    ggtitle("B)") +
    labs(x="Item response", y="P(pain | Energy)") + 
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(size=10), 
          axis.text.y = element_text(size=10)) 
}

# view each chains est + interval
ggarrange(energyPlots[[1]], energyPlots[[2]],
          energyPlots[[3]], energyPlots[[4]], nrow = 2, ncol = 2)

# estimates and intervals look pretty stable across runs so go with that
# number of PPD samples and just plot first run of age and energy probs


# save figure
pdf("Figure 2.2.pdf", width=8, height=3.5)

ggarrange(agePlots[[1]],energyPlots[[1]],  
          ncol=2, nrow=1, 
          align="h", 
          common.legend = FALSE)

dev.off()


##########################################################################

# finished, change wd & save workspace
setwd("/Users/tomrowland/Documents/PhD/Thesis write up/Chapter 2")

save.image(file="Chapter 2 Analysis.RData")

#########################################################################
