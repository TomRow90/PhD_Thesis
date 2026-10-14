##############################################
#
#    Chapter 4: Using network
#.   analysis to gain mechanistic insights
#    in animal welfare science
#
#       ***R code for network analysis***
#
##############################################

#load packages
library(welfareNet)
library(BMGM)
library(qgraph)
library(forcats)
library(ggplot2)
library(dplyr)
library(tidyr)
library(NetworkToolbox)


#--------------------------------------------#
#------------ functions used ----------------#
#--------------------------------------------#

#transformation for beta distributed variables with any precise 0s or 1s
beta.transform <- function(y, n) (y * (n - 1) + 0.5) / n


###############################################################################
#                               load data                                     #
###############################################################################

# set wd
setwd("/Users/tomrowland/Documents/PhD/Thesis write up/Chapter 4/data")

# load data
network.data <- read.table("network.data.txt", header=TRUE)

# change some column names
colnames(network.data) <- c("id", "condition", "diving bottom","opercular beat rate",
                            "scototaxis","bw stationary","thigmotaxis","oft stationary",
                            "not stationary","not latency","nor discrimination","startle",
                            "social density","home feed","activity","aggression",
                            "stationary bottom","foraging","wall swimming",
                            "vertical distribution")


# transform relevant variables to same scale that means were assessed
network.data$`diving bottom` <- network.data$`diving bottom` / 300

#squeeze transformation within [0,1]
n <- length(na.omit(network.data$activity))
network.data$activity <- beta.transform(y=network.data$activity, n=n)

n <- length(na.omit(network.data$aggression))
network.data$aggression <- beta.transform(y=network.data$aggression, n=n)

n <- length(na.omit(network.data$`wall swimming`))
network.data$`wall swimming` <- beta.transform(y=network.data$`wall swimming`, n=n)

# now separate dataset to create data sets for both conditions with only the
# reduced set of variables
ee.ind <- which(network.data$condition == "ee")
std.ind <- which(network.data$condition == "std")

ee.data <- network.data[ee.ind, c("diving bottom","opercular beat rate",
                                  "social density","activity","aggression","foraging",
                                  "wall swimming","vertical distribution")]
std.data <- network.data[std.ind,c("diving bottom","opercular beat rate",
                                   "social density","activity","aggression","foraging",
                                   "wall swimming","vertical distribution")]

# quick check
head(ee.data)
head(std.data)

################################################
#                                              #
#     Fuzzy network comparison analysis        #
#                                              #  
################################################

# estimate networks with no p value adjustments
# this is to get the P values to use for fuzzy network estimation
ee.net <- zero.order.network(data = ee.data,
                             select = "sig",
                             dcor.permutations = 100000,
                             threshold = 0.05,
                             adjust = "none")

std.net <- zero.order.network(data = std.data,
                              select = "sig",
                              dcor.permutations = 100000,
                              threshold = 0.05,
                              adjust = "none")


# estimate individual fuzzy networks from p values
ee.fuzzy <- estimate.fuzzy.network(p.value.matrix = ee.net$p.values, 
                                   prior = 0.5)
std.fuzzy <- estimate.fuzzy.network(p.value.matrix = std.net$p.values, 
                                    prior=0.5)

# plot figure 4.3 in main text
pdf("Figure 4.3.pdf", width = 8, height=3.65)

par(mfrow = c(1,2))

par(pty="s")
qgraph::qgraph(std.fuzzy$fuzzy.network, 
               layout='circle', 
               theme = 'colorblind',
               labels=colnames(std.data), 
               edge.width = 0.5,
               maximum=1, 
               colFactor = 2, 
               edge.labels=TRUE, 
               edge.label.cex=1, 
               label.scale=FALSE, 
               label.cex=0.7,
               node.label.offset=c(0.5,0.9), 
               edge.label.position=0.4)
text(x=-0.25, y=1.2, labels="A) Unenriched housing condition", font=2, cex=0.8)

par(pty="s")
qgraph::qgraph(ee.fuzzy$fuzzy.network, 
               layout='circle', 
               theme = 'colorblind',
               labels=colnames(ee.data), 
               edge.width = 0.5,
               maximum=1, 
               colFactor = 2, 
               edge.labels=TRUE, 
               edge.label.cex=1,
               label.scale=FALSE, 
               label.cex=0.7,
               node.label.offset=c(0.5,0.9), 
               edge.label.position=0.4)
text(x=-0.25, y=1.2, labels="B) Enriched housing condition", font=2, cex=0.8)

dev.off()


# construct distribution for network density of each group
ee.density <- fuzzy.density.distribution(network = ee.fuzzy$fuzzy.network)
std.density <- fuzzy.density.distribution(network = std.fuzzy$fuzzy.network)

# total number of edges
ee.density.summary <- fuzzy.summarise.density(fuzzy.net = ee.fuzzy$fuzzy.network)
ee.density.summary
std.density.summary <- fuzzy.summarise.density(fuzzy.net = std.fuzzy$fuzzy.network)
std.density.summary

# density JSDs
P <- ee.density$density.pmf$probability
Q <- std.density$density.pmf$probability
m <- as.matrix(rbind(P,Q))
density.jsd <- round(philentropy::JSD(m),2)
density.jsd

# not shown in chapter but useful visualisation
nodes <- ncol(ee.data)
edges <- nodes*(nodes - 1) / 2
K <- 0:edges

par(mfrow = c(1,1))

plot(K, ee.density$density.pmf$probability, pch = 19, 
     ylab = "Probability",
     xlab = "Density",
     cex.main = 0.5,
     cex = 0.25,
     cex.lab = 0.5,
     cex.axis = 0.5,
     ylim = c(0, 0.2),
     xlim = c(0, edges),
     col="green3")
title( main = "A) Network density PMF",
       adj = 0,
       cex.main = 0.5)
lines(K, ee.density$density.pmf$probability, col = "green3")

# add std pmf
points(K, std.density$density.pmf$probability, pch=19, cex=0.25, col = "red3")
lines(K, std.density$density.pmf$probability, col = "red3")

# add expected values
abline(v = ee.density.summary$expected.value, lty = 2, col = "green3")
abline(v = std.density.summary$expected.value, lty = 2, col = "red3")

# add text
text(5,0.15, labels=paste0("JSD = ", density.jsd), cex = 0.5, font = 2)

# and view CDF
plot(K, ee.density$density.cdf$probability, pch = 19, 
     ylab = "Probability",
     xlab = "Density",
     cex = 0.25,
     cex.lab = 0.5,
     cex.axis = 0.5,     ylim = c(0, 1),
     xlim = c(0, edges),
     col="green3")
title(main = "B) Network density CDF",
      adj = 0,      
      cex.main = 0.5)
lines(K, ee.density$density.cdf$probability, col = "green3")

points(K, std.density$density.cdf$probability, pch=19, cex=0.25, col = "red3")
lines(K, std.density$density.cdf$probability, col = "red3")


# lets get the prior implied density distributions as well so we can add
# priors to the plot
ee.density.prior <- fuzzy.density.distribution(network = ee.fuzzy$priors)
std.density.prior <- fuzzy.density.distribution(network = std.fuzzy$priors)

# construct delta distributions for prior and posterior
density.delta.prior <- fuzzy.difference.distribution(prob.x = ee.density.prior$density.pmf$probability,
                                                     prob.y = std.density.prior$density.pmf$probability)

density.delta.dist <- fuzzy.difference.distribution(prob.x = ee.density$density.pmf$probability,
                                                    prob.y = std.density$density.pmf$probability)


# summary probabilities
density.delta.probs <- fuzzy.delta.probabilities(delta = density.delta.dist)
density.delta.probs

# ---------------- prior sensitivity for network density ----------------- #

# list of prior probabilities of edge inclusion 
# for prior sensitivity analysis
priors <- seq(0.1,0.5,0.05)

# perform prior sensitivity analysis
density.prior.sensitivity <- fuzzy.density.prior.sensitivity(p.value.matrix1 = ee.net$p.values,
                                                             p.value.matrix2 = std.net$p.values,
                                                             priors = priors,
                                                             epsilon = 0)
# view results    
density.prior.sensitivity


# plot for figure 4.4
pdf("Figure 4.4.pdf", width=8, height=3.65)

par(mfrow = c(1,2),
  mar = c(2.5, 3, 2, 1),      
  mgp = c(1.5, 0.2, 0),     
  tcl = -0.2)

plot(density.delta.dist$delta, density.delta.dist$prob, 
     ylim = c(0,0.20), 
     pch = 19, 
     xlab = "Density difference", 
     ylab = "Probability", 
     cex = 0.5, 
     cex.axis = 0.6, 
     cex.lab = 0.8)
title(main = "A) Network density delta distribution", 
      cex.main = 0.8, 
      adj = 0, 
      font.main = 2,
      line = 0.5)
lines(density.delta.dist$delta, density.delta.dist$prob)
points(density.delta.prior$delta, density.delta.prior$prob, 
       pch = 1, 
       cex = 0.5)
lines(density.delta.prior$delta, density.delta.prior$prob, 
      lty = 2)
abline(v = 0, lty = 2)
text(x = -15, y = 0.17, labels = paste0("P(delta > 0) = ",round(density.delta.probs$P_greater,2)), cex = 0.55)
text(x = -15, y = 0.15, labels = paste0("P(delta = 0) = ",round(density.delta.probs$P_equal,2)), cex = 0.55)
text(x = -15, y = 0.13, labels = paste0("P(delta < 0) = ",round(density.delta.probs$P_less,2)), cex = 0.55)

plot(density.prior.sensitivity$prior, density.prior.sensitivity$`P(delta > epsilon)`, 
     ylim = c(0,1),
     xlab = "Prior probability of edge inclusion",
     ylab = "P(delta > 0)",
     pch = 19,
     cex = 0.5,
     cex.axis = 0.6, 
     cex.lab = 0.8)
points(density.prior.sensitivity$prior[9],
       density.prior.sensitivity$`P(delta > epsilon)`[9],
       col = "red3")
title(main = "B) Network density prior sensitivity", 
      cex.main = 0.8, 
      adj = 0, 
      font.main = 2,
      line = 0.5)
lines(density.prior.sensitivity$prior, density.prior.sensitivity$`P(delta > epsilon)`)

dev.off()



# ----------------- now investigation of node degrees ---------------------- #


# estimate degree distributions
ee.degree <- fuzzy.degree.distributions(network = ee.fuzzy$fuzzy.net)
std.degree <- fuzzy.degree.distributions(network = std.fuzzy$fuzzy.net)


# list of priors for node degree prior sensitivity analysis
priors <- seq(0.1,0.5,0.05)

# perform prior sensitivity analysis for each node
# the prior of 0.5 will correspond to the main networks shown
degree.prior.sensitivity <- fuzzy.degree.prior.sensitivity(p.value.matrix1 = ee.net$p.values,
                                                           p.value.matrix2 = std.net$p.values,
                                                           priors = priors,
                                                           epsilon = 0)
# collate into single data frame
degree.priors <- rbind(degree.prior.sensitivity[[1]],
                       degree.prior.sensitivity[[2]],
                       degree.prior.sensitivity[[3]],
                       degree.prior.sensitivity[[4]],
                       degree.prior.sensitivity[[5]],
                       degree.prior.sensitivity[[6]],
                       degree.prior.sensitivity[[7]],
                       degree.prior.sensitivity[[8]],
                       degree.prior.sensitivity[[9]])

# order data frame so plots shown in same order
# as nodes are displayed in all the networks
node_order <- c("diving bottom",
                "opercular beat rate",
                "social density",
                "activity",
                "aggression",
                "foraging",
                "wall swimming",
                "vertical distribution")

degree.priors$node <- factor(degree.priors$node,
                             levels = node_order)

degree.priors <- degree.priors[order(degree.priors$node), ]


# use this index to put red circle around prior that corresponds
# to network visualised in main text
idx <- which(degree.priors$prior == 0.5)


# construct figure for node degrees over range of priors
figure4.5 <- ggplot(degree.priors, aes(x = prior, y = `P(delta > epsilon)`)) +
                facet_wrap(~ node) +
                geom_point() +
                geom_line() +
                scale_y_continuous(limits = c(0, 1),
                  breaks = seq(0, 1, by = 0.2)) +
                xlab(label = "Prior probability of edge inclusion") +
                ylab(label = "P(delta > 0)") +
                geom_point(data = degree.priors[idx,], aes(x = prior, y = `P(delta > epsilon)`), 
                  shape = 1,
                  size = 3,
                  colour = "red3") + 
              theme_classic() 

figure4.5

# save figure 4.5
ggsave("Figure 4.5.pdf",  
       plot = figure4.5,
       width = 8, 
       height = 8, 
       units = "in",
       dpi = 300)


# in 0.5 model, order by P, these Ps are reported in main text
idx <- which(degree.priors$prior == 0.5)
degree.0.5 <- degree.priors[idx,]
degree.0.5[order(degree.0.5$`P(delta > epsilon)`, decreasing = TRUE),]


#---------------- fuzzy differential network --------------------#


# permutation test of network characteristics - this will use the saturated dcor 
# network in the ee.net and std.net objects. 
# With 100,000 perms takes ~40minutes or so on my laptop. Run it for less to get very similar results
# used to construct the fuzzy differential network
net.compare <- compare.zero.order.networks(data1 = ee.data, 
                                           data2 = std.data, 
                                           network1 = ee.net, 
                                           network2 = std.net,
                                           permutations = 100000,
                                           adjust = "none")

# construct fuzzy differential network
differential.fuzzy <- estimate.fuzzy.network(p.value.matrix = net.compare$differential.p.values, 
                                             prior = 0.25)

# construct network density distribution for differential network
differential.density <- fuzzy.density.distribution(network = differential.fuzzy$fuzzy.network)

# construct network density distribution for differential network under prior
differential.density.prior <- fuzzy.density.distribution(network = differential.fuzzy$priors)

# list of priors
priors <- seq(0.1,0.5,0.05)

# perform differential network density prior sensitivity analysis
differential.density.prior.sensitivity <- fuzzy.differential.density.prior.sensitivity(p.value.matrix = net.compare$differential.p.values,
                                                                                       priors = priors)

differential.density.prior.sensitivity


# construct degree distributions for differential network
differential.degrees <- fuzzy.degree.distributions(network = differential.fuzzy$fuzzy.network)

# perform sensitivity analysis for node degree distributions
# of fuzzy differential network
differential.degree.priors <- fuzzy.differential.degree.prior.sensitivity(p.value.matrix = net.compare$differential.p.values,
                                                                          priors = priors)


differential.degree.priors

# colours for differential network - we use a different colour scheme
# to highlight that the network is showing something different
node.number <- ncol(ee.data)
node.names <- colnames(ee.data)
differential.params <- net.compare$edge.differences$Difference
differential.weighted.saturated <- matrix(0, nrow = node.number, ncol = node.number)
dimnames(differential.weighted.saturated) <- list(node.names, colnames(node.names))
differential.weighted.saturated[lower.tri(differential.weighted.saturated)] <- differential.params
differential.weighted.saturated[upper.tri(differential.weighted.saturated)] <- t(differential.weighted.saturated)[upper.tri(differential.weighted.saturated)]

purple.col <- which(differential.weighted.saturated > 0)
organge.col <- which(differential.weighted.saturated < 0)

diff.edge.col <- matrix("white", nrow=ncol(ee.data), ncol=ncol(ee.data))
diff.edge.col[purple.col] <- "purple3"
diff.edge.col[organge.col] <- "darkorange2"

# colour set for plotting of node degree distributions
cols <- RColorBrewer::brewer.pal(name = "Set3", n = 12)


# plot for figure 4.6 - kind of a horrible hacky base R graphics way of producing this figure!
pdf("Figure 4.6.pdf", width = 7.5, height = 10)

par(mfrow = c(3,2),
    mar = c(2.5, 3, 2, 1),      # bottom, left, top, right
    mgp = c(1.5, 0.2, 0),     # axis title, labels, axis line
    tcl = -0.2 )


qgraph(differential.fuzzy$fuzzy.network, 
       theme = 'colorblind',
       labels = colnames(ee.data),
       edge.width = 0.5,
       maximum=1, 
       colFactor = 2, 
       edge.labels=TRUE, 
       edge.label.cex=1.5, 
       label.scale=FALSE, 
       label.cex=1,
       edge.color = diff.edge.col)
title(main ="A) Fuzzy differential network", font=2, cex.main=1.2, adj = 0, line = 1)


plot(differential.density$density.pmf$density, 
     differential.density$density.pmf$probability, 
     xlab = "Network density", 
     ylab = "Probability",
     xlim = c(0,28),
     ylim = c(0,0.2),
     cex = 0.8,
     pch = 19,
     cex.axis = 1, 
     cex.lab = 1)
lines(differential.density$density.pmf$density, 
      differential.density$density.pmf$probability, 
      lwd = 1)
points(differential.density.prior$density.pmf$density, 
       differential.density.prior$density.pmf$probability, 
       pch = 1, 
       cex = 0.8)
lines(differential.density.prior$density.pmf$density, 
      differential.density.prior$density.pmf$probability, 
      lty = 2)
title(main ="B) Differential network density", font=2, cex.main = 1.2, adj = 0, line = 1)
abline(v = sum(differential.fuzzy$fuzzy.network)/2, lty = 2)
abline(v = sum(differential.fuzzy$priors)/2, lty = 2)


plot(0:7, differential.degrees$degree.cdf[1,-1], 
     col = cols[1],
     xlab = "Node Degree", 
     ylab = "Probability",
     xlim = c(0,7),
     ylim = c(0,1),
     cex = 0.8,
     pch = 19,
     cex.axis = 1, 
     cex.lab = 1)
lines(0:7, differential.degrees$degree.cdf[1,-1], lwd = 1, col = cols[1])
title(main ="C) Differential network node degrees", font=2, cex.main = 1.2, adj = 0, line = 1)

points(0:7, differential.degrees$degree.cdf[2,-1], col = cols[3], pch = 19, cex = 0.8)
lines(0:7, differential.degrees$degree.cdf[2,-1], lwd = 1, col = cols[3])

points(0:7, differential.degrees$degree.cdf[3,-1], col = cols[4], pch = 19, cex = 0.8)
lines(0:7, differential.degrees$degree.cdf[3,-1], lwd = 1, col = cols[4])

points(0:7, differential.degrees$degree.cdf[4,-1], col = cols[5], pch = 19, cex = 0.8)
lines(0:7, differential.degrees$degree.cdf[4,-1], lwd = 1, col = cols[5])

points(0:7, differential.degrees$degree.cdf[5,-1], col = cols[6], pch = 19, cex = 0.8)
lines(0:7, differential.degrees$degree.cdf[5,-1], lwd = 1, col = cols[6])

points(0:7, differential.degrees$degree.cdf[6,-1], col = cols[7], pch = 19, cex = )
lines(0:7, differential.degrees$degree.cdf[6,-1], lwd = 1, col = cols[7])

points(0:7, differential.degrees$degree.cdf[7,-1], col = cols[8], pch = 19, cex = 0.8)
lines(0:7, differential.degrees$degree.cdf[7,-1], lwd = 1, col = cols[8])

points(0:7, differential.degrees$degree.cdf[8,-1], col = cols[9], pch = 19, cex = 0.8)
lines(0:7, differential.degrees$degree.cdf[8,-1], lwd = 1, col = cols[9])

legend("bottomright", 
       legend = c("Diving bottom",
                  "Opercular beat rate",
                  "Social density",
                  "Activity",
                  "Aggressive behaviour",
                  "Foraging",
                  "Wall swimming",
                  "Vertical distribution"), 
       col = cols[-2],
       lty = 1,
       lwd = 2,
       cex = 0.9)



plot(differential.density.prior.sensitivity$prior, 
     differential.density.prior.sensitivity$mean, 
     xlab = "Prior probability of edge inclusion", 
     ylab = "Expected value of network density",
     ylim = c(0,20),
     cex = 0.8,
     pch = 19,
     cex.axis = 1, 
     cex.lab = 1)
points(differential.density.prior.sensitivity$prior[4],
       differential.density.prior.sensitivity$mean[4],
       col = "red3", cex = 1.5)
lines(differential.density.prior.sensitivity$prior, 
      differential.density.prior.sensitivity$mean,
      lwd = 1)
title(main ="D) Expected network density prior sensitivity", font=2, cex.main = 1.2, adj = 0, line = 1)



plot(differential.degree.priors$prior, 
     differential.degree.priors$`diving bottom`, 
     col = cols[1],
     xlab = "Prior probability of edge inclusion", 
     ylab = "Expected value of node degree",
     xlim = c(0.1,0.5),
     ylim = c(0,6),
     cex = 0.8,
     pch = 19,
     cex.axis = 1, 
     cex.lab = 1)
lines(differential.degree.priors$prior, 
      differential.degree.priors$`diving bottom`, 
      col = cols[1],, lwd = 1)
points(differential.degree.priors$prior[4],
       differential.degree.priors$`diving bottom`[4],
       col = "red3", cex = 1.5)
title(main ="E) Expected node degrees prior sensitivity", font=2, cex.main = 1.2, adj = 0, line = 1)

points(differential.degree.priors$prior, 
       differential.degree.priors$`opercular beat rate`,
       col = cols[3], pch = 19, cex = 0.8)
lines(differential.degree.priors$prior, 
      differential.degree.priors$`opercular beat rate`,
      lwd = 1, col = cols[3])
points(differential.degree.priors$prior[4],
       differential.degree.priors$`opercular beat rate`[4],
       col = "red3", cex = 1.2)

points(differential.degree.priors$prior, 
       differential.degree.priors$`social density`,
       col = cols[4], pch = 19, cex = 0.8)
lines(differential.degree.priors$prior, 
      differential.degree.priors$`social density`,
      lwd = 1, col = cols[4])
points(differential.degree.priors$prior[4],
       differential.degree.priors$`social density`[4],
       col = "red3", cex = 1.5)


points(differential.degree.priors$prior, 
       differential.degree.priors$activity,
       col = cols[5], pch = 19, cex = 0.8)
lines(differential.degree.priors$prior, 
      differential.degree.priors$activity,
      lwd = 1, col = cols[5])
points(differential.degree.priors$prior[4],
       differential.degree.priors$activity[4],
       col = "red3", cex = 1.5)


points(differential.degree.priors$prior, 
       differential.degree.priors$aggression,
       col = cols[6], pch = 19, cex = 0.8)
lines(differential.degree.priors$prior, 
      differential.degree.priors$aggression,
      lwd = 1, col = cols[6])
points(differential.degree.priors$prior[4],
       differential.degree.priors$aggression[4],
       col = "red3", cex = 1.5)


points(differential.degree.priors$prior, 
       differential.degree.priors$foraging,
       col = cols[7], pch = 19, cex = 0.8)
lines(differential.degree.priors$prior, 
      differential.degree.priors$foraging,
      lwd = 1, col = cols[7])
points(differential.degree.priors$prior[4],
       differential.degree.priors$aggression[4],
       col = "red3", cex = 1.5)

points(differential.degree.priors$prior, 
       differential.degree.priors$`wall swimming`,
       col = cols[8], pch = 19, cex = 0.8)
lines(differential.degree.priors$prior, 
      differential.degree.priors$`wall swimming`,
      lwd = 1, col = cols[8])
points(differential.degree.priors$prior[4],
       differential.degree.priors$`wall swimming`[4],
       col = "red3", cex = 1.5)

points(differential.degree.priors$prior, 
       differential.degree.priors$`vertical distribution`,
       col = cols[9], pch = 19, cex = 0.8)
lines(differential.degree.priors$prior, 
      differential.degree.priors$`vertical distribution`,
      lwd = 1, col = cols[9])
points(differential.degree.priors$prior[4],
       differential.degree.priors$`vertical distribution`[4],
       col = "red3", cex = 1.5)

legend("bottomright", 
       legend = c("Diving bottom",
                  "Opercular beat rate",
                  "Social density",
                  "Activity",
                  "Aggressive behaviour",
                  "Foraging",
                  "Wall swimming",
                  "Vertical distribution"), 
       col = cols[-2],
       lty = 1,
       lwd = 2,
       cex = 0.9)


dev.off()

# Here we compute some useful probabilities.....

# P(activity degree 4 or greater)
1- differential.degrees$degree.cdf[4,5]

# compared to other nodes
1- differential.degrees$degree.cdf[1,5] # diving
1- differential.degrees$degree.cdf[2,5] # obr
1- differential.degrees$degree.cdf[3,5] # social
1- differential.degrees$degree.cdf[5,5] # agg
1- differential.degrees$degree.cdf[6,5] # foraging
1- differential.degrees$degree.cdf[7,5] # wall swimming
1- differential.degrees$degree.cdf[8,5] # vertical



##########################################
#                                        #
#    Network intervention analysis       #
#                                        #  
##########################################


# data preparation, combine housing condition into one dataset 
# where condition is binary variable
# remove activity and social density as these create problems for
# the NIA model (not shown but estimated probabilities of inclusion across run
# did not lead to extremely high correlations)
nia.data <- rbind(ee.data[,-c(3,4)], std.data[,-c(3,4)])   

# due to use of rbind we can now code housing condition with the following line
nia.data$`housing condition` <- c(rep(2,20), rep(1,20)) #so 2 = enriched and 1 = unenriched

# fit model using values selected based on 
# probability calibration simulation
nia.net <- BMGM::bmgm(X = na.omit(nia.data), 
                      type = c(rep("c",6), "m"), 
                      nburn = 10000, 
                      nsample = 50000, 
                      v_0 = 0.16, 
                      v_1 = 1, 
                      pi_beta = 2 / (7 - 1), 
                      context_spec = FALSE, 
                      cont = FALSE)

# this is a convergence check that the posterior probabilities
# are stable across separate runs
nia.net2 <- BMGM::bmgm(X = na.omit(nia.data), 
                       type = c(rep("c",6), "m"), 
                       nburn = 10000, 
                       nsample = 50000, 
                       v_0 = 0.16, 
                       v_1 = 1, 
                       pi_beta = 2 / (7 - 1), 
                       context_spec = FALSE, 
                       cont = FALSE)

# looks very good
par(mfrow = c(1,1))
plot(nia.net$inclusion_probs[lower.tri(nia.net$inclusion_probs)], 
     nia.net2$inclusion_probs[lower.tri(nia.net2$inclusion_probs)],
     xlim = c(0,1), 
     ylim = c(0,1),
     xlab = "Chain 1",
     ylab = "Chain 2")

# near perfect
cor(nia.net$inclusion_probs[lower.tri(nia.net$inclusion_probs)], 
     nia.net2$inclusion_probs[lower.tri(nia.net2$inclusion_probs)])

# identifies median probability model and sets colours
# for visualisation such that blue edges = median probability model
blue.idx <- which(nia.net$inclusion_probs > 0.5)
edge.col <- matrix("grey",nrow = 7, ncol = 7)
edge.col[blue.idx] <- "blue"

# for plotting
shapes <- c(rep("circle",6), "square")
colour <- c(rep("white",6), "lightblue")

# save this later but just view it quickly here
qgraph::qgraph(nia.net$inclusion_probs, 
               layout='circle', 
               edge.labels = TRUE, 
               edge.color = edge.col, 
               labels=colnames(nia.data), 
               shape = shapes, 
               color = colour, 
               vsize=5, 
               maximum = 1,
               colFactor = 1,
               label.font=c(rep(1,8),4),
               label.scale=FALSE, 
               label.cex=0.7, 
               edge.label.cex = 0.8)



#------------- BMGM prior sensitivity analysis ------------#

# get grid of values we want to include in 
# prior sensitivity analysis
grid <- expand.grid(
  v0 = seq(0.13, 0.20, 0.01),
  v1 = 1,
  pi = 2 / (7 - 1)
)

# create vector to store results
results <- vector("list", nrow(grid))

# open progress bar
pb <- txtProgressBar(min = 0, max = nrow(grid), style = 3)

# iterate through grid of priors and save edge inclusion probs
for (i in seq_len(nrow(grid))) {
  
  mod <- bmgm(
    X = na.omit(nia.data),
    type = c(rep("c", 6),"m"),
    nburn = 10000,
    nsample = 50000,
    v_0 = grid$v0[i],
    v_1 = grid$v1[i],
    pi_beta = grid$pi[i],
    context_spec = FALSE,
    bfdr = 0.05,
    cont = FALSE,
    verbose = FALSE
  )
  
  results[[i]] <- mod$inclusion_probs[lower.tri(mod$inclusion_probs)]
  setTxtProgressBar(pb, i)
}

close(pb)

# collate results of NIA prior sensitivity analysis
g <- do.call(rbind, results)
edge.pairs <- as.character(combn(colnames(nia.data), 2, FUN = function(x) paste(x, collapse = "-"), simplify = FALSE))
colnames(g) <- edge.pairs

# summarise across prior settings
edge_summary <- data.frame(edge = edge.pairs,
                           median = apply(g, 2, median, na.rm = TRUE),
                           min = apply(g, 2, min, na.rm = TRUE),
                           max = apply(g, 2, max, na.rm = TRUE))

# NIA BMGM prior sensitivity plot
bmgm.prior.plot <- ggplot(edge_summary,
                          aes(x = fct_reorder(edge, median, .desc = FALSE), y = median)) +
  geom_linerange(aes(ymin = min,ymax = max),
                 linewidth = 0.4) +
  geom_point(size = 0.6) +
  geom_hline(yintercept = 0.50, 
             linetype = "dotted") +
  scale_y_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.1)) +
  coord_flip() +
  labs(x = "",
       y = "Posterior Probability of Edge Inclusion",
       title = "B) Prior sensitivity analysis") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.title.position = "plot",
        plot.title = element_text(
          hjust = 0.25,
          size = 8,
          face = "bold"))

bmgm.prior.plot

#---------- figure 7 plot ---------------#

# plot figure 3 in main text
# its possible that if re-run the selected edges on the border of 0.5
# may be slightly different e.g. the housing condition - opercular beat rate
# edge may be 0.49, 0.50, or 0.51. So it could show up as grey or blue
# without using a set seed, but there are no huge differences in any
# estimated probability so the substantive conclusions remain the same as
# in the main text

pdf("Figure 4.7A.pdf", width=3.5, height=3)

par(mar=c(0, 0, 0, 0), pty = "s")


qgraph::qgraph(nia.net$inclusion_probs, 
               layout='circle', 
               edge.labels = TRUE, 
               edge.color = edge.col, 
               labels=colnames(nia.data),
               shape = shapes, 
               color = colour, 
               vsize=5, 
               maximum = 1,
               colFactor = 1,
               label.font=c(rep(3,6),4),
               label.scale=FALSE, 
               label.cex=0.6, 
               node.label.offset = c(0.5,1),
               edge.label.cex = 1)
text(x=-0.25, y=1.2, labels="A) Network Intervention Analysis", font=2, cex=0.7)

dev.off()

ggsave("Figure 4.7B.pdf", 
       plot = bmgm.prior.plot, 
       width = 3.5, 
       height = 3, 
       units = "in")



#-----------------------------------------------------#
#
#                 NIA robustness analysis
#
#-----------------------------------------------------#

#------------ BDgraph bayesian copula ---------------#

# the BDgraph algorithm seems to give quite unstable results
# e.g. some EIPs change by ~ 0.2 so many edges are quite unstable
# when selecting the median probability model. So probably not
# very reliable inference
# e.g. ......

library(BDgraph)

mods <- list()

for(i in 1:10) {
  
  mod <- bdgraph(data = nia.data,
                       method = "gcgm",
                       algorithm = "bdmcmc",
                       iter = 150000,
                       burnin = 50000,
                       g.prior = 0.5,
                       df.prior = 3,
                       save = TRUE)
  
  p <- summary(mod, vis = FALSE)$p_links
  p[lower.tri(p)] <- t(p)[lower.tri(p)]
  mods[[i]] <- p
  
}

par(mfrow = c(2,5))

for(i in 1:10) {
  
  qgraph::qgraph(mods[[i]], theme = 'colorblind', layout='circle', edge.labels = TRUE,
                 labels=colnames(nia.data), shape = shapes, color = colour, vsize=5, maximum = 1,
                 colFactor = 1, threshold = 0.5,
                 label.font=c(rep(1,8),4),label.scale=FALSE, label.cex=0.7, edge.label.cex = 0.8)
  
}




#------------ BGGM copula ---------------#

library(BGGM)

# assume relatively small effect sizes and set prior to 0.25
# implies this distribution of pcors
BGGM::plot_prior(prior_sd = 0.25, iter = 5000)

bggm.gcgm <- BGGM::explore(Y = as.matrix(nia.data), iter = 50000, type = "mixed", prior_sd = 0.25, impute = TRUE)
bggm.select <- BGGM::select(bggm.gcgm, BF_cut = 1, alternative = "two.sided")

# function to convert the BF matrix from BGGM into posterior probs
posterior_prob <- function(BF10, prior_H1 = 0.5) {
  
  prior_H0 <- 1 - prior_H1
  (BF10 * prior_H1) / (BF10 * prior_H1 + prior_H0)
  
}

# fuzzy network where edge weights are the posterior probs of inclusion
bgcgm.probs <- posterior_prob(bggm.select$BF_10, prior_H1 = 0.5)

# get median prob model structures
bmgm.mpm <- matrix(0,nrow = 7, ncol = 7)
bmgm.mpm[which(nia.net$inclusion_probs > 0.5)] <- 1

bgcgm.mpm <- matrix(0,nrow = 7, ncol = 7)
bgcgm.mpm[which(bgcgm.probs > 0.5)] <- 1

unweighted.density(bmgm.mpm) # 11 edges in bmgm
compare.sim.networks(bmgm.mpm, bgcgm.mpm, metric = "true_positives") # 9 of the same edges
compare.sim.networks(bmgm.mpm, bgcgm.mpm, metric = "false_positives") # 1 edge in bgcgm not in bmgm (condition- diving)
compare.sim.networks(bmgm.mpm, bgcgm.mpm, metric = "false_negatives") # 2 edges in bmgm not in bgcgm (diving - obr & diving - wall swimming)


# plot for appendix
pdf("Figure C2.pdf", width = 8, height = 7)

  # View BMGM model alongside BGGM_GCGM model
  par(mfrow = c(1,2))
  qgraph::qgraph(nia.net$inclusion_probs, 
                 layout='circle', 
                 edge.labels = TRUE, 
                 edge.color = edge.col, 
                 labels=colnames(nia.data), 
                 shape = shapes, 
                 color = colour, 
                 vsize=5, 
                maximum = 1,
                colFactor = 1,
                label.font=c(rep(1,8),4),
                label.scale=FALSE, 
                label.cex=0.7, 
                edge.label.cex = 1.2,
                title = "A) BMGM")

  # for visualisation such that blue edges = median probability model
  bggm.blue.idx <- which(bgcgm.probs > 0.5)
  bggm.edge.col <- matrix("grey",nrow = 7, ncol = 7)
  bggm.edge.col[bggm.blue.idx] <- "blue"

  qgraph::qgraph(bgcgm.probs, 
                 layout='circle', 
                 edge.labels = TRUE, 
                 edge.color = bggm.edge.col, 
                 labels=colnames(nia.data), 
                 shape = shapes, 
                 color = colour, 
                 vsize=5, 
                 maximum = 1,
                 colFactor = 1,
                 label.font=c(rep(1,8),4),
                 label.scale=FALSE, 
                 label.cex=0.7, 
                 edge.label.cex = 1.2,
                 title = "B) BGCGM")
  
dev.off()




# reasonably similar agreement in structure between BMGM
# and BGCGM, but BMGM tends to have more exagerated probabilities
# in either direction

# what is the correlation?
par(mfrow = c(1,1))
plot(nia.net$inclusion_probs[lower.tri(nia.net$inclusion_probs)],
     bgcgm.probs[lower.tri(bgcgm.probs)],
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "BMGM",
     ylab = "Copula")

cor(nia.net$inclusion_probs[lower.tri(nia.net$inclusion_probs)],
     bgcgm.probs[lower.tri(bgcgm.probs)], method = "pearson")


# posterior predictive distribution - which we cant yet get for model checking the BMGM     
bgcgm.ppd <- BGGM::posterior_predict(bggm.gcgm, iter = 5000)


#due to missing data, create empirical set with listwise deletion, and delete corresponding indices from PPD
ind <- which(is.na(nia.data), arr.ind=TRUE)
ind <- ind[,1]
ind <- unique(ind)
nia.data2 <- nia.data[-ind, ]
yrep <- bgcgm.ppd[-ind, , ]  

library(bayesplot)

#create ppc_bar plots for each variable for qualitative PPC
ppc_dens_overlay(y = nia.data2[,1], t(yrep[,1,1:500]), size=0.5) + ggtitle("diving") + theme(plot.title = element_text(hjust = 0.5, size=10))
ppc_dens_overlay(y = nia.data2[,2], t(yrep[,2,1:500]), size=0.5) + ggtitle("opercular beat rate") + theme(plot.title = element_text(hjust = 0.5, size=10))
ppc_dens_overlay(y = nia.data2[,3], t(yrep[,3,1:500]), size=0.5) + ggtitle("aggression") + theme(plot.title = element_text(hjust = 0.5, size=10))
ppc_dens_overlay(y = nia.data2[,4], t(yrep[,4,1:500]), size=0.5) + ggtitle("foraging") + theme(plot.title = element_text(hjust = 0.5, size=10))
ppc_dens_overlay(y = nia.data2[,5], t(yrep[,5,1:500]), size=0.5) + ggtitle("wall swimming") + theme(plot.title = element_text(hjust = 0.5, size=10))
ppc_dens_overlay(y = nia.data2[,6], t(yrep[,6,1:500]), size=0.5) + ggtitle("vertical distribution") + theme(plot.title = element_text(hjust = 0.5, size=10))
ppc_bars(y = nia.data2[,7], t(yrep[,7,1:500]), size=0.5) + ggtitle("housing condition") + theme(plot.title = element_text(hjust = 0.5, size=10))


# ppc for skew
library(e1071) # for skewness function

ppc.skew.list <- list()

for(i in 1:ncol(nia.data2)) {
  
  ppc.skew.list[[i]] <- ppc_stat(y = nia.data2[,i], yrep=t(yrep[,i,1:5000]), stat=skewness)
  ppc.skew.list[[i]]  <- ppc.skew.list[[i]] + 
                       ggtitle(colnames(nia.data2)[i]) + 
                       theme(plot.title = element_text(hjust = 0.5, size=8))
}

for(i in 1:ncol(nia.data2)) {
  
  print(ppc.skew.list[[i]])
}


# as one plot
ggarrange(ppc.skew.list[[1]],ppc.skew.list[[2]],ppc.skew.list[[3]],ppc.skew.list[[4]],ppc.skew.list[[5]],
          ppc.skew.list[[6]],ppc.skew.list[[7]], ncol=2, nrow=4, align="h", common.legend = TRUE)


#-------- prior sensitivity for BGCGM ---------------#

priors = seq(0.1,0.5,0.05)

bgcgm.results <- vector("list", length(priors))

pb <- txtProgressBar(min = 0, max = length(priors), style = 3)

for (i in 1:length(priors)) {
  
  mod <- BGGM::explore(Y = as.matrix(nia.data), 
                       iter = 50000, 
                       type = "mixed", 
                       prior_sd = priors[i], 
                       impute = TRUE,
                       progress = FALSE)
  g <- BGGM::select(mod, BF_cut = 1, alternative = "two.sided")
  p <- posterior_prob(g$BF_10, prior_H1 = 0.5)
  
  bgcgm.results[[i]] <- p[lower.tri(p)]
  
  setTxtProgressBar(pb, i)
  
}

close(pb)

bgcgm.g <- do.call(rbind, bgcgm.results)
colnames(bgcgm.g) <- edge.pairs

# summarise across prior settings
bgcgm.edge.summary <- data.frame(
  edge   = edge.pairs,
  median = apply(bgcgm.g, 2, median, na.rm = TRUE),
  min    = apply(bgcgm.g, 2, min, na.rm = TRUE),
  max    = apply(bgcgm.g, 2, max, na.rm = TRUE)
)


# NIA BGCGM prior sensitivity plot
bgcgm.prior.plot <- ggplot(bgcgm.edge.summary,
                          aes(x = fct_reorder(edge, median, .desc = FALSE), y = median)) +
  geom_linerange(aes(ymin = min,ymax = max),
                 linewidth = 0.4) +
  geom_point(size = 0.6) +
  geom_hline(yintercept = 0.50, 
             linetype = "dotted") +
  scale_y_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.1)) +
  coord_flip() +
  labs(x = "",
       y = "Posterior Probability of Edge Inclusion",
       title = "") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.title = element_text(size = 6),
        plot.title.position = "plot",
        plot.title = element_text(
          hjust = 0.25,
          size = 8,
          face = "bold"))

bgcgm.prior.plot

ggsave(filename ="Figure C3.pdf", 
       plot = bgcgm.prior.plot, 
       width = 4.5, 
       height = 4)


#----------- frequentist regularized MGM -----------------------#

# first use npn transformation on continuous variables to relax gaussian assumption
npn.data <- npn.transform(nia.data[,1:6])
npn.data$`housing condition` <- nia.data$`housing condition`
colnames(npn.data) <- c("diving.bottom","opercular.beat.rate",
                        "aggression","foraging",
                        "wall.swimming","vertical.distribution","housing.condition")


# MGM with npn transform on continuous variables
type <- c(rep("g", 6), "c")
levels <- c(rep(1,6),2)

nia.network.mgm <- mgm::mgm(as.matrix(na.omit(nia.data)), type=type, level=levels, lambdaSel="EBIC", 
                            ruleReg="OR", lambdaGam = 0)
nia.network.mgm$structure <- abs(sign(nia.network.mgm$pairwise$wadj))

par(mar=c(0, 0, 0, 0))

qgraph::qgraph(nia.network.mgm$structure, 
               layout='circle',  
               labels=colnames(nia.data), 
               shape = shapes, 
               color = colour, 
               vsize=5, 
               label.font=c(rep(1,6),4),
               label.scale=FALSE, 
               label.cex=0.5, 
               title = "FMGM")

# bootstrap stability
mgm.boot <- mgm::resample(data = as.matrix(na.omit(npn.data)), nia.network.mgm, nB=1000)


mgm::plotRes(object = mgm.boot, 
             labels = colnames(npn.data),
             quantiles = c(0.05,0.95),
             axis.ticks = seq(-0.5,2,0.5),
             decreasing = TRUE,
             cex.label = 0.5,
             layout.width.labels = 0.5,
             cex.mean = 0.3,
             cex.bg = 2)


#------------- additional figures -------------------#

# these figures basically highlight the 'threshold' or
# lower 'degrees of freedom' effect mentioned in the
# discussion e.g. if the unenriched fish were already highly 
# active, doing lots of wall swimming, and regularly entering 
# the top third of the tank, then any aggressive behaviour that 
# occurred had little effect on these behaviours, 
# and thus dependencies diminished

network.data$condition[which(network.data$condition == "ee")] <- "enriched"
network.data$condition[which(network.data$condition == "std")] <- "unenriched"
network.data$condition <- as.factor(network.data$condition)

#this plot is the one shown in supplementary material figure S14
ggplot(network.data, aes(x=aggression, y=activity, color=condition)) +
  geom_point() +
  scale_color_manual(values = c("enriched" = "green", "unenriched" = "red")) +
  labs(x = "Agressive behaviour", y = "Activity") +
  coord_cartesian(xlim=c(0,0.9), ylim=c(0,1)) +
  theme_classic() + 
  theme(axis.title = element_text(size=10), 
        plot.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10),
        legend.key.size = unit(5, "lines"),
        axis.line.x = element_line(linewidth=0.5),
        axis.line.y = element_line(linewidth=0.5),
        axis.ticks = element_line(linewidth=0.5),
        plot.margin = margin(1, 1, 1, 1),
        legend.margin = margin(t = 1, r = 1, b = 1, l = 1))

ggplot(network.data, aes(x=aggression, y=`vertical distribution`, color=condition)) +
  geom_point() +
  scale_color_manual(values = c("enriched" = "green", "unenriched" = "red")) +
  labs(x = "aggression", y = "vertical distribution") +
  theme_classic() + 
  theme(axis.title = element_text(size=10), 
        plot.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10),
        legend.key.size = unit(5, "lines"),
        axis.line.x = element_line(linewidth=0.5),
        axis.line.y = element_line(linewidth=0.5),
        axis.ticks = element_line(linewidth=0.5),
        plot.margin = margin(1, 1, 1, 1),
        legend.margin = margin(t = 1, r = 1, b = 1, l = 1))


ggplot(network.data, aes(x=activity, y=`vertical distribution`, color=condition)) +
  geom_point() +
  scale_color_manual(values = c("enriched" = "green", "unenriched" = "red")) +
  labs(x = "activity", y = "vertical distribution") +
  theme_classic() + 
  theme(axis.title = element_text(size=10), 
        plot.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10),
        legend.key.size = unit(5, "lines"),
        axis.line.x = element_line(linewidth=0.5),
        axis.line.y = element_line(linewidth=0.5),
        axis.ticks = element_line(linewidth=0.5),
        plot.margin = margin(1, 1, 1, 1),
        legend.margin = margin(t = 1, r = 1, b = 1, l = 1))

ggplot(network.data, aes(x=`wall swimming`, y=`vertical distribution`, color=condition)) +
  geom_point() +
  scale_color_manual(values = c("enriched" = "green", "unenriched" = "red")) +
  labs(x = "wall swimming", y = "vertical distribution") +
  theme_classic() + 
  theme(axis.title = element_text(size=10), 
        plot.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10),
        legend.key.size = unit(5, "lines"),
        axis.line.x = element_line(linewidth=0.5),
        axis.line.y = element_line(linewidth=0.5),
        axis.ticks = element_line(linewidth=0.5),
        plot.margin = margin(1, 1, 1, 1),
        legend.margin = margin(t = 1, r = 1, b = 1, l = 1))



ggplot(network.data, aes(x=activity, y=`wall swimming`, color=condition)) +
  geom_point() +
  scale_color_manual(values = c("enriched" = "green", "unenriched" = "red")) +
  labs(x = "activity", y = "wall swimming") +
  theme_classic() + 
  theme(axis.title = element_text(size=10), 
        plot.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10),
        legend.key.size = unit(5, "lines"),
        axis.line.x = element_line(linewidth=0.5),
        axis.line.y = element_line(linewidth=0.5),
        axis.ticks = element_line(linewidth=0.5),
        plot.margin = margin(1, 1, 1, 1),
        legend.margin = margin(t = 1, r = 1, b = 1, l = 1))


######################################
#                                    #
#                END                 #
#                                    #  
######################################