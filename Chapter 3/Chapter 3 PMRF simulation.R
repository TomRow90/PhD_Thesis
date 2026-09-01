################################################################################
#
#         Section 3.4.3: PRMF simulation in Chapter 3 of thesis
#
################################################################################

library(qgraph)
library(copula)
library(BDgraph)
library(mgm)
library(ggplot2)

setwd("/Users/tomrowland/Documents/PhD/Thesis write up/Chapter 3")


#marginal probability distributions to transform MVN into
node_margins <- c("beta","norm","exp","pois","gamma",
                     "nbinom")

beta_param_margins <- list(
  list(shape1 = 6, shape2 = 2),
  list(shape1 = 2, shape2 = 6),
  list(shape1 = 5, shape2 = 5)
)

norm_param_margins <- list(
  list(mean = 0, sd = 1),
  list(mean = 0, sd = 1),
  list(mean = 0, sd = 1)
)

exp_param_margins <- list(
  list(rate = 0.1),
  list(rate = 0.2),
  list(rate = 0.3)
)

pois_param_margins <- list(
  list(lambda = 20),
  list(lambda = 10),
  list(lambda = 50)
)

gamma_param_margins <- list(
  list(shape = 5, rate = 1),
  list(shape = 2, rate = 0.5),
  list(shape = 1, rate = 2)
)

nbinom_param_margins <- list(
  list(size = 1, prob = 0.1),
  list(size = 1, prob = 0.2),
  list(size = 1, prob = 0.3)
)


# simulation conditions 
sample_size <- c(20,50,80,100)
node_number <- c(6,10,20)
density <- c(0.3,0.4,0.6)
replications <- 100

# Due to computational burden taking my computer up I run these
# separately for each model so I dont have to run it for so long 
# thereby tying my laptop up!

# random network sim for BGCGM from BDgraph
bgcgm_sim <- parSim(
  
  sample_size = sample_size,
  node_number = node_number,
  density = density,
  replications = replications,
  nCores = 10,
  export = c("node_margins", "beta_param_margins",
             "norm_param_margins", "exp_param_margins",
             "pois_param_margins", "gamma_param_margins",
             "nbinom_param_margins"),
  
  expression = {
    
  # random network  
  true_net <- bootnet::genGGM(Nvar = node_number, 
                              p = density, 
                              parRange = c(0.5,1), 
                              propPositive = 0.5,
                              constant = 1.2, 
                              graph = "random") 
  
  lower_tri <- true_net[lower.tri(true_net)]
  pcor_idx <- which(lower_tri != 0)
  mean_pcor <- mean(abs(lower_tri[pcor_idx]))
  mean_pcor
  
  # add 1 to diagonal for pcor2cor function
  diag(true_net) <- 1
  
  # randomly define distributional form of nodes
  node_distributions <- sample(node_margins, 
                               size = node_number, 
                               replace = TRUE,
                               prob = rep(1/6,6))
  
  # sample from node distributions possible parameterisations
  param_margins <- list()
  
  for(i in 1:length(node_distributions)) {
    
    if(node_distributions[[i]] == "beta") param_margins[i] <- sample(beta_param_margins, 1, prob = rep(1/3,3))
    if(node_distributions[[i]] == "norm") param_margins[i] <- sample(norm_param_margins, 1, prob = rep(1/3,3))  
    if(node_distributions[[i]] == "exp") param_margins[i] <- sample(exp_param_margins, 1, prob = rep(1/3,3))  
    if(node_distributions[[i]] == "pois") param_margins[i] <- sample(pois_param_margins, 1, prob = rep(1/3,3))  
    if(node_distributions[[i]] == "gamma") param_margins[i] <- sample(gamma_param_margins, 1, prob = rep(1/3,3))  
    if(node_distributions[[i]] == "nbinom") param_margins[i] <- sample(nbinom_param_margins, 1, prob = rep(1/3,3))  
    
  }
  
  # transform GGM to zero order cors
  cors <- corpcor::pcor2cor(true_net)
  
  # in vector form for copula package
  cop_params <- cors[lower.tri(cors)]
  
  # create the gaussian copula.
  cop <- copula::normalCopula(param = cop_params, dim = node_number, dispstr = "un")
  
  # simulate from the copula
  u <- copula::rCopula(sample_size, cop)

  # match margin names to quantile function names
  qfuns <- paste0("q", node_distributions)
  
  # generate data
  data <- as.data.frame(mapply(function(qfun, col, params) {
    do.call(qfun, c(list(p = col), params))
  }, qfuns, as.data.frame(u), param_margins))
  
  colnames(data) <- paste0("V", seq_len(ncol(data)))
  
  # set not.cont argument for BDgraph
  not_cont <- rep(NA,node_number)

  for(i in 1:length(node_distributions)) {
    
    if(node_distributions[i] == "beta") not_cont[i] <- 0
    if(node_distributions[i] == "norm") not_cont[i] <- 0
    if(node_distributions[i] == "exp") not_cont[i] <- 0
    if(node_distributions[i] == "pois") not_cont[i] <- 1
    if(node_distributions[i] == "gamma") not_cont[i] <- 0
    if(node_distributions[i] == "nbinom") not_cont[i] <- 1

  }
  
 # estimate network
 est_net <- BDgraph::bdgraph(data = data, 
                             method = "gcgm", 
                             algorithm = "bdmcmc",
                             not.cont = not_cont,
                             iter = 50000, 
                             burnin = 10000, 
                             g.prior = 0.5, 
                             df.prior = 3)
 
  # median probability model edge selection
  est_net <- summary(est_net, vis = FALSE)$selected_g 
  
  # compare estimated to true
  comp <- welfareNet::compare.sim.networks(true_net, est_net, 
                                           metric = c("sensitivity","specificity","precision"))

  # add mean pcor to results object
  comp$mean_pcor <- mean_pcor
  
  return(comp)
  
  }
)


#saveRDS(object = bgcgm_sim, file = "bgcgm_sim.rds")
bgcgm_sim <- readRDS("bgcgm_sim.rds")

bgcgm_sim$density <- as.factor(bgcgm_sim$density)

# summaries for table
aggregate(mean_pcor ~ density + node_number, data = bgcgm_sim, FUN = mean)
aggregate(mean_pcor ~ density + node_number, data = bgcgm_sim, FUN = sd)


aggregate(sensitivity ~ sample_size + density + node_number, data = bgcgm_sim, FUN = mean)
aggregate(sensitivity ~ sample_size + density + node_number, data = bgcgm_sim, FUN = sd)


aggregate(specificity ~ sample_size + density + node_number, data = bgcgm_sim, FUN = mean)
aggregate(specificity ~ sample_size + density + node_number, data = bgcgm_sim, FUN = sd)


aggregate(precision ~ + sample_size + density + node_number, data = bgcgm_sim, FUN = mean)
aggregate(precision ~ + sample_size + density + node_number, data = bgcgm_sim, FUN = sd)

# plot - not shown in chapter
ggplot(bgcgm_sim, aes(x = factor(sample_size), y = sensitivity)) +
  facet_wrap(~ density + node_number) +
  geom_violin(scale = "width") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size=0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "red", size = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "blue", size = 0.5) +
  labs(x = "Sample Size", y = "Sensitivity", 
       title = "") +
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,0.1)) +
  theme_classic()


ggplot(bgcgm_sim, aes(x = factor(sample_size), y = specificity)) +
  facet_wrap(~ density + node_number) +
  geom_violin(scale = "width") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size=0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "red", size = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "blue", size = 0.5) +
  labs(x = "Sample Size", y = "Specificity", 
       title = "") +
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,0.1)) +
  theme_classic()



ggplot(bgcgm_sim, aes(x = factor(sample_size), y = precision)) +
  facet_wrap(~ density + node_number) +
  geom_violin(scale = "width") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size=0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "red", size = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "blue", size = 0.5) +
  labs(x = "Sample Size", y = "Precision", 
       title = "") +
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,0.1)) +
  theme_classic()



#------- now Frequentist MGM simulation --------#

# same as before so redundant
sample_size <- c(20,50,80,100)
node_number <- c(6,10,20)
density <- c(0.3,0.4,0.6)
replications <- 100

# random network sim using FMGM from mgm
fmgm_sim <- parSim(
  
  sample_size = sample_size,
  node_number = node_number,
  density = density,
  replications = replications,
  nCores = 10,
  export = c("node_margins", "beta_param_margins",
             "norm_param_margins", "exp_param_margins",
             "pois_param_margins", "gamma_param_margins",
             "nbinom_param_margins"),
  
  expression = {
    
    # random network
    true_net <- bootnet::genGGM(Nvar = node_number, 
                                p = density, 
                                parRange = c(0.5,1), 
                                propPositive = 0.5,
                                constant = 1.2, 
                                graph = "random") 
    
    lower_tri <- true_net[lower.tri(true_net)]
    pcor_idx <- which(lower_tri != 0)
    mean_pcor <- mean(abs(lower_tri[pcor_idx]))
    mean_pcor
    
    # add 1 to diagonal for pcor2cor function
    diag(true_net) <- 1
    
    # randomly define distributional form of nodes
    node_distributions <- sample(node_margins, 
                                 size = node_number, 
                                 replace = TRUE,
                                 prob = rep(1/6,6))
    
    # sample from node distributions possible parameterisations
    param_margins <- list()
    
    for(i in 1:length(node_distributions)) {
      
      if(node_distributions[[i]] == "beta") param_margins[i] <- sample(beta_param_margins, 1, prob = rep(1/3,3))
      if(node_distributions[[i]] == "norm") param_margins[i] <- sample(norm_param_margins, 1, prob = rep(1/3,3))  
      if(node_distributions[[i]] == "exp") param_margins[i] <- sample(exp_param_margins, 1, prob = rep(1/3,3))  
      if(node_distributions[[i]] == "pois") param_margins[i] <- sample(pois_param_margins, 1, prob = rep(1/3,3))  
      if(node_distributions[[i]] == "gamma") param_margins[i] <- sample(gamma_param_margins, 1, prob = rep(1/3,3))  
      if(node_distributions[[i]] == "nbinom") param_margins[i] <- sample(nbinom_param_margins, 1, prob = rep(1/3,3))  
      
    }
    
    # transform GGM to zero order cors
    cors <- corpcor::pcor2cor(true_net)
    
    # in vector form for copula package
    cop_params <- cors[lower.tri(cors)]
    
    # create gaussian copula.
    cop <- copula::normalCopula(param = cop_params, dim = node_number, dispstr = "un")
    
    # simulate from copula 
    u <- copula::rCopula(sample_size, cop)
    
    # match margin names to quantile function names
    qfuns <- paste0("q", node_distributions)
    
    # generate data
    data <- as.data.frame(mapply(function(qfun, col, params) {
      do.call(qfun, c(list(p = col), params))
    }, qfuns, as.data.frame(u), param_margins))
    
    colnames(data) <- paste0("V", seq_len(ncol(data)))
    
    # set type and levels for fmgm
    type <- rep(NA,node_number)
    level <- rep(NA, node_number)
    
    
    for(i in 1:length(node_distributions)) {
      
      if(node_distributions[i] == "beta") {
        
        type[i] <- "g"
        level[i] <- 1
      }

      if(node_distributions[i] == "norm") {
        
        type[i] <- "g"
        level[i] <- 1
      }

      if(node_distributions[i] == "exp") {
        
        type[i] <- "g" 
        level[i] <- 1
        data[,i] <- log(data[,i]) # transformation we might apply on such data
      }
      
      if(node_distributions[i] == "pois") {
        
        type[i] <- "p"
        level[i] <- 1
      }
      
      if(node_distributions[i] == "gamma") {
        
        type[i] <- "g"
        level[i] <- 1
      }
      
      if(node_distributions[i] == "nbinom") {
        
        type[i] <- "p"
        level[i] <- 1
      }
  
    }
    
    # estimate network
    est_net <- mgm::mgm(data = as.matrix(data), 
                        type = type, 
                        level = level,
                        lambdaSel = "EBIC",
                        ruleReg = "OR",
                        lambdaGam = 0)
    
    # selected structure
    est_net <- sign(abs(est_net$pairwise$wadj))
    
    # compare estimated to true
    comp <- welfareNet::compare.sim.networks(true_net, est_net, 
                                             metric = c("sensitivity","specificity","precision"))
    
    #add mean pcor to output results
    comp$mean_pcor <- mean_pcor
    
    return(comp)
    
  }
)

#saveRDS(object = fmgm_sim, file = "fmgm_sim.rds")
fmgm_sim <- readRDS("fmgm_sim.rds")

# summaries for table
aggregate(mean_pcor ~ density + node_number, data = fmgm_sim, mean)
aggregate(mean_pcor ~ density + node_number, data = fmgm_sim, sd)

aggregate(sensitivity ~ sample_size + density + node_number, data = fmgm_sim, mean)
aggregate(sensitivity ~ sample_size + density + node_number, data = fmgm_sim, sd)

aggregate(specificity ~ sample_size + density + node_number, data = fmgm_sim, mean)
aggregate(specificity ~ sample_size + density + node_number, data = fmgm_sim, sd)


aggregate(precision ~ sample_size + density + node_number, data = fmgm_sim, mean)
aggregate(precision ~ sample_size + density + node_number, data = fmgm_sim, sd)

fmgm_sim$density <- as.factor(fmgm_sim$density)

ggplot(fmgm_sim, aes(x = factor(sample_size), y = mean_pcor)) +
  facet_wrap(~ density + node_number) +
  geom_violin(scale = "width") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size=0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "red", size = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "blue", size = 0.5) +
  labs(x = "Sample Size", y = "mean_pcor", 
       title = "") +
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,0.1)) +
  theme_classic()

ggplot(fmgm_sim, aes(x = factor(sample_size), y = sensitivity)) +
  facet_wrap(~ density + node_number) +
  geom_violin(scale = "width") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size=0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "red", size = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "blue", size = 0.5) +
  labs(x = "Sample Size", y = "Sensitivity", 
       title = "") +
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,0.1)) +
  theme_classic()


ggplot(fmgm_sim, aes(x = factor(sample_size), y = specificity)) +
  facet_wrap(~ density + node_number) +
  geom_violin(scale = "width") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size=0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "red", size = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "blue", size = 0.5) +
  labs(x = "Sample Size", y = "Specificity", 
       title = "") +
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,0.1)) +
  theme_classic()



ggplot(fmgm_sim, aes(x = factor(sample_size), y = precision)) +
  facet_wrap(~ density + node_number) +
  geom_violin(scale = "width") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size=0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "red", size = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "blue", size = 0.5) +
  labs(x = "Sample Size", y = "Precision", 
       title = "") +
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,0.1)) +
  theme_classic()



#----------------------------------------------#
#---- NOW SMALL WORLD NETWORK STRUCTURES ------#
#----------------------------------------------#


sample_size <- c(20,50,80,100)
node_number <- c(6,10,20)
rewire <- 0.5
nei <- c(1,2)
replications <- 100

# smallworld network sim of FMGM using mgm 
fmgm_smallworld_sim <- parSim(
  
  sample_size = sample_size,
  node_number = node_number,
  rewire = rewire,
  nei = nei,
  replications = replications,
  nCores = 10,
  export = c("node_margins", "beta_param_margins",
             "norm_param_margins", "exp_param_margins",
             "pois_param_margins", "gamma_param_margins",
             "nbinom_param_margins"),
  
  expression = {
    
    true_net <- bootnet::genGGM(Nvar = node_number, 
                                p = rewire,
                                nei = nei,
                                parRange = c(0.5,1), 
                                propPositive = 0.5,
                                constant = 1.2, 
                                graph = "smallworld") 
    
    lower_tri <- true_net[lower.tri(true_net)]
    pcor_idx <- which(lower_tri != 0)
    mean_pcor <- mean(abs(lower_tri[pcor_idx]))
    mean_pcor
    
    density <- NetworkToolbox::conn(true_net)$density
    
    # add 1 to diagonal for pcor2cor function
    diag(true_net) <- 1
    
    # randomly define distributional form of nodes
    node_distributions <- sample(node_margins, 
                                 size = node_number, 
                                 replace = TRUE,
                                 prob = rep(1/6,6))
    
    # sample from node distributions possible parameterisations
    param_margins <- list()
    
    for(i in 1:length(node_distributions)) {
      
      if(node_distributions[[i]] == "beta") param_margins[i] <- sample(beta_param_margins, 1, prob = rep(1/3,3))
      if(node_distributions[[i]] == "norm") param_margins[i] <- sample(norm_param_margins, 1, prob = rep(1/3,3))  
      if(node_distributions[[i]] == "exp") param_margins[i] <- sample(exp_param_margins, 1, prob = rep(1/3,3))  
      if(node_distributions[[i]] == "pois") param_margins[i] <- sample(pois_param_margins, 1, prob = rep(1/3,3))  
      if(node_distributions[[i]] == "gamma") param_margins[i] <- sample(gamma_param_margins, 1, prob = rep(1/3,3))  
      if(node_distributions[[i]] == "nbinom") param_margins[i] <- sample(nbinom_param_margins, 1, prob = rep(1/3,3))  
      
    }
    
    # transform GGM to zero order cors
    cors <- corpcor::pcor2cor(true_net)
    
    # in vector form for copula package
    cop_params <- cors[lower.tri(cors)]
    
    # create the gaussian copula.
    cop <- copula::normalCopula(param = cop_params, dim = node_number, dispstr = "un")
    
    # simulate from the copula 
    u <- copula::rCopula(sample_size, cop)
    
    # match margin names to quantile function names
    qfuns <- paste0("q", node_distributions)
    
    # generate data
    data <- as.data.frame(mapply(function(qfun, col, params) {
      do.call(qfun, c(list(p = col), params))
    }, qfuns, as.data.frame(u), param_margins))
    
    colnames(data) <- paste0("V", seq_len(ncol(data)))
    
    # set type and levels for fmgm
    type <- rep(NA,node_number)
    level <- rep(NA, node_number)
    
    
    for(i in 1:length(node_distributions)) {
      
      if(node_distributions[i] == "beta") {
        
        type[i] <- "g"
        level[i] <- 1
      }
      
      if(node_distributions[i] == "norm") {
        
        type[i] <- "g"
        level[i] <- 1
      }
      
      if(node_distributions[i] == "exp") {
        
        type[i] <- "g" 
        level[i] <- 1
        data[,i] <- log(data[,i]) # transformation we might apply on such data
      }
      
      if(node_distributions[i] == "pois") {
        
        type[i] <- "p"
        level[i] <- 1
      }
      
      if(node_distributions[i] == "gamma") {
        
        type[i] <- "g"
        level[i] <- 1
      }
      
      if(node_distributions[i] == "nbinom") {
        
        type[i] <- "p"
        level[i] <- 1
      }
      
    }
    
    # now estimate network
    est_net <- mgm::mgm(data = as.matrix(data), 
                        type = type, 
                        level = level,
                        lambdaSel = "EBIC",
                        ruleReg = "OR",
                        lambdaGam = 0)
    
    # selected structure
    est_net <- sign(abs(est_net$pairwise$wadj))
    
    # compare estimated to true
    comp <- welfareNet::compare.sim.networks(true_net, est_net, 
                                             metric = c("sensitivity","specificity","precision"))
    
    #add mean pcor and density to output results
    comp$mean_pcor <- mean_pcor
    comp$density <- density
    
    return(comp)
    
  }
)

#saveRDS(object = fmgm_smallworld_sim, file = "fmgm_smallworld_sim.rds")
#readRDS("fmgm_smallworld_sim.rds")

# summaries for table
aggregate(density ~ nei + node_number, data = fmgm_smallworld_sim, mean)
aggregate(density ~ nei +  node_number, data = fmgm_smallworld_sim, sd)

aggregate(mean_pcor ~ sample_size + nei + node_number, data = fmgm_smallworld_sim, mean)
round(aggregate(mean_pcor ~ sample_size + nei + node_number, data = fmgm_smallworld_sim, sd),2)

aggregate(sensitivity ~ sample_size + nei + node_number, data = fmgm_smallworld_sim, mean)
aggregate(sensitivity ~ sample_size + nei + node_number, data = fmgm_smallworld_sim, sd)

aggregate(specificity ~ sample_size + nei + node_number, data = fmgm_smallworld_sim, mean)
aggregate(specificity ~ sample_size + nei + node_number, data = fmgm_smallworld_sim, sd)

aggregate(precision ~ sample_size + nei + node_number, data = fmgm_smallworld_sim, mean)
aggregate(precision ~ sample_size + nei + node_number, data = fmgm_smallworld_sim, sd)


# plots - not shown in chapter 
ggplot(fmgm_smallworld_sim, aes(x = factor(sample_size), y = mean_pcor)) +
  facet_wrap(~ nei + node_number) +
  geom_violin(scale = "width") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size=0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "red", size = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "blue", size = 0.5) +
  labs(x = "Sample Size", y = "mean_pcor", 
       title = "") +
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,0.1)) +
  theme_classic()

ggplot(fmgm_smallworld_sim, aes(x = factor(sample_size), y = sensitivity)) +
  facet_wrap(~ nei + node_number) +
  geom_violin(scale = "width") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size=0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "red", size = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "blue", size = 0.5) +
  labs(x = "Sample Size", y = "Sensitivity", 
       title = "") +
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,0.1)) +
  theme_classic()


ggplot(fmgm_smallworld_sim, aes(x = factor(sample_size), y = specificity)) +
  facet_wrap(~ density + node_number) +
  geom_violin(scale = "width") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size=0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "red", size = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "blue", size = 0.5) +
  labs(x = "Sample Size", y = "Specificity", 
       title = "") +
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,0.1)) +
  theme_classic()



#------- now BDgraph copula model --------#


sample_size <- c(20,50,80,100)
node_number <- c(6,10,20)
rewire <- 0.5
nei <- c(1,2)
replications <- 100

# smallworld network sim of BGCGM from bdgraph
bgcgm_smallworld_sim <- parSim(
  
  sample_size = sample_size,
  node_number = node_number,
  rewire = rewire,
  nei = nei,
  replications = replications,
  nCores = 10,
  export = c("node_margins", "beta_param_margins",
             "norm_param_margins", "exp_param_margins",
             "pois_param_margins", "gamma_param_margins",
             "nbinom_param_margins"),
  
  expression = {
    
    true_net <- bootnet::genGGM(Nvar = node_number, 
                                p = rewire,
                                nei = nei,
                                parRange = c(0.5,1), 
                                propPositive = 0.5,
                                constant = 1.2, 
                                graph = "smallworld") 
    
    lower_tri <- true_net[lower.tri(true_net)]
    pcor_idx <- which(lower_tri != 0)
    mean_pcor <- mean(abs(lower_tri[pcor_idx]))
    mean_pcor
    
    dens <- NetworkToolbox::conn(true_net)$density
    
    # add 1 to diagonal for pcor2cor function
    diag(true_net) <- 1
    
    # randomly define distributional form of nodes
    node_distributions <- sample(node_margins, 
                                 size = node_number, 
                                 replace = TRUE,
                                 prob = rep(1/6,6))
    
    # sample from node distributions possible parameterisations
    param_margins <- list()
    
    for(i in 1:length(node_distributions)) {
      
      if(node_distributions[[i]] == "beta") param_margins[i] <- sample(beta_param_margins, 1, prob = rep(1/3,3))
      if(node_distributions[[i]] == "norm") param_margins[i] <- sample(norm_param_margins, 1, prob = rep(1/3,3))  
      if(node_distributions[[i]] == "exp") param_margins[i] <- sample(exp_param_margins, 1, prob = rep(1/3,3))  
      if(node_distributions[[i]] == "pois") param_margins[i] <- sample(pois_param_margins, 1, prob = rep(1/3,3))  
      if(node_distributions[[i]] == "gamma") param_margins[i] <- sample(gamma_param_margins, 1, prob = rep(1/3,3))  
      if(node_distributions[[i]] == "nbinom") param_margins[i] <- sample(nbinom_param_margins, 1, prob = rep(1/3,3))  
      
    }
    
    # transform GGM to zero order cors
    cors <- corpcor::pcor2cor(true_net)
    
    # in vector form for copula package
    cop_params <- cors[lower.tri(cors)]
    
    # create the gaussian copula.
    cop <- copula::normalCopula(param = cop_params, dim = node_number, dispstr = "un")
    
    # simulate from the copula 
    u <- copula::rCopula(sample_size, cop)
    
    # match margin names to quantile function names
    qfuns <- paste0("q", node_distributions)
    
    # generate data
    data <- as.data.frame(mapply(function(qfun, col, params) {
      do.call(qfun, c(list(p = col), params))
    }, qfuns, as.data.frame(u), param_margins))
    
    colnames(data) <- paste0("V", seq_len(ncol(data)))
    
    # set not.cont argument for bdgraph
    not_cont <- rep(NA,node_number)
    
    for(i in 1:length(node_distributions)) {
      
      if(node_distributions[i] == "beta") not_cont[i] <- 0
      if(node_distributions[i] == "norm") not_cont[i] <- 0
      if(node_distributions[i] == "exp") not_cont[i] <- 0
      if(node_distributions[i] == "pois") not_cont[i] <- 1
      if(node_distributions[i] == "gamma") not_cont[i] <- 0
      if(node_distributions[i] == "nbinom") not_cont[i] <- 1
      
    }
    
    # now estimate network
    est_net <- BDgraph::bdgraph(data = data, 
                                method = "gcgm", 
                                algorithm = "bdmcmc",
                                not.cont = not_cont,
                                iter = 50000, 
                                burnin = 10000, 
                                g.prior = 0.5, 
                                df.prior = 3)
    
    # selected structure - median probability model
    est_net <- summary(est_net, vis = FALSE)$selected_g 
    
    # compare estimated to true
    comp <- welfareNet::compare.sim.networks(true_net, est_net, 
                                             metric = c("sensitivity","specificity","precision"))
    
    #add mean pcor and density to output
    comp$mean_pcor <- mean_pcor
    comp$density <- dens
    
    return(comp)
  }
)

#saveRDS(object = fmgm_smallworld_sim, file = "fmgm_smallworld_sim.rds")
#readRDS("fmgm_smallworld_sim.rds")

# summaries for table
aggregate(density ~ nei + node_number, data = bgcgm_smallworld_sim, mean)
aggregate(density ~ nei +  node_number, data = bgcgm_smallworld_sim, sd)

aggregate(mean_pcor ~ nei + node_number, data = bgcgm_smallworld_sim, mean)
round(aggregate(mean_pcor ~ nei + node_number, data = bgcgm_smallworld_sim, sd),2)

aggregate(sensitivity ~ sample_size + nei + node_number, data = bgcgm_smallworld_sim, mean)
aggregate(sensitivity ~ sample_size + nei + node_number, data = bgcgm_smallworld_sim, sd)

aggregate(specificity ~ sample_size + nei + node_number, data = bgcgm_smallworld_sim, mean)
aggregate(specificity ~ sample_size + nei + node_number, data = bgcgm_smallworld_sim, sd)

aggregate(precision ~ sample_size + nei + node_number, data = bgcgm_smallworld_sim, mean)
aggregate(precision ~ sample_size + nei + node_number, data = bgcgm_smallworld_sim, sd)

# plots - not shown in text
ggplot(bgcgm_smallworld_sim, aes(x = factor(sample_size), y = mean_pcor)) +
  facet_wrap(~ nei + node_number) +
  geom_violin(scale = "width") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size=0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "red", size = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "blue", size = 0.5) +
  labs(x = "Sample Size", y = "mean_pcor", 
       title = "") +
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,0.1)) +
  theme_classic()

ggplot(bgcgm_smallworld_sim, aes(x = factor(sample_size), y = sensitivity)) +
  facet_wrap(~ nei + node_number) +
  geom_violin(scale = "width") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size=0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "red", size = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "blue", size = 0.5) +
  labs(x = "Sample Size", y = "Sensitivity", 
       title = "") +
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,0.1)) +
  theme_classic()


ggplot(bgcgm_smallworld_sim, aes(x = factor(sample_size), y = specificity)) +
  facet_wrap(~ nei + node_number) +
  geom_violin(scale = "width") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size=0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "red", size = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "blue", size = 0.5) +
  labs(x = "Sample Size", y = "Specificity", 
       title = "") +
  scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,0.1)) +
  theme_classic()

###################
#
#     END
#
###################
