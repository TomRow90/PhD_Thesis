##############################################
#
#    Chapter 4: R code for conventional/
#    average treatment effect models
#
##############################################

#-----load required packages-------#

library(brms)
library(bayestestR)
library(bayesplot)
library(ggplot2)
library(ggpubr)

#----- load data -------#

# set wd
setwd("/Users/tomrowland/Documents/PhD/Thesis write up/Chapter 4/data")

network.data <- read.table("network.data.txt", header=TRUE)

# change factor level order so coefficients represent going from unenriched to enriched
network.data$condition <- factor(network.data$condition, levels=c("std","ee"))


#--------------------------------------------#
#------------ functions used ----------------#
#--------------------------------------------#

#transformation for beta distributed variables with any precise 0s or 1s
beta.transform <- function(y, n) (y * (n - 1) + 0.5) / n


# function to perform prior sensitivity analysis
prior.sensitivity.function <- function(model, priors) {
  
  results <- as.data.frame(matrix(0, nrow = length(priors), ncol = 6))
  colnames(results) <- c("prior","estimate","lower_ci","upper_ci","p_direction","BF")
  results$prior <- c("original","wider","tighter","student t", "default uninformative")
  results$BF <- NA
  
  for(i in 1:length(priors)) {
    
    updated.mod <- update(model, prior = priors[[i]])
    
    results$estimate[i] <- point_estimate(updated.mod, centrality = "median")$Median[2]
    interval <- hdi(updated.mod)
    results$lower_ci[i] <- interval$CI_low[2]
    results$upper_ci[i] <- interval$CI_high[2]
    results$p_direction[i] <- p_direction(updated.mod)$pd[2]
    
    if(i < length(priors)) {
    
    bf <- bayesfactor(updated.mod, null = 0)
    results$BF[i] <- round(exp(bf$log_BF)[2],2)
    
    }
    
  }
  
  return(results)
  
}


#--------------------------------------------#
#---------novel tank diving test-------------#
#--------------------------------------------#

# transform bottom to proportions to fit beta reg
network.data$diving.bottom.t <- network.data$diving.bottom / 300

# means and sd
div.bot.means <- tapply(network.data$diving.bottom.t, network.data$condition, mean)
div.bot.sd <- tapply(network.data$diving.bottom.t, network.data$condition, sd)

div.bot.means; div.bot.sd
div.bot.means[2] - div.bot.means[1]

diving.data <- na.omit(network.data[,c("diving.bottom.t","condition")])

#---------- prior predictive sim ----------------------#

# priors
diving.prior <- c(
  prior(normal(1, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(4, 0.1), class = "phi")
)

#fit model
div.prior.sim <- brm(diving.bottom.t ~ condition,
                     data = diving.data,
                    family = Beta(link = "logit"),
                    prior = diving.prior,
                    sample_prior = "only",
                    iter = 10000,
                    cores = 10,
                    backend = "cmdstanr")

# check prior pred dist.
pp_check(div.prior.sim, type = "dens_overlay", ndraws = 250, prefix = "ppd")


#------------- posterior sampling -------------------#

# fit model
div.bot.beta <- brm(diving.bottom.t ~ condition,
                    data = diving.data,
                    family = Beta(link = "logit"),
                    prior = diving.prior,
                    chains = 4,
                    iter = 12000,
                    warmup = 2000,
                    cores = 10,
                    backend = "cmdstanr")

# check mcmc diagnostics
plot(div.bot.beta)
summary(div.bot.beta, robust = TRUE)

# inference
hdi(div.bot.beta, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(div.bot.beta)
div.bf <- bayesfactor(div.bot.beta, null = 0) # note you may get very slightly different value but it shouldnt be substantively different
div.bf
plot(div.bf)

# plot not shown in thesis
conditional_effects(div.bot.beta, prob = 0.95)

# post pred check
pp_check(div.bot.beta, type = "dens_overlay", ndraws = 250, prefix = "ppc")

# pp check of sd in each group
diving.yrep <- posterior_predict(div.bot.beta)

ppc_stat_grouped(
  y = diving.data$diving.bottom.t,
  yrep = diving.yrep,
  group = diving.data$condition,
  stat = "sd"
)


#--------- diving bottom prior sensitivity analysis --------------#

diving.priors.list <- list(
  
  # original
  c(prior(normal(1, 1), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
 # higher uncertainty in all parameters
  c(prior(normal(1, 2), class = "Intercept"),
    prior(normal(0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),

  # tighter uncertainty in all parameters
  c(prior(normal(1, 0.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),

  # fatter tailed student t
 c(prior(student_t(3, 1, 2), class = "Intercept"),
   prior(student_t(3, 0, 2), class = "b"),
   prior(gamma(4, 0.1), class = "phi")),
 
  # default uninformative prior
 get_prior(diving.bottom.t ~ condition,
           data = diving.data,
           family = Beta(link = "logit"))
)



diving.bottom.prior.sensitivity <- prior.sensitivity.function(model = div.bot.beta,
                                                              priors = diving.priors.list)

diving.bottom.prior.sensitivity

#---------opercular beat rate------------#

obr.means <- tapply(network.data$diving.obr, network.data$condition, mean, na.rm=TRUE)
obr.sd <- tapply(network.data$diving.obr, network.data$condition, sd, na.rm=TRUE)


obr.means; obr.sd
obr.means[2] - obr.means[1]

obr.data <- na.omit(network.data[,c("condition","diving.obr")])

#---------- prior predictive sim ----------------------#

obr.prior <- c(
  prior(normal(150, 25), class = "Intercept"),
  prior(normal(0, 10), class = "b")
)

#fit model
obr.prior.sim <- brm(diving.obr ~ condition,
                     data = obr.data,
                     family = gaussian,
                     prior = obr.prior,
                     sample_prior = "only",
                     iter = 12000,
                     warmup = 2000,
                     cores = 10,
                     backend = "cmdstanr")

pp_check(obr.prior.sim, type = "dens_overlay", ndraws = 100, prefix = "ppd")


#------------- posterior sampling -------------------#

# fit model
obr.lm <- brm(diving.obr ~ condition,
              data = obr.data,
                    family = gaussian,
                    prior = obr.prior,
                    iter = 12000,
                    warmup = 2000,
                    cores = 10,
                    backend = "cmdstanr")

# check mcmc diagnostics
plot(obr.lm)
summary(obr.lm, robust = TRUE)

# inference
hdi(obr.lm, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(obr.lm)
obr.bf <- bayesfactor(obr.lm, null = 0) # this is somewhat unstable because very little mass at 0 but should always give a very large BF_10
obr.bf
plot(obr.bf)

# PPD check
pp_check(obr.lm, type = "dens_overlay", ndraws = 250, prefix = "ppc")

# pp check of sd in each group
obr.yrep <- posterior_predict(obr.lm)

ppc_stat_grouped(
  y = obr.data$diving.obr,
  yrep = obr.yrep,
  group = obr.data$condition,
  stat = "sd"
)


#--------- opercular beat rate prior sensitivity analysis --------------#

obr.priors.list <- list(
  
  # original
  c(prior(normal(150, 25), class = "Intercept"),
  prior(normal(0, 10), class = "b")),

  # higher uncertainty in all parameters
  c(prior(normal(150, 50), class = "Intercept"),
    prior(normal(0, 20), class = "b")),
  
  # tighter uncertainty in all parameters
  c(prior(normal(150, 12.5), class = "Intercept"),
    prior(normal(0, 5), class = "b")),
  
  # fatter tailed student t
  c(prior(student_t(3, 150, 30), class = "Intercept"),
    prior(student_t(3, 0, 15), class = "b")),
  
  # default uninformative prior
  get_prior(diving.obr ~ condition,
          data = obr.data,
          family = gaussian)
)


obr.prior.sensitivity <- prior.sensitivity.function(model = obr.lm,
                                                    priors = obr.priors.list)

obr.prior.sensitivity

# plots for diving test variables - not shown in chapter text where focus is on
# network analysis

div.plot <- ggplot(network.data, aes(x = condition, y = diving.bottom.t)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(0,1)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=7)) +
  labs(y="Proportion of time in bottom zone in diving test", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))


#plot
obr.plot <- ggplot(network.data, aes(x = condition, y = diving.obr)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(130,200)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=7)) +
  labs(y="Opercular beat rate (bpm) in diving test", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))


# plot not shown in thesis
ggarrange(div.plot, obr.plot,
          nrow=1, ncol=2,
          labels=c("A)","B)"),
          font.label=list(size=7))


#--------------------------------------------#
#-------------black white test---------------#
#--------------------------------------------#

#------- black duration / scototaxis ---------#

# transform to proportion
network.data$scototaxis.t <- network.data$bw.black / 300

bw.black.means <- tapply(network.data$scototaxis.t, network.data$condition, mean, na.rm=TRUE)
bw.black.sd <- tapply(network.data$scototaxis.t, network.data$condition, sd, na.rm=TRUE)

bw.black.means; bw.black.sd
bw.black.means[2] - bw.black.means[1]

# two tanks with precise zeros
ind0 <- which(network.data$scototaxis.t == 0) # 2
ind1 <- which(network.data$scototaxis.t == 1) # 0

# squeeze transformation within [0,1]
n <- length(na.omit(network.data$scototaxis.t))
network.data$scototaxis.t2 <- beta.transform(y=network.data$scototaxis.t, n=n)

scototaxis.data <- na.omit(network.data[,c("scototaxis.t2","condition")])

#---------- prior predictive sim ----------------------#

scototaxis.prior <- c(
  prior(normal(1, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(4, 0.1), class = "phi")
)

# fit model
scototaxis.prior.sim <- brm(scototaxis.t2 ~ condition,
                            data = scototaxis.data,
                            family = Beta(link = "logit"),
                            prior = scototaxis.prior,
                            sample_prior = "only",
                            iter = 10000,
                            warmup = 2000,
                            cores = 10,
                            backend = "cmdstanr")

pp_check(scototaxis.prior.sim, type = "dens_overlay", ndraws = 500, prefix = "ppd")

#------------- posterior sampling -------------------#

# fit model
scototaxis.beta <- brm(scototaxis.t2 ~ condition,
                       data = scototaxis.data,
                       family = Beta(link = "logit"),
                       prior = scototaxis.prior,
                       iter = 12000,
                       warmup = 2000,
                       cores = 10,
                       backend = "cmdstanr")

# check mcmc diagnostics
plot(scototaxis.beta)
summary(scototaxis.beta, robust = TRUE)

# inference
hdi(scototaxis.beta, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(scototaxis.beta)
scototaxis.bf <- bayesfactor(scototaxis.beta, null = 0)
scototaxis.bf
plot(scototaxis.bf)

# PPD
pp_check(scototaxis.beta, type = "dens_overlay", ndraws =250, prefix = "ppc")

# pp check of sd in each group
scototaxis.yrep <- posterior_predict(scototaxis.beta)

ppc_stat_grouped(
  y = scototaxis.data$scototaxis.t,
  yrep = scototaxis.yrep,
  group = scototaxis.data$condition,
  stat = "sd"
)


#--------- scototaxis prior sensitivity analysis --------------#

scototaxis.priors.list <- list(
  
  # original
  c(prior(normal(1, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(4, 0.1), class = "phi")),
  
  # higher uncertainty in all parameters
  c(prior(normal(1, 2), class = "Intercept"),
    prior(normal(0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),

  # tighter uncertainty in all parameters
  c(prior(normal(1, 0.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # fatter tailed student t
  c(prior(student_t(3, 1, 2), class = "Intercept"),
    prior(student_t(3, 0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # default uninformative prior
  get_prior(scototaxis.t2 ~ condition,
            data = scototaxis.data,
            family = Beta(link = "logit"))
)


scototaxis.prior.sensitivity <- prior.sensitivity.function(model = scototaxis.beta,
                                                              priors = scototaxis.priors.list)

scototaxis.prior.sensitivity


#---------stationary-------------#

# transform stationary to proportion
network.data$bw.stationary.t <- network.data$bw.stationary / 300

bw.st.means <- tapply(network.data$bw.stationary.t, network.data$condition, mean, na.rm=TRUE)
bw.st.sd <- tapply(network.data$bw.stationary.t, network.data$condition, sd, na.rm=TRUE)

bw.st.means; bw.st.sd
bw.st.means[2] - bw.st.means[1]

# any exact 0s and 1s? no
ind0 <- which(network.data$bw.stationary.t == 0) # 0
ind1 <- which(network.data$bw.stationary.t == 1) # 0

bw.stationary.data <- na.omit(network.data[,c("bw.stationary.t","condition")])


#---------- prior predictive sim ----------------------#

bw.stationary.prior <- c(
  prior(normal(1, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(4, 0.1), class = "phi")
)

# fit model
bw.stationary.prior.sim <- brm(bw.stationary.t ~ condition,
                            data = bw.stationary.data,
                            family = Beta(link = "logit"),
                            prior = bw.stationary.prior,
                            sample_prior = "only",
                            iter = 12000,
                            warmup = 2000,
                            cores = 10,
                            backend = "cmdstanr")

pp_check(bw.stationary.prior.sim, type = "dens_overlay", ndraws = 100, prefix = "ppd")

#------------- posterior sampling -------------------#

# fit model
bw.stationary.beta <- brm(bw.stationary.t ~ condition,
                       data = bw.stationary.data,
                       family = Beta(link = "logit"),
                       prior = bw.stationary.prior,
                       iter = 12000,
                       warmup = 2000,
                       cores = 10,
                       backend = "cmdstanr")

# check mcmc diagnostics
plot(bw.stationary.beta)
summary(bw.stationary.beta, robust = TRUE)

# inference
hdi(bw.stationary.beta, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(bw.stationary.beta)
bw.stationary.bf <- bayesfactor(bw.stationary.beta, null = 0)
bw.stationary.bf
plot(bw.stationary.bf)

# ppd 
pp_check(bw.stationary.beta, type = "dens_overlay", ndraws = 250, prefix = "ppc")

# pp check of sd in each group
bw.stationary.yrep <- posterior_predict(bw.stationary.beta)

ppc_stat_grouped(
  y = bw.stationary.data$bw.stationary.t,
  yrep = bw.stationary.yrep,
  group = bw.stationary.data$condition,
  stat = "sd"
)


#--------- bw stationary prior sensitivity analysis --------------#

bw.stationary.priors.list <- list(
  
  # original
  c(prior(normal(1, 1), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # higher uncertainty in all parameters
  c(prior(normal(1, 2), class = "Intercept"),
    prior(normal(0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # tighter uncertainty in all parameters
  c(prior(normal(1, 0.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # fatter tailed student t
  c(prior(student_t(3, 1, 2), class = "Intercept"),
    prior(student_t(3, 0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # default uninformative prior
  get_prior(bw.stationary.t ~ condition,
            data = bw.stationary.data,
            family = Beta(link = "logit"))
)

bw.stationary.prior.sensitivity <- prior.sensitivity.function(model = bw.stationary.beta,
                                                           priors = bw.stationary.priors.list)

bw.stationary.prior.sensitivity


# plots
black.plot <- ggplot(network.data, aes(x = condition, y = scototaxis.t)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(0,1)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=7)) +
  labs(y="Proportion of time in black zone in black-white test", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))

#plot
black.stat.plot <- ggplot(network.data, aes(x = condition, y = bw.stationary.t)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(0,1)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=7)) +
  labs(y="Proportion of time stationary in black-white test", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))


ggarrange(black.plot, black.stat.plot, 
          labels=c("A)","B)"),
          font.label=list(size=7),
          label.x = 0.025,
          label.y = 1)


#--------------------------------------------#
#-------------- Open Field Test -------------#
#--------------------------------------------#

# transform to proportion
network.data$oft.stationary.t <- network.data$oft.stationary / 300

oft.st.means <- tapply(network.data$oft.stationary.t, network.data$condition, mean, na.rm=TRUE)
oft.st.sd <- tapply(network.data$oft.stationary.t, network.data$condition, sd, na.rm=TRUE)

oft.st.means; oft.st.sd
oft.st.means[2] - oft.st.means[1]

# any precise 0 or 1? no
ind0 <- which(network.data$oft.stationary.t == 0)#0
ind1 <- which(network.data$oft.stationary.t == 1) #0

oft.stationary.data <- na.omit(network.data[,c("oft.stationary.t","condition")])

#---------- prior predictive sim ----------------------#

oft.stationary.prior <- priors <- c(
  prior(normal(1, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(4, 0.1), class = "phi")
)

# fit model
oft.stationary.prior.sim <- brm(oft.stationary.t ~ condition,
                            data = oft.stationary.data,
                            family = Beta(link = "logit"),
                            prior = oft.stationary.prior,
                            sample_prior = "only",
                            iter = 10000,
                            warmup = 2000,
                            cores = 10,
                            backend = "cmdstanr")

pp_check(oft.stationary.prior.sim, type = "dens_overlay", ndraws = 250, prefix = "ppd")

#------------- posterior sampling -------------------#

# fit model
oft.stationary.beta <- brm(oft.stationary.t ~ condition,
                       data = oft.stationary.data,
                       family = Beta(link = "logit"),
                       prior = oft.stationary.prior,
                       chains = 4,
                       iter = 12000,
                       warmup = 2000,
                       cores = 10,
                       backend = "cmdstanr")

# check mcmc diagnostics
plot(oft.stationary.beta)
summary(oft.stationary.beta, robust = TRUE)

# inference
hdi(oft.stationary.beta, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(oft.stationary.beta)
oft.stationary.bf <- bayesfactor(oft.stationary.beta, null = 0)
oft.stationary.bf
plot(oft.stationary.bf)

pp_check(oft.stationary.beta, type = "dens_overlay", ndraws = 250, prefix = "ppc")

# pp check of sd in each group
oft.stationary.yrep <- posterior_predict(oft.stationary.beta)

ppc_stat_grouped(
  y = oft.stationary.data$oft.stationary.t,
  yrep = oft.stationary.yrep,
  group = oft.stationary.data$condition,
  stat = "sd"
)

#--------- oft stationary prior sensitivity analysis --------------#

oft.stationary.priors.list <- list(
  
  # original
  c(prior(normal(1, 1), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # higher uncertainty in all parameters
  c(prior(normal(1, 2), class = "Intercept"),
    prior(normal(0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # tighter uncertainty in all parameters
  c(prior(normal(1, 0.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # fatter tailed student t
  c(prior(student_t(3, 1, 2), class = "Intercept"),
    prior(student_t(3, 0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # default uninformative prior
  get_prior(oft.stationary.t ~ condition,
            data = oft.stationary.data,
            family = Beta(link = "logit"))
  
)

oft.stationary.prior.sensitivity <- prior.sensitivity.function(model = oft.stationary.beta,
                                                              priors = oft.stationary.priors.list)

oft.stationary.prior.sensitivity


#--------external/thigmotaxis---------#


# transform to proportion
network.data$thigmotaxis.t <- network.data$oft.external / 300

oft.thig.means <- tapply(network.data$thigmotaxis.t, network.data$condition, mean, na.rm=TRUE)
oft.thig.sd <- tapply(network.data$thigmotaxis.t, network.data$condition, sd, na.rm=TRUE)

oft.thig.means; oft.thig.sd

# any precise 0 or 1? no
ind0 <- which(network.data$thigmotaxis.t == 0) #0
ind1 <- which(network.data$thigmotaxis.t == 1) #0

oft.thigmotaxis.data <- na.omit(network.data[,c("thigmotaxis.t","condition")])


#---------- prior predictive sim ----------------------#


oft.thigmotaxis.prior <- c(
  prior(normal(1, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(4, 0.1), class = "phi")
)

# fit model
oft.thigmotaxis.prior.sim <- brm(thigmotaxis.t ~ condition,
                                data = oft.thigmotaxis.data,
                                family = Beta(link = "logit"),
                                prior = oft.thigmotaxis.prior,
                                sample_prior = "only",
                                iter = 10000,
                                warmup = 2000,
                                cores = 10,
                                backend = "cmdstanr")

pp_check(oft.stationary.prior.sim, type = "dens_overlay", ndraws = 250, prefix = "ppd")



#------------- posterior sampling -------------------#

# fit model
oft.thigmotaxis.beta <- brm(thigmotaxis.t ~ condition,
                           data = oft.thigmotaxis.data,
                           family = Beta(link = "logit"),
                           prior = oft.thigmotaxis.prior,
                           chains = 4,
                           iter = 12000,
                           warmup = 2000,
                           cores = 10,
                           backend = "cmdstanr")

# check mcmc diagnostics
plot(oft.thigmotaxis.beta)
summary(oft.thigmotaxis.beta, robust = TRUE)

# inference
hdi(oft.thigmotaxis.beta, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(oft.thigmotaxis.beta)
oft.thigmotaxis.bf <- bayesfactor(oft.thigmotaxis.beta)
oft.thigmotaxis.bf
plot(oft.thigmotaxis.bf)

pp_check(oft.thigmotaxis.beta, type = "dens_overlay", ndraws = 250, prefix = "ppc")

# pp check of sd in each group
thigmotaxis.yrep <- posterior_predict(oft.thigmotaxis.beta)

ppc_stat_grouped(
  y = oft.thigmotaxis.data$thigmotaxis.t,
  yrep = thigmotaxis.yrep,
  group = oft.thigmotaxis.data$condition,
  stat = "sd"
)


#--------- oft thigmotaxis prior sensitivity analysis --------------#

oft.thigmotaxis.priors.list <- list(
  
  # original
  c(prior(normal(1, 1), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # higher uncertainty in all parameters
  c(prior(normal(1, 2), class = "Intercept"),
    prior(normal(0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # tighter uncertainty in all parameters
  c(prior(normal(1, 0.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # fatter tailed student t
  c(prior(student_t(3, 1, 2), class = "Intercept"),
    prior(student_t(3, 0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # default uninformative prior
  get_prior(thigmotaxis.t ~ condition,
            data = oft.thigmotaxis.data,
            family = Beta(link = "logit"))
  
)

oft.thigmotaxis.prior.sensitivity <- prior.sensitivity.function(model = oft.thigmotaxis.beta,
                                                               priors = oft.thigmotaxis.priors.list)

oft.thigmotaxis.prior.sensitivity


# plots
oft.stat.plot <-  ggplot(network.data, aes(x = condition, y = oft.stationary.t)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(0,1)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=7)) +
  labs(y="Proportion of time stationary in open field test", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))

thig.plot <-  ggplot(network.data, aes(x = condition, y = thigmotaxis.t)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(0,1)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=7)) +
  labs(y="Proportion of time in external zone in open field test", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))


ggarrange(thig.plot, oft.stat.plot, 
          nrow=1, ncol=2,
          labels=c("A)","B)"),
          font.label=list(size=7),
          label.x = 0.025,
          label.y = 1)


#--------------------------------------------#
#----------- Novel Object Test --------------#
#--------------------------------------------#

#-------not stationary variable-------#

# transform to proportion
network.data$not.stationary.t <- network.data$not.stationary / 300

not.st.means <- tapply(network.data$not.stationary.t, network.data$condition, mean, na.rm=TRUE)
not.st.sd <- tapply(network.data$not.stationary.t, network.data$condition, sd, na.rm=TRUE)

not.st.means; not.st.sd

# any precise 0 or 1? no
ind0 <- which(network.data$not.stationary.t == 0) #0
ind1 <- which(network.data$not.stationary.t == 1) #0

not.stationary.data <- na.omit(network.data[,c("not.stationary.t","condition")])

#---------- prior predictive sim ----------------------#


not.stationary.prior <- priors <- c(
  prior(normal(1, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(4, 0.1), class = "phi")
)

# fit model
not.stationary.prior.sim <- brm(not.stationary.t ~ condition,
                           data = not.stationary.data,
                           family = Beta(link = "logit"),
                           prior = not.stationary.prior,
                           sample_prior = "only",
                           chains = 4,
                           iter = 12000,
                           warmup = 2000,
                           cores = 10,
                           backend = "cmdstanr")

pp_check(not.stationary.prior.sim, type = "dens_overlay", ndraws = 250, prefix = "ppd")


#------------- posterior sampling -------------------#


# fit model
not.stationary.beta <- brm(not.stationary.t ~ condition,
                            data = not.stationary.data,
                            family = Beta(link = "logit"),
                            prior = not.stationary.prior,
                            chains = 4,
                            iter = 12000,
                            warmup = 2000,
                            cores = 10,
                            backend = "cmdstanr")

# check mcmc diagnostics
plot(not.stationary.beta)
summary(not.stationary.beta, robust = TRUE)

# inference
hdi(not.stationary.beta, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(not.stationary.beta)
not.stationary.bf <- bayesfactor(not.stationary.beta)
not.stationary.bf
plot(not.stationary.bf)

# ppd
pp_check(not.stationary.beta, type = "dens_overlay", ndraws = 250, prefix = "ppc")

# pp check of sd in each group
not.stationary.yrep <- posterior_predict(not.stationary.beta)

ppc_stat_grouped(
  y = not.stationary.data$not.stationary.t,
  yrep = not.stationary.yrep,
  group = not.stationary.data$condition,
  stat = "sd"
)


#--------- not stationary prior sensitivity analysis --------------#

not.stationary.priors.list <- list(
  
  # original
  c(prior(normal(1, 1), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # higher uncertainty in all parameters
  c(prior(normal(1, 2), class = "Intercept"),
    prior(normal(0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # tighter uncertainty in all parameters
  c(prior(normal(1, 0.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # fatter tailed student t
  c(prior(student_t(3, 1, 2), class = "Intercept"),
    prior(student_t(3, 0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # default uninformative prior
  get_prior(not.stationary.t ~ condition,
            data = not.stationary.data,
            family = Beta(link = "logit"))
  
)

not.stationary.prior.sensitivity <- prior.sensitivity.function(model = not.stationary.beta,
                                                                priors = not.stationary.priors.list)

not.stationary.prior.sensitivity


#----------not latency variable-------------#


# transform to prop
network.data$not.latency.t <- network.data$not.latency / 300

not.lat.means <- tapply(network.data$not.latency.t, network.data$condition, mean, na.rm=TRUE)
not.lat.sd <- tapply(network.data$not.latency.t, network.data$condition, sd, na.rm=TRUE)

not.lat.means; not.lat.sd

# any precise 0 or 1? no
which(network.data$not.latency.t == 0)
which(network.data$not.latency.t == 1)

not.latency.data <- network.data[,c("not.latency.t","condition")]


#---------- prior predictive sim ----------------------#


not.latency.prior <- priors <- c(
  prior(normal(0,1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(4, 0.1), class = "phi")
)


# fit model
not.latency.prior.sim <- brm(not.latency.t ~ condition,
                           data = not.latency.data,
                           family = Beta(link = "logit"),
                           prior = not.latency.prior,
                           sample_prior = "only",
                           chains = 4,
                           iter = 12000,
                           warmup = 2000,
                           cores = 10,
                           backend = "cmdstanr")


pp_check(not.latency.prior.sim, type = "dens_overlay", ndraws = 250, prefix = "ppd")

#------------- posterior sampling -------------------#

# fit model
not.latency.beta <- brm(not.latency.t ~ condition,
                           data = not.latency.data,
                           family = Beta(link = "logit"),
                           prior = not.latency.prior,
                           chains = 4,
                           iter = 12000,
                           warmup = 2000,
                           cores = 10,
                           backend = "cmdstanr")

# check mcmc diagnostics
plot(not.latency.beta)
summary(not.latency.beta, robust = TRUE)

# inference
hdi(not.latency.beta, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(not.latency.beta)
not.latency.bf <- bayesfactor(not.latency.beta)
not.latency.bf
plot(not.latency.bf)


#--------- not latency prior sensitivity analysis --------------#

not.latency.priors.list <- list(
  
  # original
  c(prior(normal(1, 1), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # higher uncertainty in all parameters
  c(prior(normal(1, 2), class = "Intercept"),
    prior(normal(0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # tighter uncertainty in all parameters
  c(prior(normal(1, 0.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # fatter tailed student t
  c(prior(student_t(3, 0, 2), class = "Intercept"),
    prior(student_t(3, 0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # default uninformative prior
  get_prior(not.latency.t ~ condition,
            data = not.latency.data,
            family = Beta(link = "logit"))
  
)

not.latency.prior.sensitivity <- prior.sensitivity.function(model = not.latency.beta,
                                                               priors = not.latency.priors.list)

not.latency.prior.sensitivity

# plot
not.stat.plot <- ggplot(network.data, aes(x = condition, y = not.stationary.t)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(0,1)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=7)) +
  labs(y="Proportion of time stationary in novel object test", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))


# plot
lat.plot <- ggplot(network.data, aes(x = condition, y = not.latency)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(0,300)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=7)) +
  labs(y="Latency (s) to enter object zone in novel object test", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))


ggarrange(not.stat.plot, lat.plot, 
          nrow=1, ncol=2,
          labels=c("A)","B)"),
          font.label=list(size=7),
          label.x = 0.025,
          label.y = 1)


#--------------------------------------------#
#--------- Novel Object Recognition ---------#
#--------------------------------------------#

nor.means <- tapply(network.data$nor.discrimination, network.data$condition, mean, na.rm=TRUE)
nor.sd <- tapply(network.data$nor.discrimination, network.data$condition, sd, na.rm=TRUE)

nor.means; nor.sd

# transform to beta 0,1 interval

network.data$nor.discrimination.t <- (network.data$nor.discrimination + 1) / 2

# check all within 0,1
any(network.data$nor.discrimination.t == 1)
any(network.data$nor.discrimination.t == 0)

nor.data <- na.omit(network.data[,c("nor.discrimination.t","condition")])

#---------- prior predictive sim ----------------------#


nor.prior <- priors <- c(
  prior(normal(0,1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(4, 0.1), class = "phi")
)


# fit model
nor.prior.sim <- brm(nor.discrimination.t ~ condition,
                             data = nor.data,
                             family = Beta(link = "logit"),
                             prior = nor.prior,
                             sample_prior = "only",
                             chains = 4,
                             iter = 12000,
                             warmup = 2000,
                             cores = 10,
                             backend = "cmdstanr")


pp_check(nor.prior.sim, type = "dens_overlay", ndraws = 100, prefix = "ppd")

#------------- posterior sampling -------------------#

# fit model
nor.beta <- brm(nor.discrimination.t ~ condition,
                        data = nor.data,
                        family = Beta(link = "logit"),
                        prior = nor.prior,
                        chains = 4,
                        iter = 12000,
                        warmup = 2000,
                        cores = 10,
                        backend = "cmdstanr")

# check mcmc diagnostics
plot(nor.beta)
summary(nor.beta, robust = TRUE)

# inference
hdi(nor.beta, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(nor.beta)
nor.bf <- bayesfactor(nor.beta)
nor.bf
plot(nor.bf)

# visual
plot(conditional_effects(nor.beta))

# pp check
pp_check(nor.beta, type = "dens_overlay", ndraws = 250, prefix = "ppc")

# pp check of sd in each group
nor.yrep <- posterior_predict(nor.beta)

ppc_stat_grouped(
  y = nor.data$nor.discrimination.t,
  yrep = nor.yrep,
  group = nor.data$condition,
  stat = "sd"
)

#--------- nor prior sensitivity analysis --------------#

nor.priors.list <- list(
  
  # original
  c(prior(normal(1, 1), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # higher uncertainty in all parameters
  c(prior(normal(1, 2), class = "Intercept"),
    prior(normal(0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # tighter uncertainty in all parameters
  c(prior(normal(1, 0.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # fatter tailed student t
  c(prior(student_t(3, 0, 2), class = "Intercept"),
    prior(student_t(3, 0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # default uninformative prior
  get_prior(nor.discrimination.t ~ condition,
            data = nor.data,
            family = Beta(link = "logit"))
  
)

nor.prior.sensitivity <- prior.sensitivity.function(model = nor.beta,
                                                            priors = nor.priors.list)

nor.prior.sensitivity

# plot
nor.plot <- ggplot(network.data, aes(x = condition, y = nor.discrimination)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(-1,1)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=7)) +
  labs(y="Discrimination index in object recognition", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))


nor.plot


#--------------------------------------------#
#-------- response to mild startle ----------#
#--------------------------------------------#

# which ate or not
ind0 <- which(network.data$startle == 180)
length(ind0) #10

# transform to proportion and apply squeeze transformation
network.data$startle.t <- network.data$startle / 180
n <- length(na.omit(network.data$startle.t))
network.data$startle.t2 <- beta.transform(y=network.data$startle.t, n=n)

startle.data <- na.omit(network.data[,c("startle.t2","condition")])

#---------- prior predictive sim ----------------------#


startle.prior <- priors <- c(
  prior(normal(0,1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(4, 0.1), class = "phi")
)


# fit model
startle.prior.sim <- brm(startle.t2 ~ condition,
                     data = startle.data,
                     family = Beta(link = "logit"),
                     prior = startle.prior,
                     sample_prior = "only",
                     chains = 4,
                     iter = 12000,
                     warmup = 2000,
                     cores = 10,
                     backend = "cmdstanr")


pp_check(startle.prior.sim, type = "dens_overlay", ndraws = 100, prefix = "ppd")

#------------- posterior sampling -------------------#

# fit model
startle.beta <- brm(startle.t2 ~ condition,
                data = startle.data,
                family = Beta(link = "logit"),
                prior = startle.prior,
                chains = 4,
                iter = 12000,
                warmup = 2000,
                cores = 10,
                backend = "cmdstanr")

# check mcmc diagnostics
plot(startle.beta)
summary(startle.beta, robust = TRUE)

# inference
hdi(startle.beta, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(startle.beta)
startle.bf <- bayesfactor(startle.beta)
startle.bf
plot(startle.bf)

# pp check
pp_check(startle.beta, type = "dens_overlay", ndraws = 250, prefix = "ppc")

# pp check of sd in each group
startle.yrep <- posterior_predict(startle.beta)

ppc_stat_grouped(
  y = startle.data$startle.t2,
  yrep = startle.yrep,
  group = startle.data$condition,
  stat = "sd"
)

# how many in std did not eat
std.ind <- which(network.data$condition == "std")
length(which(network.data$startle.t[std.ind] == 1)) #6

# how many in ee did not eat?
ee.ind <- which(network.data$condition == "ee")
length(which(network.data$startle.t[ee.ind] == 1)) #4


#--------- startle sensitivity analysis --------------#

startle.priors.list <- list(
  
  # original
  c(prior(normal(0, 1), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # higher uncertainty in all parameters
  c(prior(normal(0, 2), class = "Intercept"),
    prior(normal(0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # tighter uncertainty in all parameters
  c(prior(normal(0, 0.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # fatter tailed student t
  c(prior(student_t(3, 0, 2), class = "Intercept"),
    prior(student_t(3, 0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # default uninformative prior
  get_prior(startle.t2 ~ condition,
            data = startle.data,
            family = Beta(link = "logit"))
  
)

startle.prior.sensitivity <- prior.sensitivity.function(model = startle.beta,
                                                    priors = startle.priors.list)

startle.prior.sensitivity

#plot
startle.plot <- ggplot(network.data, aes(x = condition, y = startle)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(0,180),
                     breaks = seq(0, 180, by = 20)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=6)) +
  labs(y="Latency to resume eating", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))

startle.plot


#--------------------------------------------#
#--------- social network assay -------------#
#--------------------------------------------#

soci.means <- tapply(network.data$social.dens, network.data$condition, mean, na.rm=TRUE)
soci.sd <- tapply(network.data$social.dens, network.data$condition, sd, na.rm=TRUE)

soci.means; soci.sd

social.data <- na.omit(network.data[,c("social.dens","condition")])

#---------- prior predictive sim ----------------------#


social.prior <- priors <- c(
  prior(normal(0,1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(4, 0.1), class = "phi")
)


# fit model
social.prior.sim <- brm(social.dens ~ condition,
                     data = social.data,
                     family = Beta(link = "logit"),
                     prior = social.prior,
                     sample_prior = "only",
                     chains = 4,
                     iter = 12000,
                     warmup = 2000,
                     cores = 10,
                     backend = "cmdstanr")


pp_check(social.prior.sim, type = "dens_overlay", ndraws = 100, prefix = "ppd")

#------------- posterior sampling -------------------#

# fit model
social.beta <- brm(social.dens ~ condition,
                data = social.data,
                family = Beta(link = "logit"),
                prior = social.prior,
                chains = 4,
                iter = 12000,
                warmup = 2000,
                cores = 10,
                backend = "cmdstanr")

# check mcmc diagnostics
plot(social.beta)
summary(social.beta, robust = TRUE)

# inference
hdi(social.beta, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(social.beta)
social.bf <- bayesfactor(social.beta)
social.bf
plot(social.bf)

# pp check
pp_check(social.beta, type = "dens_overlay", ndraws = 250, prefix = "ppc")

# pp check of sd in each group
social.yrep <- posterior_predict(social.beta)

ppc_stat_grouped(
  y = social.data$social.dens,
  yrep = social.yrep,
  group = social.data$condition,
  stat = "sd"
)


#--------- social sensitivity analysis --------------#

social.priors.list <- list(
  
  # original
  c(prior(normal(0, 1), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # higher uncertainty in all parameters
  c(prior(normal(0, 2), class = "Intercept"),
    prior(normal(0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # tighter uncertainty in all parameters
  c(prior(normal(0, 0.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # fatter tailed student t
  c(prior(student_t(3, 0, 2), class = "Intercept"),
    prior(student_t(3, 0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # default uninformative prior
  get_prior(social.dens ~ condition,
            data = social.data,
            family = Beta(link = "logit"))
  
)

# perform prior sensitivity analysis
social.prior.sensitivity <- prior.sensitivity.function(model = social.beta,
                                                        priors = social.priors.list)

social.prior.sensitivity

#plots
soc.plot <- ggplot(network.data, aes(x = condition, y = social.dens)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(0,1)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=8)) +
  labs(y="Mean social network density", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))


soc.plot


#--------- repeat analysis on tanks which had all 4 fish ------------#

# tanks that had 3 fish includes:
# IDrep1tank4,5,6,8,9,10
# IDrep2tank1,5
# IDrep3tank4,8
# IDrep4tank1,5,9,10

idx <- which(network.data$id %in% c("IDrep1tank4",
                                    "IDrep1tank5",
                                    "IDrep1tank6",
                                    "IDrep1tank8",
                                    "IDrep1tank10",
                                    "IDrep2tank1",
                                    "IDrep2tank5",
                                    "IDrep3tank4",
                                    "IDrep3tank8",
                                    "IDrep4tank1",
                                    "IDrep4tank5",
                                    "IDrep4tank9",
                                    "IDrep4tank10"))

# is that equal numbers? very close = 7 std and 6 ee
# so densities wont be heavily biased lower for one condition over the other
summary(network.data[idx,"condition"])

social.data2 <- na.omit(network.data[-idx,c("social.dens","condition")])

# fit model
social.beta2 <- brm(social.dens ~ condition,
                   data = social.data2,
                   family = Beta(link = "logit"),
                   prior = social.prior,
                   chains = 4,
                   iter = 12000,
                   warmup = 2000,
                   cores = 10,
                   backend = "cmdstanr")

# check mcmc diagnostics
plot(social.beta2)
summary(social.beta2, robust = TRUE)

# inference - very similar results to using all data points
hdi(social.beta2, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(social.beta2)
social.bf2 <- bayesfactor(social.beta2)
social.bf2
plot(social.bf2)


# plot
soc.plot2 <- ggplot(social.data2, aes(x = condition, y = social.dens)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(0,1)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=8)) +
  labs(y="Mean social network density", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))

soc.plot2

#--------------------------------------------#
#-------- home tank routine feeding ---------#
#--------------------------------------------#

# transform to proportion
network.data$home.feed.t <- network.data$appetite / 600

# mean and sds
home.feed.means <- tapply(network.data$home.feed.t, network.data$condition, mean, na.rm=TRUE)
home.feed.sd <- tapply(network.data$home.feed.t, network.data$condition, sd, na.rm=TRUE)

home.feed.means; home.feed.sd

home.feed.data <- na.omit(network.data[,c("home.feed.t","condition")])


#---------- prior predictive sim ----------------------#


home.feed.prior <- priors <- c(
  prior(normal(-1,1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(4, 0.1), class = "phi")
)


# fit model
home.feed.prior.sim <- brm(home.feed.t ~ condition,
                     data = home.feed.data,
                     family = Beta(link = "logit"),
                     prior = home.feed.prior,
                     sample_prior = "only",
                     chains = 4,
                     iter = 12000,
                     warmup = 2000,
                     cores = 10,
                     backend = "cmdstanr")


pp_check(home.feed.prior.sim, type = "dens_overlay", ndraws = 100, prefix = "ppd")

#------------- posterior sampling -------------------#

# fit model
home.feed.beta <- brm(home.feed.t ~ condition,
                data = home.feed.data,
                family = Beta(link = "logit"),
                prior = home.feed.prior,
                chains = 4,
                iter = 12000,
                warmup = 2000,
                cores = 10,
                backend = "cmdstanr")

# check mcmc diagnostics
plot(home.feed.beta)
summary(home.feed.beta, robust = TRUE)

# inference
hdi(home.feed.beta, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(home.feed.beta)
home.feed.bf <- bayesfactor(home.feed.beta)
home.feed.bf
plot(home.feed.bf)

# pp check
pp_check(home.feed.beta, type = "dens_overlay", ndraws = 250, prefix = "ppc")

# pp check of sd in each group
home.feed.yrep <- posterior_predict(home.feed.beta)

ppc_stat_grouped(
  y = home.feed.data$home.feed.t,
  yrep = home.feed.yrep,
  group = home.feed.data$condition,
  stat = "sd"
)


#--------- home tank feed sensitivity analysis --------------#

home.feed.priors.list <- list(
  
  # original
  c(prior(normal(-1, 1), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # higher uncertainty in all parameters
  c(prior(normal(-1, 2), class = "Intercept"),
    prior(normal(0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # tighter uncertainty in all parameters
  c(prior(normal(-1, 0.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # fatter tailed student t
  c(prior(student_t(3, -1, 2), class = "Intercept"),
    prior(student_t(3, 0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # default uninformative prior
  get_prior(home.feed.t ~ condition,
            data = home.feed.data,
            family = Beta(link = "logit"))
  
)

home.feed.prior.sensitivity <- prior.sensitivity.function(model = home.feed.beta,
                                                        priors = home.feed.priors.list)

home.feed.prior.sensitivity

# plot of untransformed data
home.feed.plot <- ggplot(network.data, aes(x = condition, y = appetite)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(0,600),
                     breaks = seq(0, 600, by = 100)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=8)) +
  labs(y="Routine feeding: Latency to eat (s)", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))


home.feed.plot


#--------------------------------------------#
#--------- home tank behaviour---------------#
#--------------------------------------------#


#--------activity variable--------#

tapply(network.data$activity, network.data$condition, mean, na.rm=TRUE)
tapply(network.data$activity, network.data$condition, sd, na.rm=TRUE)

#any precise 0 or 1? yes
which(network.data$activity == 0)
length(which(network.data$activity == 1)) #1
ind1 <- which(network.data$activity == 1)

#squeeze transform to within [0,1]
n <- length(na.omit(network.data$activity))
network.data$activity.t <- beta.transform(y=network.data$activity, n=n)

activity.data <- na.omit(network.data[,c("activity.t","condition")])

#---------- prior predictive sim ----------------------#


act.prior <-c(prior(normal(0,1), class = "Intercept"),
              prior(normal(0, 1), class = "b"),
              prior(gamma(4, 0.1), class = "phi"))

#fit model
act.prior.sim <- brm(activity.t ~ condition,
                     data = activity.data,
                     family = Beta(link = "logit"),
                     prior = act.prior,
                     sample_prior = "only",
                     chains = 4,
                     iter = 12000,
                     warmup = 2000,
                     cores = 10,
                     backend = "cmdstanr")


pp_check(act.prior.sim, type = "dens_overlay", ndraws = 100, prefix = "ppd")


#------------- posterior sampling -------------------#

#fit model
act.beta <- brm(activity.t ~ condition,
                data = activity.data,
                family = Beta(link = "logit"),
                prior = act.prior,
                chains = 4,
                iter = 12000,
                warmup = 2000,
                cores = 10,
                backend = "cmdstanr")

#check mcmc diagnostics
plot(act.beta)
summary(act.beta, robust = TRUE)

# inference
hdi(act.beta, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(act.beta)
act.bf <- bayesfactor(act.beta)
act.bf
plot(act.bf)

# pp check
pp_check(act.beta, type = "dens_overlay", ndraws = 250, prefix = "ppc")

act.yrep <- posterior_predict(act.beta)

ppc_stat_grouped(
  y = activity.data$activity.t,
  yrep = act.yrep,
  group = activity.data$condition,
  stat = "sd"
)



#--------- activity prior sensitivity analysis --------------#

activity.priors.list <- list(
  
  # original
  c(prior(normal(0,1), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # higher uncertainty in all parameters
  c(prior(normal(0,2), class = "Intercept"),
    prior(normal(0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # tighter uncertainty in all parameters
  c(prior(normal(0,0.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # student t fatter tail version of original
  c(prior(student_t(3, 0, 2), class = "Intercept"),
    prior(student_t(3, 0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # default uninformative prior
  get_prior(activity.t ~ condition,
            data = activity.data,
            family = Beta(link = "logit"))
)



activity.prior.sensitivity <- prior.sensitivity.function(model = act.beta,
                                                         priors = activity.priors.list)

activity.prior.sensitivity


#------aggressive behaviour variable-----------#

# mean & sd
tapply(network.data$agg.beh, network.data$condition, mean, na.rm=TRUE)
tapply(network.data$agg.beh, network.data$condition, sd, na.rm=TRUE)


#any precise 0 or 1? yes
ind0 <- which(network.data$agg.beh == 0) #3
length(which(network.data$agg.beh == 1)) #0

#squeeze transform to within [0,1]
n <- length(na.omit(network.data$agg.beh))
network.data$aggression.t <- beta.transform(y=network.data$agg.beh, n=n)

# define new dataframe- makes ppc later easier
agg.data <- network.data[,c("aggression.t","condition")]
agg.data <- na.omit(agg.data)


#---------- prior predictive sim ----------------------#


agg.prior <- priors <- c(
  prior(normal(-1,1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(4, 0.1), class = "phi")
)


#fit model
agg.prior.sim <- brm(aggression.t ~ condition,
                     data = agg.data,
                     family = Beta(link = "logit"),
                     prior = agg.prior,
                     sample_prior = "only",
                     chains = 4,
                     iter = 12000,
                     warmup = 2000,
                     cores = 10,
                     backend = "cmdstanr")


pp_check(agg.prior.sim, type = "dens_overlay", ndraws = 100, prefix = "ppd")

#------------- posterior sampling -------------------#

#fit model
agg.beta <- brm(aggression.t ~ condition,
                data = agg.data,
                family = Beta(link = "logit"),
                prior = agg.prior,
                chains = 4,
                iter = 12000,
                warmup = 2000,
                cores = 10,
                backend = "cmdstanr")

#check mcmc diagnostics
plot(agg.beta)
summary(agg.beta)

# inference
hdi(agg.beta, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(agg.beta)
agg.bf <- bayesfactor(agg.beta)
agg.bf
plot(agg.bf)

# pp check
pp_check(agg.beta, type = "dens_overlay", ndraws = 250, prefix = "ppc")

# pp check of sd in each group
agg.yrep <- posterior_predict(agg.beta)

ppc_stat_grouped(
  y = agg.data$aggression.t,
  yrep = agg.yrep,
  group = agg.data$condition,
  stat = "sd"
)

# refit modelling precision parameter given SD not well modelled in model above

# what are default priors for this model
get_prior(bf(
  aggression.t ~ condition,
  phi ~ condition),
  data = agg.data,
  family = Beta(link = "logit"))

# reset priors
agg.prior2 <- c(
  prior(normal(-1, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b", coef = "conditionee"),
  prior(normal(log(40), 1), class = "Intercept", dpar = "phi"),
  prior(normal(0, 1), class = "b", coef = "conditionee", dpar = "phi"))


#fit model
agg.prior.sim2 <- brm(bf(
                     aggression.t ~ condition,
                     phi ~ condition),
                     data = agg.data,
                     family = Beta(link = "logit"),
                     prior = agg.prior2,
                     sample_prior = "only",
                     chains = 4,
                     iter = 12000,
                     warmup = 2000,
                     cores = 10,
                     backend = "cmdstanr")

# prior sim check
pp_check(agg.prior.sim2, type = "dens_overlay", ndraws = 100, prefix = "ppd")

# get prior samples for latter BF
agg.prior.samples <- as_draws_df(agg.prior.sim2)

#fit model
agg.beta2 <- brm(bf(
                aggression.t ~ condition,
                phi ~ condition),
                data = agg.data,
                family = Beta(link = "logit"),
                prior = agg.prior2,
                chains = 4,
                iter = 12000,
                warmup = 2000,
                cores = 10,
                backend = "cmdstanr")

#check mcmc diagnostics
plot(agg.beta2)
summary(agg.beta2, robust = TRUE)

# get posterior samples
agg.post <- as_draws_df(agg.beta2)

hdi(agg.post$b_conditionee, ci = c(0.5, 0.68, 0.95, 0.99))
hdi(agg.post$b_phi_conditionee, ci = c(0.5, 0.68, 0.95, 0.99))

p_direction(agg.beta2)
p_direction(agg.post$b_phi_conditionee)

agg.bf <- bayesfactor_parameters(posterior = agg.post$b_conditionee,
                                 prior = agg.prior.samples$b_conditionee)
agg.bf
plot(agg.bf)

agg.bf2 <- bayesfactor_parameters(posterior = agg.post$b_phi_conditionee,
                                  prior = agg.prior.samples$b_phi_conditionee)
agg.bf2
plot(agg.bf2)

# pp check
pp_check(agg.beta2, type = "dens_overlay", ndraws = 250, prefix = "ppc")

# pp check of sd in each group
agg.yrep2 <- posterior_predict(agg.beta2)

# not perfect but better than before 
ppc_stat_grouped(
  y = agg.data$aggression.t,
  yrep = agg.yrep2,
  group = agg.data$condition,
  stat = "sd"
)

#--------- aggression prior sensitivity analysis --------------#

aggression.priors.list <- list(
  
  # original
  c(prior(normal(-1, 1), class = "Intercept"),
    prior(normal(0, 1), class = "b", coef = "conditionee"),
    prior(normal(log(40), 1), class = "Intercept", dpar = "phi"),
    prior(normal(0, 1), class = "b", coef = "conditionee", dpar = "phi")),
  
  
  # higher uncertainty in all parameters
  c(prior(normal(-1,2), class = "Intercept"),
    prior(normal(0, 2), class = "b", coef = "conditionee"),
    prior(normal(log(40), 2), class = "Intercept", dpar = "phi"),
    prior(normal(0,2), class="b", coef = "conditionee", dpar="phi")),

  # tighter uncertainty in all parameters
  c(prior(normal(-1,0.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "b", coef = "conditionee"),
    prior(normal(log(40), 0.5), class = "Intercept", dpar = "phi"),
    prior(normal(0,0.5), class="b", coef = "conditionee", dpar="phi")),
  
  # student t fatter tail version of original
  c(prior(student_t(3,-1.1,2), class = "Intercept"),
    prior(student_t(3, 0, 2), class = "b", coef = "conditionee"),
    prior(student_t(3, 0, 2), class = "Intercept", dpar = "phi"),
    prior(student_t(3, 0, 2), class = "b", coef = "conditionee", dpar = "phi")),
  
  # default uninformative priors
  get_prior(bf(
    aggression.t ~ condition,
    phi ~ condition),
    data = agg.data,
    family = Beta(link = "logit"))
)

  
aggression.prior.sensitivity <- prior.sensitivity.function(model = agg.beta2,
                                                           priors = aggression.priors.list)

aggression.prior.sensitivity


#--------stationary on bottom substrate variable-------------#

# for convenience change variable name here (currently freezing, but better named differently)
colnames(network.data)[17] <- "stationary.bottom"

#mean and sd
tapply(network.data$stationary.bottom, network.data$condition, mean, na.rm=TRUE)
tapply(network.data$stationary.bottom, network.data$condition, sd, na.rm=TRUE)

#any precise 0 or 1? yes
ind0 <- which(network.data$stationary.bottom == 0) #3
ind1 <- which(network.data$stationary.bottom == 1) #5

#squeeze transform to within [0,1]
n <- length(na.omit(network.data$stationary.bottom))
network.data$stationary.bottom.t <- beta.transform(y=network.data$stationary.bottom, n = n)

stationary.bottom.data <- na.omit(network.data[,c("stationary.bottom.t","condition")])

#---------- prior predictive sim ----------------------#


stationary.bottom.prior <- priors <- c(
  prior(normal(0,1.5), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(4, 0.1), class = "phi")
)


#fit model
stationary.bottom.prior.sim <- brm(stationary.bottom.t ~ condition,
                     data = stationary.bottom.data,
                     family = Beta(link = "logit"),
                     prior = stationary.bottom.prior,
                     sample_prior = "only",
                     chains = 4,
                     iter = 12000,
                     warmup = 2000,
                     cores = 10,
                     backend = "cmdstanr")


pp_check(stationary.bottom.prior.sim, type = "dens_overlay", ndraws = 250, prefix = "ppd")

#------------- posterior sampling -------------------#

#fit model
stationary.bottom.beta <- brm(stationary.bottom.t ~ condition,
                data = stationary.bottom.data,
                family = Beta(link = "logit"),
                prior = stationary.bottom.prior,
                chains = 4,
                iter = 12000,
                warmup = 2000,
                cores = 10,
                backend = "cmdstanr")

#check mcmc diagnostics
plot(stationary.bottom.beta)
summary(stationary.bottom.beta, robust = TRUE)

# inference
hdi(stationary.bottom.beta, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(stationary.bottom.beta)
stationary.bottom.bf <- bayesfactor(stationary.bottom.beta)
stationary.bottom.bf
plot(stationary.bottom.bf)

# pp check
pp_check(stationary.bottom.beta, type = "dens_overlay", ndraws = 250, prefix = "ppc")

# pp check of sd in each group
stationary.bottom.yrep <- posterior_predict(stationary.bottom.beta)

ppc_stat_grouped(
  y = stationary.bottom.data$stationary.bottom.t,
  yrep = stationary.bottom.yrep,
  group = stationary.bottom.data$condition,
  stat = "sd"
)

# hmm its not bad, but modelling precision may improve it

stationary.bottom.prior2 <- priors <- c(
  prior(normal(0,1.5), class = "Intercept"),
  prior(normal(0, 1), class = "b", coef = "conditionee"),
  prior(normal(log(40), 1), class = "Intercept", dpar = "phi"),
  prior(normal(0,1), class="b", coef = "conditionee", dpar="phi"))


# fit model
stationary.bottom.prior.sim2 <- brm(bf(
                     stationary.bottom.t ~ condition,
                     phi ~ condition),
                     data = stationary.bottom.data,
                     family = Beta(link = "logit"),
                     prior = stationary.bottom.prior2,
                     sample_prior = "only",
                     chains = 4,
                     iter = 12000,
                     warmup = 2000,
                     cores = 10,
                     backend = "cmdstanr")

# Prior check
pp_check(stationary.bottom.prior.sim2, type = "dens_overlay", ndraws = 100, prefix = "ppd")

stationary.bottom.prior.samples <- as_draws_df(stationary.bottom.prior.sim2)

#fit model
stationary.bottom.beta2 <- brm(bf(
                stationary.bottom.t ~ condition,
                phi ~ condition),
                data = stationary.bottom.data,
                family = Beta(link = "logit"),
                prior = stationary.bottom.prior2,
                chains = 4,
                iter = 12000,
                warmup = 2000,
                cores = 10,
                backend = "cmdstanr")

# check mcmc diagnostics
plot(stationary.bottom.beta2)
summary(stationary.bottom.beta2, robust = TRUE)

# inference
stationary.bottom.post <- as_draws_df(stationary.bottom.beta2)

hdi(stationary.bottom.post$b_conditionee, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(stationary.bottom.post$b_conditionee)

hdi(stationary.bottom.post$b_phi_conditionee, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(stationary.bottom.post$b_phi_conditionee)

stationary.bottom.bf <- bayesfactor_parameters(posterior = stationary.bottom.post$b_conditionee,
                                                prior = stationary.bottom.prior.samples$b_conditionee)
stationary.bottom.bf
plot(stationary.bottom.bf)

stationary.bottom.bf2 <- bayesfactor_parameters(posterior = stationary.bottom.post$b_phi_conditionee,
                                  prior = stationary.bottom.prior.samples$b_phi_conditionee)
stationary.bottom.bf2
plot(stationary.bottom.bf2)

# pp check
pp_check(stationary.bottom.beta2, type = "dens_overlay", ndraws = 250, prefix = "ppc")

# pp check of sd in each group
stationary.bottom.yrep2 <- posterior_predict(stationary.bottom.beta2)

# looks a bit better
ppc_stat_grouped(
  y = stationary.bottom.data$stationary.bottom.t,
  yrep = stationary.bottom.yrep2,
  group = stationary.bottom.data$condition,
  stat = "sd"
)


#--------- stationary bottom prior sensitivity analysis --------------#

stationary.bottom.priors.list <- list(
  
  # original
  c(prior(normal(0,1.5), class = "Intercept"),
    prior(normal(0, 1), class = "b", coef = "conditionee"),
    prior(normal(log(40), 1), class = "Intercept", dpar = "phi"),
    prior(normal(0,1), class="b", coef = "conditionee", dpar="phi")),
  
  # higher uncertainty in all parameters
  c(prior(normal(0,3), class = "Intercept"),
    prior(normal(0, 2), class = "b", coef = "conditionee"),
    prior(normal(log(40), 2), class = "Intercept", dpar = "phi"),
    prior(normal(0,2), class="b", coef = "conditionee", dpar="phi")),
  
  # tighter uncertainty in all parameters
  c(prior(normal(0,0.75), class = "Intercept"),
    prior(normal(0, 0.5), class = "b", coef = "conditionee"),
    prior(normal(log(40), 0.5), class = "Intercept", dpar = "phi"),
    prior(normal(0,0.5), class="b", coef = "conditionee", dpar="phi")),
  
  # student t fatter tail version of original
  c(prior(student_t(3,0, 2), class = "Intercept"),
    prior(student_t(3, 0, 2), class = "b", coef = "conditionee"),
    prior(student_t(3, 0, 2), class = "Intercept", dpar = "phi"),
    prior(student_t(3, 0, 2), class = "b", coef = "conditionee", dpar = "phi")),
  
  # default uninformative priors
  get_prior(bf(
    stationary.bottom.t ~ condition,
    phi ~ condition),
    data = stationary.bottom.data,
    family = Beta(link = "logit"))
)


stationary.bottom.prior.sensitivity <- prior.sensitivity.function(model = stationary.bottom.beta2,
                                                           priors = stationary.bottom.priors.list)

stationary.bottom.prior.sensitivity


#--------foraging variable-----------#

#mean and sd
tapply(network.data$foraging, network.data$condition, mean, na.rm=TRUE)
tapply(network.data$foraging, network.data$condition, sd, na.rm=TRUE)

#any precise 0 or 1? no
which(network.data$foraging == 0)
which(network.data$foraging == 1)

foraging.data <- na.omit(network.data[,c("foraging","condition")])

#---------- prior predictive sim ----------------------#


foraging.prior <- c(prior(normal(0,1.5), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(4, 0.1), class = "phi"))


#fit model
foraging.prior.sim <- brm(foraging ~ condition,
                                   data = foraging.data,
                                   family = Beta(link = "logit"),
                                   prior = foraging.prior,
                                   sample_prior = "only",
                                   chains = 4,
                                   iter = 12000,
                                   warmup = 2000,
                                   cores = 10,
                                   backend = "cmdstanr")


pp_check(foraging.prior.sim, type = "dens_overlay", ndraws = 250, prefix = "ppd")

#------------- posterior sampling -------------------#

#fit model
foraging.beta <- brm(foraging ~ condition,
                              data = foraging.data,
                              family = Beta(link = "logit"),
                              prior = foraging.prior,
                              chains = 4,
                              iter = 12000,
                              warmup = 2000,
                              cores = 10,
                              backend = "cmdstanr")

#check mcmc diagnostics
plot(foraging.beta)
summary(foraging.beta, robust = TRUE)

# inference
hdi(foraging.beta, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(foraging.beta)
foraging.bf <- bayesfactor(foraging.beta)
foraging.bf
plot(foraging.bf)

# pp check
pp_check(foraging.beta, type = "dens_overlay", ndraws = 250, prefix = "ppc")

# pp check of sd in each group
foraging.yrep <- posterior_predict(foraging.beta)

ppc_stat_grouped(
  y = foraging.data$foraging,
  yrep = foraging.yrep,
  group = foraging.data$condition,
  stat = "sd"
)


#--------- foraging prior sensitivity analysis --------------#

foraging.priors.list <- list(
  
  # original
  c(prior(normal(0,1.5), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # higher uncertainty in all parameters
  c(prior(normal(0,3), class = "Intercept"),
    prior(normal(0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # tighter uncertainty in all parameters
  c(prior(normal(0,0.75), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # student t fatter tail version of original
  c(prior(student_t(3, 0, 2), class = "Intercept"),
    prior(student_t(3, 0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # default uninformative priors
  get_prior(bf(
    foraging ~ condition,
    phi ~ condition),
    data = foraging.data,
    family = Beta(link = "logit"))
)


foraging.prior.sensitivity <- prior.sensitivity.function(model = foraging.beta,
                                                                  priors = foraging.priors.list)

foraging.prior.sensitivity



#------wall swimming--------#

#mean and sd
tapply(network.data$wall.swimming, network.data$condition, mean, na.rm=TRUE)
tapply(network.data$wall.swimming, network.data$condition, sd, na.rm=TRUE)

#any precise 0 or 1? yes
ind0 <- which(network.data$wall.swimming == 0) # 5
which(network.data$wall.swimming == 1) #0

#squeeze transform to within [0,1]
n <- length(na.omit(network.data$wall.swimming))
network.data$wall.swimming.t <- beta.transform(y=network.data$wall.swimming, n=n)

wall.swimming.data <- na.omit(network.data[,c("wall.swimming.t","condition")])

#---------- prior predictive sim ----------------------#


wall.swimming.prior <- priors <- c(
  prior(normal(0,1.5), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(4, 0.1), class = "phi")
)


#fit model
wall.swimming.prior.sim <- brm(wall.swimming.t ~ condition,
                                   data = wall.swimming.data,
                                   family = Beta(link = "logit"),
                                   prior = wall.swimming.prior,
                                   sample_prior = "only",
                                   chains = 4,
                                   iter = 12000,
                                   warmup = 2000,
                                   cores = 10,
                                   backend = "cmdstanr")


pp_check(wall.swimming.prior.sim, type = "dens_overlay", ndraws = 250, prefix = "ppd")

#------------- posterior sampling -------------------#

#fit model
wall.swimming.beta <- brm(wall.swimming.t ~ condition,
                              data = wall.swimming.data,
                              family = Beta(link = "logit"),
                              prior = wall.swimming.prior,
                              chains = 4,
                              iter = 12000,
                              warmup = 2000,
                              cores = 10,
                              backend = "cmdstanr")

#check mcmc diagnostics
plot(wall.swimming.beta)
summary(wall.swimming.beta, robust = TRUE)

# inference
hdi(wall.swimming.beta, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(wall.swimming.beta)
wall.swimming.bf <- bayesfactor(wall.swimming.beta)
wall.swimming.bf
plot(wall.swimming.bf)

# pp check
pp_check(wall.swimming.beta, type = "dens_overlay", ndraws = 250, prefix = "ppc")

# pp check of sd in each group
wall.swimming.yrep <- posterior_predict(wall.swimming.beta)

ppc_stat_grouped(
  y = wall.swimming.data$wall.swimming.t,
  yrep = wall.swimming.yrep,
  group = wall.swimming.data$condition,
  stat = "sd"
)

#--------- foraging prior sensitivity analysis --------------#

wall.swimming.priors.list <- list(
  
  # original
  c(prior(normal(0,1.5), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # higher uncertainty in all parameters
  c(prior(normal(0,3), class = "Intercept"),
    prior(normal(0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # tighter uncertainty in all parameters
  c(prior(normal(0,0.75), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # student t fatter tail version of original
  c(prior(student_t(3, 0, 2), class = "Intercept"),
    prior(student_t(3, 0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # default uninformative priors
  get_prior(bf(
    wall.swimming.t ~ condition,
    phi ~ condition),
    data = wall.swimming.data,
    family = Beta(link = "logit"))
)


wall.swimming.prior.sensitivity <- prior.sensitivity.function(model = wall.swimming.beta,
                                                              priors = wall.swimming.priors.list)

wall.swimming.prior.sensitivity

#-------- vertical distribution ---------#

colnames(network.data)[20] <- "vertical.distribution"

tapply(network.data$vertical.distribution, network.data$condition, mean, na.rm=TRUE)
tapply(network.data$vertical.distribution, network.data$condition, sd, na.rm=TRUE)

ind0 <- which(network.data$vertical.distribution == 0) #1
which(network.data$vertical.distribution == 1) #0

n <- length(na.omit(network.data$vertical.distribution))
network.data$vertical.distribution.t <- beta.transform(y=network.data$vertical.distribution, n=n)

vertical.distribution.data <- na.omit(network.data[,c("vertical.distribution.t","condition")])

#---------- prior predictive sim ----------------------#


vertical.distribution.prior <- priors <- c(
  prior(normal(-1,1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(4, 0.1), class = "phi")
)


#fit model
vertical.distribution.prior.sim <- brm(vertical.distribution.t ~ condition,
                               data = vertical.distribution.data,
                               family = Beta(link = "logit"),
                               prior = vertical.distribution.prior,
                               sample_prior = "only",
                               chains = 4,
                               iter = 12000,
                               warmup = 2000,
                               cores = 10,
                               backend = "cmdstanr")


pp_check(vertical.distribution.prior.sim, type = "dens_overlay", ndraws = 250, prefix = "ppd")

#------------- posterior sampling -------------------#

#fit model
vertical.distribution.beta <- brm(vertical.distribution.t ~ condition,
                          data = vertical.distribution.data,
                          family = Beta(link = "logit"),
                          prior = vertical.distribution.prior,
                          chains = 4,
                          iter = 12000,
                          warmup = 2000,
                          cores = 10,
                          backend = "cmdstanr")

#check mcmc diagnostics
plot(vertical.distribution.beta)
summary(vertical.distribution.beta, robust = TRUE)

# inference
hdi(vertical.distribution.beta, ci = c(0.5, 0.68, 0.95, 0.99))
p_direction(vertical.distribution.beta)
vertical.distribution.bf <- bayesfactor(vertical.distribution.beta)
vertical.distribution.bf
plot(vertical.distribution.bf)

# pp check
pp_check(vertical.distribution.beta, type = "dens_overlay", ndraws = 250, prefix = "ppc")

# pp check of sd in each group
vertical.distribution.yrep <- posterior_predict(vertical.distribution.beta)

ppc_stat_grouped(
  y = vertical.distribution.data$vertical.distribution.t,
  yrep = vertical.distribution.yrep,
  group = vertical.distribution.data$condition,
  stat = "sd"
)

#--------- vertical distribution sensitivity analysis --------------#

vertical.distribution.priors.list <- list(
  
  # original
  c(prior(normal(-1, 1), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # higher uncertainty in all parameters
  c(prior(normal(-1, 2), class = "Intercept"),
    prior(normal(0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # tighter uncertainty in all parameters
  c(prior(normal(-1, 0.5), class = "Intercept"),
    prior(normal(0, 0.5), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # fatter tailed student t
  c(prior(student_t(3, -1, 2), class = "Intercept"),
    prior(student_t(3, 0, 2), class = "b"),
    prior(gamma(4, 0.1), class = "phi")),
  
  # default uninformative prior
  get_prior(vertical.distribution.t ~ condition,
            data = vertical.distribution.data,
            family = Beta(link = "logit"))
  
)

vertical.distribution.prior.sensitivity <- prior.sensitivity.function(model = vertical.distribution.beta,
                                                          priors = vertical.distribution.priors.list)

vertical.distribution.prior.sensitivity


#---------home tank plots--------------#

#plot
act.plot <- ggplot(network.data, aes(x = condition, y = activity)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(0,1)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=7)) +
  labs(y="Proportion of videos activity occurred", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))

#plot
agg.plot <- ggplot(network.data, aes(x = condition, y = agg.beh)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(0,1)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=8)) +
  labs(y="Proportion of videos aggression occurred", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))

#plot
for.plot <- ggplot(network.data, aes(x = condition, y = foraging)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(0,1)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=8)) +
  labs(y="Proportion of videos foraging occurred", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))

#plot
sub.plot <- ggplot(network.data, aes(x = condition, y = stationary.bottom)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(0,1)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=8)) +
  labs(y="Proportion of videos fish stationary on floor", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))

#plot
top.plot <- ggplot(network.data, aes(x = condition, y = vertical.distribution)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(0,1)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=8)) +
  labs(y="Mean proportion of fish in top third of tank", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))


#plot
wall.plot <- ggplot(network.data, aes(x = condition, y = wall.swimming)) +
  geom_violin(width = 0.5) +  
  stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black", size = 0.2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "gray", size = 0.2) +
  scale_fill_manual(values = c("ee" = "grey", "std" = "white"), guide = "none") + 
  geom_jitter(width = 0.1, height = 0, size = 1, alpha = 0.7) +
  scale_y_continuous(limits=c(0,1)) +
  theme_classic(base_size=8) + 
  theme(axis.title.y = element_text(size=8)) +
  labs(y="Proportion of videos wall swimming occurred", x="Condition") +
  scale_x_discrete(labels=c("unenriched","enrichment"))


ggarrange(act.plot, agg.plot,
          for.plot, sub.plot,
          wall.plot, top.plot,
          labels = c("A)","B)","C)","D)",
                     "E)","F)"),
          font.label=list(size=7),
          label.x = 0.025,
          label.y = 1,
          nrow=3, ncol=2)

######################################
#                                    #
#                END                 #
#                                    #  
######################################
