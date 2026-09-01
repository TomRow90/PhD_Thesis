#####################################################################
#
#        Section 3.5.2: Power simulations for marginal dependence
#
####################################################################

# load required packages
library(simDAG)
library(ggplot2)
library(ggpubr)
library(wesanderson)
library(welfareNet)
library(parSim)


setwd("/Users/tomrowland/Documents/PhD/Thesis write up/Chapter 3")


#------------------------ functions used -----------------------------#

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


# estimate mixed pair mutual information
estimate.mi <- function(x,y,type,permutations) {
  
  if(type[1] == "c" & type[2] == "c") {
    
    empirical.mi <- mpmi::cmi.pw(x, y)$mi
    
  }
  
  if(type[1] == "d" & type[2] == "d") {
    
    empirical.mi  <- mpmi::dmi.pw(x, y)$mi
    
  }
  
  if(type[1] == "c" & type[2] == "d") {
    
    empirical.mi  <- mpmi::mmi.pw(x, y)$mi
    
  }
  
  if(type[1] == "d" & type[2] == "c") {
    
    empirical.mi  <- mpmi::mmi.pw(y, x)$mi
    
  }
  
  
  return(empirical.mi)
}

# permutation function for mixed pair p value
mi.perm <- function(x,y, type) {
  
  data <- as.data.frame(cbind(x,y))
  colnames(data) <- c("x","y")
  data <- na.omit(data)
  
  n.obs <- nrow(data)
  
  # generate permuted dataset
  s <- sample(1:n.obs,n.obs,replace=FALSE)
  perm.x <- data[s,1]
  perm.y <- data[,2]
  
  # test on permuted data
  if(type[1] == "c" & type[2] == "c") {
    
    test <- mpmi::cmi.pw(perm.x, perm.y)$mi
    
    return(test)
    
  }
  
  if(type[1] == "d" & type[2] == "d") {
    
    test <- mpmi::dmi.pw(perm.x, perm.y)$mi
    
    return(test)
    
  }
  
  if(type[1] == "c" & type[2] == "d") {
    
    test <- mpmi::mmi.pw(perm.x, perm.y)$mi
    
    return(test)
    
  }
  
  if(type[1] == "d" & type[2] == "c") {
    
    test <- mpmi::mmi.pw(perm.y, perm.x)$mi
    
    return(test)
    
  }
  
}


# function that takes data, calculates mixed pair MI and performs test if requested
mi.test <- function(x,y,type,test=TRUE,permutations) {
  
  #empirical estimate
  if(type[1] == "c" & type[2] == "c") {
    
    empirical.mi <- mpmi::cmi.pw(x, y)$mi
    
  }
  
  if(type[1] == "d" & type[2] == "d") {
    
    empirical.mi  <- mpmi::dmi.pw(x, y)$mi
    
  }
  
  if(type[1] == "c" & type[2] == "d") {
    
    empirical.mi  <- mpmi::mmi.pw(x, y)$mi
    
  }
  
  if(type[1] == "d" & type[2] == "c") {
    
    empirical.mi  <- mpmi::mmi.pw(y, x)$mi
    
  }
  
  if(test == TRUE) {
    
    # run permutation function n times
    perms <- replicate(permutations, mi.perm(x,y,type), simplify = TRUE)
    
    b <- sum(perms >= empirical.mi)
    
    m <- length(perms)
    
    p <- (b + 1) / (m + 1)
    
    return(list(mutual.information = empirical.mi,
                p.value = p))
    
  } else {
    
    return(list(mutual.information = empirical.mi))
    
  }
  
  
}


# helper function to get power for each simulated sample size
get.power <- function(data, sample.size, p, alpha, type) {
  
  ind <- which(data$sample.size == sample.size & data$p == p)
  
  dep.ind <- which(colnames(data) == type)
  
  na.ind <- which(is.na(data[ind, dep.ind]) == TRUE)
  
  power <- sum(data[ind,dep.ind] < alpha, na.rm=TRUE) / (length(ind) - length(na.ind))
  
  return(power)
  
}


# I spline function for fitting power curves
# set up with help from info from this stackexchange post:
# https://stats.stackexchange.com/questions/519465/how-to-correctly-use-i-splines-for-monotone-non-decreasing-increasing-regressio

spline.fun <- function(sample.size, power, minSampleSize, maxSampleSize, correlation=NULL, method) {
  
  # degree 3 = cubic spline
  degree <- 3
  
  n <- length(sample.size) - 2 # - 2 because we will be using LOOCV
  
  # max knots
  min_df <- degree + 1
  max_df <- n
  
  # MSE vector
  MSE <- rep(NA, length(min_df:max_df))
  
  # prediction data sets to store so can return optimal one
  pred_d <- list()
  
  # look through different degrees of freedom / knot numbers
  for (k in min_df:max_df) {
    
    errors <- rep(NA, n)
    
    # find optimal 
    for (i in 1:n) {
      
      # Leave-one-out split
      x_train <- sample.size[-i]
      y_train <- power[-i]
      x_test <- sample.size[i]
      y_test <- power[i]
      
      # create the I-spline
      ispline <- splines2::iSpline(x_train, df = k, degree = degree, intercept = TRUE)
      
      ispline <- cbind(1, ispline)
      
      # creating `D` matrix and `d` vector.
      d_mat <- crossprod(ispline, ispline)
      d_vec <- crossprod(ispline, y_train)
      a_mat <- diag(1, ncol(ispline))
      
      b_vec <- rep(0, ncol(ispline))
      
      # freeing the first parameter/intercept.
      a_mat[1, 1] <- 0
      
      # solve
      qp.solve <- quadprog::solve.QP(Dmat = d_mat, dvec = d_vec, Amat = t(a_mat), bvec = b_vec)$solution
      
      # full range of sample sizes
      Ns <- minSampleSize:maxSampleSize
      
      # create I-spline basis for the full sample size range
      full_spline <- splines2::iSpline(Ns, df=k, degree = degree, intercept = TRUE)
      
      # add the same intercept column you used before
      full_spline <- cbind(1, full_spline)
      
      # interpolated statistics
      g_hat <- full_spline %*% qp.solve
      
      d <- as.data.frame(cbind(Ns,g_hat))
      colnames(d) <- c("Ns","g_hat")
      
      ind <- which(d[,1] == x_test)
      
      y_pred <- d[ind,2]
      
      # prediction error
      loocv_error <- (y_pred - y_test)^2
      
      # store
      errors[i] <- loocv_error
      
    }
    
    MSE[k-3] <- mean(errors)
    
  }
  
  # find df which minimise mean squared error in LOOCV 
  best_df <- which.min(MSE)
  
  # refit based on best df
  df_range <- min_df:max_df
  df <- df_range[best_df]
  
  # create the I-spline
  ispline <- splines2::iSpline(sample.size, df = df, degree = degree, intercept = TRUE)
  
  ispline <- cbind(1, ispline)
  
  # creating `D` matrix and `d` vector.
  d_mat <- crossprod(ispline, ispline)
  d_vec <- crossprod(ispline, power)
  
  a_mat <- diag(1, ncol(ispline))
  
  b_vec <- rep(0, ncol(ispline))
  
  # freeing the first parameter/intercept
  a_mat[1, 1] <- 0
  
  # solve
  qp.solve <- quadprog::solve.QP(Dmat = d_mat, dvec = d_vec, Amat = t(a_mat), bvec = b_vec)$solution
  
  # full range of sample sizes
  Ns <- minSampleSize:maxSampleSize
  
  # create I-spline basis for the full sample size range
  full_spline <- splines2::iSpline(Ns, df = df, degree = degree, intercept = TRUE)
  
  # add the same intercept column used before
  full_spline <- cbind(1, full_spline)
  
  # interpolated statistics
  g_hat <- full_spline %*% qp.solve
  
  # how to return results
  if(is.numeric(correlation)) {
    
    d <- as.data.frame(cbind(Ns,g_hat))
    d$correlation = rep(correlation, nrow(d))
    d$method <- rep(method, nrow(d))
    
    colnames(d) <- c("sample size","power","correlation","method")
    
  } else {
    
    d <- as.data.frame(cbind(Ns,g_hat))
    d$method <- rep(method, nrow(d))
    
    colnames(d) <- c("sample size","power","method")
  }
  
  
  return(d)
  
}

# function for plotting power curves 
plot_power <- function(df, title, coefs) {
  
  g <- ggplot(transform(df, power_plot = pmin(power, 1)),
              aes(x = `sample size`,
                  y = power_plot,
                  color = method,
                  linetype = correlation,
                  group = interaction(method, correlation))) +
    geom_line(size = 1) +
    geom_hline(yintercept=0.05, color="black") + 
    scale_color_manual(
      values = c(
        "Distance" = wes_palette("FantasticFox1")[1],
        "Pearson" = wes_palette("FantasticFox1")[2],
        "Spearman" = wes_palette("FantasticFox1")[3],
        "TIC" = wes_palette("AsteroidCity1")[4],
        "Mutual" = wes_palette("FantasticFox1")[5]
      )
    ) +
    scale_linetype_manual(
      values = setNames(
        c("dotted", "longdash","solid"),
        as.character(coefs)
      )
    ) +
    guides(
      color = guide_legend(order = 1),
      linetype = guide_legend(order = 2)
    ) +
    scale_y_continuous(limits = c(0, 1), breaks=seq(0,1,0.1)) +
    scale_x_continuous(limits = c(10, 100), breaks=seq(0,100,10)) +
    theme_classic() +
    labs(title = title,
         x = "Sample size",
         y = "Power",
         color = "Method",
         linetype = "Coefficient")
  
  return(g)
  
}



#------------------------------------------------------------------------------#
#                                simulations                                   #
#------------------------------------------------------------------------------#


#---------------- beta to beta ---------------------------#

sample.size <- seq(10,100,10)
p <- c(0,1,2.5)
reps=500

#run sim using parSim
beta.beta.sim <- parSim(
  
  sample.size = sample.size,
  p = p,
  reps = reps,
  nCores = 10,
  export = c("node_beta","mi.perm","mi.test","estimate.mi"),
  
  expression = {
    
    library(simDAG)
    
    dag <- empty_dag() +
      node("x", type="rbeta", shape1=1, shape2=6) +
      node("y", type=node_beta, parents="x", betas=p, intercept=1, phi=40)
    
    data <- as.data.frame(sim_from_dag(dag, sample.size))
    
    x <- data[,1]
    y <- data[,2]
    
    dc <- energy::dcor.test(x, y, R=500)
    dc.p <- dc$p.value
    
    pc <- cor.test(x, y, method="pearson")
    pc.p <- pc$p.value
    
    sc <- coin::spearman_test(y ~ x, distribution=coin::approximate(nresample=500)) 
    sc.p <- coin::pvalue(sc)[1]
    
    mi <- mi.test(x, y, type=c("c","c"), permutations=500)
    mi.p <- mi$p.value
    
    tic <- minerva::mictools(as.matrix(data), nperm=500, seed = as.integer(runif(1, min = 0, max = 2^31-1)), p.adjust.method="none")
    tic.p <- tic$pval$pval
    
    data.frame(
      dc.p = dc.p,
      pc.p = pc.p,
      sc.p = sc.p,
      mi.p = mi.p,
      tic.p = tic.p)
  }
)


# create dataframe to store power at each sample size
beta.beta.power <- as.data.frame(matrix(0, nrow=length(sample.size)*length(p), ncol = 7))
colnames(beta.beta.power) <- c("sample.size","correlation","pearson.power","spearman.power","distance.power","tic.power","mutual.power")

beta.beta.power$sample.size <- rep(sample.size, length(p))
beta.beta.power$correlation <- rep(p, 1, each=length(sample.size))

# loop to estimate power at each sample size for each dependence measure
for(i in 1:nrow(beta.beta.power)) {
  
  beta.beta.power$pearson.power[i] <- get.power(data=beta.beta.sim, sample.size=beta.beta.power$sample.size[i],
                                                p=beta.beta.power$correlation[i], type="pc.p", alpha=0.05)
  
  beta.beta.power$spearman.power[i] <- get.power(data=beta.beta.sim, sample.size=beta.beta.power$sample.size[i],
                                                 p=beta.beta.power$correlation[i], type="sc.p", alpha=0.05)
  
  beta.beta.power$distance.power[i] <- get.power(data=beta.beta.sim, sample.size=beta.beta.power$sample.size[i],
                                                 p=beta.beta.power$correlation[i], type="dc.p",alpha=0.05)
  
  beta.beta.power$mutual.power[i] <- get.power(data=beta.beta.sim, sample.size=beta.beta.power$sample.size[i],
                                               p=beta.beta.power$correlation[i], type="mi.p",alpha=0.05)
  
  beta.beta.power$tic.power[i] <- get.power(data=beta.beta.sim, sample.size=beta.beta.power$sample.size[i],
                                            p=beta.beta.power$correlation[i], type="tic.p",alpha=0.05)
  
}


#fit power curves to estimated power for each dependence strength for each dependence measure
ind <- which(beta.beta.power$correlation == 0)
pearson.beta.beta.type1.error <- spline.fun(sample.size=sample.size, power=beta.beta.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Pearson")
spearman.beta.beta.type1.error <- spline.fun(sample.size=sample.size, power=beta.beta.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Spearman")
distance.beta.beta.type1.error <- spline.fun(sample.size=sample.size, power=beta.beta.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Distance")
tic.beta.beta.type1.error <- spline.fun(sample.size=sample.size, power=beta.beta.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "TIC")
mutual.beta.beta.type1.error <- spline.fun(sample.size=sample.size, power=beta.beta.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Mutual")


ind <- which(beta.beta.power$correlation == 1)
pearson.beta.beta.low.power <- spline.fun(sample.size=sample.size, power=beta.beta.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1, method = "Pearson")
spearman.beta.beta.low.power <- spline.fun(sample.size=sample.size, power=beta.beta.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1, method = "Spearman")
distance.beta.beta.low.power <- spline.fun(sample.size=sample.size, power=beta.beta.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1, method = "Distance")
tic.beta.beta.low.power <- spline.fun(sample.size=sample.size, power=beta.beta.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1, method = "TIC")
mutual.beta.beta.low.power <- spline.fun(sample.size=sample.size, power=beta.beta.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1, method = "Mutual")

ind <- which(beta.beta.power$correlation == 2.5)
pearson.beta.beta.high.power <- spline.fun(sample.size=sample.size, power=beta.beta.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 2.5, method = "Pearson")
spearman.beta.beta.high.power <- spline.fun(sample.size=sample.size, power=beta.beta.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 2.5, method = "Spearman")
distance.beta.beta.high.power <- spline.fun(sample.size=sample.size, power=beta.beta.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 2.5, method = "Distance")
tic.beta.beta.high.power <- spline.fun(sample.size=sample.size, power=beta.beta.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 2.5, method = "TIC")
mutual.beta.beta.high.power <- spline.fun(sample.size=sample.size, power=beta.beta.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 2.5, method = "Mutual")


# combined for plotting
beta.beta.df <- rbind(pearson.beta.beta.type1.error,pearson.beta.beta.low.power,pearson.beta.beta.high.power,
                      spearman.beta.beta.type1.error,spearman.beta.beta.low.power,spearman.beta.beta.high.power,
                      distance.beta.beta.type1.error,distance.beta.beta.low.power,distance.beta.beta.high.power,
                      tic.beta.beta.type1.error,tic.beta.beta.low.power,tic.beta.beta.high.power,
                      mutual.beta.beta.type1.error, mutual.beta.beta.low.power, mutual.beta.beta.high.power)

beta.beta.df$correlation <- as.factor(beta.beta.df$correlation)
beta.beta.df$method <- as.factor(beta.beta.df$method)

# figure
beta.beta.plot <- plot_power(df = beta.beta.df, title = "(Beta,Beta)", coefs = c(0,1,2.5))
beta.beta.plot


#------------ gamma to beta ------------------#

sample.size <- seq(10,100,10)
p <- c(0,0.4,0.8)
reps=500

#run sim using parSim
gamma.beta.sim <- parSim(
  
  sample.size = sample.size,
  p = p,
  reps = reps,
  nCores = 10,
  export = c("node_beta","mi.perm","mi.test","estimate.mi"),
  
  expression = {
    
    library(simDAG)
    
    dag <- empty_dag() +
      node("x", type="rgamma", shape=1, rate=2) +
      node("y", type=node_beta, parents="x", betas=p, intercept=1, phi=20)
    
    data <- as.data.frame(sim_from_dag(dag, sample.size))
    
    x <- data[,1]
    y <- data[,2]
    
    dc <- energy::dcor.test(x, y, R=500)
    dc.p <- dc$p.value
    
    pc <- cor.test(x, y, method="pearson")
    pc.p <- pc$p.value
    
    sc <- coin::spearman_test(y ~ x, distribution=coin::approximate(nresample=500)) 
    sc.p <- coin::pvalue(sc)[1]
    
    mi <- mi.test(x, y, type=c("c","c"), permutations=500)
    mi.p <- mi$p.value
    
    tic <- minerva::mictools(as.matrix(data), nperm=500, seed = as.integer(runif(1, min = 0, max = 2^31-1)), p.adjust.method="none")
    tic.p <- tic$pval$pval
    
    data.frame(
      dc.p = dc.p,
      pc.p = pc.p,
      sc.p = sc.p,
      mi.p = mi.p,
      tic.p = tic.p)
    
  }
)



gamma.beta.power <- as.data.frame(matrix(0, nrow=length(sample.size)*length(p), ncol = 7))
colnames(gamma.beta.power) <- c("sample.size","correlation","pearson.power","spearman.power","distance.power","tic.power","mutual.power")

gamma.beta.power$sample.size <- rep(sample.size, length(p))
gamma.beta.power$correlation <- rep(p, 1, each=length(sample.size))

for(i in 1:nrow(gamma.beta.power)) {
  
  gamma.beta.power$pearson.power[i] <- get.power(data=gamma.beta.sim, sample.size=gamma.beta.power$sample.size[i],
                                                 p=gamma.beta.power$correlation[i], type="pc.p", alpha=0.05)
  
  gamma.beta.power$spearman.power[i] <- get.power(data=gamma.beta.sim, sample.size=gamma.beta.power$sample.size[i],
                                                  p=gamma.beta.power$correlation[i], type="sc.p", alpha=0.05)
  
  gamma.beta.power$distance.power[i] <- get.power(data=gamma.beta.sim, sample.size=gamma.beta.power$sample.size[i],
                                                  p=gamma.beta.power$correlation[i], type="dc.p", alpha=0.05)
  
  gamma.beta.power$tic.power[i] <- get.power(data=gamma.beta.sim, sample.size=gamma.beta.power$sample.size[i],
                                             p=gamma.beta.power$correlation[i], type="tic.p", alpha=0.05)
  
  gamma.beta.power$mutual.power[i] <- get.power(data=gamma.beta.sim, sample.size=gamma.beta.power$sample.size[i],
                                                p=gamma.beta.power$correlation[i], type="mi.p", alpha=0.05)
  
}


#fit power curves
ind <- which(gamma.beta.power$correlation == 0)
pearson.gamma.beta.type1.error <- spline.fun(sample.size=sample.size, power=gamma.beta.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Pearson")
spearman.gamma.beta.type1.error <- spline.fun(sample.size=sample.size, power=gamma.beta.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Spearman")
distance.gamma.beta.type1.error <- spline.fun(sample.size=sample.size, power=gamma.beta.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Distance")
tic.gamma.beta.type1.error <- spline.fun(sample.size=sample.size, power=gamma.beta.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "TIC")
mutual.gamma.beta.type1.error <- spline.fun(sample.size=sample.size, power=gamma.beta.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Mutual")

ind <- which(gamma.beta.power$correlation == 0.4)
pearson.gamma.beta.low.power <- spline.fun(sample.size=sample.size, power=gamma.beta.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.4, method = "Pearson")
spearman.gamma.beta.low.power <- spline.fun(sample.size=sample.size, power=gamma.beta.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.4, method = "Spearman")
distance.gamma.beta.low.power <- spline.fun(sample.size=sample.size, power=gamma.beta.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.4, method = "Distance")
tic.gamma.beta.low.power <- spline.fun(sample.size=sample.size, power=gamma.beta.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.4, method = "TIC")
mutual.gamma.beta.low.power <- spline.fun(sample.size=sample.size, power=gamma.beta.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.4, method = "Mutual")


ind <- which(gamma.beta.power$correlation == 0.8)
pearson.gamma.beta.high.power <- spline.fun(sample.size=sample.size, power=gamma.beta.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.8, method = "Pearson")
spearman.gamma.beta.high.power <- spline.fun(sample.size=sample.size, power=gamma.beta.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.8, method = "Spearman")
distance.gamma.beta.high.power <- spline.fun(sample.size=sample.size, power=gamma.beta.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.8, method = "Distance")
tic.gamma.beta.high.power <- spline.fun(sample.size=sample.size, power=gamma.beta.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.8, method = "TIC")
mutual.gamma.beta.high.power <- spline.fun(sample.size=sample.size, power=gamma.beta.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.8, method = "Mutual")

gamma.beta.df <- rbind(pearson.gamma.beta.type1.error,pearson.gamma.beta.low.power,pearson.gamma.beta.high.power,
                       spearman.gamma.beta.type1.error,spearman.gamma.beta.low.power,spearman.gamma.beta.high.power,
                       distance.gamma.beta.type1.error,distance.gamma.beta.low.power,distance.gamma.beta.high.power,
                       tic.gamma.beta.type1.error,tic.gamma.beta.low.power,tic.gamma.beta.high.power,
                       mutual.gamma.beta.type1.error, mutual.gamma.beta.low.power, mutual.gamma.beta.high.power)

gamma.beta.df$correlation <- as.factor(gamma.beta.df$correlation)


gamma.beta.plot <- plot_power(df = gamma.beta.df, title = "(Gamma,Beta)", coefs = c(0,0.4,0.8))
gamma.beta.plot


#---------------- binary to beta  ---------------------------#

sample.size <- seq(10,100,10)
p <- c(0,0.4,0.8)
reps=500

#run sim using parSim
binary.beta.sim <- parSim(
  
  sample.size = sample.size,
  p = p,
  reps = reps,
  nCores = 10,
  export = c("node_beta","mi.perm","mi.test","estimate.mi"),
  
  expression = {
    
    library(simDAG)
    
    dag <- empty_dag() +
      node("x", type="rbinom", size=1, prob=0.5) +
      node("y", type=node_beta, parents="x", betas=p, intercept=1, phi=10)
    
    data <- as.data.frame(sim_from_dag(dag, sample.size))
    
    x <- data[,1]
    y <- data[,2]
    
    dc <- energy::dcor.test(x, y, R=500)
    dc.p <- dc$p.value
    
    pc <- cor.test(x, y, method="pearson")
    pc.p <- pc$p.value
    
    sc <- coin::spearman_test(y ~ x, distribution=coin::approximate(nresample=500)) 
    sc.p <- coin::pvalue(sc)[1]
    
    mi <- mi.test(x, y, type=c("d","c"), permutations=500)
    mi.p <- mi$p.value
    
    tic <- minerva::mictools(as.matrix(cbind(x,y)), seed = as.integer(runif(1, min = 0, max = 2^31-1)), nperm=500, p.adjust.method="none")
    tic.p <- tic$pval$pval
    
    data.frame(
      dc.p = dc.p,
      pc.p = pc.p,
      sc.p = sc.p,
      mi.p = mi.p,
      tic.p = tic.p)
    
  }
)



binary.beta.power <- as.data.frame(matrix(0, nrow=length(sample.size)*length(p), ncol = 7))
colnames(binary.beta.power) <- c("sample.size","correlation","pearson.power","spearman.power","distance.power","mutual.power","tic.power")

binary.beta.power$sample.size <- rep(sample.size, length(p))
binary.beta.power$correlation <- rep(p, 1, each=length(sample.size))

for(i in 1:nrow(binary.beta.power)) {
  
  binary.beta.power$pearson.power[i] <- get.power(data=binary.beta.sim, sample.size=binary.beta.power$sample.size[i],
                                                  p=binary.beta.power$correlation[i], type="pc.p", alpha=0.05)
  
  binary.beta.power$spearman.power[i] <- get.power(data=binary.beta.sim, sample.size=binary.beta.power$sample.size[i],
                                                   p=binary.beta.power$correlation[i], type="sc.p", alpha=0.05)
  
  binary.beta.power$distance.power[i] <- get.power(data=binary.beta.sim, sample.size=binary.beta.power$sample.size[i],
                                                   p=binary.beta.power$correlation[i], type = "dc.p", alpha=0.05)
  
  binary.beta.power$tic.power[i] <- get.power(data=binary.beta.sim, sample.size=binary.beta.power$sample.size[i],
                                              p=binary.beta.power$correlation[i], type = "tic.p", alpha=0.05)
  
  binary.beta.power$mutual.power[i] <- get.power(data=binary.beta.sim, sample.size=binary.beta.power$sample.size[i],
                                                 p=binary.beta.power$correlation[i], type = "mi.p", alpha=0.05)
}


#fit power curves
ind <- which(binary.beta.power$correlation == 0)
pearson.binary.beta.type1.error <- spline.fun(sample.size=sample.size, power=binary.beta.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Pearson")
spearman.binary.beta.type1.error <- spline.fun(sample.size=sample.size, power=binary.beta.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Spearman")
distance.binary.beta.type1.error <- spline.fun(sample.size=sample.size, power=binary.beta.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Distance")
tic.binary.beta.type1.error <- spline.fun(sample.size=sample.size, power=binary.beta.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "TIC")
mutual.binary.beta.type1.error <- spline.fun(sample.size=sample.size, power=binary.beta.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Mutual")


ind <- which(binary.beta.power$correlation == 0.4)

pearson.binary.beta.low.power <- spline.fun(sample.size=sample.size, power=binary.beta.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.4, method = "Pearson")
spearman.binary.beta.low.power <- spline.fun(sample.size=sample.size, power=binary.beta.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.4, method = "Spearman")
distance.binary.beta.low.power <- spline.fun(sample.size=sample.size, power=binary.beta.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.4, method = "Distance")
mutual.binary.beta.low.power <- spline.fun(sample.size=sample.size, power=binary.beta.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.4, method = "Mutual")
tic.binary.beta.low.power <- spline.fun(sample.size=sample.size, power=binary.beta.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.4, method = "TIC")

ind <- which(binary.beta.power$correlation == 0.8)
pearson.binary.beta.high.power <- spline.fun(sample.size=sample.size, power=binary.beta.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.8, method = "Pearson")
spearman.binary.beta.high.power <- spline.fun(sample.size=sample.size, power=binary.beta.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.8, method = "Spearman")
distance.binary.beta.high.power <- spline.fun(sample.size=sample.size, power=binary.beta.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.8, method = "Distance")
mutual.binary.beta.high.power <- spline.fun(sample.size=sample.size, power=binary.beta.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.8, method = "Mutual")
tic.binary.beta.high.power <- spline.fun(sample.size=sample.size, power=binary.beta.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.8, method = "TIC")


binary.beta.df <- rbind(pearson.binary.beta.type1.error,pearson.binary.beta.low.power,pearson.binary.beta.high.power,
                        spearman.binary.beta.type1.error, spearman.binary.beta.low.power, spearman.binary.beta.high.power,
                        distance.binary.beta.type1.error,distance.binary.beta.low.power,distance.binary.beta.high.power,
                        tic.binary.beta.type1.error, tic.binary.beta.low.power, tic.binary.beta.high.power,
                        mutual.binary.beta.type1.error,mutual.binary.beta.low.power,mutual.binary.beta.high.power)

binary.beta.df$correlation <- as.factor(binary.beta.df$correlation)

binary.beta.plot <- plot_power(df = binary.beta.df, 
                               title = "(Binary,Beta)", 
                               coefs = c(0,0.4,0.8))
binary.beta.plot


#---------------- negative binomial to beta ---------------------------#

sample.size <- seq(10,100,10)
p <- c(0,0.05,0.1)
reps=500

#run sim using parSim
nbinom.beta.sim <- parSim(
  
  sample.size = sample.size,
  p = p,
  reps = reps,
  nCores = 10,
  export = c("node_beta","mi.perm","mi.test","estimate.mi"),
  
  expression = {
    
    library(simDAG)
    
    dag <- empty_dag() +
      node("x", type="rnbinom", size=0.9, mu=5) +
      node("y", type=node_beta, parents="x", betas=p, intercept=-2, phi=10)
    
    data <- as.data.frame(sim_from_dag(dag, sample.size))
    
    x <- data[,1]
    y <- data[,2]
    
    dc <- energy::dcor.test(x, y, R=500)
    dc.p <- dc$p.value
    
    pc <- cor.test(x, y, method="pearson")
    pc.p <- pc$p.value
    
    sc <- coin::spearman_test(y ~ x, distribution=coin::approximate(nresample=500)) 
    sc.p <- coin::pvalue(sc)[1]
    
    mi <- mi.test(x, y, type=c("d","c"), permutations=500)
    mi.p <- mi$p.value
    
    tic <- minerva::mictools(as.matrix(cbind(x,y)), seed = as.integer(runif(1, min = 0, max = 2^31-1)), nperm=500, p.adjust.method="none")
    tic.p <- tic$pval$pval
    
    data.frame(
      dc.p = dc.p,
      pc.p = pc.p,
      sc.p = sc.p,
      mi.p = mi.p,
      tic.p = tic.p)
    
  }
)



nbinom.beta.power <- as.data.frame(matrix(0, nrow=length(sample.size)*length(p), ncol = 7))
colnames(nbinom.beta.power) <- c("sample.size","correlation","pearson.power","spearman.power","distance.power","mutual.power","tic.power")

nbinom.beta.power$sample.size <- rep(sample.size, length(p))
nbinom.beta.power$correlation <- rep(p, 1, each=length(sample.size))

for(i in 1:nrow(nbinom.beta.power)) {
  
  nbinom.beta.power$pearson.power[i] <- get.power(data=nbinom.beta.sim, sample.size=nbinom.beta.power$sample.size[i],
                                                  p=nbinom.beta.power$correlation[i], type="pc.p", alpha=0.05)
  
  nbinom.beta.power$spearman.power[i] <- get.power(data=nbinom.beta.sim, sample.size=nbinom.beta.power$sample.size[i],
                                                   p=nbinom.beta.power$correlation[i], type="sc.p", alpha=0.05)
  
  nbinom.beta.power$distance.power[i] <- get.power(data=nbinom.beta.sim, sample.size=nbinom.beta.power$sample.size[i],
                                                   p=nbinom.beta.power$correlation[i], type="dc.p", alpha=0.05)
  
  nbinom.beta.power$tic.info.power[i] <- get.power(data=nbinom.beta.sim, sample.size=nbinom.beta.power$sample.size[i],
                                                   p=nbinom.beta.power$correlation[i], type="tic.p", alpha=0.05)
  
  
  nbinom.beta.power$mutual.info.power[i] <- get.power(data=nbinom.beta.sim, sample.size=nbinom.beta.power$sample.size[i],
                                                      p=nbinom.beta.power$correlation[i], type="mi.p", alpha=0.05)
}


#fit power curves
ind <- which(nbinom.beta.power$correlation == 0)
pearson.nbinom.beta.type1.error <- spline.fun(sample.size=sample.size, power=nbinom.beta.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Pearson")
spearman.nbinom.beta.type1.error <- spline.fun(sample.size=sample.size, power=nbinom.beta.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Spearman")
distance.nbinom.beta.type1.error <- spline.fun(sample.size=sample.size, power=nbinom.beta.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Distance")
mutual.nbinom.beta.type1.error <- spline.fun(sample.size=sample.size, power=nbinom.beta.power$mutual.info.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Mutual")
tic.nbinom.beta.type1.error <- spline.fun(sample.size=sample.size, power=nbinom.beta.power$tic.info.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "TIC")


ind <- which(nbinom.beta.power$correlation == 0.05)
pearson.nbinom.beta.low.power <- spline.fun(sample.size=sample.size, power=nbinom.beta.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.05, method = "Pearson")
spearman.nbinom.beta.low.power <- spline.fun(sample.size=sample.size, power=nbinom.beta.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.05, method = "Spearman")
distance.nbinom.beta.low.power <- spline.fun(sample.size=sample.size, power=nbinom.beta.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.05, method = "Distance")
mutual.nbinom.beta.low.power <- spline.fun(sample.size=sample.size, power=nbinom.beta.power$mutual.info.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.05, method = "Mutual")
tic.nbinom.beta.low.power <- spline.fun(sample.size=sample.size, power=nbinom.beta.power$tic.info.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.05, method = "TIC")


ind <- which(nbinom.beta.power$correlation == 0.1)
pearson.nbinom.beta.high.power <- spline.fun(sample.size=sample.size, power=nbinom.beta.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.1, method = "Pearson")
spearman.nbinom.beta.high.power <- spline.fun(sample.size=sample.size, power=nbinom.beta.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.1, method = "Spearman")
distance.nbinom.beta.high.power <- spline.fun(sample.size=sample.size, power=nbinom.beta.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.1, method = "Distance")
mutual.nbinom.beta.high.power <- spline.fun(sample.size=sample.size, power=nbinom.beta.power$mutual.info.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.1, method = "Mutual")
tic.nbinom.beta.high.power <- spline.fun(sample.size=sample.size, power=nbinom.beta.power$tic.info.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.1, method = "TIC")

nbinom.beta.df <- rbind(pearson.nbinom.beta.type1.error,pearson.nbinom.beta.low.power,pearson.nbinom.beta.high.power,
                        spearman.nbinom.beta.type1.error,spearman.nbinom.beta.low.power,spearman.nbinom.beta.high.power,
                        distance.nbinom.beta.type1.error,distance.nbinom.beta.low.power,distance.nbinom.beta.high.power,
                        mutual.nbinom.beta.type1.error,mutual.nbinom.beta.low.power,mutual.nbinom.beta.high.power,
                        tic.nbinom.beta.type1.error,tic.nbinom.beta.low.power,tic.nbinom.beta.high.power)
nbinom.beta.df$correlation <- as.factor(nbinom.beta.df$correlation)


nbinom.beta.plot <- plot_power(df = nbinom.beta.df, 
                               title = "(NegBinomial,Beta)",
                               coefs = c(0,0.05,0.1))
nbinom.beta.plot


#---------------- negative binomial to negative binomial ---------------------------#


sample.size <- seq(10,100,10)
p <- c(0, 0.04, 0.08)
reps=500

#run sim using parSim
nbinom.nbinom.sim <- parSim(
  
  sample.size = sample.size,
  p = p,
  reps = reps,
  nCores = 10,
  export = c("mi.perm","mi.test","estimate.mi"),
  
  expression = {
    
    library(simDAG)
    
    dag <- empty_dag() +
      node("x", type="rnbinom", size=0.9, mu=5) +
      node("y", type="negative_binomial", parents="x", betas=p, intercept=1, theta=10)
    data <- as.data.frame(sim_from_dag(dag, sample.size))
    
    x <- data[,1]
    y <- data[,2]
    
    dc <- energy::dcor.test(x, y, R=500)
    dc.p <- dc$p.value
    
    pc <- cor.test(x, y, method="pearson")
    pc.p <- pc$p.value
    
    sc <- coin::spearman_test(y ~ x, distribution=coin::approximate(nresample=500)) 
    sc.p <- coin::pvalue(sc)[1]
    
    mi <- mi.test(x, y, type=c("d","d"), permutations=500)
    mi.p <- mi$p.value
    
    tic <- minerva::mictools(as.matrix(cbind(x,y)), seed = as.integer(runif(1, min = 0, max = 2^31-1)), nperm=500, p.adjust.method="none")
    tic.p <- tic$pval$pval
    
    data.frame(
      dc.p = dc.p,
      pc.p = pc.p,
      sc.p = sc.p,
      mi.p = mi.p,
      tic.p = tic.p)
  }
)



nbinom.nbinom.power <- as.data.frame(matrix(0, nrow=length(sample.size)*length(p), ncol = 7))
colnames(nbinom.nbinom.power) <- c("sample.size","correlation","pearson.power","spearman.power","distance.power","mutual.power","tic.power")

nbinom.nbinom.power$sample.size <- rep(sample.size, length(p))
nbinom.nbinom.power$correlation <- rep(p, 1, each=length(sample.size))

for(i in 1:nrow(nbinom.nbinom.power)) {
  
  nbinom.nbinom.power$pearson.power[i] <- get.power(data=nbinom.nbinom.sim, sample.size=nbinom.nbinom.power$sample.size[i],
                                                    p=nbinom.nbinom.power$correlation[i], type="pc.p", alpha=0.05)
  
  nbinom.nbinom.power$spearman.power[i] <- get.power(data=nbinom.nbinom.sim, sample.size=nbinom.nbinom.power$sample.size[i],
                                                     p=nbinom.nbinom.power$correlation[i], type="sc.p", alpha=0.05)
  
  
  nbinom.nbinom.power$distance.power[i] <- get.power(data=nbinom.nbinom.sim, sample.size=nbinom.nbinom.power$sample.size[i],
                                                     p=nbinom.nbinom.power$correlation[i], type="dc.p", alpha=0.05)
  
  nbinom.nbinom.power$mutual.power[i] <- get.power(data=nbinom.nbinom.sim, sample.size=nbinom.nbinom.power$sample.size[i],
                                                   p=nbinom.nbinom.power$correlation[i], type="mi.p", alpha=0.05)
  
  nbinom.nbinom.power$tic.power[i] <- get.power(data=nbinom.nbinom.sim, sample.size=nbinom.nbinom.power$sample.size[i],
                                                p=nbinom.nbinom.power$correlation[i], type="tic.p", alpha=0.05)
  
}


#fit power curves
ind <- which(nbinom.nbinom.power$correlation == 0)
pearson.nbinom.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=nbinom.nbinom.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Pearson")
spearman.nbinom.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=nbinom.nbinom.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Spearman")
distance.nbinom.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=nbinom.nbinom.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Distance")
mutual.nbinom.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=nbinom.nbinom.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Mutual")
tic.nbinom.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=nbinom.nbinom.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "TIC")


ind <- which(nbinom.nbinom.power$correlation == 0.04)
pearson.nbinom.nbinom.low.power <- spline.fun(sample.size=sample.size, power=nbinom.nbinom.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.04, method = "Pearson")
spearman.nbinom.nbinom.low.power <- spline.fun(sample.size=sample.size, power=nbinom.nbinom.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.04, method = "Spearman")
distance.nbinom.nbinom.low.power <- spline.fun(sample.size=sample.size, power=nbinom.nbinom.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.04, method = "Distance")
mutual.nbinom.nbinom.low.power <- spline.fun(sample.size=sample.size, power=nbinom.nbinom.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.04, method = "Mutual")
tic.nbinom.nbinom.low.power <- spline.fun(sample.size=sample.size, power=nbinom.nbinom.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.04, method = "TIC")


ind <- which(nbinom.nbinom.power$correlation == 0.08)
pearson.nbinom.nbinom.high.power <- spline.fun(sample.size=sample.size, power=nbinom.nbinom.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.08, method = "Pearson")
spearman.nbinom.nbinom.high.power <- spline.fun(sample.size=sample.size, power=nbinom.nbinom.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.08, method = "Spearman")
distance.nbinom.nbinom.high.power <- spline.fun(sample.size=sample.size, power=nbinom.nbinom.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.08, method = "Distance")
mutual.nbinom.nbinom.high.power <- spline.fun(sample.size=sample.size, power=nbinom.nbinom.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.08, method = "Mutual")
tic.nbinom.nbinom.high.power <- spline.fun(sample.size=sample.size, power=nbinom.nbinom.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.08, method = "TIC")

nbinom.nbinom.df <- rbind(pearson.nbinom.nbinom.type1.error,pearson.nbinom.nbinom.low.power,pearson.nbinom.nbinom.high.power,
                          spearman.nbinom.nbinom.type1.error,spearman.nbinom.nbinom.low.power,spearman.nbinom.nbinom.high.power,
                          distance.nbinom.nbinom.type1.error,distance.nbinom.nbinom.low.power,distance.nbinom.nbinom.high.power,
                          mutual.nbinom.nbinom.type1.error,mutual.nbinom.nbinom.low.power,mutual.nbinom.nbinom.high.power,
                          tic.nbinom.nbinom.type1.error,tic.nbinom.nbinom.low.power,tic.nbinom.nbinom.high.power)
nbinom.nbinom.df$correlation <- as.factor(nbinom.nbinom.df$correlation)

nbinom.nbinom.plot <- plot_power(df = nbinom.nbinom.df,
                                 title = c("(NegBinomial,NegBinomial"),
                                 coefs = c(0,0.04,0.08))
nbinom.nbinom.plot


#---------------- gamma to negative binomial ---------------------------#


sample.size <- seq(10,100,10)
p <- c(0,0.1,0.2)
reps=500

#run sim using parSim
gamma.nbinom.sim <- parSim(
  
  sample.size = sample.size,
  p = p,
  reps = reps,
  nCores = 10,
  export = c("mi.perm","mi.test","estimate.mi"),
  
  expression = {
    
    library(simDAG)
    
    dag <- empty_dag() +
      node("x", type="rgamma", shape=3, rate=1) +
      node("y", type="negative_binomial", parents="x", betas=p, intercept=1, theta=50)
    
    data <- as.data.frame(sim_from_dag(dag, sample.size))
    
    x <- data[,1]
    y <- data[,2]
    
    dc <- energy::dcor.test(x, y, R=500)
    dc.p <- dc$p.value
    
    pc <- cor.test(x, y, method="pearson")
    pc.p <- pc$p.value
    
    sc <- coin::spearman_test(y ~ x, distribution=coin::approximate(nresample=500)) 
    sc.p <- coin::pvalue(sc)[1]
    
    mi <- mi.test(x, y, type=c("c","d"), permutations=500)
    mi.p <- mi$p.value
    
    tic <- minerva::mictools(as.matrix(cbind(x,y)), seed = as.integer(runif(1, min = 0, max = 2^31-1)), nperm=500, p.adjust.method="none")
    tic.p <- tic$pval$pval
    
    data.frame(
      dc.p = dc.p,
      pc.p = pc.p,
      sc.p = sc.p,
      mi.p = mi.p,
      tic.p = tic.p)
    
  }
)



gamma.nbinom.power <- as.data.frame(matrix(0, nrow=length(sample.size)*length(p), ncol = 7))
colnames(gamma.nbinom.power) <- c("sample.size","correlation","pearson.power","spearman.power","distance.power","mutual.power","tic.power")

gamma.nbinom.power$sample.size <- rep(sample.size, length(p))
gamma.nbinom.power$correlation <- rep(p, 1, each=length(sample.size))

for(i in 1:nrow(gamma.nbinom.power)) {
  
  gamma.nbinom.power$pearson.power[i] <- get.power(data=gamma.nbinom.sim, sample.size=gamma.nbinom.power$sample.size[i],
                                                   p=gamma.nbinom.power$correlation[i], type="pc.p", alpha=0.05)
  
  gamma.nbinom.power$spearman.power[i] <- get.power(data=gamma.nbinom.sim, sample.size=gamma.nbinom.power$sample.size[i],
                                                    p=gamma.nbinom.power$correlation[i], type="sc.p", alpha=0.05)
  
  gamma.nbinom.power$distance.power[i] <- get.power(data=gamma.nbinom.sim, sample.size=gamma.nbinom.power$sample.size[i],
                                                    p=gamma.nbinom.power$correlation[i], type="dc.p", alpha=0.05)
  
  gamma.nbinom.power$mutual.power[i] <- get.power(data=gamma.nbinom.sim, sample.size=gamma.nbinom.power$sample.size[i],
                                                  p=gamma.nbinom.power$correlation[i], type="mi.p", alpha=0.05)
  
  gamma.nbinom.power$tic.power[i] <- get.power(data=gamma.nbinom.sim, sample.size=gamma.nbinom.power$sample.size[i],
                                               p=gamma.nbinom.power$correlation[i], type="tic.p", alpha=0.05)
  
}


#fit power curves
ind <- which(gamma.nbinom.power$correlation == 0)
pearson.gamma.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=gamma.nbinom.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Pearson")
spearman.gamma.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=gamma.nbinom.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Spearman")
distance.gamma.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=gamma.nbinom.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Distance")
mutual.gamma.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=gamma.nbinom.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Mutual")
tic.gamma.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=gamma.nbinom.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "TIC")


ind <- which(gamma.nbinom.power$correlation == 0.1)
pearson.gamma.nbinom.low.power <- spline.fun(sample.size=sample.size, power=gamma.nbinom.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.1, method = "Pearson")
spearman.gamma.nbinom.low.power <- spline.fun(sample.size=sample.size, power=gamma.nbinom.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.1, method = "Spearman")
distance.gamma.nbinom.low.power <- spline.fun(sample.size=sample.size, power=gamma.nbinom.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.1, method = "Distance")
mutual.gamma.nbinom.low.power <- spline.fun(sample.size=sample.size, power=gamma.nbinom.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.1, method = "Mutual")
tic.gamma.nbinom.low.power <- spline.fun(sample.size=sample.size, power=gamma.nbinom.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.1, method = "TIC")


ind <- which(gamma.nbinom.power$correlation == 0.2)
pearson.gamma.nbinom.high.power <- spline.fun(sample.size=sample.size, power=gamma.nbinom.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.2, method = "Pearson")
spearman.gamma.nbinom.high.power <- spline.fun(sample.size=sample.size, power=gamma.nbinom.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.2, method = "Spearman")
distance.gamma.nbinom.high.power <- spline.fun(sample.size=sample.size, power=gamma.nbinom.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.2, method = "Distance")
mutual.gamma.nbinom.high.power <- spline.fun(sample.size=sample.size, power=gamma.nbinom.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.2, method = "Mutual")
tic.gamma.nbinom.high.power <- spline.fun(sample.size=sample.size, power=gamma.nbinom.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.2, method = "TIC")

gamma.nbinom.df <- rbind(pearson.gamma.nbinom.type1.error,pearson.gamma.nbinom.low.power,pearson.gamma.nbinom.high.power,
                         spearman.gamma.nbinom.type1.error,spearman.gamma.nbinom.low.power,spearman.gamma.nbinom.high.power,
                         distance.gamma.nbinom.type1.error,distance.gamma.nbinom.low.power,distance.gamma.nbinom.high.power,
                         mutual.gamma.nbinom.type1.error,mutual.gamma.nbinom.low.power,mutual.gamma.nbinom.high.power,
                         tic.gamma.nbinom.type1.error,tic.gamma.nbinom.low.power,tic.gamma.nbinom.high.power)

gamma.nbinom.df$correlation <- as.factor(gamma.nbinom.df$correlation)

gamma.nbinom.plot <- plot_power(df = gamma.nbinom.df, 
                                title = "(Gamma,NegBinomial)",
                                coefs = c(0,0.1,0.2))
gamma.nbinom.plot


#---------------- Binary to negative binomial ---------------------------#

sample.size <- seq(10,100,10)
p <- c(0,0.3,0.6)
reps=500

#run sim using parSim
binary.nbinom.sim <- parSim(
  
  sample.size = sample.size,
  p = p,
  reps = reps,
  nCores = 10,
  export = c("mi.perm","mi.test","estimate.mi"),
  
  expression = {
    
    library(simDAG)
    
    dag <- empty_dag() +
      node("x", type="rbinom", size=1, prob=0.5) +
      node("y", type="negative_binomial", parents="x", betas=p, intercept=2, theta=10)
    
    data <- as.data.frame(sim_from_dag(dag, sample.size))
    
    x <- data[,1]
    y <- data[,2]
    
    dc <- energy::dcor.test(x, y, R=500)
    dc.p <- dc$p.value
    
    pc <- cor.test(x, y, method="pearson")
    pc.p <- pc$p.value
    
    sc <- coin::spearman_test(y ~ x, distribution=coin::approximate(nresample=500)) 
    sc.p <- coin::pvalue(sc)[1]
    
    mi <- mi.test(x, y, type=c("d","d"), permutations=500)
    mi.p <- mi$p.value
    
    tic <- minerva::mictools(as.matrix(cbind(x,y)), seed = as.integer(runif(1, min = 0, max = 2^31-1)), nperm=500, p.adjust.method="none")
    tic.p <- tic$pval$pval
    
    data.frame(
      dc.p = dc.p,
      pc.p = pc.p,
      sc.p = sc.p,
      mi.p = mi.p,
      tic.p = tic.p)
  }
)



binary.nbinom.power <- as.data.frame(matrix(0, nrow=length(sample.size)*length(p), ncol = 7))
colnames(binary.nbinom.power) <- c("sample.size","correlation","pearson.power","spearman.power","distance.power","mutual.power","tic.power")

binary.nbinom.power$sample.size <- rep(sample.size, length(p))
binary.nbinom.power$correlation <- rep(p, 1, each=length(sample.size))

for(i in 1:nrow(binary.nbinom.power)) {
  
  binary.nbinom.power$pearson.power[i] <- get.power(data=binary.nbinom.sim, sample.size=binary.nbinom.power$sample.size[i],
                                                    p=binary.nbinom.power$correlation[i], type="pc.p", alpha=0.05)
  
  binary.nbinom.power$spearman.power[i] <- get.power(data=binary.nbinom.sim, sample.size=binary.nbinom.power$sample.size[i],
                                                     p=binary.nbinom.power$correlation[i], type="sc.p", alpha=0.05)
  
  binary.nbinom.power$distance.power[i] <- get.power(data=binary.nbinom.sim, sample.size=binary.nbinom.power$sample.size[i],
                                                     p=binary.nbinom.power$correlation[i], type="dc.p", alpha=0.05)
  
  binary.nbinom.power$mutual.power[i] <- get.power(data=binary.nbinom.sim, sample.size=binary.nbinom.power$sample.size[i],
                                                   p=binary.nbinom.power$correlation[i], type="mi.p", alpha=0.05)
  
  binary.nbinom.power$tic.power[i] <- get.power(data=binary.nbinom.sim, sample.size=binary.nbinom.power$sample.size[i],
                                                p=binary.nbinom.power$correlation[i], type="tic.p", alpha=0.05)
  
}


#fit power curves
ind <- which(binary.nbinom.power$correlation == 0)
pearson.binary.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=binary.nbinom.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Pearson")
spearman.binary.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=binary.nbinom.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Spearman")
distance.binary.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=binary.nbinom.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Distance")
mutual.binary.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=binary.nbinom.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Mutual")
tic.binary.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=binary.nbinom.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "TIC")


ind <- which(binary.nbinom.power$correlation == 0.3)
pearson.binary.nbinom.low.power <- spline.fun(sample.size=sample.size, power=binary.nbinom.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.3, method = "Pearson")
spearman.binary.nbinom.low.power <- spline.fun(sample.size=sample.size, power=binary.nbinom.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.3, method = "Spearman")
distance.binary.nbinom.low.power <- spline.fun(sample.size=sample.size, power=binary.nbinom.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.3, method = "Distance")
mutual.binary.nbinom.low.power <- spline.fun(sample.size=sample.size, power=binary.nbinom.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.3, method = "Mutual")
tic.binary.nbinom.low.power <- spline.fun(sample.size=sample.size, power=binary.nbinom.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.3, method = "TIC")


ind <- which(binary.nbinom.power$correlation == 0.6)
pearson.binary.nbinom.high.power <- spline.fun(sample.size=sample.size, power=binary.nbinom.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.6, method = "Pearson")
spearman.binary.nbinom.high.power <- spline.fun(sample.size=sample.size, power=binary.nbinom.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.6, method = "Spearman")
distance.binary.nbinom.high.power <- spline.fun(sample.size=sample.size, power=binary.nbinom.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.6, method = "Distance")
mutual.binary.nbinom.high.power <- spline.fun(sample.size=sample.size, power=binary.nbinom.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.6, method = "Mutual")
tic.binary.nbinom.high.power <- spline.fun(sample.size=sample.size, power=binary.nbinom.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.6, method = "TIC")

binary.nbinom.df <- rbind(pearson.binary.nbinom.type1.error,pearson.binary.nbinom.low.power,pearson.binary.nbinom.high.power,
                          spearman.binary.nbinom.type1.error,spearman.binary.nbinom.low.power,spearman.binary.nbinom.high.power,
                          distance.binary.nbinom.type1.error,distance.binary.nbinom.low.power,distance.binary.nbinom.high.power,
                          mutual.binary.nbinom.type1.error, mutual.binary.nbinom.low.power, mutual.binary.nbinom.high.power,
                          tic.binary.nbinom.type1.error, tic.binary.nbinom.low.power, tic.binary.nbinom.high.power)

binary.nbinom.df$correlation <- as.factor(binary.nbinom.df$correlation)

binary.nbinom.plot <- plot_power(df = binary.nbinom.df,
                                 title = "(Binary,NegBinomial)",
                                 coefs = c(0,0.3,0.6))
binary.nbinom.plot

ind <- which(binary.nbinom.sim$sample.size == 40 & binary.nbinom.sim$p == 0)
length(which(binary.nbinom.sim$dc.p[ind] < 0.05)) / reps
length(which(binary.nbinom.sim$pc.p[ind] < 0.05)) / reps
length(which(binary.nbinom.sim$mi.p[ind] < 0.05)) / reps
length(which(binary.nbinom.sim$tic.p[ind] < 0.05)) / reps


ind <- which(binary.nbinom.df$correlation == 0.4)
plot(binary.nbinom.df$sample.size[ind], binary.nbinom.df$distance.power[ind], ylim=c(0,1))
lines(binary.nbinom.df$sample.size[ind], binary.nbinom.df$distance.power[ind])

points(binary.nbinom.df$sample.size[ind], binary.nbinom.df$pearson.power[ind], pch=16)
lines(binary.nbinom.df$sample.size[ind], binary.nbinom.df$pearson.power[ind], pch=10)


#---------------- Beta to negative binomial  ---------------------------#


sample.size <- seq(10,100,10)
p <- c(0,1.5,3)
reps=500

#run sim using parSim
beta.nbinom.sim <- parSim(
  
  sample.size = sample.size,
  p = p,
  reps = reps,
  nCores = 10,
  export = c("mi.perm","mi.test","estimate.mi"),
  
  expression = {
    
    library(simDAG)
    
    dag <- empty_dag() +
      node("x", type="rbeta", shape1=1, shape2=6) +
      node("y", type="negative_binomial", parents="x", betas=p, intercept=2, theta=3)
    
    data <- as.data.frame(sim_from_dag(dag, sample.size))
    
    x <- data[,1]
    y <- data[,2]
    
    dc <- energy::dcor.test(x, y, R=500)
    dc.p <- dc$p.value
    
    pc <- cor.test(x, y, method="pearson")
    pc.p <- pc$p.value
    
    sc <- coin::spearman_test(y ~ x, distribution=coin::approximate(nresample=500)) 
    sc.p <- coin::pvalue(sc)[1]
    
    mi <- mi.test(x, y, type=c("c","d"), permutations=500)
    mi.p <- mi$p.value
    
    tic <- minerva::mictools(as.matrix(cbind(x,y)), seed = as.integer(runif(1, min = 0, max = 2^31-1)), nperm=500, p.adjust.method="none")
    tic.p <- tic$pval$pval
    
    data.frame(
      dc.p = dc.p,
      pc.p = pc.p,
      sc.p = sc.p,
      mi.p = mi.p,
      tic.p = tic.p)
    
  }
)


beta.nbinom.power <- as.data.frame(matrix(0, nrow=length(sample.size)*length(p), ncol = 7))
colnames(beta.nbinom.power) <- c("sample.size","correlation","pearson.power","spearman.power","distance.power","mutual.power","tic.power")

beta.nbinom.power$sample.size <- rep(sample.size, length(p))
beta.nbinom.power$correlation <- rep(p, 1, each=length(sample.size))

for(i in 1:nrow(beta.nbinom.power)) {
  
  beta.nbinom.power$pearson.power[i] <- get.power(data=beta.nbinom.sim, sample.size=beta.nbinom.power$sample.size[i],
                                                  p=beta.nbinom.power$correlation[i], type="pc.p", alpha=0.05)
  
  beta.nbinom.power$spearman.power[i] <- get.power(data=beta.nbinom.sim, sample.size=beta.nbinom.power$sample.size[i],
                                                   p=beta.nbinom.power$correlation[i], type="sc.p", alpha=0.05)
  
  beta.nbinom.power$distance.power[i] <- get.power(data=beta.nbinom.sim, sample.size=beta.nbinom.power$sample.size[i],
                                                   p=beta.nbinom.power$correlation[i], type="dc.p", alpha=0.05)
  
  beta.nbinom.power$mutual.power[i] <- get.power(data=beta.nbinom.sim, sample.size=beta.nbinom.power$sample.size[i],
                                                 p=beta.nbinom.power$correlation[i], type="mi.p", alpha=0.05)
  
  beta.nbinom.power$tic.power[i] <- get.power(data=beta.nbinom.sim, sample.size=beta.nbinom.power$sample.size[i],
                                              p=beta.nbinom.power$correlation[i], type="tic.p", alpha=0.05)
  
}


#fit power curves
ind <- which(beta.nbinom.power$correlation == 0)
pearson.beta.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=beta.nbinom.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Pearson")
spearman.beta.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=beta.nbinom.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Spearman")
distance.beta.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=beta.nbinom.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Distance")
mutual.beta.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=beta.nbinom.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Mutual")
tic.beta.nbinom.type1.error <- spline.fun(sample.size=sample.size, power=beta.nbinom.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "TIC")


ind <- which(beta.nbinom.power$correlation == 1.5)
pearson.beta.nbinom.low.power <- spline.fun(sample.size=sample.size, power=beta.nbinom.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1.5, method = "Pearson")
spearman.beta.nbinom.low.power <- spline.fun(sample.size=sample.size, power=beta.nbinom.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1.5, method = "Spearman")
distance.beta.nbinom.low.power <- spline.fun(sample.size=sample.size, power=beta.nbinom.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1.5, method = "Distance")
mutual.beta.nbinom.low.power <- spline.fun(sample.size=sample.size, power=beta.nbinom.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1.5, method = "Mutual")
tic.beta.nbinom.low.power <- spline.fun(sample.size=sample.size, power=beta.nbinom.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1.5, method = "TIC")


ind <- which(beta.nbinom.power$correlation == 3)
pearson.beta.nbinom.high.power <- spline.fun(sample.size=sample.size, power=beta.nbinom.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 3, method = "Pearson")
spearman.beta.nbinom.high.power <- spline.fun(sample.size=sample.size, power=beta.nbinom.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 3, method = "Spearman")
distance.beta.nbinom.high.power <- spline.fun(sample.size=sample.size, power=beta.nbinom.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 3, method = "Distance")
mutual.beta.nbinom.high.power <- spline.fun(sample.size=sample.size, power=beta.nbinom.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 3, method = "Mutual")
tic.beta.nbinom.high.power <- spline.fun(sample.size=sample.size, power=beta.nbinom.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 3, method = "TIC")

beta.nbinom.df <- rbind(pearson.beta.nbinom.type1.error,pearson.beta.nbinom.low.power,pearson.beta.nbinom.high.power,
                        spearman.beta.nbinom.type1.error,spearman.beta.nbinom.low.power,spearman.beta.nbinom.high.power,
                        distance.beta.nbinom.type1.error,distance.beta.nbinom.low.power,distance.beta.nbinom.high.power,
                        mutual.beta.nbinom.type1.error, mutual.beta.nbinom.low.power, mutual.beta.nbinom.high.power,
                        tic.beta.nbinom.type1.error, tic.beta.nbinom.low.power, tic.beta.nbinom.high.power)

beta.nbinom.df$correlation <- as.factor(beta.nbinom.df$correlation)

beta.nbinom.plot <- plot_power(df = beta.nbinom.df,
                               title = "(Beta,NegBinomial)",
                               coefs = c(0,1.5,3))
beta.nbinom.plot



#---------------- gamma to gamma ---------------------------#


sample.size <- seq(10,100,10)
p <- c(0,0.05,0.1)
reps=500

#run sim using parSim
gamma.gamma.sim <- parSim(
  
  sample.size = sample.size,
  p = p,
  reps = reps,
  nCores = 10,
  export = c("node_gamma","mi.perm","mi.test","estimate.mi"),
  
  expression = {
    
    library(simDAG)
    
    dag <- empty_dag() +
      node("x", type="rgamma", shape=4, rate=1) +
      node("y", type=node_gamma, parents="x", betas=p, intercept=3.2, phi=0.08)
    
    data <- as.data.frame(sim_from_dag(dag, sample.size))
    
    x <- data[,1]
    y <- data[,2]
    
    dc <- energy::dcor.test(x, y, R=500)
    dc.p <- dc$p.value
    
    pc <- cor.test(x, y, method="pearson")
    pc.p <- pc$p.value
    
    sc <- coin::spearman_test(y ~ x, distribution=coin::approximate(nresample=500)) 
    sc.p <- coin::pvalue(sc)[1]
    
    mi <- mi.test(x, y, type=c("c","c"), permutations=500)
    mi.p <- mi$p.value
    
    tic <- minerva::mictools(as.matrix(cbind(x,y)), seed = as.integer(runif(1, min = 0, max = 2^31-1)), nperm=500, p.adjust.method="none")
    tic.p <- tic$pval$pval
    
    data.frame(
      dc.p = dc.p,
      pc.p = pc.p,
      sc.p = sc.p,
      mi.p = mi.p,
      tic.p = tic.p)
    
  }
)



gamma.gamma.power <- as.data.frame(matrix(0, nrow=length(sample.size)*length(p), ncol = 7))
colnames(gamma.gamma.power) <- c("sample.size","correlation","pearson.power","spearman.power","distance.power","mutual.power","tic.power")

gamma.gamma.power$sample.size <- rep(sample.size, length(p))
gamma.gamma.power$correlation <- rep(p, 1, each=length(sample.size))

for(i in 1:nrow(gamma.gamma.power)) {
  
  gamma.gamma.power$pearson.power[i] <- get.power(data=gamma.gamma.sim, sample.size=gamma.gamma.power$sample.size[i],
                                                  p=gamma.gamma.power$correlation[i], type="pc.p", alpha=0.05)
  
  gamma.gamma.power$spearman.power[i] <- get.power(data=gamma.gamma.sim, sample.size=gamma.gamma.power$sample.size[i],
                                                   p=gamma.gamma.power$correlation[i], type="sc.p", alpha=0.05)
  
  gamma.gamma.power$distance.power[i] <- get.power(data=gamma.gamma.sim, sample.size=gamma.gamma.power$sample.size[i],
                                                   p=gamma.gamma.power$correlation[i], type="dc.p", alpha=0.05)
  
  gamma.gamma.power$mutual.power[i] <- get.power(data=gamma.gamma.sim, sample.size=gamma.gamma.power$sample.size[i],
                                                 p=gamma.gamma.power$correlation[i], type="mi.p", alpha=0.05)
  
  gamma.gamma.power$tic.power[i] <- get.power(data=gamma.gamma.sim, sample.size=gamma.gamma.power$sample.size[i],
                                              p=gamma.gamma.power$correlation[i], type="tic.p", alpha=0.05)
  
}


#fit power curves
ind <- which(gamma.gamma.power$correlation == 0)
pearson.gamma.gamma.type1.error <- spline.fun(sample.size=sample.size, power=gamma.gamma.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Pearson")
spearman.gamma.gamma.type1.error <- spline.fun(sample.size=sample.size, power=gamma.gamma.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Spearman")
distance.gamma.gamma.type1.error <- spline.fun(sample.size=sample.size, power=gamma.gamma.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Distance")
mutual.gamma.gamma.type1.error <- spline.fun(sample.size=sample.size, power=gamma.gamma.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Mutual")
tic.gamma.gamma.type1.error <- spline.fun(sample.size=sample.size, power=gamma.gamma.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "TIC")


ind <- which(gamma.gamma.power$correlation == 0.05)
pearson.gamma.gamma.low.power <- spline.fun(sample.size=sample.size, power=gamma.gamma.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.05, method = "Pearson")
spearman.gamma.gamma.low.power <- spline.fun(sample.size=sample.size, power=gamma.gamma.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.05, method = "Spearman")
distance.gamma.gamma.low.power <- spline.fun(sample.size=sample.size, power=gamma.gamma.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.05, method = "Distance")
mutual.gamma.gamma.low.power <- spline.fun(sample.size=sample.size, power=gamma.gamma.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.05, method = "Mutual")
tic.gamma.gamma.low.power <- spline.fun(sample.size=sample.size, power=gamma.gamma.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.05, method = "TIC")


ind <- which(gamma.gamma.power$correlation == 0.1)
pearson.gamma.gamma.high.power <- spline.fun(sample.size=sample.size, power=gamma.gamma.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.1, method = "Pearson")
spearman.gamma.gamma.high.power <- spline.fun(sample.size=sample.size, power=gamma.gamma.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.1, method = "Spearman")
distance.gamma.gamma.high.power <- spline.fun(sample.size=sample.size, power=gamma.gamma.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.1, method = "Distance")
mutual.gamma.gamma.high.power <- spline.fun(sample.size=sample.size, power=gamma.gamma.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.1, method = "Mutual")
tic.gamma.gamma.high.power <- spline.fun(sample.size=sample.size, power=gamma.gamma.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.1, method = "TIC")

gamma.gamma.df <- rbind(pearson.gamma.gamma.type1.error,pearson.gamma.gamma.low.power,pearson.gamma.gamma.high.power,
                        spearman.gamma.gamma.type1.error,spearman.gamma.gamma.low.power,spearman.gamma.gamma.high.power,
                        distance.gamma.gamma.type1.error, distance.gamma.gamma.low.power, distance.gamma.gamma.high.power,
                        mutual.gamma.gamma.type1.error, mutual.gamma.gamma.low.power, mutual.gamma.gamma.high.power,
                        tic.gamma.gamma.type1.error, tic.gamma.gamma.low.power, tic.gamma.gamma.high.power)

gamma.gamma.df$correlation <- as.factor(gamma.gamma.df$correlation)

gamma.gamma.plot <- plot_power(df = gamma.gamma.df,
                               title = "(Gamma,Gamma)",
                               coefs = c(0,0.05,0.1))
gamma.gamma.plot



#---------------- negative binomial to gamma ---------------------------#


sample.size <- seq(10,100,10)
p <- c(0,0.02,0.04)
reps=500

#run sim using parSim
nbinom.gamma.sim <- parSim(
  
  sample.size = sample.size,
  p = p,
  reps = reps,
  nCores = 10,
  export = c("node_gamma","mi.perm","mi.test","estimate.mi"),
  
  expression = {
    
    library(simDAG)
    
    dag <- empty_dag() +
      node("x", type="rnbinom", size=1, mu=5) +
      node("y", type=node_gamma, parents="x", betas=p, intercept=1, phi=0.1)
    
    data <- as.data.frame(sim_from_dag(dag, sample.size))
    
    x <- data[,1]
    y <- data[,2]
    
    dc <- energy::dcor.test(x, y, R=500)
    dc.p <- dc$p.value
    
    pc <- cor.test(x, y, method="pearson")
    pc.p <- pc$p.value
    
    sc <- coin::spearman_test(y ~ x, distribution=coin::approximate(nresample=500)) 
    sc.p <- coin::pvalue(sc)[1]
    
    mi <- mi.test(x, y, type=c("d","c"), permutations=500)
    mi.p <- mi$p.value
    
    tic <- minerva::mictools(as.matrix(cbind(x,y)), seed = as.integer(runif(1, min = 0, max = 2^31-1)), nperm=500, p.adjust.method="none")
    tic.p <- tic$pval$pval
    
    data.frame(
      dc.p = dc.p,
      pc.p = pc.p,
      sc.p = sc.p,
      mi.p = mi.p,
      tic.p = tic.p)
    
  }
)



nbinom.gamma.power <- as.data.frame(matrix(0, nrow=length(sample.size)*length(p), ncol = 6))
colnames(nbinom.gamma.power) <- c("sample.size","correlation","pearson.power","distance.power","mutual.power","tic.power")

nbinom.gamma.power$sample.size <- rep(sample.size, length(p))
nbinom.gamma.power$correlation <- rep(p, 1, each=length(sample.size))

for(i in 1:nrow(nbinom.gamma.power)) {
  
  nbinom.gamma.power$pearson.power[i] <- get.power(data=nbinom.gamma.sim, sample.size=nbinom.gamma.power$sample.size[i],
                                                   p=nbinom.gamma.power$correlation[i], type="pc.p", alpha=0.05)
  
  nbinom.gamma.power$spearman.power[i] <- get.power(data=nbinom.gamma.sim, sample.size=nbinom.gamma.power$sample.size[i],
                                                    p=nbinom.gamma.power$correlation[i], type="sc.p", alpha=0.05)
  
  nbinom.gamma.power$distance.power[i] <- get.power(data=nbinom.gamma.sim, sample.size=nbinom.gamma.power$sample.size[i],
                                                    p=nbinom.gamma.power$correlation[i], type="dc.p", alpha=0.05)
  
  nbinom.gamma.power$mutual.power[i] <- get.power(data=nbinom.gamma.sim, sample.size=nbinom.gamma.power$sample.size[i],
                                                  p=nbinom.gamma.power$correlation[i], type="mi.p", alpha=0.05)
  
  nbinom.gamma.power$tic.power[i] <- get.power(data=nbinom.gamma.sim, sample.size=nbinom.gamma.power$sample.size[i],
                                               p=nbinom.gamma.power$correlation[i], type="tic.p", alpha=0.05)
  
}


#fit power curves
ind <- which(nbinom.gamma.power$correlation == 0)
pearson.nbinom.gamma.type1.error <- spline.fun(sample.size=sample.size, power=nbinom.gamma.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Pearson")
spearman.nbinom.gamma.type1.error <- spline.fun(sample.size=sample.size, power=nbinom.gamma.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Spearman")
distance.nbinom.gamma.type1.error <- spline.fun(sample.size=sample.size, power=nbinom.gamma.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Distance")
mutual.nbinom.gamma.type1.error <- spline.fun(sample.size=sample.size, power=nbinom.gamma.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Mutual")
tic.nbinom.gamma.type1.error <- spline.fun(sample.size=sample.size, power=nbinom.gamma.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "TIC")


ind <- which(nbinom.gamma.power$correlation == 0.02)
pearson.nbinom.gamma.low.power <- spline.fun(sample.size=sample.size, power=nbinom.gamma.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.02, method = "Pearson")
spearman.nbinom.gamma.low.power <- spline.fun(sample.size=sample.size, power=nbinom.gamma.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.02, method = "Spearman")
distance.nbinom.gamma.low.power <- spline.fun(sample.size=sample.size, power=nbinom.gamma.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.02, method = "Distance")
mutual.nbinom.gamma.low.power <- spline.fun(sample.size=sample.size, power=nbinom.gamma.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.02, method = "Mutual")
tic.nbinom.gamma.low.power <- spline.fun(sample.size=sample.size, power=nbinom.gamma.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.02, method = "TIC")


ind <- which(nbinom.gamma.power$correlation == 0.04)
pearson.nbinom.gamma.high.power <- spline.fun(sample.size=sample.size, power=nbinom.gamma.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.04, method = "Pearson")
spearman.nbinom.gamma.high.power <- spline.fun(sample.size=sample.size, power=nbinom.gamma.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.04, method = "Spearman")
distance.nbinom.gamma.high.power <- spline.fun(sample.size=sample.size, power=nbinom.gamma.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.04, method = "Distance")
mutual.nbinom.gamma.high.power <- spline.fun(sample.size=sample.size, power=nbinom.gamma.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.04, method = "Mutual")
tic.nbinom.gamma.high.power <- spline.fun(sample.size=sample.size, power=nbinom.gamma.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.04, method = "TIC")

nbinom.gamma.df <- rbind(pearson.nbinom.gamma.type1.error, pearson.nbinom.gamma.low.power, pearson.nbinom.gamma.high.power,
                         spearman.nbinom.gamma.type1.error, spearman.nbinom.gamma.low.power, spearman.nbinom.gamma.high.power,
                         distance.nbinom.gamma.type1.error, distance.nbinom.gamma.low.power, distance.nbinom.gamma.high.power,
                         mutual.nbinom.gamma.type1.error, mutual.nbinom.gamma.low.power, mutual.nbinom.gamma.high.power,
                         tic.nbinom.gamma.type1.error, tic.nbinom.gamma.low.power, tic.nbinom.gamma.high.power)

nbinom.gamma.df$correlation <- as.factor(nbinom.gamma.df$correlation)

nbinom.gamma.plot <- plot_power(df = nbinom.gamma.df,
                                title = "(NegBinomial,Gamma",
                                coefs = c(0,0.02,0.04))
nbinom.gamma.plot



#---------------- beta to gamma ---------------------------#


sample.size <- seq(10,100,10)
p <- c(0,0.6,1.5)
reps=500

#run sim using parSim
beta.gamma.sim <- parSim(
  
  sample.size = sample.size,
  p = p,
  reps = reps,
  nCores = 10,
  export = c("node_gamma","mi.perm","mi.test","estimate.mi"),
  
  expression = {
    
    library(simDAG)
    
    dag <- empty_dag() +
      node("x", type="rbeta", shape1=2, shape2=6) +
      node("y", type=node_gamma, parents="x", betas=p, intercept=2, phi=0.1)
    
    data <- as.data.frame(sim_from_dag(dag, sample.size))
    
    x <- data[,1]
    y <- data[,2]
    
    dc <- energy::dcor.test(x, y, R=500)
    dc.p <- dc$p.value
    
    pc <- cor.test(x, y, method="pearson")
    pc.p <- pc$p.value
    
    sc <- coin::spearman_test(y ~ x, distribution=coin::approximate(nresample=500)) 
    sc.p <- coin::pvalue(sc)[1]
    
    mi <- mi.test(x, y, type=c("c","c"), permutations=500)
    mi.p <- mi$p.value
    
    tic <- minerva::mictools(as.matrix(cbind(x,y)), seed = as.integer(runif(1, min = 0, max = 2^31-1)), nperm=500, p.adjust.method="none")
    tic.p <- tic$pval$pval
    
    data.frame(
      dc.p = dc.p,
      pc.p = pc.p,
      sc.p = sc.p,
      mi.p = mi.p,
      tic.p = tic.p)
    
  }
)



beta.gamma.power <- as.data.frame(matrix(0, nrow=length(sample.size)*length(p), ncol = 7))
colnames(beta.gamma.power) <- c("sample.size","correlation","pearson.power","spearman.power","distance.power","mutual.info.power","TIC.power")

beta.gamma.power$sample.size <- rep(sample.size, length(p))
beta.gamma.power$correlation <- rep(p, 1, each=length(sample.size))

for(i in 1:nrow(beta.gamma.power)) {
  
  beta.gamma.power$pearson.power[i] <- get.power(data=beta.gamma.sim, sample.size=beta.gamma.power$sample.size[i],
                                                 p=beta.gamma.power$correlation[i], type="pc.p", alpha=0.05)
  
  beta.gamma.power$spearman.power[i] <- get.power(data=beta.gamma.sim, sample.size=beta.gamma.power$sample.size[i],
                                                  p=beta.gamma.power$correlation[i], type="sc.p", alpha=0.05)
  
  beta.gamma.power$distance.power[i] <- get.power(data=beta.gamma.sim, sample.size=beta.gamma.power$sample.size[i],
                                                  p=beta.gamma.power$correlation[i], type="dc.p", alpha=0.05)
  
  beta.gamma.power$mutual.power[i] <- get.power(data=beta.gamma.sim, sample.size=beta.gamma.power$sample.size[i],
                                                p=beta.gamma.power$correlation[i], type="mi.p", alpha=0.05)
  
  beta.gamma.power$tic.power[i] <- get.power(data=beta.gamma.sim, sample.size=beta.gamma.power$sample.size[i],
                                             p=beta.gamma.power$correlation[i], type="tic.p", alpha=0.05)
  
}


#fit power curves
ind <- which(beta.gamma.power$correlation == 0)
pearson.beta.gamma.type1.error <- spline.fun(sample.size=sample.size, power=beta.gamma.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Pearson")
spearman.beta.gamma.type1.error <- spline.fun(sample.size=sample.size, power=beta.gamma.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Spearman")
distance.beta.gamma.type1.error <- spline.fun(sample.size=sample.size, power=beta.gamma.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Distance")
mutual.beta.gamma.type1.error <- spline.fun(sample.size=sample.size, power=beta.gamma.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Mutual")
tic.beta.gamma.type1.error <- spline.fun(sample.size=sample.size, power=beta.gamma.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "TIC")


ind <- which(beta.gamma.power$correlation == 0.6)
pearson.beta.gamma.low.power <- spline.fun(sample.size=sample.size, power=beta.gamma.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.6, method = "Pearson")
spearman.beta.gamma.low.power <- spline.fun(sample.size=sample.size, power=beta.gamma.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.6, method = "Spearman")
distance.beta.gamma.low.power <- spline.fun(sample.size=sample.size, power=beta.gamma.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.6, method = "Distance")
mutual.beta.gamma.low.power <- spline.fun(sample.size=sample.size, power=beta.gamma.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.6, method = "Mutual")
tic.beta.gamma.low.power <- spline.fun(sample.size=sample.size, power=beta.gamma.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.6, method = "TIC")


ind <- which(beta.gamma.power$correlation == 1.5)
pearson.beta.gamma.high.power <- spline.fun(sample.size=sample.size, power=beta.gamma.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1.5, method = "Pearson")
spearman.beta.gamma.high.power <- spline.fun(sample.size=sample.size, power=beta.gamma.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1.5, method = "Spearman")
distance.beta.gamma.high.power <- spline.fun(sample.size=sample.size, power=beta.gamma.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1.5, method = "Distance")
mutual.beta.gamma.high.power <- spline.fun(sample.size=sample.size, power=beta.gamma.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1.5, method = "Mutual")
tic.beta.gamma.high.power <- spline.fun(sample.size=sample.size, power=beta.gamma.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1.5, method = "TIC")

beta.gamma.df <- rbind(pearson.beta.gamma.type1.error, pearson.beta.gamma.low.power, pearson.beta.gamma.high.power,
                       spearman.beta.gamma.type1.error, spearman.beta.gamma.low.power, spearman.beta.gamma.high.power,
                       distance.beta.gamma.type1.error,distance.beta.gamma.low.power, distance.beta.gamma.high.power,
                       mutual.beta.gamma.type1.error, mutual.beta.gamma.low.power, mutual.beta.gamma.high.power,
                       tic.beta.gamma.type1.error, tic.beta.gamma.low.power, tic.beta.gamma.high.power)

beta.gamma.df$correlation <- as.factor(beta.gamma.df$correlation)

beta.gamma.plot <- plot_power(df = beta.gamma.df,
                              title = "(Beta,Gamma)",
                              coefs = c(0,0.6,1.5))
beta.gamma.plot


#---------------- binary to gamma ---------------------------#

sample.size <- seq(10,100,10)
p <- c(0,0.5,1)
reps=500

#run sim using parSim
binary.gamma.sim <- parSim(
  
  sample.size = sample.size,
  p = p,
  reps = reps,
  nCores = 10,
  export = c("node_gamma","mi.perm","mi.test","estimate.mi"),
  
  expression = {
    
    library(simDAG)
    
    dag <- empty_dag() +
      node("x", type="rbinom", size=1, prob=0.5) +
      node("y", type=node_gamma, parents="x", betas=p, intercept=2, phi=0.6)
    
    data <- as.data.frame(sim_from_dag(dag, sample.size))
    
    x <- data[,1]
    y <- data[,2]
    
    dc <- energy::dcor.test(x, y, R=500)
    dc.p <- dc$p.value
    
    pc <- cor.test(x, y, method="pearson")
    pc.p <- pc$p.value
    
    sc <- coin::spearman_test(y ~ x, distribution=coin::approximate(nresample=500)) 
    sc.p <- coin::pvalue(sc)[1]
    
    mi <- mi.test(x, y, type=c("d","c"), permutations=500)
    mi.p <- mi$p.value
    
    tic <- minerva::mictools(as.matrix(cbind(x,y)), seed = as.integer(runif(1, min = 0, max = 2^31-1)), nperm=500, p.adjust.method="none")
    tic.p <- tic$pval$pval
    
    data.frame(
      dc.p = dc.p,
      pc.p = pc.p,
      sc.p = sc.p,
      mi.p = mi.p,
      tic.p = tic.p)
    
  }
)



binary.gamma.power <- as.data.frame(matrix(0, nrow=length(sample.size)*length(p), ncol = 7))
colnames(binary.gamma.power) <- c("sample.size","correlation","pearson.power","spearman.power","distance.power","mutual.power","tic.power")

binary.gamma.power$sample.size <- rep(sample.size, length(p))
binary.gamma.power$correlation <- rep(p, 1, each=length(sample.size))

for(i in 1:nrow(binary.gamma.power)) {
  
  binary.gamma.power$pearson.power[i] <- get.power(data=binary.gamma.sim, sample.size=binary.gamma.power$sample.size[i],
                                                   p=binary.gamma.power$correlation[i], type="pc.p", alpha=0.05)
  
  binary.gamma.power$spearman.power[i] <- get.power(data=binary.gamma.sim, sample.size=binary.gamma.power$sample.size[i],
                                                    p=binary.gamma.power$correlation[i], type = "sc.p", alpha=0.05)
  
  binary.gamma.power$distance.power[i] <- get.power(data=binary.gamma.sim, sample.size=binary.gamma.power$sample.size[i],
                                                    p=binary.gamma.power$correlation[i], type = "dc.p", alpha=0.05)
  
  binary.gamma.power$mutual.power[i] <- get.power(data=binary.gamma.sim, sample.size=binary.gamma.power$sample.size[i],
                                                  p=binary.gamma.power$correlation[i], type = "mi.p", alpha=0.05)
  
  binary.gamma.power$tic.power[i] <- get.power(data=binary.gamma.sim, sample.size=binary.gamma.power$sample.size[i],
                                               p=binary.gamma.power$correlation[i], type = "tic.p", alpha=0.05)
  
}


#fit power curves
ind <- which(binary.gamma.power$correlation == 0)
pearson.binary.gamma.type1.error <- spline.fun(sample.size=sample.size, power=binary.gamma.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Pearson")
spearman.binary.gamma.type1.error <- spline.fun(sample.size=sample.size, power=binary.gamma.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Spearman")
distance.binary.gamma.type1.error <- spline.fun(sample.size=sample.size, power=binary.gamma.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Distance")
mutual.binary.gamma.type1.error <- spline.fun(sample.size=sample.size, power=binary.gamma.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "Mutual")
tic.binary.gamma.type1.error <- spline.fun(sample.size=sample.size, power=binary.gamma.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0, method = "TIC")


ind <- which(binary.gamma.power$correlation == 0.5)
pearson.binary.gamma.low.power <- spline.fun(sample.size=sample.size, power=binary.gamma.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.5, method = "Pearson")
spearman.binary.gamma.low.power <- spline.fun(sample.size=sample.size, power=binary.gamma.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.5, method = "Spearman")
distance.binary.gamma.low.power <- spline.fun(sample.size=sample.size, power=binary.gamma.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.5, method = "Distance")
mutual.binary.gamma.low.power <- spline.fun(sample.size=sample.size, power=binary.gamma.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.5, method = "Mutual")
tic.binary.gamma.low.power <- spline.fun(sample.size=sample.size, power=binary.gamma.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 0.5, method = "TIC")


ind <- which(binary.gamma.power$correlation == 1)
pearson.binary.gamma.high.power <- spline.fun(sample.size=sample.size, power=binary.gamma.power$pearson.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1, method = "Pearson")
spearman.binary.gamma.high.power <- spline.fun(sample.size=sample.size, power=binary.gamma.power$spearman.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1, method = "Spearman")
distance.binary.gamma.high.power <- spline.fun(sample.size=sample.size, power=binary.gamma.power$distance.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1, method = "Distance")
mutual.binary.gamma.high.power <- spline.fun(sample.size=sample.size, power=binary.gamma.power$mutual.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1, method = "Mutual")
tic.binary.gamma.high.power <- spline.fun(sample.size=sample.size, power=binary.gamma.power$tic.power[ind], minSampleSize = 10, maxSampleSize = 100, correlation = 1, method = "TIC")

binary.gamma.df <- rbind(pearson.binary.gamma.type1.error, pearson.binary.gamma.low.power, pearson.binary.gamma.high.power,
                         spearman.binary.gamma.type1.error, spearman.binary.gamma.low.power, spearman.binary.gamma.high.power,
                         distance.binary.gamma.type1.error, distance.binary.gamma.low.power, distance.binary.gamma.high.power,
                         mutual.binary.gamma.type1.error, mutual.binary.gamma.low.power, mutual.binary.gamma.high.power,
                         tic.binary.gamma.type1.error, tic.binary.gamma.low.power, tic.binary.gamma.high.power)

binary.gamma.df$correlation <- as.factor(binary.gamma.df$correlation)

binary.gamma.plot <- plot_power(df = binary.gamma.df,
                                title = "(Binary,Gamma)",
                                coefs = c(0,0.5,1))
binary.gamma.plot



#------------------------------------#
#--------- SAVE PLOTS ---------------#
#------------------------------------#

setwd("/Users/tomrowland/Documents/PhD/Thesis write up/Chapter 3")

pdf("Figure 3.1.pdf", width=32, height=18)

ggarrange(beta.beta.plot, binary.beta.plot, gamma.beta.plot, nbinom.beta.plot,
          nbinom.nbinom.plot, binary.nbinom.plot, gamma.nbinom.plot, beta.nbinom.plot,
          gamma.gamma.plot, binary.gamma.plot, nbinom.gamma.plot, beta.nbinom.plot,
          nrow=3, ncol=4)

dev.off()


###########
#
#   END
#
##########
