#### Adv Meas Spring 2020 Level Test Simulation ####
#### Author: Joseph DeWeese ####
#### 2020-05-05 ####

# load packages and functions
  library(tidyverse)
  library(irtplay)
  library(irtoys)
  source("lab.irf.r")
  source("gen.resp.r")
  source("level_scatter_shape_test_jnd_20200312.R")

# sim people
  set.seed(456)
  
  # parameters
    n_people <- 10000
    n_dims <- 3
    mu <- rep(0, n_dims)
    sigma <- array(0, c(n_dims, n_dims))
    diag(sigma) <- 1
    change <- c(0, 0, .25, .5, .75, 1)
  
  # initial thetas
    tbase <- MASS::mvrnorm(n_people, mu, sigma)
  
  # add change amounts
    thetas <- vector("list", length = length(change)+2)
    names(thetas) <- c("t_base", "t_0", "t_25", "t_50", "t_75", "t_1", "t_n75", "t_n150")
    
    for(i in seq_along(change)){
      thetas[[i]] <- tbase + change[i]
    }
    
    cn1 <- matrix(c(rep(.75, n_people), rep(.75, n_people), rep(0, n_people)), ncol = 3)
    cn2 <- matrix(c(rep(1.5, n_people), rep(0, n_people), rep(0, n_people)), ncol = 3)
    
    thetas$t_n75 <- tbase + cn1
    thetas$t_n150 <- tbase + cn2

# sim items
  n_items <- 40
  
  # a
    a_mu <- rep(1, n_dims)
    a_sigma <- array(0, c(n_dims, n_dims))
    diag(a_sigma) <- .15^2
    
    as <- MASS::mvrnorm(n_items, a_mu, a_sigma)
  
  # b
    b_mu <- rep(0, n_dims)
    b_sigma <- array(0, c(n_dims, n_dims))
    diag(b_sigma) <- 1
    
    bs <- MASS::mvrnorm(n_items, b_mu, b_sigma)
  
  # c
    c_min <- .0
    c_max <- .2
    cs <- cbind(runif(n_items, c_min, c_max),
                runif(n_items, c_min, c_max),
                runif(n_items, c_min, c_max))
  
  # combine
    ips <- vector("list", length = n_dims)
    names(ips) <- paste0("dim", 1:3)
    for(i in seq_len(n_dims)) {
      ips[[i]] <- cbind(as[,i], bs[,i], cs[,i])
      colnames(ips[[i]]) <- paste0(c("a", "b", "c"), i)
    }

# generate responses
  d.thetas <- map(thetas, as.data.frame) %>%
      map(~rename_all(., ~str_replace(., "V", "dim")))
  
  
  resps <- map(d.thetas,
               ~map2(., ips,
                     ~gen.resp(.y, .x, D= 1.7)
               )
  )

# estimate thetas with irtoys
  # irtoys uses D = 1, so need to multiply a by 1.7 
  ips.toys <- map(ips, ~ . %*% diag(c(1.7,1,1)))
  
  ests <- map(resps, 
               ~map2(., ips.toys,
                     ~mlebme(.x, .y, method = "ML")))
  
  # combine into dfs within condition (so all dims in one df)
    ests.d <- map(ests, as.data.frame)
    
# calculate level stats
    
  # function, specific to this sim's data organiziation
    comp_lvl <- function(data.list, base, comp, est.names, sem.names) {
      avg1 <- apply(data.list[[base]][, est.names], 1, mean)
      avg2 <- apply(data.list[[comp]][, est.names], 1, mean)
      var1 <- apply(data.list[[base]][, sem.names], 1, function(x) sum(x^2)/n_dims^2)
      var2 <- apply(data.list[[comp]][, sem.names], 1, function(x) sum(x^2)/n_dims^2)
      Z <- (avg2 - avg1) / sqrt(var2+var1)
      pvals <- 2 * (1 - pnorm(abs(Z)))
      return(list(avg = data.frame(avg1 = avg1, avg2 = avg2), 
                  var = data.frame(var1 = var1, var2 = var2),
                  Z = Z, p_val = pvals))
    }
    
  # use fn to calc level
    est.names <- c("dim1.est", "dim2.est", "dim3.est")
    sem.names <- c("dim1.sem", "dim2.sem", "dim3.sem")
    cond.names <- c(t_0 = "t_0", t_25 = "t_25", t_50 = "t_50", 
                    t_75 = "t_75", t_1 = "t_1", t_n75 = "t_n75", t_n150 = "t_n150")
    
    lvls <- map(cond.names, ~comp_lvl(ests.d, "t_base", .,
                                      est.names, sem.names))
    
    
# corr .5 ####
  sigma2 <- array(.5, c(n_dims, n_dims))
  diag(sigma2) <- 1   
  
  # initial thetas
    tbase2 <- MASS::mvrnorm(n_people, mu, sigma2)
  
  # add change amounts
    thetas2 <- vector("list", length = length(change)+2)
    names(thetas2) <- c("t_base", "t_0", "t_25", "t_50", "t_75", "t_1", "t_n75", "t_n150")
    
    for(i in seq_along(change)){
      thetas2[[i]] <- tbase2 + change[i]
    }
    
    cn1 <- matrix(c(rep(.75, n_people), rep(.75, n_people), rep(0, n_people)), ncol = 3)
    cn2 <- matrix(c(rep(1.5, n_people), rep(0, n_people), rep(0, n_people)), ncol = 3)
    
    thetas2$t_n75 <- tbase2 + cn1
    thetas2$t_n150 <- tbase2 + cn2
    
  # generate responses
    d.thetas2 <- map(thetas2, as.data.frame) %>%
      map(~rename_all(., ~str_replace(., "V", "dim")))
    
    
    resps2 <- map(d.thetas2,
                 ~map2(., ips,
                       ~gen.resp(.y, .x, D= 1.7)
                 )
    )
  
  # estimate thetas with irtoys
  # irtoys uses D = 1, so need to multiply a by 1.7 
  ests2 <- map(resps2, 
              ~map2(., ips.toys,
                    ~mlebme(.x, .y, method = "ML")))
  
  # combine into dfs within condition (so all dims in one df)
  ests2.d <- map(ests2, as.data.frame)
  
  # calculate level stats
  
  # use fn to calc level
  lvls2 <- map(cond.names, ~comp_lvl(ests2.d, "t_base", .,
                                    est.names, sem.names))  
    
    
# create list of everything and output it ####
  out10k2 <- list(thetas = d.thetas,
              resps = resps,
              ips = ips,
              ips.toys = ips.toys,
              ests = ests.d,
              lvls = lvls,
              thetas2 = d.thetas2,
              resps2 = resps2,
              ests2 = ests2.d,
              lvls2 = lvls2)

  save(out10k2, file = "level_sim_results_10k2.Rdata")
  
#### end ####