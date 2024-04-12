# What do we want to do here?
# See if singular values in orthogonal design have selection bias
# See if categorical variables with more categories are preferred over the ones with less in mboost

library(mboost)
library(tidyverse)
library(pracma)
library(gglasso)
library(mvtnorm)

run_experiment <- function(iteration, intercept = TRUE, scaled = TRUE, type = 'ortho'){
  if(type == 'ortho'){
  dat<- pracma::randortho(100, type = 'orthonormal')[,1:5,drop=F] %>%
    as.data.frame() 
  } else if(type == 'independent'){
    dat <- matrix( rnorm(100*5,mean=0,sd=1), 100, 5) %>% as.data.frame() 
  } else if(type == 'singular'){
    dat <- pracma::randortho(100, type = 'orthonormal')[,1:5,drop=F]%*%diag(c(1,2,1,5,1))%*%
      t(pracma::randortho(5, type = 'orthonormal')[,1:5,drop=F]) %>% as.data.frame()
  } else if(type == 'cor'){
    dat <- mvtnorm::rmvnorm(n=100, mean = rep(0,5), sigma = matrix(c(1,  0.9,0,0,0.5, 
                                                              0.9,1,  0,0,0.5,
                                                              0,  0,  1,0,0,
                                                              0,  0,  0,1,0,
                                                              0.5,0.5,0,0,1
                                                              ),nrow = 5)) %>%
      as.data.frame()
  }
  dat <- dat %>%
    mutate(V2 = V2*2, V4 = V4*5) %>%
    mutate(y = rnorm(n=100))
  if(scaled == TRUE){
    dat <- dat %>%
      mutate_all(function(x){as.numeric(scale(x))})
  }
  mb_model_unreg <- mboost(y ~bols(V1, V2)+bols(V3, V4) + bols(V5),
                     data = dat, control = boost_control(mstop = 1, nu = 1))
  mb_model_0.5 <- mboost(y ~bols(V1, V2, df = 0.5)+
                           bols(V3, V4, df = 0.5) +
                           bols(V5, df = 0.5),
                         data = dat, control = boost_control(mstop = 1, nu = 1))
  mb_model_0.1 <- mboost(y ~bols(V1, V2, df = 0.1)+
                           bols(V3, V4, df = 0.1) +
                           bols(V5, df = 0.1),
                         data = dat, control = boost_control(mstop = 1, nu = 1))
  mb_model_1 <- mboost(y ~bols(V1, V2, df = 1)+
                           bols(V3, V4, df = 1) +
                           bols(V5, df = 1),
                         data = dat, control = boost_control(mstop = 1, nu = 1))
  mb_model_lam <- mboost(y ~bols(V1, V2, lambda = 10)+
                           bols(V3, V4, lambda = 10) +
                           bols(V5, lambda = 10),
                         data = dat, control = boost_control(mstop = 1, nu = 1))
  gr_lasso <- gglasso(x = as.matrix(dat %>% select(-y)), y = dat$y,
                       group = c(1,1,2,2,3))
  gr_lasso_coef <- gr_lasso$beta %>% 
    as.data.frame() %>% 
    rownames_to_column() %>%
    gather(key = 'lambda', value = 'coef', -rowname) %>%
    filter(coef != 0) %>%
    slice(1:5) 
  gr_lasso_coef <- gr_lasso_coef %>%
   filter(lambda == gr_lasso_coef$lambda[1])
  if(dim(gr_lasso_coef)[1] > 2){gr_sel <- NA} else{
    gr_sel <- case_when(any(gr_lasso_coef$rowname %in% c('V1','V2'))~ 1, 
                        any(gr_lasso_coef$rowname %in% c('V3'))~3, T ~ 2)
  }
  ret <- data.frame(model = c('mb', 'mb_0.5', 'mb_1', 'mb_0.1', 'mb_lambda', 'grla'), 
                    sel = c(mb_model_unreg$xselect(), mb_model_0.5$xselect(),
                            mb_model_1$xselect(), mb_model_0.1$xselect(),
                            mb_model_lam$xselect(), gr_sel))
  return(ret)
}

library(furrr)
plan(multisession)
iterations <- 100000
# Define the 12 data scenarios
params <- tidyr::expand_grid(iteration = 1:iterations, scaled = c(TRUE,FALSE),
                             type = c('ortho', 'independent', 'singular', 'cor'))
# run simulation in parallel
set.seed(5)
tictoc::tic()
results_df <- params %>%
  mutate( 
    res = future_pmap(., .f = run_experiment, .progress = T, .options=furrr_options(seed = TRUE))
  ) %>%
  unnest(cols = res)
tictoc::toc()
saveRDS(file = 'results/small_group_1_raw.RDS', results_df)
results_df %>%
  group_by(scaled, model, type, sel) %>%
  tally() %>%
  filter(!is.na(sel)) %>%
  mutate(prop = round(n/sum(n,na.rm = T), 3)) %>%
saveRDS(file = 'results/small_group_1.RDS')
plan(sequential)

################################################################################
################################################################################
# group of one vs group of 15

run_experiment <- function(iteration, intercept = TRUE, scaled = TRUE, type = 'ortho'){
  if(type == 'ortho'){
    dat<- pracma::randortho(100, type = 'orthonormal')[,1:16,drop=F] %>%
      as.data.frame() 
  } else if(type == 'independent'){
    dat <- matrix( rnorm(100*16,mean=0,sd=1), 100, 16) %>% as.data.frame() 
  } else if(type == 'singular'){
    dat <- pracma::randortho(100, type = 'orthonormal')[,1:16,drop=F]%*%diag(c(1,rep(15,5), rep(1,5), rep(2,5)))%*%
      t(pracma::randortho(16, type = 'orthonormal')[,1:16,drop=F]) %>% as.data.frame()
  } else if(type == 'cor'){
    dat <- mvtnorm::rmvnorm(n=100, mean = 1:16, sigma = ones(16,16)-0.2 + diag(c(rep(0.2,16)))) %>%
      as.data.frame()
  }
  dat <- dat %>%
    mutate(y = rnorm(n=100))
  if(scaled == TRUE){
    dat <- dat %>%
      mutate_all(function(x){as.numeric(scale(x))})
  }
  mb_model_unreg <- mboost(y ~bols(V1)+bols(V2, V3, V4,V5, V6, V7,V8, V9, V10,V11, V12, V13,V14,V15),
                           data = dat, control = boost_control(mstop = 1, nu = 1))
  mb_model_0.5 <- mboost(y ~bols(V1,df =0.5)+bols(V2, V3, V4,V5, V6, V7,V8, V9, V10,V11, V12, V13,V14,V15, df=0.5),
                         data = dat, control = boost_control(mstop = 1, nu = 1))
  mb_model_0.1 <- mboost(y ~bols(V1,df =0.1)+bols(V2, V3, V4,V5, V6, V7,V8, V9, V10,V11, V12, V13,V14,V15, df=0.1),
                         data = dat, control = boost_control(mstop = 1, nu = 1))
  mb_model_1 <- mboost(y ~bols(V1,df =1)+bols(V2, V3, V4,V5, V6, V7,V8, V9, V10,V11, V12, V13,V14,V15, df=1),
                       data = dat, control = boost_control(mstop = 1, nu = 1))
  mb_model_lam <- mboost(y ~bols(V1,lambda =10)+bols(V2, V3, V4,V5, V6, V7,V8, V9, V10,V11, V12, V13,V14,V15, lambda=10),
                         data = dat, control = boost_control(mstop = 1, nu = 1))
  gr_lasso <- gglasso(x = as.matrix(dat %>% select(-y)), y = dat$y,
                      group = c(1,rep(2,15)))
  gr_lasso_coef <- gr_lasso$beta %>% 
    as.data.frame() %>% 
    rownames_to_column() %>%
    gather(key = 'lambda', value = 'coef', -rowname) %>%
    filter(coef != 0) %>%
    slice(1:16) 
  gr_lasso_coef <- gr_lasso_coef %>%
    filter(lambda == gr_lasso_coef$lambda[1])
  if(dim(gr_lasso_coef)[1] > 16){gr_sel <- NA} else{
    gr_sel <- case_when(any(gr_lasso_coef$rowname %in% c('V1'))~ 1, T ~ 2)
  }
  ret <- data.frame(model = c('mb', 'mb_0.5', 'mb_1', 'mb_0.1', 'mb_lambda', 'grla'), 
                    sel = c(mb_model_unreg$xselect(), mb_model_0.5$xselect(),
                            mb_model_1$xselect(), mb_model_0.1$xselect(),
                            mb_model_lam$xselect(), gr_sel))
  return(ret)
}
# group size 1 vs 15
library(furrr)
plan(multisession)
iterations <- 100000
# Define the 12 data scenarios
params <- tidyr::expand_grid(iteration = 1:iterations, scaled = c(TRUE,FALSE),
                             type = c('ortho', 'independent', 'singular', 'cor'))
# run simulation in parallel
set.seed(5)
tictoc::tic()
results_df <- params %>%
  mutate( 
    res = future_pmap(., .f = run_experiment, .progress = T, .options=furrr_options(seed = TRUE))
  ) %>%
  unnest(cols = res)
tictoc::toc()
saveRDS(file = 'results/big_group_1_raw.RDS', results_df)
results_df %>%
  group_by(scaled, model, type, sel) %>%
  tally() %>%
  filter(!is.na(sel)) %>%
  mutate(prop = round(n/sum(n,na.rm = T), 3)) %>%
  saveRDS(file = 'results/big_group_1.RDS')
plan(sequential)


# Categorical data 
x <- sample( LETTERS[1:4], 10000, replace=TRUE, prob=c(0.1, 0.2, 0.65, 0.05) )
