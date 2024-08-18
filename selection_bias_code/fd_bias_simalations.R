library(mboost)
library(tidyverse)
library(pracma)
library(gglasso)
library(mvtnorm)
library(furrr)

get_parameters <- function(df, group_df, block_fit, block_size = 1000, max_steps = 50, tol = 0.001, rate = 0.1, initialization = F){
  # Initialization
  scaling <- as.data.frame(t(rep(1,length(unique(group_df$var_name)))))
  colnames(scaling) <- group_df$var_name
  standard_df <- df %>% mutate_all(function(x){as.numeric(scale(x))})
  big_results_df <- data.frame()
  opt_perc <- 1/length(unique(group_df$group_name))
  best_scale <- rep(1,dim(group_df)[1])
  if(initialization == T){
    scaling <- as.data.frame(t(c(1,rep(0.1,length(unique(group_df$var_name))-1))))
    best_scale <- c(1,rep(0.1,length(unique(group_df$var_name))-1))
  }
  # steps
  for(i in 1:max_steps){
    param_df <- expand_grid(iter = 1:block_size)
    model_df <- as.data.frame(as.matrix(standard_df)%*%diag(tail(scaling,1)))
    colnames(model_df) <- colnames(standard_df)
    results_df <- param_df %>%
      mutate( 
        sel = future_pmap(., .f = ~block_fit(.iter, model_df), .progress = T, .options=furrr_options(seed = TRUE))
      ) %>%
      unnest(cols = sel) %>% 
      mutate(step = i)
    results_df <-
      results_df %>%
      group_by(step, sel) %>%
      tally() %>%
      right_join(
        group_df %>%
          select(group_name) %>%
          group_by(group_name) %>% slice(1) %>%
          mutate(step = i),
        by = c('step', c('sel' = 'group_name'))
      ) %>%
      arrange(sel) %>%
      group_by(step) %>%
      mutate(n = case_when(is.na(n) ~ 1, T ~ n),
             perc = n / sum(n)) %>%
      right_join(group_df, by = c('sel' = 'group_name')) %>%
      mutate(
        current_scale = as.vector(t(tail(scaling,1))),
        best_scale = best_scale,
        scale_factor = opt_perc / (perc),
        new_scale = scale_factor * current_scale,
        new_scale = best_scale * (1 - rate) + new_scale * rate,
        loss = sum((opt_perc - perc) ^ 2)
      ) %>%
      arrange(sel)
   scaling <- scaling %>%
     rbind(results_df$new_scale)
   big_results_df <- big_results_df %>%
     bind_rows(results_df)

   best_step <- filter(big_results_df,loss == min(loss))$step[1]
   best_scale <- filter(big_results_df, loss == min(loss))$current_scale[1:dim(group_df)[1]]
   if(best_step == 1){best_step <- 2}
   if(all(abs((results_df$perc - opt_perc)/opt_perc) < tol)){
     print(paste0('stopped after', i, ' iterations'))
     break
   }
  }
  return(list('results' = big_results_df, 'scaling' = tail(scaling,1), 'best' = best_scale))
}
## Helper function Tuning grlasso
fit_lasso_01 <- function(iter, df){
  
  sim_df <- df %>% dplyr::mutate(y = rnorm(dim(df)[1]))
 
  gr_lasso <- gglasso(x = as.matrix(sim_df %>% select(-y)), y = sim_df$y,
                      group = c(1,1,2,2,3), lambda = 0.1)
  gr_lasso_coef <- gr_lasso$beta %>% 
    as.data.frame() %>% 
    rownames_to_column() %>%
    gather(key = 'lambda', value = 'coef', -rowname) %>%
    filter(coef != 0) %>%
    slice(1:5) 
  gr_lasso_coef <- gr_lasso_coef %>%
    filter(lambda == gr_lasso_coef$lambda[1])
  gr_lasso_coef
  gr_lasso_coef <- gr_lasso_coef %>% 
    mutate(sel = case_when(rowname %in% c('V1','V2')~ 1, 
                           rowname %in% c('V3','V4')~2, 
                           rowname %in% c('V5')~ 3))
  unique(gr_lasso_coef$sel)
}
fit_lasso_005 <- function(iter, df){
  
  sim_df <- df %>% dplyr::mutate(y = rnorm(dim(df)[1]))
  
  gr_lasso <- gglasso(x = as.matrix(sim_df %>% select(-y)), y = sim_df$y,
                      group = c(1,1,2,2,3), lambda = 0.05)
  gr_lasso_coef <- gr_lasso$beta %>% 
    as.data.frame() %>% 
    rownames_to_column() %>%
    gather(key = 'lambda', value = 'coef', -rowname) %>%
    filter(coef != 0) %>%
    slice(1:5) 
  gr_lasso_coef <- gr_lasso_coef %>%
    filter(lambda == gr_lasso_coef$lambda[1])
  gr_lasso_coef
  gr_lasso_coef <- gr_lasso_coef %>% 
    mutate(sel = case_when(rowname %in% c('V1','V2')~ 1, 
                           rowname %in% c('V3','V4')~2, 
                           rowname %in% c('V5')~ 3))
  unique(gr_lasso_coef$sel)
}

group_df <- data.frame(var_name = c('V1', 'V2', 'V3', 'V4', 'V5'),
                       group_name = as.integer(c(1,1,2,2,3)))

fit_mb_lam_10 <- function(iter, df){
  formula <- as.formula(y ~bols(V1, V2, lambda = 10)+
                          bols(V3, V4, lambda = 10) +
                          bols(V5, lambda = 10))
  sim_df <- df %>% dplyr::mutate(y = rnorm(dim(df)[1]))
  mb_model <- mboost(formula = formula, data = sim_df, control = boost_control(mstop = 1))
  mb_model$xselect()
}
fit_mb_df <- function(iter, df){
  formula <- as.formula(y ~bols(V1, V2, df = 0.5)+
                          bols(V3, V4,df = 0.5) +
                          bols(V5,df = 0.5))
  sim_df <- df %>% dplyr::mutate(y = rnorm(dim(df)[1]))
  mb_model <- mboost(formula = formula, data = sim_df, control = boost_control(mstop = 1))
  mb_model$xselect()
}


### Simulations to show bias in methods
set.seed(1)
dat_ortho <- pracma::randortho(100, type = 'orthonormal')[,1:5,drop=F] %>%
  as.data.frame() 
dat_singular <- pracma::randortho(100, type = 'orthonormal')[,1:5,drop=F]%*%diag(c(1,2,1,5,1))%*%
  t(pracma::randortho(5, type = 'orthonormal')[,1:5,drop=F]) %>% as.data.frame()
dat_independent <- matrix( rnorm(100*5,mean=0,sd=1), 100, 5) %>% as.data.frame() 
dat_cor <- mvtnorm::rmvnorm(n=100, mean = rep(0,5), sigma = matrix(c(1,  0.9,0,0,0.5, 
                                                                     0.9,1,  0,0,0.5,
                                                                     0,  0,  1,0,0,
                                                                     0,  0,  0,1,0,
                                                                     0.5,0.5,0,0,1),nrow = 5)) %>%
  as.data.frame()
dat_list <- list('ortho' = dat_ortho,'singular' = dat_singular,
                 'ind' = dat_independent, 'col' = dat_cor)

params_list_1 <- dat_list %>% 
  lapply(
    function(dat){
      plan(multisession)
      mb_params_df <- get_parameters(df = dat, group_df = group_df, block_size = 10000, block_fit = fit_mb_df, rate = 0.3)
      plan(sequential)
      plan(multisession)
      lasso_params_01 <- get_parameters(df = dat, group_df = group_df, block_size = 10000, block_fit = fit_lasso_01, rate = 0.3)
      plan(sequential)
      print('0.5')
      plan(multisession)
      lasso_params_005 <- get_parameters(df = dat, group_df = group_df, block_size = 10000, block_fit = fit_lasso_005, rate = 0.3)
      plan(sequential)
      plan(multisession)
      mb_params_lam_10 <- get_parameters(df = dat, group_df = group_df, block_size = 10000, block_fit = fit_mb_lam_10, rate = 0.8)
      plan(sequential)
      print('1')
      return(list('mb_df' = mb_params_df, 'lasso_01' = lasso_params_01, 'lasso_005' = lasso_params_005,
                  'mb_lam_10' = mb_params_lam_10))
    }
  )
#saveRDS(params_list_1, 'results_fixed/params_list_1.RDS')
params_list_1 <- readRDS('results_fixed/params_list_1.RDS')

run_experiment <- function(iteration){
  ret <- data.frame()
  for(beta_type in c('zero', 'high', 'low', 'equal')){
  for(type in c('ortho', 'ind', 'singular', 'col')){
  if(type == 'ortho'){
    dat <- dat_ortho
  } else if(type == 'ind'){
    dat <- dat_independent
  } else if(type == 'singular'){
    dat <- dat_singular
  } else if(type == 'col'){
    dat <- dat_cor
  }
    dat <- dat %>%
      dplyr::mutate_all(function(x){as.numeric(scale(x))})
    if(beta_type == 'zero'){
      beta <- rep(0,5)
    }
    if(beta_type == 'high'){
      beta <- c(0,0,0.1,0,0)
    } 
    if(beta_type == 'low'){
      beta <- c(0,0,0,0,0.1)
    } 
    if(beta_type == 'equal'){
      beta <- c(0.1/sqrt(2), 0.1/sqrt(2), 0.1/sqrt(2), 0.1/sqrt(2), 0.1)
    }
    if(beta_type == 'one'){
      beta <- c(0.1,0,0.1,0,0.1)
    }
    dat <- dat %>%
     dplyr::mutate(y = V1*beta[1]+ V2*beta[2]+ V3*beta[3]+ V4*beta[4]+ V5*beta[5]+ rnorm(n=1))
  mb_model_unreg <- mboost::mboost(y ~bols(V1, V2)+bols(V3, V4) + bols(V5),
                     data = dat, control = boost_control(mstop = 1, nu = 1))
  mb_model_0.5 <- mboost::mboost(y ~bols(V1, V2, df = 0.5) +
                           bols(V3, V4, df = 0.5) +
                           bols(V5, df = 0.5),
                         data = dat, control = boost_control(mstop = 1, nu = 1))
  mb_model_0.5_sc <- mboost::mboost(dat$y ~bols(V1, V2, df = 0.5) +
                                   bols(V3, V4, df = 0.5) +
                                   bols(V5, df = 0.5),
                                 data = as.data.frame(as.matrix(dat %>% select(-y))%*%diag(params_list_1[[type]]$mb_df$scaling)),
                            control = boost_control(mstop = 1, nu = 1))
  mb_model_lam <- mboost::mboost(y ~bols(V1, V2, lambda = 10) +
                           bols(V3, V4, lambda = 10) +
                           bols(V5, lambda = 10),
                         data = dat, control = boost_control(mstop = 1, nu = 1))
  mb_model_lam_sc <- mboost::mboost(dat$y ~bols(V1, V2, lambda = 10) +
                           bols(V3, V4, lambda = 10) +
                           bols(V5, lambda = 10),
                         data =  as.data.frame(as.matrix(dat %>% select(-y))%*%diag(params_list_1[[type]]$mb_lam_10$scaling)),
                         control = boost_control(mstop = 1, nu = 1))
  gr_lasso_01_sc <- gglasso(x = as.matrix(dat %>% select(-y))%*%diag(params_list_1[[type]]$lasso_01$scaling),
                         y = dat$y, group = c(1,1,2,2,3), lambda = 0.1)
  gr_lasso_01_sc_coef <- gr_lasso_01_sc$beta %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    tidyr::gather(key = 'lambda', value = 'coef', -rowname) %>%
    dplyr::filter(coef != 0) %>%
    dplyr::slice(1:5)  %>%
    #dplyr::filter(lambda == gr_lasso_01_sc$lambda[1]) %>%
    dplyr::mutate(sel = case_when(rowname %in% c('V1','V2')~ 1,
                           rowname %in% c('V3','V4')~2,
                           rowname %in% c('V5')~ 3))
  if(length(gr_lasso_01_sc_coef$sel) == 0){gr_lasso_01_sc_coef_sel <- NA} else{gr_lasso_01_sc_coef_sel <- gr_lasso_01_sc_coef$sel}
  gr_lasso_01_sc_sel <- data.frame(model = 'grla01_sc', sel = unique(gr_lasso_01_sc_coef_sel), type = type, beta_type = beta_type)
  gr_lasso_005_sc <- gglasso(x = as.matrix(dat %>% select(-y))%*%diag(params_list_1[[type]]$lasso_005$scaling),
                         y = dat$y, group = c(1,1,2,2,3), lambda = 0.05)
  gr_lasso_005_sc_coef <- gr_lasso_005_sc$beta %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    tidyr::gather(key = 'lambda', value = 'coef', -rowname) %>%
    dplyr::filter(coef != 0) %>%
    dplyr::slice(1:5)  %>%
    dplyr::mutate(sel = case_when(rowname %in% c('V1','V2')~ 1,
                           rowname %in% c('V3','V4')~2,
                           rowname %in% c('V5')~ 3))
  if(length(gr_lasso_005_sc_coef$sel) == 0){gr_lasso_005_sc_coef_sel <- NA} else{gr_lasso_005_sc_coef_sel <- gr_lasso_005_sc_coef$sel}
  gr_lasso_005_sc_sel <- data.frame(model = 'grla005_sc', sel = unique(gr_lasso_005_sc_coef_sel), type = type, beta_type = beta_type)
  # without scaling
   gr_lasso_01 <- gglasso(x = as.matrix(dat %>% select(-y)),
                         y = dat$y, group = c(1,1,2,2,3), lambda = 0.1)
  gr_lasso_01_coef <- gr_lasso_01$beta %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    tidyr::gather(key = 'lambda', value = 'coef', -rowname) %>%
    dplyr::filter(coef != 0) %>%
    dplyr::slice(1:5)
  gr_lasso_01_coef <- gr_lasso_01_coef %>%
    dplyr::mutate(sel = case_when(rowname %in% c('V1','V2')~ 1,
                           rowname %in% c('V3','V4')~2,
                           rowname %in% c('V5')~ 3))
  if(length(gr_lasso_01_coef$sel) == 0){gr_lasso_01_coef_sel <- NA} else{gr_lasso_01_coef_sel <- gr_lasso_01_coef$sel}
  gr_lasso_01_sel <- data.frame(model = 'grla01', sel = unique(gr_lasso_01_coef_sel), type = type, beta_type = beta_type)
  gr_lasso_005 <- gglasso(x = as.matrix(dat %>% select(-y)),
                         y = dat$y, group = c(1,1,2,2,3), lambda = 0.05)
  gr_lasso_005_coef <- gr_lasso_005$beta %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    tidyr::gather(key = 'lambda', value = 'coef', -rowname) %>%
    dplyr::filter(coef != 0) %>%
    dplyr::slice(1:5)
  gr_lasso_005_coef <- gr_lasso_005_coef %>%
    dplyr::mutate(sel = case_when(rowname %in% c('V1','V2')~ 1,
                           rowname %in% c('V3','V4')~2,
                           rowname %in% c('V5')~ 3))
  if(length(gr_lasso_005_coef$sel) == 0){gr_lasso_005_coef_sel <- NA} else{gr_lasso_005_coef_sel <- gr_lasso_005_coef$sel}
  gr_lasso_005_sel <- data.frame(model = 'grla005', sel = unique(gr_lasso_005_coef_sel), type = type, beta_type = beta_type)
  ret <- ret %>% 
    bind_rows(data.frame(model = c('mb', 'mb_0.5', 'mb_0.5_sc','mb_lam', 'mb_lam_sc'),
              sel = c(mb_model_unreg$xselect(), mb_model_0.5$xselect(),
                      mb_model_0.5_sc$xselect(),
                      mb_model_lam$xselect(), mb_model_lam_sc$xselect()), type = type, beta_type = beta_type) %>%
    rbind(gr_lasso_005_sc_sel) %>%
    rbind(gr_lasso_005_sel) %>%
    rbind(gr_lasso_01_sc_sel) %>%
    rbind(gr_lasso_01_sel) %>%
    filter(!is.na(sel)))
  }
  }
  return(ret)
}

library(furrr)
plan(multisession)
iterations <- 50000
# Define the 12 data scenarios
params <- tidyr::expand_grid(iteration = 1:iterations)
# run simulation in parallel
set.seed(5)
tictoc::tic()
results_df <- params %>%
  dplyr::mutate( 
    res = future_pmap(., .f = ~run_experiment(.iteration), .progress = T,
                      .options=furrr_options(seed = TRUE, globals = c('params_list_1', 'run_experiment', 'dat_ortho', 'dat_independent', 'dat_singular', 'dat_cor'),
                                             packages = c('mboost', 'tidyverse', 'gglasso')))
  ) %>%
  unnest(cols = res)
tictoc::toc()
saveRDS(file = 'results_fixed/small_group_1_raw.RDS', results_df)
results_df %>%
  dplyr::group_by(model, type, beta_type, sel) %>%
  dplyr::tally() %>%
  dplyr::filter(!is.na(sel)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(model, type, beta_type) %>%
  dplyr::mutate(prop = round(n/sum(n,na.rm = T), 3)) %>%
saveRDS(file = 'results_fixed/small_group_1.RDS')
plan(sequential)
print('juhu 1')
################################################################################
################################################################################
# group of one vs group of 15

fit_lasso_01 <- function(iter, df){
  sim_df <- df %>% dplyr::mutate(y = rnorm(dim(df)[1]))
  gr_lasso <- gglasso(x = as.matrix(sim_df %>% select(-y)), y = sim_df$y,
                      group = c(1, rep(2,15)), lambda = 0.1)
  gr_lasso_coef <- gr_lasso$beta %>% 
    as.data.frame() %>% 
    rownames_to_column() %>%
    gather(key = 'lambda', value = 'coef', -rowname) %>%
    filter(coef != 0) %>%
    slice(1:16) 
  gr_lasso_coef <- gr_lasso_coef %>%
    filter(lambda == gr_lasso_coef$lambda[1])
  gr_lasso_coef
  gr_lasso_coef <- gr_lasso_coef %>% 
    mutate(sel = case_when(rowname %in% c('V1')~ 1, 
                           rowname %in% c( 'V2','V3','V4','V5','V6','V7','V8','V9','V10',
                                           'V11','V12','V13','V14','V15','V16')~2))
  unique(gr_lasso_coef$sel)
}
fit_lasso_005 <- function(iter, df){
  sim_df <- df %>% dplyr::mutate(y = rnorm(dim(df)[1]))
  gr_lasso <- gglasso(x = as.matrix(sim_df %>% select(-y)), y = sim_df$y,
                      group = c(1, rep(2,15)), lambda = 0.05)
  gr_lasso_coef <- gr_lasso$beta %>% 
    as.data.frame() %>% 
    rownames_to_column() %>%
    gather(key = 'lambda', value = 'coef', -rowname) %>%
    filter(coef != 0) %>%
    slice(1:16) 
  gr_lasso_coef <- gr_lasso_coef %>%
    filter(lambda == gr_lasso_coef$lambda[1])
  gr_lasso_coef
  gr_lasso_coef <- gr_lasso_coef %>% 
    mutate(sel = case_when(rowname %in% c('V1')~ 1, 
                           rowname %in% c( 'V2','V3','V4','V5','V6','V7','V8','V9','V10',
                                           'V11','V12','V13','V14','V15','V16')~2))
  unique(gr_lasso_coef$sel)
}

group_df <- data.frame(var_name = c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8',
                                    'V9', 'V10', 'V11', 'V12', 'V13', 'V14', 'V15', 'V16'),
                       group_name = as.integer(c(1, rep(2,15))))

fit_mb_lam_10 <- function(iter, df){
  formula <- as.formula(y ~bols(V1,lambda =10)+bols(V2, V3, V4,V5, V6, V7,V8, V9, V10,V11, V12, V13,V14,V15,V16, lambda=10))
  sim_df <- df %>% dplyr::mutate(y = rnorm(dim(df)[1]))
  mb_model <- mboost(formula = formula, data = sim_df, control = boost_control(mstop = 1))
  mb_model$xselect()
}
fit_mb_df <- function(iter, df){
  formula <- as.formula(y ~bols(V1,df =0.5)+bols(V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16, df=0.5))
  sim_df <- df %>% dplyr::mutate(y = rnorm(dim(df)[1]))
  mb_model <- mboost(formula = formula, data = sim_df, control = boost_control(mstop = 1))
  mb_model$xselect()
}


### Simulations to show bias in methods
set.seed(1)
dat_ortho <-  pracma::randortho(100, type = 'orthonormal')[,1:16,drop=F] %>%
  as.data.frame() 
dat_independent <-  matrix( rnorm(100*16,mean=0,sd=1), 100, 16) %>% as.data.frame() 
dat_singular <-  pracma::randortho(100, type = 'orthonormal')[,1:16,drop=F]%*%diag(c(1,rep(15,5), rep(1,5), rep(2,5)))%*%
  t(pracma::randortho(16, type = 'orthonormal')[,1:16,drop=F]) %>% as.data.frame()
dat_cor <- mvtnorm::rmvnorm(n=100, mean = 1:16, sigma = ones(16,16)-0.2 + diag(c(rep(0.2,16)))) %>%
  as.data.frame()
dat_list <- list('ortho' = dat_ortho,'singular' = dat_singular,
                 'ind' = dat_independent, 'col' = dat_cor)

params_list_2 <- dat_list %>% 
  lapply(
    function(dat){
      plan(multisession)
      mb_params_df <- get_parameters(df = dat, group_df = group_df, block_size = 10000, block_fit = fit_mb_df, rate = 0.3, initialization = T)
      plan(sequential)
      plan(multisession)
      lasso_params_01 <- get_parameters(df = dat, group_df = group_df, block_size = 10000, block_fit = fit_lasso_01, rate = 0.3)
      plan(sequential)
      print('0.5')
      plan(multisession)
      lasso_params_005 <- get_parameters(df = dat, group_df = group_df, block_size = 10000, block_fit = fit_lasso_005, rate = 0.3)
      plan(sequential)
      plan(multisession)
      mb_params_lam_10 <- get_parameters(df = dat, group_df = group_df, block_size = 10000, block_fit = fit_mb_lam_10, rate = 0.3, initialization = T)
      plan(sequential)
      print('1')
      return(list('mb_df' = mb_params_df, 'lasso_01' = lasso_params_01, 'lasso_005' = lasso_params_005,
                  'mb_lam_10' = mb_params_lam_10))
    }
  )
saveRDS(params_list_2, 'results_fixed/params_list_2.RDS')
params_list_2 <- readRDS('results_fixed/params_list_2.RDS')

run_experiment <- function(iteration){
  ret <- data.frame()
  for(beta_type in c('zero', 'high', 'low', 'equal')){
    for(type in c('ortho', 'ind', 'singular', 'col')){
      if(type == 'ortho'){
        dat <- dat_ortho
      } else if(type == 'ind'){
        dat <- dat_independent
      } else if(type == 'singular'){
        dat <- dat_singular
      } else if(type == 'col'){
        dat <- dat_cor
      }
      dat <- dat %>%
        dplyr::mutate_all(function(x){as.numeric(scale(x))})
      if(beta_type == 'zero'){
        beta <- rep(0,16)
      }
      if(beta_type == 'high'){
        beta <- c(0, rep(0.1,15))
      } 
      if(beta_type == 'low'){
        beta <- c(0.1,rep(0,15))
      } 
      if(beta_type == 'equal'){
        beta <- c(0, rep(0.1/sqrt(15),15))
      }
      if(beta_type == 'one'){
        beta <- c(0.1, 0.1, rep(0, 14))
      }
      dat <- dat %>%
        dplyr::mutate(y = V1*beta[1]+V2*beta[2]+V3*beta[3]+V4*beta[4]+V5*beta[5]+
                          V6*beta[6]+V7*beta[7]+V8*beta[8]+V9*beta[9]+V10*beta[10]+
                          V11*beta[11]+V12*beta[12]+V13*beta[13]+V14*beta[14]+V15*beta[15]+rnorm(n=1))
      mb_model_unreg <- mboost::mboost(y ~bols(V1,df=0.5)+bols(V2, V3, V4,V5, V6, V7,V8, V9, V10,V11, V12, V13,V14,V15,V16, df=0.5),
                                       data = dat, control = boost_control(mstop = 1, nu = 1))
      mb_model_0.5 <- mboost::mboost(y ~bols(V1,df=0.5)+bols(V2, V3, V4,V5, V6, V7,V8, V9, V10,V11, V12, V13,V14,V15,V16, df=0.5),
                                     data = dat, control = boost_control(mstop = 1, nu = 1))
      mb_model_0.5_sc <- mboost::mboost(dat$y ~bols(V1,df=0.5)+bols(V2, V3, V4,V5, V6, V7,V8, V9, V10,V11, V12, V13,V14,V15,V16, df=0.5),
                                        data = as.data.frame(as.matrix(dat %>% select(-y))%*%diag(params_list_2[[type]]$mb_df$scaling)),
                                        control = boost_control(mstop = 1, nu = 1))
      mb_model_lam <- mboost::mboost(dat$y ~bols(V1,lambda =10)+bols(V2, V3, V4,V5, V6, V7,V8, V9, V10,V11, V12, V13,V14,V15,V16, lambda=10),
                                     data = dat, control = boost_control(mstop = 1, nu = 1))
      mb_model_lam_sc <- mboost::mboost(dat$y ~bols(V1,lambda =10)+bols(V2, V3, V4,V5, V6, V7,V8, V9, V10,V11, V12, V13,V14,V15,V16, lambda=10),
                                        data =  as.data.frame(as.matrix(dat %>% select(-y))%*%diag(params_list_2[[type]]$mb_lam_10$scaling)),
                                        control = boost_control(mstop = 1, nu = 1))
      gr_lasso_01_sc <- gglasso(x = as.matrix(dat %>% select(-y))%*%diag(params_list_2[[type]]$lasso_01$scaling),
                                y = dat$y, group = c(1, rep(2,15)), lambda = 0.1)
      gr_lasso_01_sc_coef <- gr_lasso_01_sc$beta %>%
        as.data.frame() %>%
        tibble::rownames_to_column() %>%
        tidyr::gather(key = 'lambda', value = 'coef', -rowname) %>%
        dplyr::filter(coef != 0) %>%
        dplyr::slice(1:5)  %>%
        #dplyr::filter(lambda == gr_lasso_01_sc$lambda[1]) %>%
        dplyr::mutate(sel = case_when(rowname %in% c('V1')~ 1, 
                                      rowname %in% c( 'V2','V3','V4','V5','V6','V7','V8','V9','V10',
                                                      'V11','V12','V13','V14','V15','V16')~2))
      if(length(gr_lasso_01_sc_coef$sel) == 0){gr_lasso_01_sc_coef_sel <- NA} else{gr_lasso_01_sc_coef_sel <- gr_lasso_01_sc_coef$sel}
      gr_lasso_01_sc_sel <- data.frame(model = 'grla01_sc', sel = unique(gr_lasso_01_sc_coef_sel), type = type, beta_type = beta_type)
      gr_lasso_005_sc <- gglasso(x = as.matrix(dat %>% select(-y))%*%diag(params_list_2[[type]]$lasso_005$scaling),
                                 y = dat$y, group = c(1, rep(2,15)), lambda = 0.05)
      gr_lasso_005_sc_coef <- gr_lasso_005_sc$beta %>%
        as.data.frame() %>%
        tibble::rownames_to_column() %>%
        tidyr::gather(key = 'lambda', value = 'coef', -rowname) %>%
        dplyr::filter(coef != 0) %>%
        dplyr::slice(1:5)  %>%
        dplyr::mutate(sel = case_when(rowname %in% c('V1')~ 1, 
                                      rowname %in% c( 'V2','V3','V4','V5','V6','V7','V8','V9','V10',
                                                      'V11','V12','V13','V14','V15','V16')~2))
      if(length(gr_lasso_005_sc_coef$sel) == 0){gr_lasso_005_sc_coef_sel <- NA} else{gr_lasso_005_sc_coef_sel <- gr_lasso_005_sc_coef$sel}
      gr_lasso_005_sc_sel <- data.frame(model = 'grla005_sc', sel = unique(gr_lasso_005_sc_coef_sel), type = type, beta_type = beta_type)
      # without scaling
      gr_lasso_01 <- gglasso(x = as.matrix(dat %>% select(-y)),
                             y = dat$y, group = c(1, rep(2,15)), lambda = 0.1)
      gr_lasso_01_coef <- gr_lasso_01$beta %>%
        as.data.frame() %>%
        tibble::rownames_to_column() %>%
        tidyr::gather(key = 'lambda', value = 'coef', -rowname) %>%
        dplyr::filter(coef != 0) %>%
        dplyr::slice(1:5)
      gr_lasso_01_coef <- gr_lasso_01_coef %>%
        dplyr::mutate(sel = case_when(rowname %in% c('V1')~ 1, 
                                      rowname %in% c( 'V2','V3','V4','V5','V6','V7','V8','V9','V10',
                                                      'V11','V12','V13','V14','V15','V16')~2))
      if(length(gr_lasso_01_coef$sel) == 0){gr_lasso_01_coef_sel <- NA} else{gr_lasso_01_coef_sel <- gr_lasso_01_coef$sel}
      gr_lasso_01_sel <- data.frame(model = 'grla01', sel = unique(gr_lasso_01_coef_sel), type = type, beta_type = beta_type)
      gr_lasso_005 <- gglasso(x = as.matrix(dat %>% select(-y)),
                              y = dat$y, group = c(1, rep(2,15)), lambda = 0.05)
      gr_lasso_005_coef <- gr_lasso_005$beta %>%
        as.data.frame() %>%
        tibble::rownames_to_column() %>%
        tidyr::gather(key = 'lambda', value = 'coef', -rowname) %>%
        dplyr::filter(coef != 0) %>%
        dplyr::slice(1:5)
      gr_lasso_005_coef <- gr_lasso_005_coef %>%
        dplyr::mutate(sel = case_when(rowname %in% c('V1')~ 1, 
                                      rowname %in% c( 'V2','V3','V4','V5','V6','V7','V8','V9','V10',
                                                      'V11','V12','V13','V14','V15','V16')~2))
      if(length(gr_lasso_005_coef$sel) == 0){gr_lasso_005_coef_sel <- NA} else{gr_lasso_005_coef_sel <- gr_lasso_005_coef$sel}
      gr_lasso_005_sel <- data.frame(model = 'grla005', sel = unique(gr_lasso_005_coef_sel), type = type, beta_type = beta_type)
      ret <- ret %>% 
        bind_rows(data.frame(model = c('mb', 'mb_0.5', 'mb_0.5_sc','mb_lam', 'mb_lam_sc'),
                             sel = c(mb_model_unreg$xselect(), mb_model_0.5$xselect(),
                                     mb_model_0.5_sc$xselect(),
                                     mb_model_lam$xselect(), mb_model_lam_sc$xselect()), type = type, beta_type = beta_type) %>%
                    rbind(gr_lasso_005_sc_sel) %>%
                    rbind(gr_lasso_005_sel) %>%
                    rbind(gr_lasso_01_sc_sel) %>%
                    rbind(gr_lasso_01_sel) %>%
                    filter(!is.na(sel)))
    }
  }
  return(ret)
}

library(furrr)
plan(multisession)
iterations <- 50000
# Define the 12 data scenarios
params <- tidyr::expand_grid(iteration = 1:iterations)
# run simulation in parallel
set.seed(5)
tictoc::tic()
results_df <- params %>%
  dplyr::mutate( 
    res = future_pmap(., .f = ~run_experiment(.iteration), .progress = T,
                      .options=furrr_options(seed = TRUE, globals = c('params_list_2', 'run_experiment', 'dat_ortho', 'dat_independent', 'dat_singular', 'dat_cor'),
                                             packages = c('mboost', 'tidyverse', 'gglasso')))
  ) %>%
  unnest(cols = res)
tictoc::toc()
saveRDS(file = 'results_fixed/big_group_1_raw.RDS', results_df)
results_df %>%
  dplyr::group_by(model, type, beta_type, sel) %>%
  dplyr::tally() %>%
  dplyr::filter(!is.na(sel)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(model, type, beta_type) %>%
  dplyr::mutate(prop = round(n/sum(n,na.rm = T), 3)) %>%
  saveRDS(file = 'results_fixed/big_group_1.RDS')
plan(sequential)


print('juhu 2')
################################################################################
################################################################################
# Categorical data 



fit_lasso_01 <- function(iter, df){
  
  sim_df <- df %>% dplyr::mutate(y = rnorm(dim(df)[1]))
  
  gr_lasso <- gglasso(x = as.matrix(sim_df %>% select(-y)), y = sim_df$y,
                      group = c(1,1,2,2,3,4), lambda = 0.1)
  gr_lasso_coef <- gr_lasso$beta %>% 
    as.data.frame() %>% 
    rownames_to_column() %>%
    gather(key = 'lambda', value = 'coef', -rowname) %>%
    filter(coef != 0) %>%
    slice(1:6) 
  gr_lasso_coef <- gr_lasso_coef %>%
    filter(lambda == gr_lasso_coef$lambda[1])
  gr_lasso_coef
  gr_lasso_coef <- gr_lasso_coef %>% 
    mutate(sel = case_when(rowname %in% c('V1B','V1C')~ 1, 
                           rowname %in% c('V2B','V2C')~2, 
                           rowname %in% c('V3B')~ 3, rowname %in% c('V4')~ 4))
  unique(gr_lasso_coef$sel)
}
fit_lasso_005 <- function(iter, df){
  
  sim_df <- df %>% dplyr::mutate(y = rnorm(dim(df)[1]))
  
  gr_lasso <- gglasso(x = as.matrix(sim_df %>% select(-y)), y = sim_df$y,
                      group = c(1,1,2,2,3,4), lambda = 0.05)
  gr_lasso_coef <- gr_lasso$beta %>% 
    as.data.frame() %>% 
    rownames_to_column() %>%
    gather(key = 'lambda', value = 'coef', -rowname) %>%
    filter(coef != 0) %>%
    slice(1:6) 
  gr_lasso_coef <- gr_lasso_coef %>%
    filter(lambda == gr_lasso_coef$lambda[1])
  gr_lasso_coef
  gr_lasso_coef <- gr_lasso_coef %>% 
    mutate(sel = case_when(rowname %in% c('V1B','V1C')~ 1, 
                           rowname %in% c('V2B','V2C')~2, 
                           rowname %in% c('V3B')~ 3, rowname %in% c('V4')~ 4))
  unique(gr_lasso_coef$sel)
}

group_df <- data.frame(var_name = c('V1B','V1C','V2C','V2B','V3B','V4'),
                       group_name = as.integer(c(1,1,2,2,3,4)))

fit_mb_lam_10 <- function(iter, df){
  formula <- as.formula(y ~bols(V1B,V1C,df=0.5)+bols(V2C,V2B,df=0.5)+bols(V3B,df=0.5)+bols(V4,df=0.5))
  sim_df <- df %>% dplyr::mutate(y = rnorm(dim(df)[1]))
  mb_model <- mboost(formula = formula, data = sim_df, control = boost_control(mstop = 1))
  mb_model$xselect()
}
fit_mb_df <- function(iter, df){
  formula <- as.formula(y~bols(V1B,V1C,df=0.5)+bols(V2C,V2B,df=0.5)+bols(V3B,df=0.5)+bols(V4,df=0.5))
  sim_df <- df %>% dplyr::mutate(y = rnorm(dim(df)[1]))
  mb_model <- mboost(formula = formula, data = sim_df, control = boost_control(mstop = 1))
  mb_model$xselect()
}


### Simulations to show bias in methods
set.seed(1)
dat <- data.frame(V1 = factor(sample(LETTERS[1:3], 120, replace=TRUE, prob=c(0.7, 0.15,0.15))), 
                  V2 = factor(sample(LETTERS[1:3], 120, replace=TRUE, prob=c(1/3, 1/3,1/3))),
                  V3 = factor(sample(LETTERS[1:2], 120, replace= TRUE, prob=c(0.5,0.5))),
                  V4 = rnorm(120))
dat_model <- as.data.frame(model.matrix(~., dat)[,-1])

      plan(multisession)
      mb_params_df <- get_parameters(df = dat_model, group_df = group_df, block_size = 10000, block_fit = fit_mb_df, rate = 0.3)
      plan(sequential)
      plan(multisession)
      lasso_params_01 <- get_parameters(df = dat_model, group_df = group_df, block_size = 10000, block_fit = fit_lasso_01, rate = 0.3)
      plan(sequential)
      print('0.5')
      plan(multisession)
      lasso_params_005 <- get_parameters(df = dat_model, group_df = group_df, block_size = 10000, block_fit = fit_lasso_005, rate = 0.3)
      plan(sequential)
      plan(multisession)
      mb_params_lam_10 <- get_parameters(df = dat_model, group_df = group_df, block_size = 10000, block_fit = fit_mb_lam_10, rate = 0.8)
      plan(sequential)
#saveRDS(params_list_1, 'results_fixed/params_list_1.RDS')
#params_list_1 <- readRDS('results_fixed/params_list_1.RDS')

run_experiment <- function(iteration){
  ret <- data.frame()
  for(beta_type in c('zero', 'high', 'low', 'equal')){
      dat_model <- dat_model %>%
        dplyr::mutate_all(function(x){as.numeric(scale(x))})
      if(beta_type == 'zero'){
        beta <- rep(0,6)
      }
      if(beta_type == 'high'){
        beta <- c(0,0,0.1,0,0,0)
      } 
      if(beta_type == 'low'){
        beta <- c(0,0,0,0,0,0.1)
      } 
      if(beta_type == 'equal'){
        beta <- c(0.1/sqrt(2), 0.1/sqrt(2), 0.1/sqrt(2), 0.1/sqrt(2), 0.1, 0.1)
      }
      if(beta_type == 'one'){
        beta <- c(0.1,0,0.1,0,0.1,0.1)
      }
      dat_model <- dat_model %>%
        dplyr::mutate(y = V1B*beta[1]+ V1C*beta[2]+ V2C*beta[3]+ V2B*beta[4]+ V3B*beta[5]+ V4*beta[6]+rnorm(1))
      #dat$y <- dat_model$y
      mb_model_unreg <- mboost::mboost(y ~bols(V1B,V1C)+bols(V2C,V2B)+bols(V3B)+bols(V4),
                                       data = dat_model, control = boost_control(mstop = 1, nu = 1))
      mb_model_0.5 <- mboost::mboost(y ~bols(V1B,V1C,df=0.5)+bols(V2C,V2B,df=0.5)+bols(V3B,df=0.5)+bols(V4,df=0.5),
                                     data = dat_model, control = boost_control(mstop = 1, nu = 1))
      dat_mb <- as.data.frame(as.matrix(dat_model %>% select(-y))%*%diag(mb_params_df$scaling))
      colnames(dat_mb) <- colnames(dat_model)[1:6]
      mb_model_0.5_sc <- mboost::mboost(dat_model$y ~bols(V1B,V1C,df=0.5)+bols(V2C,V2B,df=0.5)+bols(V3B,df=0.5)+bols(V4,df=0.5),
                                        data = dat_mb,
                                        control = boost_control(mstop = 1, nu = 1))
      mb_model_lam <- mboost::mboost(y ~bols(V1B,V1C,lambda=10)+bols(V2C,V2B,lambda=10)+bols(V3B,lambda=10)+bols(V4,lambda=10),
                                     data = dat_model, control = boost_control(mstop = 1, nu = 1))
      dat_lam_mb <- as.data.frame(as.matrix(dat_model %>% select(-y))%*%diag(mb_params_lam_10$scaling))
      colnames(dat_lam_mb) <- colnames(dat_model)[1:6]
      mb_model_lam_sc <- mboost::mboost(dat_model$y ~bols(V1B,V1C,lambda=10)+bols(V2C,V2B,lambda=10)+bols(V3B,lambda=10)+bols(V4,lambda=10),
                                        data =  dat_lam_mb,
                                        control = boost_control(mstop = 1, nu = 1))
      gr_lasso_01_sc <- gglasso(x = as.matrix(dat_model %>% select(-y))%*%diag(lasso_params_01$scaling),
                                y = dat_model$y, group = c(1,1,2,2,3,4), lambda = 0.1)
      gr_lasso_01_sc_coef <- gr_lasso_01_sc$beta %>%
        as.data.frame() %>%
        tibble::rownames_to_column() %>%
        tidyr::gather(key = 'lambda', value = 'coef', -rowname) %>%
        dplyr::filter(coef != 0) %>%
        dplyr::slice(1:6)  %>%
        #dplyr::filter(lambda == gr_lasso_01_sc$lambda[1]) %>%
        dplyr::mutate(sel = case_when(rowname %in% c('V1B','V1C')~ 1, 
                                      rowname %in% c('V2B','V2C')~2, 
                                      rowname %in% c('V3B')~ 3, rowname %in% c('V4')~ 4))
      if(length(gr_lasso_01_sc_coef$sel) == 0){gr_lasso_01_sc_coef_sel <- NA} else{gr_lasso_01_sc_coef_sel <- gr_lasso_01_sc_coef$sel}
      gr_lasso_01_sc_sel <- data.frame(model = 'grla01_sc', sel = unique(gr_lasso_01_sc_coef_sel), beta_type = beta_type)
      gr_lasso_005_sc <- gglasso(x = as.matrix(dat_model %>% select(-y))%*%diag(lasso_params_005$scaling),
                                 y = dat_model$y, group = c(1,1,2,2,3,4), lambda = 0.05)
      gr_lasso_005_sc_coef <- gr_lasso_005_sc$beta %>%
        as.data.frame() %>%
        tibble::rownames_to_column() %>%
        tidyr::gather(key = 'lambda', value = 'coef', -rowname) %>%
        dplyr::filter(coef != 0) %>%
        dplyr::slice(1:6)  %>%
        dplyr::mutate(sel = case_when(rowname %in% c('V1B','V1C')~ 1, 
                                      rowname %in% c('V2B','V2C')~2, 
                                      rowname %in% c('V3B')~ 3, rowname %in% c('V4')~ 4))
      if(length(gr_lasso_005_sc_coef$sel) == 0){gr_lasso_005_sc_coef_sel <- NA} else{gr_lasso_005_sc_coef_sel <- gr_lasso_005_sc_coef$sel}
      gr_lasso_005_sc_sel <- data.frame(model = 'grla005_sc', sel = unique(gr_lasso_005_sc_coef_sel), beta_type = beta_type)
      # without scaling
      gr_lasso_01 <- gglasso(x = as.matrix(dat_model %>% select(-y)),
                             y = dat_model$y, group = c(1,1,2,2,3,4), lambda = 0.1)
      gr_lasso_01_coef <- gr_lasso_01$beta %>%
        as.data.frame() %>%
        tibble::rownames_to_column() %>%
        tidyr::gather(key = 'lambda', value = 'coef', -rowname) %>%
        dplyr::filter(coef != 0) %>%
        dplyr::slice(1:6)
      gr_lasso_01_coef <- gr_lasso_01_coef %>%
        dplyr::mutate(sel = case_when(rowname %in% c('V1B','V1C')~ 1, 
                                      rowname %in% c('V2B','V2C')~2, 
                                      rowname %in% c('V3B')~ 3, rowname %in% c('V4')~ 4))
      if(length(gr_lasso_01_coef$sel) == 0){gr_lasso_01_coef_sel <- NA} else{gr_lasso_01_coef_sel <- gr_lasso_01_coef$sel}
      gr_lasso_01_sel <- data.frame(model = 'grla01', sel = unique(gr_lasso_01_coef_sel), beta_type = beta_type)
      gr_lasso_005 <- gglasso(x = as.matrix(dat_model %>% select(-y)),
                              y = dat_model$y, group = c(1,1,2,2,3,4), lambda = 0.05)
      gr_lasso_005_coef <- gr_lasso_005$beta %>%
        as.data.frame() %>%
        tibble::rownames_to_column() %>%
        tidyr::gather(key = 'lambda', value = 'coef', -rowname) %>%
        dplyr::filter(coef != 0) %>%
        dplyr::slice(1:6)
      gr_lasso_005_coef <- gr_lasso_005_coef %>%
        dplyr::mutate(sel = case_when(rowname %in% c('V1B','V1C')~ 1, 
                                      rowname %in% c('V2B','V2C')~2, 
                                      rowname %in% c('V3B')~ 3, rowname %in% c('V4')~ 4))
      if(length(gr_lasso_005_coef$sel) == 0){gr_lasso_005_coef_sel <- NA} else{gr_lasso_005_coef_sel <- gr_lasso_005_coef$sel}
      gr_lasso_005_sel <- data.frame(model = 'grla005', sel = unique(gr_lasso_005_coef_sel),  beta_type = beta_type)
      ret <- ret %>% 
        bind_rows(data.frame(model = c('mb', 'mb_0.5', 'mb_0.5_sc','mb_lam', 'mb_lam_sc'),
                             sel = c(mb_model_unreg$xselect(), mb_model_0.5$xselect(),
                                     mb_model_0.5_sc$xselect(),
                                     mb_model_lam$xselect(), mb_model_lam_sc$xselect()), beta_type = beta_type) %>%
                    rbind(gr_lasso_005_sc_sel) %>%
                    rbind(gr_lasso_005_sel) %>%
                    rbind(gr_lasso_01_sc_sel) %>%
                    rbind(gr_lasso_01_sel) %>%
                    filter(!is.na(sel)))
  }
  return(ret)
}

library(furrr)
plan(multisession)
iterations <- 50000
# Define the 12 data scenarios
params <- tidyr::expand_grid(iteration = 1:iterations)
# run simulation in parallel
set.seed(5)
tictoc::tic()
results_df <- params %>%
  dplyr::mutate( 
    res = future_pmap(., .f = ~run_experiment(.iteration), .progress = T,
                      .options=furrr_options(seed = TRUE, globals = c('params_list_1', 'run_experiment', 
                                                                      'mb_params_df','mb_params_lam_10', 'lasso_params_005', 'lasso_params_01', 'dat_model'),
                                             packages = c('mboost', 'tidyverse', 'gglasso')))
  ) %>%
  unnest(cols = res)
tictoc::toc()
saveRDS(file = 'results_fixed/cat_1_raw.RDS', results_df)
results_df %>%
  dplyr::group_by(model, beta_type, sel) %>%
  dplyr::tally() %>%
  dplyr::filter(!is.na(sel)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(model,  beta_type) %>%
  dplyr::mutate(prop = round(n/sum(n,na.rm = T), 3)) %>%
  saveRDS(file = 'results_fixed/cat_1.RDS')
plan(sequential)

