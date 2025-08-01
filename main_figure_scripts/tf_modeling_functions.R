suppressPackageStartupMessages(library(PseudotimeDE))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(slingshot))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(irlba))
library(Seurat)
library(RColorBrewer)
library(ggplot2)

### FUNCTIONS ###

# Sum of square of differences (SSD)
ssd <- function(array1, array2) {
  diff = array1 - array2
  return(sum(diff * diff))
}

# Excitatory Hill Function
hill_ex <- function(X,X0,n) {
  a = (X/X0)**n
  return(a/(1+a))
}

hill_inh <- function(X,X0,n) {
  a = (X/X0)**n
  return(1/(1+a))
}

target_model <- function(t, Xs, param) { # only target expression
  
  tf_1 = as.numeric(Xs[1])
  target = as.numeric(Xs[2])
  
  g0 = param[1]
  g1 = param[2] 
  k = param[3] # k=1/t between max and min exprs
  n_1 = param[4]
  th_1 = param[5]
  
  dydt = g0 * (g1 + (1-g1) * hill_inh(tf_1, th_1, n_1)) - k * target
  
  dydt = as.numeric(dydt)
  return(dydt) # dydt
}

# 4th order Runge-Kutta (RK4) for a generic multi-variable system, see Part 3A
RK4_generic <- function(derivs, Xn, X, t.total, dt, param) {
  # derivs: the function of the derivatives 
  # X0: initial condition, a vector of multiple variables
  # t.total: total simulation time, assuming t starts from 0 at the beginning
  # dt: time step size 
  
  t_all = length(t.total) 
  X_all = rep(0, t_all)
  X_all[1] = unlist(Xn)[2] 
  
  for (i in 1:(t_all-1)) {
    
    dtx = dt[i]
    
    t_0 = t.total[i]
    t_0.5 = t_0 + 0.5*dtx
    t_1 = t_0 + dtx
    
    tf_0 = X[i] # corresp TF val
    tf_0.5 = mean(X[i], X[i+1])
    tf_1 = X[i+1]
    
    k1 = dtx * derivs(t_0,   Xs = c(tf_0, X_all[i]),          param=param)
    k2 = dtx * derivs(t_0.5, Xs = c(tf_0.5, X_all[i] + k1/2), param=param)
    k3 = dtx * derivs(t_0.5, Xs = c(tf_0.5, X_all[i] + k2/2), param=param)
    k4 = dtx * derivs(t_1,   Xs = c(tf_1, X_all[i] + k3),     param=param)
    
    X_all[i+1] = X_all[i] + (k1+2*k2+2*k3+k4)/6
    
  } 
  return(X_all) # outputting predicted gene expression, X(t)
}

cal_error <- function(derivs, data_exp, param){  
  
  # t.total <- 1:nrow(data_exp)
  t = data_exp[,1]
  t.total = t
  dt = diff(t)
  
  Xs = data_exp[,2] # regulator expression
  Xn = unlist(data_exp[1,2:3])  # the first data point as the initial condition
  
  sim = RK4_generic(derivs = derivs, Xn = Xn, X = Xs,
                    t.total = t.total, dt = dt, param)
  
  ind = seq(1, nrow(data_exp), 5)
  
  err_x = ssd(data_exp[ind,1], t.total[ind])
  err_y = ssd(data_exp[ind,3], sim[ind])
  
  error = err_x+err_y
  return(error)
}

# FUNCTIONS
model_sim <- function(param, data_exp, derivs, lower_range, upper_range) {
  
  final_error <- c()
  f <- function(param) {
    return(cal_error(derivs = derivs, data_exp = as.matrix(data_exp), 
                     param=param))
  }
  
  fitted = nlminb(param, objective = f, 
                  lower=lower_range, 
                  upper=upper_range,
                  control=list(trace=1, rel.tol=1e-5)) # updating error
  p_est <- fitted$par
  final_error <- c(final_error, f(p_est))
  
  t = data_exp[,1]
  t.total = t
  dt = diff(t)
  
  X = data_exp[,2]
  Xn = as.numeric(data_exp[1,2:3])
  
  results_simu = RK4_generic(derivs = derivs, Xn = Xn, X = X,
                             t.total = t.total, dt = dt, param = p_est)
  
  return(list(results_simu, final_error, p_est))
}
