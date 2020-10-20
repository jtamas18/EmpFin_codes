# structure:
# list X stores all relevant information
# the elements of X are other lists ("say Y-s"DGP_[n]"s), each storing information for 
# different data generating processes (dgp-s)
# the "DGP_[n]" lists store the dgp parameters ("DGP_[n]$gdp_parameters"),
# the simulated trajectories ("DGP_[n]$trajectories"), the model selection parameters 
# (in "DGP_[n]$modelling$model_[n]$model_selection_params) and the
# "chosen_orders" matrix (in "DGP_[n]$modelling$model_[n]$chosen_orders")
# in which the (i,j) element is the ratio "(times ARMA(i,j) was chosen)/
# (number of simulated trajectories)

def_dgp_parameters <- function(X, p_lag, q_lag, ar_coefs, ma_coefs, 
                               intercept, len = 100, noise = 1){
  
  # stores the data generating process parameters in list "dgp_param".
  # returns list "X" with dgp_param" as first elemenet.
  # X will store all other values related to the analysis
  
  # inputs: 
  #analysis storage list X
  # DGP parameters
  
  # store parameters:
  dgp_case <- list()      # this list contains all relevant info for this DGP
  dgp_parameters <- list()
  dgp_parameters$ts_length <- len
  dgp_parameters$noise <- noise
  dgp_parameters$p_lags <- p_lag
  dgp_parameters$q_lags <- q_lag
  dgp_parameters$intercept <- intercept
  dgp_parameters$ar_coefs <- ar_coefs
  dgp_parameters$ma_coefs <- ma_coefs
  dgp_parameters$sim_num <- 1000
  
  dgp_case$dgp_parameters <- dgp_parameters
  
  # append DGP list to storage matrix X:
  X[[length(X)+1]] <- dgp_case
  
  # give name "DGP_[n]" key to the DGP list in X: 
  names(X) <- c(names(X[1:length(names(X))-1]), paste0("DGP_", toString(length(X))))
  
  #Check whether the process determined by the ARMA coefficients is stationary
  #It is carried out with the help of the phi-polinom
  #The real and imaginary parts of the roots of the polinom are collected in the cplx_matrix
  #The lenght of the complex roots are calculated and compared to 1
  #If any of the roots is outside the unit circle the User is asked to modify the AR-coefficients
  roots <- polyroot(c(1, -dgp_parameters$ar_coefs))
  cplx_matrix <- matrix(0, nrow = length(roots), ncol = 4)
  for (i in 1:length(roots)) {
    cplx_matrix[i,1] <- Re(roots)[i]
    cplx_matrix[i,2] <- Im(roots)[i]
    cplx_matrix[i,3] <- sqrt(sum(cplx_matrix[i,1:2]^2))
    cplx_matrix[i,4] <- cplx_matrix[i,3] > 1
  }
  if(prod(cplx_matrix[,4]) == 0){
    print("The process determined by the coefficients is not stationary. Please modify the components of ar_coefs accordingly!")
  } else {
    
  }
  
  return(X)
}

simulate_dgp <- function(X, dgp_index){
  # simulates the DGP defined in "dgp_parameters" times "dgp_parameters$sim_num"
  # stores the sim_num different trajectories in a matrix (size= ts_length * sim_num)
  # where each column is a different trajecotry
  # appends the matrix to X$DGP_[n]$trajectories  
  
  # inputs: 
  # the storage list X
  # the index of the DGP in X; eg. for dgp_index = 1, values used in the simulation
  # are taken from "X$DGP_1$..."
  
  # take parameters from "X$DGP_[n]$dgp_parameters:
  ts_length <- X[[dgp_index]]$dgp_parameters$ts_length
  noise <- X[[dgp_index]]$dgp_parameters$noise
  p_lags <- X[[dgp_index]]$dgp_parameters$p_lags
  q_lags <- X[[dgp_index]]$dgp_parameters$q_lags
  intercept <- X[[dgp_index]]$dgp_parameters$intercept
  ar_coefs <- X[[dgp_index]]$dgp_parameters$ar_coefs
  ma_coefs <- X[[dgp_index]]$dgp_parameters$ma_coefs
  sim_num <- X[[dgp_index]]$dgp_parameters$sim_num
  
  
  trajectories <- matrix(nrow = ts_length, ncol = sim_num)  
  
  # calculate expected value of ARMA(p,q) process:
  expected_value <- intercept/(1-sum(ar_coefs))
  
  # "lengthening" length of the time series because of the lags:
  lengthening <- max(q_lags,p_lags)
  
  # simulate errors with variance = noise:
  Us <- rnorm(n = sim_num * (ts_length + lengthening), sd = noise)
  Us <- matrix(Us, nrow = (ts_length + lengthening), ncol = sim_num, byrow = FALSE) 
  
  # simulate the trajectories:
  trajectories <- matrix(nrow = (ts_length + lengthening), ncol = sim_num)
  
  for (t in c(1:(ts_length+lengthening))){
    if (t <= lengthening){
      trajectories[t,] <- expected_value
    } 
    else {
      trajectories[t,] <- intercept  + rev(ar_coefs) %*% trajectories[(t-p_lags):(t-1),] +
        rev(ma_coefs) %*% Us[(t-q_lags):(t-1),] + Us[t,]
    }
  }
  
  # drop data only used for generation ("lengthening"):
  trajectories <- trajectories[(lengthening+1):nrow(trajectories),]
  
  # store trajectories in storage matrix X
  X[[dgp_index]]$trajectories <- trajectories
  
  return(X)
  
}

def_model_fitting_params <- function(X, dgp_index, inf_criterion, max_p = 4, max_q = 4){
  
  # stores model selection parameters in storage matrix X
  
  # inputs:
  # storage list X
  # the index of the DGP in X; eg. for dgp_index = 1, values used in the simulation
  # inf_criterion used for model selection: either "aic" (Akaike Information Criterion) or
  # "bic" (Bayesian Information Criterion)
  # maximum number of p and q lags considered
  
  
  # store model_params in list
  model_params <- list("inf_crit" = inf_criterion,
                       "max_p_order" = max_p, "max_q_order" = max_q)
  
  # create storage place for the model in X and store it:
  
  # if "modelling" exists in X, pull it, otherwise create it
  if ( "modelling" %in% names(X[[dgp_index]])){
    modelling <- X[[dgp_index]]$modelling
  } else {
    modelling <- list()
  }
  
  new_model <- list()
  new_model$model_selection_params <- model_params
  
  modelling[[length(modelling) + 1]] <- new_model
  X[[dgp_index]]$modelling <- modelling
  
  # give name "model_[n]" key to the model list in X: 
  names(X[[dgp_index]]$modelling) <- 
    c(names(X[[dgp_index]]$modelling[1:length(names(X[[dgp_index]]$modelling))-1]), 
      paste0("model_", toString(length(X[[dgp_index]]$modelling))))
  
  return(X)
  
}

model_fitting <- function(X,dgp_index, model_index){
  # fits different ARIMA(p,q) models to the simulated trajectiories in X[[dgp_index]]$trajectories
  # p,q orders go from (0,0) to (map_p, max_q) defined in
  # X[[dgp_index]]$modelling[[model_index]]$model_selection_params
  # returns X with a "chosen_orders" matrix attached:  
  # in which the (i,j) element is the ratio "(times ARMA(i,j) was chosen)/
  # (number of simulated trajectories)
  
  # take parameters from X
  max_p  <- X[[dgp_index]]$modelling[[model_index]]$model_selection_params$max_p_order
  max_q  <- X[[dgp_index]]$modelling[[model_index]]$model_selection_params$max_q_order
  inf_crit_used <- X[[dgp_index]]$modelling[[model_index]]$model_selection_params$inf_crit
  sim_num <- X[[dgp_index]]$dgp_parameters$sim_num
  
  #initialize result storage matrix
  chosen_orders <- matrix(0, nrow = (max_p + 1), ncol = (max_q + 1))
  
  # name columns and rows like "q = 1" etc...
  colname <- rep("q = ", times = (max_q +1))
  rowname <- rep("p = ", times = (max_p + 1))
  for (col_num in 1:(max_q+1)){
    colname[col_num] <- paste0(colname[col_num], toString(col_num -1))
  }
  for (row_num in 1:(max_p+1)){
    rowname[row_num] <- paste0(rowname[row_num], toString(row_num -1))
  }
  
  colnames(chosen_orders) <- colname
  rownames(chosen_orders) <- rowname
  
  # iterate over each trajectory in trajectories matrix
  # iterate over each p,q pair
  # fit arma(p,q) model for each pair
  # choose best (p_0,q_0) order based on inf criterion
  # add +1 to cell (p_0,q_0) in the chosen_orders matrix
  
  
  for (column in 1:ncol(X[[dgp_index]]$trajectories)){ 
    
    trajectory <- X[[dgp_index]]$trajectories[,column]
    best_orders <- c(0,0) # stores the (p,q) orders of the model with the lowest inf crit so far
    
    
    #Let's use auto.arima which considers only those series being stationary based on some stationerity test
    fitted_arima <- forecast::auto.arima(trajectory, max.p = max_p, max.q = max_q, stationary = TRUE, ic = inf_crit_used, test = "adf") 
    best_orders[1] <- length(fitted_arima$model$phi)
    best_orders[2] <- length(fitted_arima$model$theta)
    
    
    # adds +1 the cell representing the orders chosen
    chosen_orders[best_orders[1]+1,best_orders[2]+1] <- 
      chosen_orders[best_orders[1] +1,best_orders[2]+1] + 1
    
    # making numbers in chosen_order proportional to sim_num
    #chosen_orders <- chosen_orders / sim_num
    
    # store chosen_orders matrix in X:
    X[[dgp_index]]$modelling[[model_index]]$chosen_orders <- chosen_orders
    
  }
  
  return(X)
  
}

create_summary_table <- function(X, inf_crit_investigated = c("aic", "bic"), p_lags = 1, q_lags = 1){
  
  true_ratio <- matrix(0, nrow = length(possible_length_values), ncol = length(noise_levels))
  # name columns and rows like "q = 1" etc...
  colname <- rep("noise = ", times = (length(noise_levels)))
  rowname <- rep("n = ", times = (length(possible_length_values)))
  
  for (i in 1:(length(possible_length_values))) {
    rowname[i] <- paste0("n = ", toString(possible_length_values[i]))
  }
  for (j in 1:(length(noise_levels))) {
    colname[j] <- paste0("noise = ", toString(noise_levels[j]))
  }
  colnames(true_ratio) <- colname
  rownames(true_ratio) <- rowname
  
  
  if (inf_crit_investigated == "aic") {
    model_index <- 1
  } else {
    model_index <- 2
  }
  
  index <- 0
  for (i in 1:length(possible_length_values)) {
    for (j in 1:length(noise_levels)) {
      index <- index + 1
      true_ratio[i,j] <- X[[index]]$modelling[[model_index]]$chosen_orders[p_lags + 1, q_lags +1]/sum(X[[index]]$modelling[[model_index]]$chosen_orders)
    }
    
  }
  
  print(paste0("The identification performance of ", inf_crit_investigated, " is given in the below table"))
  true_ratio
}


# example:
# how to use the functions:

X <- list()
X <- def_dgp_parameters(X, len = 25, noise = 0.5, p_lag = 1, q_lag = 1, 
                        ar_coefs = c(-0.5), ma_coefs = c(0.1), intercept = 0.7)
X <- simulate_dgp(X, 1)
X <- def_model_fitting_params(X, 1, inf_criterion = "aic", max_p = 2, max_q = 2)
X <- model_fitting(X,1,1)
create_summary_table(X, inf_crit_investigated = "aic")

#Let's use the above defined functions for the purposes of the analysis!
#Creating a grid based on different
#Noise levels and
#Time-series length

possible_length_values <- c(25,50,100,1000)
noise_levels <- c(0.5,1,2,4)

#define X
X <- list()
#Definition of the altogether 16 (= 4 * 4) DGP

for (i in 1:length(possible_length_values)) {
  for (j in 1:length(noise_levels)) {
    X <- def_dgp_parameters(X, len = possible_length_values[i], noise = noise_levels[j], p_lag = 1, q_lag = 1, 
                            ar_coefs = c(-0.5), ma_coefs = c(0), intercept = 0.6 )
  }
}



#Definition of the model fitting parameters

for (i in 1:(length(possible_length_values)*length(noise_levels))) {
  set.seed(i)
  X <- simulate_dgp(X, i)
  X <- def_model_fitting_params(X, i, inf_criterion = "aic", max_p = 2, max_q = 2)
  X <- def_model_fitting_params(X, i, inf_criterion = "bic", max_p = 2, max_q = 2)
}

#Fitting ARMA models to the simulated series from 16 DGP. The best model is selected either by AIC (j = 1) or by BIC (j = 2)

for (i in 1:(length(possible_length_values)*length(noise_levels))) {
  for (j in 1:2) {
    X <- model_fitting(X,i,j)
  }
}



#Summarizing the result for either "aic" or "bic" for different noise-sample size combinations
create_summary_table(X, inf_crit_investigated = "bic")