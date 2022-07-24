source("stratified_functions.R")



################## MIRLA Replication #########################

#Kalamazoo replication
# start by auditing only contest with smallest margin
# order is c(whitmer, schuette, other/invalid)
reported_votes_CVR <- c(3765, 1349, 180)
reported_votes_noCVR <- c(16934, 4220, 1218)
dm_CVR <- (reported_votes_CVR[1] - reported_votes_CVR[2]) / sum(reported_votes_CVR)
dm_noCVR <- (reported_votes_noCVR[1] - reported_votes_noCVR[2]) / sum(reported_votes_noCVR)
#stratum weights from the total number of votes in each stratum
N <- c(sum(reported_votes_CVR), sum(reported_votes_noCVR))
w <- N/sum(N)
diluted_margin <- ((reported_votes_CVR[1] + reported_votes_noCVR[1]) - (reported_votes_CVR[2] + reported_votes_noCVR[2])) / sum(reported_votes_CVR + reported_votes_noCVR)
#plurarlity contest assorter bound
u <- 1

#define tuning parameters, statistics, and assorters for CVR stratum
assorter_mean_CVR <- mean(rep(c(1,0,1/2), reported_votes_CVR)) 
v <- (reported_votes_CVR[1] - reported_votes_CVR[2]) / sum(reported_votes_CVR)
min_SD <- .05
n_CVR <- 8
#the value for a CVR assorter when omega_i = 0 (CVR tally matches MVR)
#the assorter max is 1: 1 - omega_i = 1
B_correct <- u
CVR_samples <- rep(B_correct, n_CVR)
d_1 <- 2
f_1 <- .1
epsilon_1 <- 1/(2*N[1])
eta_10 <- B_correct
u_CVR <- 2 * u

#define samples and tuning parameters for no-CVR stratum; there is some randomness in the no-CVR samples, only summaries were given
#may need to change path on other machines
sampled_noCVR <- read_csv("../mirla18/log/Kalamazoo-sampled-election-day.csv", skip = 1)
n_noCVR <- 32
eta_20 <- mean(rep(c(1,0,1/2), reported_votes_noCVR))
d_2 <- 5
f_2 <- 0
epsilon_2 <- 1/(2*N[2])

#define a function that computes stratified alpha P-value given a value of mu_01 for optimization
#this function mostly scopes its arguments to the global environment
pval_kalamazoo <- function(mu_01, noCVR_samples, martingale = "alpha", combine = "intersection"){
  raw_mu_tilde_01 <- u + mu_01 - assorter_mean_CVR
  raw_mu_02 <- (0.5 - w[1]*mu_01)/w[2]
  
  #adjust for sampling w/o replacement
  mu_tilde_01 <- (N[1] * raw_mu_tilde_01 - cumsum(CVR_samples)) / (N[1] - 1:length(CVR_samples) + 1)
  mu_02 <- (N[2] * raw_mu_02 - cumsum(noCVR_samples)) / (N[2] - 1:length(noCVR_samples) + 1)
  
  if(martingale == "alpha"){
    #compute alpha martingale for CVR stratum
    sigma_1 <- pmax(c(1/2, lag(get_sigma_hat_squared(CVR_samples))[2:length(CVR_samples)]), min_SD)
    eta_1i <- pmin(pmax(((cumsum(CVR_samples) + d_1*eta_10) / (d_1 + 1:n_CVR - 1) + f_1 * u_CVR/sigma_1) / (1 + f_1/sigma_1), mu_tilde_01 + epsilon_1), u_CVR*(1-.Machine$double.eps))
    terms_1 <- CVR_samples / mu_tilde_01 * (eta_1i - mu_tilde_01)/(u_CVR - mu_tilde_01) + (u_CVR - eta_1i)/(u_CVR - mu_tilde_01)
    
    #compute alpha martingale for no CVR stratum
    sigma_2 <- pmax(c(1/2, lag(get_sigma_hat_squared(noCVR_samples))[2:length(noCVR_samples)]), min_SD)
    eta_2i <- pmin(pmax(((cumsum(noCVR_samples) + d_2*eta_20) / (d_2 + 1:n_noCVR - 1) + f_2 * u/sigma_2) / (1 + f_2/sigma_2), mu_02 + epsilon_2), u*(1-.Machine$double.eps)) 
    terms_2 <- noCVR_samples / mu_02 * (eta_2i - mu_02)/(u - mu_02) + (u - eta_2i)/(u - mu_02)
  } else if(martingale == "eb"){
    #rescaling to [0,1]
    rescaled_CVR_samples <- CVR_samples / u_CVR
    rescaled_noCVR_samples <- noCVR_samples / u
    rescaled_mu_tilde_01 <- mu_01 / u_CVR
    rescaled_mu_02 <- mu_02 / u
    
    lambda_1 <- pmin(sqrt(log(2/.05) / (c(1/2, lag(get_sigma_hat_squared(rescaled_CVR_samples))[2:n_CVR]) * n_CVR)), 0.75)
    lambda_2 <- pmin(sqrt(log(2/.05) / (c(1/2, lag(get_sigma_hat_squared(rescaled_noCVR_samples))[2:n_noCVR]) * n_noCVR)), 0.75)
    v_1 <- get_empirical_Vn(rescaled_CVR_samples)
    v_2 <- get_empirical_Vn(rescaled_noCVR_samples)
    psi_E_1 <- get_psi_E(lambda_1)
    psi_E_2 <- get_psi_E(lambda_2)
    
    terms_1 <- exp(lambda_1 * (rescaled_CVR_samples - rescaled_mu_tilde_01) - v_1 * psi_E_1)
    terms_2 <- exp(lambda_2 * (rescaled_noCVR_samples - rescaled_mu_02) - v_2 * psi_E_2)
  } else{
    stop("input a valid method")
  }
  
  mart_1 <- prod(terms_1) 
  mart_2 <- prod(terms_2)
  
  if(combine == "intersection"){
    intersection_mart <- mart_1 * mart_2
    p_value <- min(1, 1/intersection_mart)
  } else if(combine == "fisher"){
    combined <- -2 * (log(min(1, 1/mart_1)) + log(min(1, 1/mart_2)))
    p_value <- 1 - pchisq(combined, df = 4)
  }
  p_value
}

#get the P-value for a particular order of CVR ballots
get_p_values <- function(){
  #there is randomness in the ballot-polling stratum 
  noCVR_samples <- shuffle(rep(c(1,0,1/2), c(23, 8, 1)))
  #optimize p values
  global_p_val_alpha_intersection <- optimize(
    function(x){pval_kalamazoo(mu_0 = x, noCVR_samples = noCVR_samples, martingale = "alpha", combine = "intersection")}, 
    interval = c(assorter_mean_CVR - 1, assorter_mean_CVR + 1),
    maximum = T)$objective
  global_p_val_eb_intersection <- optimize(
    function(x){pval_kalamazoo(mu_0 = x, noCVR_samples = noCVR_samples, martingale = "eb", combine = "intersection")}, 
    interval = c(assorter_mean_CVR - 1, assorter_mean_CVR + 1), 
    maximum = T)$objective
  global_p_val_alpha_fisher <- optimize(
    function(x){pval_kalamazoo(mu_0 = x, noCVR_samples = noCVR_samples, martingale = "alpha", combine = "fisher")}, 
    interval = c(assorter_mean_CVR - 1, assorter_mean_CVR + 1), 
    maximum = T)$objective
  global_p_val_eb_fisher <- optimize(
    function(x){pval_kalamazoo(mu_0 = x, noCVR_samples = noCVR_samples, martingale = "eb", combine = "fisher")}, 
    interval = c(assorter_mean_CVR - 1, assorter_mean_CVR + 1), 
    maximum = T)$objective
  c("alpha_intersection" = global_p_val_alpha_intersection,
             "eb_intersection" = global_p_val_eb_intersection,
             "alpha_fisher" = global_p_val_alpha_fisher,
             "eb_fisher" = global_p_val_eb_fisher)
}

p_value_distributions <- replicate(10000, get_p_values()) %>%  t()





