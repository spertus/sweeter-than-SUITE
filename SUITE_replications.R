source("stratified_functions.R")

#for michigan replication
load("~/Dropbox/RLAs/stratified-inference/elec.strat/data/MN_Senate_2006.rda")

############# CORLA replications ################
#replicating simulations in CORLA18: https://github.com/pbstark/CORLA18/blob/master/code/hybrid-audit-example-1.ipynb
##### example 1 ######
#upper bound on assorters for plurality contests
u <- 1

#example 1, medium sized election, close race
#define stratum 1 (CVR) population 
stratum_1_reported <- c(rep(1, 45500), rep(0, 49500), rep(1/2, 5000))
stratum_1_hand <- stratum_1_reported
stratum_1_v <- 2 * mean(stratum_1_reported) - 1
stratum_1_range <- c(0, (1 + 1/u) / (2-stratum_1_v/u))
stratum_1_omegas <- stratum_1_reported - stratum_1_hand 
stratum_1_population <- (1 - stratum_1_omegas/u) / (2 - stratum_1_v/u)

#define stratum 2 (no CVR) population
stratum_2_population <- c(rep(1, 7500), rep(0, 1500), rep(1/2, 1000))

#combine 
population <- c(stratum_1_population, stratum_2_population)
strata <- c(rep(1, length(stratum_1_population)), rep(2, length(stratum_2_population)))
sample_sizes <- matrix(c(700, 500), nrow = 1, ncol = 2)

eta_0 <- c(1 / (2 - stratum_1_v / stratum_1_range[2]), mean(stratum_2_population))


power_frame <- run_stratified_simulation(
  population = population, 
  strata = strata, 
  sample_sizes = sample_sizes, 
  n_sims = 1e4, 
  alpha = 0.1, 
  pars = list(
    d = 20,
    eta_0 = eta_0),
  bounds = c(0, max(stratum_1_range[2], 1))
  )
save(power_frame, file = "CORLA_replication_frame_1")

##### example 2 ######

#upper bound on assorters for plurality contests
u <- 1

#example 1, medium sized election, close race
#define stratum 1 (CVR) population 
stratum_1_reported <- c(rep(1, 1.102e6), rep(0, 7.03e5), rep(1/2, 9.5e4))
stratum_1_hand <- stratum_1_reported
stratum_1_v <- 2 * mean(stratum_1_reported) - 1
stratum_1_range <- c(0, (1 + 1/u) / (2-v/u))
stratum_1_omegas <- stratum_1_reported - stratum_1_hand 
stratum_1_population <- (1 - stratum_1_omegas/u) / (2 - stratum_1_v/u)

#define stratum 2 (no CVR) population
stratum_2_population <- c(rep(1, 4.25e4), rep(0, 5.25e4), rep(1/2, 5e3))

#combine 
population <- c(stratum_1_population, stratum_2_population)
strata <- c(rep(1, length(stratum_1_population)), rep(2, length(stratum_2_population)))
sample_sizes <- matrix(c(700, 500), nrow = 1, ncol = 2)

eta_0 <- c(1 / (2 - stratum_1_v / stratum_1_range[2]), mean(stratum_2_population))


power_frame <- run_stratified_simulation(
  population = population, 
  strata = strata, 
  sample_sizes = sample_sizes, 
  n_sims = 1e3, 
  alpha = .1,
  pars = list(
    d = 20,
    eta_0 = eta_0),
  bounds = c(0, max(stratum_1_range[2], 1))
  )

save(power_frame, file = "CORLA_replication_frame_2")


################## MIRLA Replications #########################

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


#the SHANGRLA assorter value of a correct reported vote in the CVR stratum
#note that all were correct in MIRLA CVR sample
assorter_max <- 1
v <- (reported_votes_CVR[1] - reported_votes_CVR[2]) / sum(reported_votes_CVR)
B_correct <- 1 / (2 - v/assorter_max)

n_CVR <- 8
CVR_samples <- rep(B_correct, n_CVR)
d_1 <- 5
epsilon_1 <- 1/(2*N[1])
eta_10 <- B_correct
u_CVR <- 2 / (2 - v/assorter_max)

u_noCVR <- 1
n_noCVR <- 32
noCVR_samples <- shuffle(rep(c(1,0,1/2), c(23, 8, 1)))
eta_20 <- mean(rep(c(1,0,1/2), reported_votes_noCVR))
d_2 <- 10
epsilon_2 <- 1/(2*N[2])

#define a function to compute stratified alpha P-value given a value of mu_01
pval_kalamazoo <- function(mu_01, martingale = "alpha", combine = "intersection"){
  mu_02 <- (0.5 - w[1]*mu_01)/w[2]
  
  #compute alpha martingale for CVR stratum
  if(martingale == "alpha"){
    eta_1i <- pmin(pmax((cumsum(CVR_samples) + d_1*eta_10) / (d_1 + 1:n_CVR - 1), mu_01 + epsilon_1), u_CVR*(1-.Machine$double.eps))
    terms_1 <- CVR_samples / mu_01 * (eta_1i - mu_01)/(u_CVR - mu_01) + (u_CVR - eta_1i)/(u_CVR - mu_01)
    #compute alpha martingale for no CVR stratum
    eta_2i <- pmin(pmax((cumsum(noCVR_samples) + d_2*eta_20) / (d_2 + 1:n_noCVR - 1), mu_02 + epsilon_2), u_noCVR*(1-.Machine$double.eps)) 
    terms_2 <- noCVR_samples / mu_02 * (eta_2i - mu_02)/(u_noCVR - mu_02) + (u_noCVR - eta_2i)/(u_noCVR - mu_02)
  } else if(martingale == "eb"){
    #rescaling to [0,1]
    rescaled_CVR_samples <- CVR_samples / u_CVR
    rescaled_noCVR_samples <- noCVR_samples / u_noCVR
    rescaled_mu_01 <- mu_01 / u_CVR
    rescaled_mu_02 <- mu_02 / u_noCVR
    
    lambda_1 <- pmin(sqrt(log(2/.05) / (c(1/2, lag(get_sigma_hat_squared(rescaled_CVR_samples))[2:n_CVR]) * n_CVR)), 0.5)
    lambda_2 <- pmin(sqrt(log(2/.05) / (c(1/2, lag(get_sigma_hat_squared(rescaled_noCVR_samples))[2:n_noCVR]) * n_noCVR)), 0.5)
    v_1 <- get_empirical_Vn(rescaled_CVR_samples)
    v_2 <- get_empirical_Vn(rescaled_noCVR_samples)
    psi_E_1 <- get_psi_E(lambda_1)
    psi_E_2 <- get_psi_E(lambda_2)
    
    terms_1 <- exp(lambda_1 * (rescaled_CVR_samples - rescaled_mu_01) - v_1 * psi_E_1)
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

global_p_val_alpha_intersection <- optimize(function(x){pval_kalamazoo(mu_0 = x, martingale = "alpha", combine = "intersection")}, interval = c(0, min(u,0.5/w[1])), maximum = T)$objective
global_p_val_eb_intersection <- optimize(function(x){pval_kalamazoo(mu_0 = x, martingale = "eb", combine = "intersection")}, interval = c(0, min(u,0.5/w[1])), maximum = T)$objective
global_p_val_alpha_fisher <- optimize(function(x){pval_kalamazoo(mu_0 = x, martingale = "alpha", combine = "fisher")}, interval = c(0, min(u,0.5/w[1])), maximum = T)$objective
global_p_val_eb_fisher <- optimize(function(x){pval_kalamazoo(mu_0 = x, martingale = "eb", combine = "fisher")}, interval = c(0, min(u,0.5/w[1])), maximum = T)$objective








