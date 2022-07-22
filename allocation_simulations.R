source("stratified_functions.R")

########## ALPHA allocation rules ############
#helper function to get sample sizes from array
get_total_sample_size <- function(x, alpha){
  if(any(x[,1] < alpha)){
    sum(x[1:min(which(x[,1] < alpha)),2:3])
  } else{
    NA
  }
}

run_allocation_simulation <- function(reported_tally, hand_tally, strata, n_sims = 300, alpha = .05){
  reported_mean <- tapply(reported_tally, strata, mean)
  assorter_bound <- 1
  omegas <- reported_tally - hand_tally
  # u <- (1 + 1/assorter_bound) / (2 - v/assorter_bound)
  # stratum_1 <- (1 - omegas[strata == 1] / assorter_bound) / (2 - v[1] / assorter_bound)
  # stratum_2 <- (1 - omegas[strata == 2] / assorter_bound) / (2 - v[2] / assorter_bound)
  stratum_1 <- assorter_bound - omegas[strata == 1]
  stratum_2 <- assorter_bound - omegas[strata == 2]
  u <- rep(2 * assorter_bound, 2)

  reported_margin <- 2*mean(reported_tally) - 1
  true_margin <- 2*mean(hand_tally) - 1
  understatements <- mean(c(stratum_1, stratum_2) == 1)
  overstatements <- mean(c(stratum_1, stratum_2) == 0)
  d <- c(20, 20)
  eta_0 <- u
  
  results_alpha_equal_product <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, d = d, eta_0 = eta_0, A_c = reported_mean, u = u, replace = FALSE, combine = "product", rule = "equal"))
  results_alpha_ucb_product <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, d = d, eta_0 = eta_0, A_c = reported_mean, u = u, replace = FALSE, combine = "product", rule = "alpha_ucb"))
  results_alpha_equal_fisher <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, d = d, eta_0 = eta_0, A_c = reported_mean, u = u, replace = FALSE, combine = "fisher", rule = "equal"))
  results_alpha_ucb_fisher <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, d = d, eta_0 = eta_0, A_c = reported_mean, u = u, replace = FALSE, combine = "fisher", rule = "alpha_ucb"))
  
  results_alpha_f01_equal_product <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, d = d, f = c(.01,.01), eta_0 = eta_0, A_c = reported_mean, u = u, replace = FALSE, combine = "product", rule = "equal"))  
  results_alpha_f01_ucb_product <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, d = d, f = c(.01,.01), eta_0 = eta_0, A_c = reported_mean, u = u, replace = FALSE, combine = "product", rule = "alpha_ucb"))
  results_alpha_f01_equal_fisher <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, d = d, f = c(.01,.01), eta_0 = eta_0, A_c = reported_mean, u = u, replace = FALSE, combine = "fisher", rule = "equal"))
  results_alpha_f01_ucb_fisher <- replicate(n_sims, get_two_strata_alpha(stratum_1, stratum_2, mu_0 = 0.5, d = d, f = c(.01,.01), eta_0 = eta_0, A_c = reported_mean, u = u, replace = FALSE, combine = "fisher", rule = "alpha_ucb"))
  
  results_eb_equal_product <- replicate(n_sims, get_two_strata_EB(stratum_1, stratum_2, mu_0 = 0.5, A_c = reported_mean, u = u, replace = FALSE, combine = "product", rule = "equal"))
  results_eb_ucb_product <- replicate(n_sims, get_two_strata_EB(stratum_1, stratum_2, mu_0 = 0.5, A_c = reported_mean, u = u, replace = FALSE, combine = "product", rule = "eb_ucb"))
  results_eb_equal_fisher <- replicate(n_sims, get_two_strata_EB(stratum_1, stratum_2, mu_0 = 0.5, A_c = reported_mean, u = u, replace = FALSE, combine = "fisher", rule = "equal"))
  results_eb_ucb_fisher <- replicate(n_sims, get_two_strata_EB(stratum_1, stratum_2, mu_0 = 0.5, A_c = reported_mean, u = u, replace = FALSE, combine = "fisher", rule = "eb_ucb"))
  
  #record stopping times
  stopping_times_alpha_equal_product <- apply(results_alpha_equal_product, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_alpha_ucb_product <- apply(results_alpha_ucb_product, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_alpha_equal_fisher <- apply(results_alpha_equal_fisher, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_alpha_ucb_fisher <- apply(results_alpha_ucb_fisher, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_alpha_f01_equal_product <- apply(results_alpha_f01_equal_product, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_alpha_f01_ucb_product <- apply(results_alpha_f01_ucb_product, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_alpha_f01_equal_fisher <- apply(results_alpha_f01_equal_fisher, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_alpha_f01_ucb_fisher <- apply(results_alpha_f01_ucb_fisher, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_eb_equal_product <- apply(results_eb_equal_product, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_eb_ucb_product <- apply(results_eb_ucb_product, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_eb_equal_fisher <- apply(results_eb_equal_fisher, 3, function(x){min(which(x[,1] < alpha))})
  stopping_times_eb_ucb_fisher <- apply(results_eb_ucb_fisher, 3, function(x){min(which(x[,1] < alpha))})
  
  #record total samples
  total_samples_alpha_equal_product <- apply(results_alpha_equal_product, 3, get_total_sample_size, alpha = alpha)
  total_samples_alpha_ucb_product <- apply(results_alpha_ucb_product, 3, get_total_sample_size, alpha = alpha)
  total_samples_alpha_equal_fisher <- apply(results_alpha_equal_fisher, 3, get_total_sample_size, alpha = alpha)
  total_samples_alpha_ucb_fisher <- apply(results_alpha_ucb_fisher, 3, get_total_sample_size, alpha = alpha)
  total_samples_alpha_f01_equal_product <- apply(results_alpha_f01_equal_product, 3, get_total_sample_size, alpha = alpha)
  total_samples_alpha_f01_ucb_product <- apply(results_alpha_f01_ucb_product, 3, get_total_sample_size, alpha = alpha)
  total_samples_alpha_f01_equal_fisher <- apply(results_alpha_f01_equal_fisher, 3, get_total_sample_size, alpha = alpha)
  total_samples_alpha_f01_ucb_fisher <- apply(results_alpha_f01_ucb_fisher, 3, get_total_sample_size, alpha = alpha)
  total_samples_eb_equal_product <- apply(results_eb_equal_product, 3, get_total_sample_size, alpha = alpha)
  total_samples_eb_ucb_product <- apply(results_eb_ucb_product, 3, get_total_sample_size, alpha = alpha)
  total_samples_eb_equal_fisher <- apply(results_eb_equal_fisher, 3, get_total_sample_size, alpha = alpha)
  total_samples_eb_ucb_fisher <- apply(results_eb_ucb_fisher, 3, get_total_sample_size, alpha = alpha)
  
  #gather into dataframe and return
  power_frame <- data.frame("stop" = stopping_times_alpha_equal_product, "total_samples" = total_samples_alpha_equal_product, "simulation" = 1:n_sims, "allocation" = "equal", "combined" = "product", "martingale" = "alpha") %>%
    bind_rows(
      data.frame("stop" = stopping_times_alpha_ucb_product, "total_samples" = total_samples_alpha_ucb_product, "simulation" = 1:n_sims, "allocation" = "ucb", "combined" = "product", "martingale" = "alpha"),
      data.frame("stop" = stopping_times_alpha_equal_fisher, "total_samples" = total_samples_alpha_equal_fisher, "simulation" = 1:n_sims, "allocation" = "equal", "combined" = "fisher", "martingale" = "alpha"),
      data.frame("stop" = stopping_times_alpha_ucb_fisher, "total_samples" = total_samples_alpha_ucb_fisher, "simulation" = 1:n_sims, "allocation" = "ucb", "combined" = "fisher", "martingale" = "alpha"),
      data.frame("stop" = stopping_times_alpha_f01_equal_product, "total_samples" = total_samples_alpha_f01_equal_product, "simulation" = 1:n_sims, "allocation" = "equal", "combined" = "product", "martingale" = "alpha_f01"),
      data.frame("stop" = stopping_times_alpha_f01_ucb_product, "total_samples" = total_samples_alpha_f01_ucb_product, "simulation" = 1:n_sims, "allocation" = "ucb", "combined" = "product", "martingale" = "alpha_f01"),
      data.frame("stop" = stopping_times_alpha_f01_equal_fisher, "total_samples" = total_samples_alpha_f01_equal_fisher, "simulation" = 1:n_sims, "allocation" = "equal", "combined" = "fisher", "martingale" = "alpha_f01"),
      data.frame("stop" = stopping_times_alpha_f01_ucb_fisher, "total_samples" = total_samples_alpha_f01_ucb_fisher, "simulation" = 1:n_sims, "allocation" = "ucb", "combined" = "fisher", "martingale" = "alpha_f01"),
      data.frame("stop" = stopping_times_eb_equal_product, "total_samples" = total_samples_eb_equal_product, "simulation" = 1:n_sims, "allocation" = "equal", "combined" = "product", "martingale" = "eb"),
      data.frame("stop" = stopping_times_eb_ucb_product, "total_samples" = total_samples_eb_ucb_product, "simulation" = 1:n_sims, "allocation" = "ucb", "combined" = "product", "martingale" = "eb"),
      data.frame("stop" = stopping_times_eb_equal_fisher, "total_samples" = total_samples_eb_equal_fisher, "simulation" = 1:n_sims, "allocation" = "equal", "combined" = "fisher", "martingale" = "eb"),
      data.frame("stop" = stopping_times_eb_ucb_fisher, "total_samples" = total_samples_eb_ucb_fisher, "simulation" = 1:n_sims, "allocation" = "ucb", "combined" = "fisher", "martingale" = "eb")
    ) %>%
    mutate(
      risk_limit = alpha,
      overstatements = overstatements, 
      understatements = understatements, 
      reported_margin = round(reported_margin, 3), 
      true_margin = round(true_margin, 3))
  power_frame
}


#allocation simulations
reported_tallies <- rbind(
  c(rep(0, 500), rep(1,500), rep(0, 400), rep(1, 600)),
  c(rep(0, 500), rep(1,500), rep(0, 450), rep(1, 550)),
  c(rep(0, 500), rep(1,500), rep(0, 490), rep(1, 510))
)

#CVR is correct in first instance and understates the margin. CVR is exactly correct in second. CVR overstates in third and 4th while outcome is correct, CVR is wrong in 5th AND the reported outcome is wrong.
#all are 2 vote overstatements (vote for loser recorded as vote for winner)
hand_tallies <- rbind(
  c(rep(0, 500), rep(1,500), rep(0, 400), rep(1, 600)),
  c(rep(0, 500), rep(1,500), rep(0, 450), rep(1, 550)),
  c(rep(0, 500), rep(1,500), rep(0, 490), rep(1, 510)),
  c(rep(0, 500), rep(1,500), rep(0, 500), rep(1, 500))
)

risk_limits <- c(.01, .05, .1)
outer_frames <- list()
for(i in 1:length(risk_limits)){
  inner_frames <- list()
  for(j in 1:nrow(reported_tallies)){
    inner_inner_frames <- list()
    for(k in 1:length(risk_limits)){
      inner_inner_frames[[k]] <- run_allocation_simulation(
        hand_tally = hand_tallies[k,],
        reported_tally = reported_tallies[j,],
        strata = c(rep(1,1000), rep(2,1000)), 
        n_sims = 100,
        alpha = risk_limits[i]
      )
      print(paste("i = ", i, ", j = ",j, " k = ", k, sep=""))
    }
    inner_frames[[j]] <- inner_inner_frames %>% reduce(bind_rows)
  }
  outer_frames[[i]] <- inner_frames %>% reduce(bind_rows)
}


allocation_frame <- outer_frames %>% 
  reduce(bind_rows) %>%
  mutate(finite_stop = ifelse(is.infinite(stop), 1000, stop)) %>%
  mutate(unconditional_total_samples = ifelse(is.na(total_samples), 2000, total_samples)) %>%
  #mutate(reported_margin = paste("Reported Margin =", reported_margin), true_margin = paste("True Margin =", true_margin)) %>%
  #mutate(allocation = recode(allocation, equal = "Proportional", threshold = "Thresholded", ucb = "UCB")) %>%
  as_tibble()


save(allocation_frame, file = "risklimits_allocation_power_frame_n100")
# load("risklimits_allocation_power_frame")
# allocation_frame_toplot <- allocation_frame %>%
#   filter(allocation %in% c("ucb", "equal")) %>%
#   filter(martingale %in% c("eb","alpha","alpha_f01")) %>%
#   filter(risk_limit == 0.05) %>%
#   filter(reported_margin == 0.10) %>%
#   filter(true_margin %in% c(.05, .01)) %>%
#   mutate(allocation = recode(allocation, equal = "Proportional", ucb = "Lower-sided test")) %>%
#   mutate(combined = ifelse(combined == "fisher", "Fisher", "Intersection")) %>%
#   mutate(martingale = recode(martingale, eb = "Empirical Bernstein", alpha = "ALPHA-ST", alpha_f01 = "ALPHA-UB")) %>%
#   mutate(martingale = factor(martingale, levels = c("ALPHA-ST", "ALPHA-UB", "Empirical Bernstein"))) %>%
#   mutate(true_margin_long = paste("True margin =", true_margin))
# 
# ggplot(allocation_frame_toplot, aes(x = unconditional_total_samples, linetype = combined, color = allocation)) +
#   stat_ecdf(size = 1.5, alpha = .75) +
#   scale_color_manual(values = c("darkorange3","steelblue","firebrick")) +
#   facet_grid(true_margin ~ martingale) +
#   xlim(0, 2000) +
#   theme_bw() +
#   theme(text = element_text(size = 18), axis.text = element_text(size = 14)) +
#   labs(x = "Total Samples", y = "Cumulative probability of stopping", color = "Allocation Rule", linetype = "Combination Rule")
# 
# 
# allocation_table <- allocation_frame %>%
#   filter(martingale %in% c("eb","alpha","alpha_f01")) %>%
#   filter(risk_limit == 0.05) %>%
#   group_by(
#     martingale,
#     allocation,
#     combined,
#     reported_margin,
#     true_margin) %>%
#   summarize(
#     expected_total_samples = mean(unconditional_total_samples),
#     percentile_total_samples = quantile(unconditional_total_samples, .9)) %>%
#   #mutate(total_samples = paste(round(expected_total_samples), " (", round(percentile_total_samples), ")", sep = "")) %>%
#   mutate(allocation = recode(allocation, equal = "Proportional", ucb = "Lower-sided test")) %>%
#   mutate(combined = ifelse(combined == "fisher", "Fisher", "Intersection")) %>%
#   mutate(martingale = recode(martingale, eb = "Empirical Bernstein", alpha = "ALPHA-ST", alpha_f01 = "ALPHA-UB")) %>%
#   pivot_longer(cols = c("expected_total_samples", "percentile_total_samples"), names_to = "total_samples") %>%
#   select(reported_margin, martingale, combined, allocation, true_margin, total_samples, value) %>%
#   mutate(value = as.character(round(value))) %>%
#   filter(true_margin != 0) %>%
#   pivot_wider(names_from = c("true_margin", "total_samples"), values_from = "value") %>%
#   arrange(reported_margin, martingale, combined, allocation)
# 
# print(xtable::xtable(allocation_table), include.rownames = FALSE)
# 
# 
# #ratios of workloads in every scenario and method, depending on the risk-limit
# risk_limit_ratios <- allocation_frame %>%
#   filter(martingale %in% c("eb","alpha","alpha_f01")) %>%
#   group_by(
#     martingale,
#     allocation,
#     combined,
#     risk_limit,
#     reported_margin,
#     true_margin) %>%
#   summarize(expected_total_samples = mean(unconditional_total_samples)) %>%
#   mutate(allocation = recode(allocation, equal = "Proportional", ucb = "Lower-sided test")) %>%
#   mutate(combined = ifelse(combined == "fisher", "Fisher", "Intersection")) %>%
#   mutate(martingale = recode(martingale, eb = "Empirical Bernstein", alpha = "ALPHA-ST", alpha_f01 = "ALPHA-UB")) %>%
#   pivot_longer(cols = c("expected_total_samples"), names_to = "total_samples") %>%
#   select(reported_margin, risk_limit, martingale, combined, allocation, true_margin, total_samples, value) %>%
#   filter(true_margin != 0) %>%
#   pivot_wider(names_from = c("risk_limit", "total_samples"), values_from = "value", names_prefix = "rl_") %>%
#   mutate(one_over_five = rl_0.01_expected_total_samples / rl_0.05_expected_total_samples, 
#          ten_over_five = rl_0.1_expected_total_samples / rl_0.05_expected_total_samples)
# 
# 
# #compute ratios of workloads in each scenario compared to the smallest, then take geometric means
# workload_ratios <- allocation_frame %>%
#   filter(risk_limit == 0.05) %>%
#   filter(martingale %in% c("eb","alpha","alpha_f01")) %>%
#   group_by(
#     martingale,
#     allocation,
#     combined,
#     reported_margin,
#     true_margin) %>%
#   summarize(expected_total_samples = mean(unconditional_total_samples)) %>%
#   ungroup() %>%
#   group_by(reported_margin, true_margin) %>%
#   mutate(workload_ratio = expected_total_samples / min(expected_total_samples)) %>%
#   group_by(martingale, allocation, combined) %>%
#   summarize(geo_mean = exp(mean(log(workload_ratio)))) %>%
#   mutate(allocation = recode(allocation, equal = "Proportional", ucb = "Lower-sided test")) %>%
#   mutate(combined = ifelse(combined == "fisher", "Fisher", "Intersection")) %>%
#   mutate(martingale = recode(martingale, eb = "Empirical Bernstein", alpha = "ALPHA-ST", alpha_f01 = "ALPHA-UB")) %>%
#   select(martingale, combined, allocation, geo_mean) %>%
#   arrange(martingale, combined, allocation, geo_mean)
# print(xtable::xtable(workload_ratios), include.rownames = FALSE)


