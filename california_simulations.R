### running a highly stratified audit on hypothetical California data
library(tidyverse)
library(readxl)
source("stratified_functions.R")
#data is from here: https://www.sos.ca.gov/elections/prior-elections/statewide-election-results/general-election-november-3-2020/statement-vote
set.seed(1337)

california_results <- read_excel("csv-candidates.xlsx")
total_votes <- california_results %>%
  filter(`Contest Name` == "President") %>%
  group_by(`County Name`) %>%
  summarize(votes = sum(`Vote Total`))
biden_trump_votes <- california_results %>%
  filter(`Candidate Name` %in% c("Joseph R. Biden", "Donald J. Trump")) %>%
  select(`County Name`, `Candidate Name`, votes = `Vote Total`) %>%
  pivot_wider(names_from = `Candidate Name`, values_from = votes)

votes <- biden_trump_votes %>% 
  inner_join(total_votes, by = "County Name") %>%
  rename(county = `County Name`, biden = `Joseph R. Biden`, trump = `Donald J. Trump`) %>%
  mutate(invalids = votes - biden - trump)

reported_votes_list <- list()
for(i in 1:nrow(votes)){
  reported_votes_list[[i]] <- c(
    rep(1, votes$biden[i]),
    rep(0, votes$trump[i]),
    rep(1/2, votes$invalids[i])
  )
}

strata <- rep(1:58, votes$votes)
#start with ballot polling; and assume reported tallies are accurate
true_votes <- reported_votes_list %>% 
  reduce(c)

weights <- prop.table(table(strata))


#it takes 4.154 seconds to run this one time on all the data.
n_strata <- round_strata_sizes(weights * 40000) + 10
p_val <- get_stratified_pvalue(
  population = true_votes, 
  strata = strata, 
  n = n_strata, 
  mu_0 = 0.5,
  method = "empirical_bernstein", 
  pool = "fisher", 
  alpha = 0.05, 
  bounds = c(0,1))


#simulations
num_sims <- 300
n_grid <- seq(30000, 40000, by = 1000)
EB_pval <- matrix(NA, nrow = num_sims, ncol = length(n_grid))
for(i in 1:length(n_grid)){
  n_strata <- round_strata_sizes(weights * n_grid[i]) + 10
  EB_pval[,i] <- replicate(num_sims, get_stratified_pvalue(
    population = true_votes, 
    strata = strata, 
    n = n_strata, 
    mu_0 = 0.5,
    method = "empirical_bernstein", 
    pool = "fisher", 
    alpha = 0.05, 
    bounds = c(0,1))
  )
}

results_list <- list(
  n = n_grid,
  pvals = EB_pval
)
save(results_list, file = "california_audit_results")

