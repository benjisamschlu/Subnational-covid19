


#######################################################################################
##      PROJECT: Sub-national estimation of covid impact in BEL at district level
##      -- by Benjamin-Samuel Schlüter --
##      UCLouvain
##      09/04/2021
#######################################################################################
#
# CODE AIM: Project life expectancy at birth at the district level
#           while accounting for various sources of uncertainty.
#
#######################################################################################
#
#
# Notes:
# 1) 



############################################################################################################################################################################################


rm(list = ls())



#######################################################################################################################################################################################################################
# ----- Load pkg ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######################################################################################################################################################################################################################


packages <- c("tidyverse", "ggplot2", "rstan", "tidybayes", "ggpubr", "viridis", "scales", "janitor",
              "splines", "rstan", "tidybayes", "HMDHFDplus")
invisible( lapply(packages, library, character.only = TRUE))


# Functions to obtain life table values (mix between GC Camarda
# and M. Alexander code)
source("./code/functions/derive_lifetable_values_by5y.R")
source("./code/functions/derive_chiang_interval.R")


####################################################################################################################
# ----- Load & tidy data --------------------------------------------------------------------------------------------------
####################################################################################################################

# Death counts to model
# Set last age at 95 yo
d <- readRDS("./data/tidy/df_ageyearsexadmins_extrapol2020exp.rda")

# Total death year*district for model
df <- d %>% 
        filter(year %in% as.character(1991:2020)) %>% 
        mutate(age.gp = factor(age,
                               levels = c(0, 1, seq(5, 100, 5)),
                               labels = c(0, 1, seq(5, 95, 5), 95)),
               sex = factor(sex,
                            levels = c("f", "m"),
                            labels = c("Female", "Male")),
               year = as.character(year),
               year = as.numeric(year),
               dist_i = factor(dist),
               dist_i = as.numeric(dist_i)) %>% 
        group_by(year, dist, age.gp) %>% 
        # sum age group >95 yo
        summarise(dth = sum(dth),
                  pop = sum(pop),
                  exp = sum(exp)) %>% 
        ungroup() 


# Load 2021 pop data from statbel
# to avoid extrapolation
pop.2021 <- openxlsx::read.xlsx("./data/raw/TF_SOC_POP_STRUCT_2021.xlsx",
                                startRow = 1,
                                colNames = TRUE)
# Tidy data st it can be merged with df
df.pop.2021 <- pop.2021 %>% 
        rename("dist" = TX_ADM_DSTR_DESCR_FR,
               "age" = CD_AGE,
               "pop" = MS_POPULATION,
               "sex" = CD_SEX) %>% 
        mutate(age.gp = cut(age,
                            breaks = c(0, 1, seq(5, 95, 5), Inf),
                            labels = c(0, 1, seq(5, 95, 5)),
                            right = FALSE)) %>% 
        group_by(dist, age.gp) %>% 
        summarise(pop = sum(pop)) %>% 
        ungroup() %>% 
        # st dist name match btw sources
        mutate(dist = gsub("’", "'", dist))
# Check that dist names match
# all(unique(df.pop.2021$dist) %in% unique(df$dist))

# Extrapolation led to tot exp of 11,521,002 in 2020
sum(df$exp[df$year == 2020])

# Careful !! Following highly dependent on order
# Compute exposure
exp.2020 = round((df$pop[df$year == 2020] + df.pop.2021$pop)/2,0)

# Replace new exp in data set
df$exp[df$year == 2020] = exp.2020

df = df %>% 
        mutate(mx = dth/exp,
               # required for life table computation
               age = as.character(age.gp) %>% 
                       as.numeric())

# Obtain life expectancy for all district
df.e0 = df %>% 
        group_by(year, dist) %>% 
        mutate(ex = derive_ex_values(mx, age)) %>% 
        ungroup() 
# Check that e0 are fine
# ggplot(df.e0 %>% filter(age == 0), 
#        aes(x = ex, y = dist)) +
#         geom_point() +
#         facet_wrap(~year)

# Obtain life exepctancy for Bel
df.e0.Bel = df %>% 
        group_by(year, age) %>% 
        summarise(dth = sum(dth),
                  exp = sum(exp),
                  mx = dth/exp) %>% 
        ungroup() %>% 
        group_by(year) %>% 
        mutate(ex = derive_ex_values(mx, age))

# Join both
df.e0 = df.e0 %>% 
        left_join(df.e0.Bel %>% select(year, age, ex),
                  by = c("year", "age")) %>% 
        # Only look at life exp at birth
        filter(age == 0) %>% 
        # compute diff btw national and 
        # district e0
        mutate(diff_e = ex.x - ex.y)

# Compute Chiang confidence interval for 2020
df.e0.chiang = df %>% 
        filter(year == 2020) %>% 
        group_by(year, dist) %>% 
        mutate(ex = derive_ex_values(mx, age),
               # extract median of ex
               ex.m.chg = CIex(age, mx, dth)[2,],
               # extract 2.5 quantile
               ex.l.chg = CIex(age, mx, dth)[1,],
               # extract 97.5 quantile
               ex.u.chg = CIex(age, mx, dth)[3,]) %>% 
        ungroup()


####################################################################################################################
# ----- 2ND MODEL:REGIONAL DIFF TREND E0 --------------------------------------------------------------------------------------------------
####################################################################################################################


# Adaptation of Sevcikova & Raftery (2021)

# ----- Prior predictive checks ----------------------------------------------------------------------------------------------

require(truncnorm)

set.seed(2)
nsims = 1000

# Vague priors led to implausible values for diff in e0
# -> Weakly informative priors
# Generate draws for hyper parameters
mu_alpha = rnorm(nsims, 0, 1)
mu_beta = rnorm(nsims, 0, 1)
mu_sigma = rnorm(nsims, 0, 1)
sigma_alpha = rtruncnorm(nsims, a=0, mean = 0, sd = 1)
sigma_beta = rtruncnorm(nsims, a=0, mean = 0, sd = 1)
sigma_sigma = rtruncnorm(nsims, a=0, mean = 0, sd = 1)
# Generate draws for parameters
alpha = rnorm(nsims, mu_alpha, sigma_alpha)
beta = rnorm(nsims, mu_beta, sigma_beta)
sigma = rtruncnorm(nsims, a=0, mean = mu_sigma, sd = sigma_sigma)

# Create data set of sims
dsims = tibble(t = 0:28)
# Prior predictive check
for (i in 1:nsims){
        this_mu = alpha[i] + beta[i] * dsims$t
        dsims[paste0(i)] = this_mu + rnorm(nrow(dsims), 0, sigma[i])
}
# Change shape for ploting
dsl <- dsims %>% 
        pivot_longer(`1`:`1000`, names_to = "sim", values_to = "sim_diff")
# 95% values belongs to +- 49:
# the distribution clearly has mass around all
# plausible values for the difference in national
# vs subnational e0
quantile(dsl$sim_diff, probs = c(0.975, 0.5, 0.025))

# Plot prior predictive check
dsl %>% 
        ggplot(aes(sim_diff)) +
        geom_density() + 
        geom_segment(aes(col="95th quantiles"), x=-49.4, xend=-49.4 ,y=0, yend=0.0025, 
                     size = 2) +
        geom_segment(aes(col="95th quantiles"), x=49.6, xend=49.6 ,y=0, yend=0.0025, 
                     size = 2) +
        geom_text(label = -49.4,
                  x = -49.4,
                  y = 0.0035) +
        geom_text(label = 49.6,
                  x = 49.6,
                  y = 0.0035) +
        theme_bw() +
        theme(legend.position = c(0.2, 0.8),
              legend.title = element_blank(),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
              plot.margin = grid::unit(c(0,0,0,2), "mm")) +
        scale_color_manual(values = c("95th quantiles" = "#9C2964FF")) +
        labs(y = "Density",
             x = "Difference in e0, years")

# ggsave(filename = "./figures/prior_predictive_check_diff_e0.pdf",
#        device = "pdf")


# ----- STAN section -----------------------------------------------------------------------------------------------------------

# Length of dimensions
p = length(unique(df.e0.Bel$year[df.e0.Bel$year != 2020]))
d = length(unique(df.e0$dist))
t = df.e0.Bel %>% 
        filter(age == 0,
               year < 2020) %>% 
        mutate(t = year - 1991) %>% 
        pull(t)
# store diff nat vs subnat e0 in array
D_e = matrix(df.e0 %>%
                     filter(year < 2020) %>% 
                     arrange(dist, year) %>% 
                     pull(diff_e),
             nrow = p,
             ncol = d)
# data for stan
stan_data = list(
        p = p,
        d = d,
        t = t,
        D_e = D_e)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)

# Weakly informative priors as in prior
# predictive check
# run MCMC
fit <- stan("./code/stan/diff_e0_subnat_trend.stan",
            data = stan_data,
            iter = 2000,
            seed = 240)


# ----- Model's outputs -------------------------------------------------------------------------------------------------------------------

# Extract estimates
pars = array(as.matrix(fit, c("alpha", "beta")),
             dim = c(4000, 43, 2),
             dimnames = list(paste0("sim", 1:4000), unique(df.e0$dist), c("alpha", "beta"))) 
pars_median = apply(pars, 2:3, quantile, probs = 0.5)
sigma = as.matrix(fit, c("sigma"))
sigma_median = apply(sigma, 2, quantile, probs = 0.5)

# National posterior draws of é0 in 2020
# estimated previously
e0_lc = as.numeric( readRDS("./data/estimates/e0_2020_draws_LC.rda") )

# Combine national posterior draws of é0 in 2020
# with posterior predictive draws for the diff 
# nat vs subnat e0 during 2020, had covid19 not happened
# --> account for uncertainty in nation Lee-Carter projection
# + uncertainty in parameters estimation for the diff nat-subnat e0
# + data variation reflecting the size of the district
proj_e0_2020 = sapply(1:43, function(d) {
        
        e0_lc + ( pars[ , d, 1] + pars[ , d, 2] * (2020-1991) )  + rnorm(4000, mean = 0, sd = sigma[ , d])
})



# ----- Diff e0 posterior ranking distribution ------------------------------------------------------------------------------------------------

# Generate all possible realisations of e0
# with Chiang's method
obs_e0_chiang = sapply(unique(df.e0$dist), function(d){
        
        temp.df = df %>% 
                filter(year == 2020,
                       dist == d)
        out = with(temp.df, CIex(age, mx, dth, ns = 4000, all.draws = TRUE))
        
})
# For each realisation and posterior draws,
# compute the difference
diff_e0_draws = obs_e0_chiang - proj_e0_2020


# posterior ranking distribution
rank_mat <- matrix(0, ncol = 43, nrow = 43)
dist_id <- 1:43
nsim = 4000
# Algorithm to obtain posterior rank probabilities
for (k in 1:nsim){
        # cbind dist_id and theta
        temp <- cbind(dist_id, diff_e0_draws[k, ])
        # order dist_id by diff in e0 and store dist_id in rank order
        rank <- temp[order(temp[,2], decreasing = FALSE), 1]
        # Add +1 to the correct position
        for (j in 1:43){
                # add 1 for the rank (col) of each district (row)
                rank_mat[rank[j], j] <- rank_mat[rank[j], j] + 1      
        }
}
# Check: all(rowSums(rank_mat) == nsim)

# Convert rank_mat into rank_prob_mat
rank_prob_mat <- rank_mat/nsim
# Obtain expected rank for each district
exp_rank <- as.numeric(rank_prob_mat %*% c(1:43))
df_exp_rank <- tibble(dist_i = 1:43,
                      exp_rank = exp_rank,
                      dist = unique(df.e0$dist)) 
# Exctract order according to expected rank for plot
order_dist <- df_exp_rank %>% 
        mutate(dist = case_when(substr(dist, 16, 17) == "de" ~ substr(dist, 19, nchar(dist)),
                                substr(dist, 16, 17) == "d'" ~ substr(dist, 18, nchar(dist))) ) %>% 
        arrange(desc(exp_rank)) %>% 
        pull(dist)

# Convert rank_prob_mat into df
df_rank_prob <- rank_prob_mat %>% 
        as_tibble() %>% 
        mutate(dist_i = 1:43) %>% 
        pivot_longer(!dist_i, values_to = "prob", names_to = "rank") %>% 
        mutate(rank = substr(rank, 2, 3),
               # Order matter for plot
               rank = factor(rank,
                             levels = 1:43,
                             labels = 1:43),
               dist = unique(df.e0$dist)[dist_i]) %>% 
        mutate(dist = case_when(substr(dist, 16, 17) == "de" ~ substr(dist, 19, nchar(dist)),
                                substr(dist, 16, 17) == "d'" ~ substr(dist, 18, nchar(dist))),
               dist = factor(dist,
                             levels = order_dist,
                             labels = order_dist))

# Plot posterior rank probability with heat map
ggplot(df_rank_prob, aes(x = rank, y = dist, fill= prob)) +
        geom_tile(col = "white") +
        scale_fill_gradient(low="white", high="blue") +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_text(face="bold", size=10),
              legend.position = c(0.9, 0.9),
              plot.margin = grid::unit(c(0.01,0.01,0.01,2), "mm")) +
        labs(x = "Rank",
             fill = "Probability")

# ggsave(filename = "./figures/posterior_rank_prob_diff_e0.pdf",
#        device = "pdf")



# ----- Posterior predictive checks ----------------------------------------------------------------------------------------------

# Obtain mu draws
mu_draws = lapply(0:28, function(t) {
        
        pars[, , "alpha"] + pars[, , "beta"] * t
})
# Obtain diff_e0 draws
diff_e0_rep = array(NA, 
                    dim = c(29, 43, 4000),
                    dimnames = list(1991:2019, unique(df.e0$dist), paste0("sim", 1:4000)))
set.seed(3)
for (t in 1:29){
        for(i in 1:4000){
                diff_e0_rep[t, ,i] = rnorm(43, 
                                           mean = mu_draws[[t]][i, ],
                                           sd = sigma[i, ])
        }
}
# store in data set
df.diff_e0_rep = diff_e0_rep %>% 
        as.data.frame.table() %>% 
        rename("year" = Var1,
               "dist" = Var2,
               "sim"= Var3,
               "diff_e" = Freq)
# Combine draws and realiSation to compare
# and performe posterior pred check
set.seed(4)
df.ppc = rbind(df.diff_e0_rep %>% 
                       # only select some simulation for easier visualization
                       filter(sim %in% paste0("sim", sample(1:4000, 11))) %>% 
                       # remove factor for rbind()
                       mutate(year = as.character(year),
                              dist = as.character(dist),
                              sim = as.character(sim)), 
               df.e0 %>% 
                       filter(year < 2020) %>% 
                       select(year, dist, diff_e) %>% 
                       mutate(sim = "observed")
               )
# Plot posterior predictive check
df.ppc %>% 
        mutate(year = as.numeric(year)) %>% 
        ggplot(aes(x = year, y = diff_e, col = ifelse(sim == "observed", "observed", "sim"))) +
        geom_point() +
        theme_bw() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle=45),
              axis.title.x = element_blank()) +
        facet_wrap(~ sim) +
        labs(y = "Diff. in life expectancy at birth, years")

# ggsave(filename = "./figures/posterior_predictive_check_diff_e0.pdf",
#        device = "pdf")


# ----- Distribution of residuals from model ----------------------------------------------------------------------------------------------

# Create df and create expected value
# to be filled out
df.res.multilevel = df.e0 %>% 
        filter(year < 2020) %>% 
        select(!c(dth, pop, exp, mx)) %>% 
        mutate(diff_e_hat = NA)
# fill with estimates
for (d in unique(df.e0$dist)) {
        
        diff.e.hat = pars_median[d, 1] + pars_median[d, 2] * (0:28)
        df.res.multilevel$diff_e_hat[df.res.multilevel$dist == d] = diff.e.hat
}
# Compute residuals
df.res.multilevel = df.res.multilevel %>% 
        mutate(res = diff_e - diff_e_hat) 
# Compute sd of residuals
sd.res.multi = sd(df.res.multilevel$res) 

# Plot residuals with +-2 sd
ggplot(df.res.multilevel, aes(x = year, y = res, col = dist, group = dist)) +
        geom_point() +
        geom_hline(yintercept = c(2*sd.res.multi, -2*sd.res.multi),
                   linetype = "dashed") +
        theme_bw() +
        theme(legend.position = "none",
              axis.title.x = element_blank(),
              panel.border = element_blank(),
              axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black")) +
        labs(y = "Residuals")



# ----- Assess change in e0 distributions over years ------------------------------------------------------------------------------------------------

# Generate all possible realisations of e0
# with Chiang's method and then, compute for each 
# year its associated standard deviation

# Container for sd(e0) iterations
sd_e0 = matrix(NA,
               nrow = 4000,
               ncol = 6)
for (y in 2015:2020) {
        
        obs_e0_chiang = sapply(unique(df.e0$dist), function(d){
                
                temp.df = df %>% 
                        filter(year == y,
                               dist == d)
                out = with(temp.df, CIex(age, mx, dth, ns = 4000, all.draws = TRUE))
                
        })
        sd_e0[, (y-2014)] <- apply(obs_e0_chiang, 1, sd)

}
# Quantiles for draws of sd for each year
apply(sd_e0, 2, quantile, probs = c(0.025, 0.5, 0.975))

# Plot
sd_e0 %>% 
        as.data.frame() %>% 
        pivot_longer(V1:V6, names_to = "year", values_to = "values") %>% 
        mutate(year = case_when(year == "V1" ~ "2015 -- 1.40 (1.30-1.51)",
                                year == "V2" ~ "2016 -- 1.56 (1.46-1.67)",
                                year == "V3" ~ "2017 -- 1.46 (1.35-1.58)",
                                year == "V4" ~ "2018 -- 1.49 (1.39-1.60)",
                                year == "V5" ~ "2019 -- 1.47 (1.37-1.58)",
                                year == "V6" ~ "2020 -- 1.72 (1.63-1.81)")) %>% 
        ggplot(aes(x = values, col = year, fill = year, group = year)) +
        geom_density(alpha=0.3) +
        scale_color_manual(values = c("2020 -- 1.72 (1.63-1.81)" = "#460B6AFF",
                                      "2019 -- 1.47 (1.37-1.58)" = "#238A8DFF",
                                      "2018 -- 1.49 (1.39-1.60)" = "#20A386FF",
                                      "2017 -- 1.46 (1.35-1.58)" = "#3CBC75FF",
                                      "2016 -- 1.56 (1.46-1.67)" = "#74D055FF",
                                      "2015 -- 1.40 (1.30-1.51)" = "#B8DE29FF")) +
        scale_fill_manual(values = c("2020 -- 1.72 (1.63-1.81)" = "#460B6AFF",
                                     "2019 -- 1.47 (1.37-1.58)" = "#238A8DFF",
                                     "2018 -- 1.49 (1.39-1.60)" = "#20A386FF",
                                     "2017 -- 1.46 (1.35-1.58)" = "#3CBC75FF",
                                     "2016 -- 1.56 (1.46-1.67)" = "#74D055FF",
                                     "2015 -- 1.40 (1.30-1.51)" = "#B8DE29FF")) +
        theme_bw() +
        theme(legend.position = c(0.2, 0.65),
              legend.title = element_blank(),
              panel.grid = element_blank(),
              plot.margin = grid::unit(c(0,0,0,2), "mm")) +
        xlim(c(1, 1.9)) +
        labs(y = expression(f(sigma[e^0])),
             x = expression(sigma[e^0])) +
        annotate("text", x = 1.17, y = 7, label = "Year -- Median (95% C.I)", 
                 size = 4, hjust = .5)

# ggsave(filename = "./figures/density_sigma_e0.pdf",
#        device = "pdf")


# ----- Draws for differences between observed and projected -------------------------------------------------------------------------------------------------------------------------

# Suppose no life table uncertainty (= no Chiang draws)
diff_e0_draws_no_obs_uncert <- matrix(rep(df.e0.chiang$ex[df.e0.chiang$age == 0], each = 4000), 
                                      ncol = 43, 
                                      byrow = FALSE) - proj_e0_2020

diff_e0_draws_no_obs_uncert <-  apply(diff_e0_draws_no_obs_uncert, 2, quantile, probs = c(0.025, 0.1, 0.5, 0.9, 0.975)) %>% 
        as.data.frame.table() %>% 
        pivot_wider(names_from = Var1, values_from = Freq) %>% 
        rename("dist" = Var2,
               l.95 = "2.5%",
               l.80 = "10%",
               median = "50%",
               u.80 = "90%",
               u.95 = "97.5%") %>% 
        mutate(dist = colnames(obs_e0_chiang),
               dist = case_when(substr(dist, 16, 17) == "de" ~ substr(dist, 19, nchar(dist)),
                                substr(dist, 16, 17) == "d'" ~ substr(dist, 18, nchar(dist))),
               dist = factor(dist) )



apply(diff_e0_draws, 2, quantile, probs = c(0.025, 0.1, 0.5, 0.9, 0.975)) %>% 
        as.data.frame.table() %>% 
        pivot_wider(names_from = Var1, values_from = Freq) %>% 
        rename("dist" = Var2,
               l.95 = "2.5%",
               l.80 = "10%",
               median = "50%",
               u.80 = "90%",
               u.95 = "97.5%") %>% 
        mutate(dist = as.character(dist),
               dist = case_when(substr(dist, 16, 17) == "de" ~ substr(dist, 19, nchar(dist)),
                                substr(dist, 16, 17) == "d'" ~ substr(dist, 18, nchar(dist))),
               dist = factor(dist) ) %>% 
        ggplot(aes(x = median, y = fct_reorder(dist, median, .desc = TRUE) )) +
        geom_linerange(aes(xmin = l.95, xmax = u.95, col = "Projection and LT uncertainties 95%"), 
                           size = 1) +
        geom_linerange(aes(xmin = l.80, xmax = u.80, col = "Projection and LT uncertainties 80%"),
                           size = 1) +
        geom_pointinterval(data = diff_e0_draws_no_obs_uncert,
                           aes(xmin = l.80, xmax = u.80, col = "Projection uncertainty 80%"),
                           size = 2,
                           point_size = 2) +
        theme_bw() +
        geom_vline(xintercept = 0, 
                   linetype = "dashed") +
        ggforestplot::geom_stripes(odd = "#33333333", even = "#00000000") +
        theme_bw() +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_text(face = "bold",
                                         size = 12),
              legend.position = c(0.855, 0.85),
              legend.text = element_text(size = 9.5),
              panel.border = element_blank(),
              panel.grid.major.y = element_blank(),
              axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
              axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
              plot.margin = grid::unit(c(0.01,0.01,0.01,2), "mm")) +
        scale_color_manual("Uncertainty scenario",
                           values = c("Projection uncertainty 80%" = "#C03A76FF",
                                      "Projection and LT uncertainties 80%" = "#3F4788FF",
                                      "Projection and LT uncertainties 95%" = "#FBB91FFF")) +
        labs(x = "Loss of life expectancy, years")

# ggsave(filename = "./figures/lle0_80_95_alt.pdf",
#        device = "pdf")


# Save data set for maps
lle0_maps <- apply(diff_e0_draws, 2, quantile, probs = c(0.025, 0.1, 0.5, 0.9, 0.975)) %>% 
        as.data.frame.table() %>% 
        pivot_wider(names_from = Var1, values_from = Freq) %>% 
        rename("dist" = Var2,
               l.95 = "2.5%",
               l.80 = "10%",
               median = "50%",
               u.80 = "90%",
               u.95 = "97.5%") %>% 
        mutate(dist = as.character(dist),
               width_u = u.95-l.95) %>% 
        select(dist, median, width_u)

# saveRDS(lle0_maps,
#         "./data/estimates/lle0_for_maps.rda")

# lle0_maps <- readRDS("./data/estimates/lle0_for_maps.rda")



# ----- Info for manuscript text -------------------------------------------------------------------------------------

required_nbr <- apply(diff_e0_draws, 2, quantile, probs = c(0.025, 0.1, 0.5, 0.9, 0.975)) %>% 
        as.data.frame.table() %>% 
        pivot_wider(names_from = Var1, values_from = Freq) %>% 
        rename("dist" = Var2,
               l.95 = "2.5%",
               l.80 = "10%",
               median = "50%",
               u.80 = "90%",
               u.95 = "97.5%") %>% 
        mutate(dist = as.character(dist),
               dist = case_when(substr(dist, 16, 17) == "de" ~ substr(dist, 19, nchar(dist)),
                                substr(dist, 16, 17) == "d'" ~ substr(dist, 18, nchar(dist))),
               dist = factor(dist) )

required_nbr %>% filter(dist %in% c("Bruxelles-Capitale", "Arlon", "Mons"))
