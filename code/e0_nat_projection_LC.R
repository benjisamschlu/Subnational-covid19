


#######################################################################################
##      PROJECT: Sub-national estimation of covid impact in BEL at district level
##      -- by Benjamin-Samuel Schlüter --
##      UCLouvain
##      09/04/2021
#######################################################################################
#
# CODE AIM: Project national life expectancy at birth with L-C model
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

####################################################################################################################
# ----- Load & tidy data --------------------------------------------------------------------------------------------------
####################################################################################################################

# Death counts to model
# Set last age at 95 yo
d <- readRDS("./data/tidy/df_ageyearsexadmins_extrapol2020exp.rda")

# Total death year*district for model
df <- d %>% 
        filter(year %in% as.character(1991:2020)) %>% 
        mutate(sex = factor(sex,
                            levels = c("f", "m"),
                            labels = c("Female", "Male")),
               year = as.character(year),
               year = as.numeric(year),
               dist_i = factor(dist),
               dist_i = as.numeric(dist_i)) %>% 
        group_by(year, dist, age) %>% 
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
                            breaks = c(0, 1, seq(5, 100, 5), Inf),
                            labels = c(0, 1, seq(5, 100, 5)),
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

# Careful !! Following highly dependent on order of data set
# Compute exposure
exp.2020 = round((df$pop[df$year == 2020] + df.pop.2021$pop)/2,0)

# Replace new exp in data set
df$exp[df$year == 2020] = exp.2020

df = df %>% 
        mutate(mx = dth/exp,
               # required for life table computation
               age = as.character(age) %>% 
                       as.numeric())

####################################################################################################################
# ----- LEE-CARTER MODEL --------------------------------------------------------------------------------------------------
####################################################################################################################


# Obtain log mx matrix
m_tx <- df %>% 
        # not using 2020 
        filter(year < 2020) %>% 
        group_by(year, age) %>% 
        summarise(dth = sum(dth),
                  exp = sum(exp)) %>% 
        ungroup() %>% 
        mutate(mx = dth/exp) %>% 
        select(-c(dth, exp)) %>% 
        pivot_wider(names_from = age, values_from = mx) %>% 
        select(-year) %>% 
        as.matrix()

logm_tx <- log(m_tx)
# Average mortality rates
ax <- apply(logm_tx, 2, mean)

# Variables
ages <- unique(df$age)
years <- seq(1991, 2019, 1)

# Perform SVD

# Subtract average rates
swept_logm_tx <- sweep(logm_tx, 2, ax)

svd_mx <- svd(swept_logm_tx)

bx <- svd_mx$v[, 1]/sum(svd_mx$v[, 1])
kt <- svd_mx$d[1] * svd_mx$u[, 1] * sum(svd_mx$v[, 1])

# Plot obtained values
lc_age_df <- tibble(age = ages, ax = ax, bx = bx)
lc_time_df <- tibble(year = years, kt = kt)

p1 <- ggplot(lc_age_df, aes(age, ax)) + 
        geom_line(lwd = 1.1) + 
        ggtitle("ax values")

p2 <- ggplot(lc_age_df, aes(age, bx)) + 
        geom_line(lwd = 1.1) + 
        ggtitle("bx values")

p3 <- ggplot(lc_time_df, aes(year, kt)) + 
        geom_line(lwd = 1.1) + 
        ggtitle("kt values")

# Assess the fit
df %>% 
        filter(year < 2020) %>% 
        group_by(year, age) %>% 
        summarise(dth = sum(dth),
                  exp = sum(exp)) %>% 
        ungroup() %>% 
        mutate(mx = dth/exp) %>% 
        select(-c(dth, exp)) %>%
        left_join(lc_age_df) %>% 
        left_join(lc_time_df) %>% 
        mutate(estimate = exp(ax+bx*kt)) %>% 
        select(year, age, mx, estimate) %>% 
        filter(year %in% seq(1991, 2019, by = 5)) %>% 
        ggplot(aes(age, mx)) + 
        geom_point(aes(color = "mx")) + 
        geom_line(aes(age, estimate, color = "fit"), lwd = 1.1)+
        facet_grid(~year) + 
        scale_y_log10() + 
        ylab("mortality rate") +
        scale_color_brewer(name = "",palette = "Set1") + 
        theme_bw()

# Create object for Stan
D <- df %>% 
        filter(year < 2020) %>% 
        group_by(year, age) %>% 
        summarise(dth = sum(dth)) %>% 
        ungroup() %>% 
        pivot_wider(names_from = age, values_from = dth) %>% 
        select(-year) %>% 
        as.matrix()

P <- df %>% 
        filter(year < 2020) %>% 
        group_by(year, age) %>% 
        summarise(exp = sum(exp)) %>% 
        ungroup() %>% 
        pivot_wider(names_from = age, values_from = exp) %>% 
        select(-year) %>% 
        as.matrix()

stan_data <- list(D = D, 
                  P = P, 
                  nyears = length(years), 
                  nages = length(ages),
                  nprojyears = 5,
                  ax = ax,
                  bx = bx,
                  kt = kt)
# Run MCMC
nchain = 4
niter = 4000

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)

# run Lee-Carter stan code obtained from Monica Alexander
# github account (assumes RW2 for k_t)
fit <- stan("./code/stan/lee_carter.stan",
            data = stan_data,
            chains = nchain,
            iter=niter, 
            seed = 123)


####################################################################################################################
# ----- E0 PROJECTION --------------------------------------------------------------------------------------------------
####################################################################################################################


# Extract mx posterior draws
proj_2020 = fit %>% 
        gather_draws(log_mu_proj[p,x]) %>% 
        filter(p == 1) %>%
        ungroup() %>% 
        select(-c(p, .chain, .iteration, .variable)) %>% 
        mutate(age = c(0, 1, seq(5, 100, 5))[x],
               mx = exp(.value)) %>% 
        select(-x)
# Compute life exptancy at all ages for each iteration
lt_2020 = proj_2020 %>% 
        group_by(.draw) %>% 
        mutate(ex = derive_ex_values(mx, age)) %>% 
        ungroup()
# extract median estimate and 95% quantiles
# from posterior e0 in 2020
lt_2020 %>% 
        group_by(age) %>% 
        summarise(e0_median = quantile(ex, probs = 0.5),
               e0_low = quantile(ex, probs = 0.025),
               e0_up = quantile(ex, probs = 0.975)) %>% 
        ungroup() %>% 
        filter(age == 0)

# Store e0 posterior draws in 2020
e0_2020_LC <- lt_2020 %>% 
        filter(age == 0) %>% 
        select(-c(.draw, .value, age, mx)) %>% 
        as.matrix()

# saveRDS(e0_2020_LC,
#         "./data/estimates/e0_2020_draws_LC.rda")

# National projection of e0 for 2020 according to Lee-Carter is 
# Median: 81.7 and 95% CI= 80.9-82.4



