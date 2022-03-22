


#######################################################################################
##      PROJECT: Sub-national estimation of covid impact in BEL at district level
##      -- by Benjamin-Samuel Schlüter --
##      UCLouvain
##      09/04/2021
#######################################################################################
#
# CODE AIM: Use SMR to assess how heterogenous the shock has been within Belgium
#
#######################################################################################
#
#
# Notes:
# 1) Instead of left join everywhere, create district vector with same order as dist_i



############################################################################################################################################################################################


rm(list = ls())





#######################################################################################################################################################################################################################
# ----- Load pkg ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#######################################################################################################################################################################################################################


packages <- c("tidyverse", "ggplot2", "rstan", "tidybayes", "ggpubr", "viridis", "scales", "janitor",
              "splines", "rstan", "tidybayes")
invisible( lapply(packages, library, character.only = TRUE))




####################################################################################################################
# ----- Load & tidy data --------------------------------------------------------------------------------------------------
####################################################################################################################


# Data from civil register and death certificates #
# ----------------------------------------------- #

d = readRDS("./data/tidy/df_ageyearsexadmins_extrapol2020exp.rda")

# Create counts at the district level
df = d %>% 
        mutate(sex = factor(sex,
                            levels = c("f", "m"),
                            labels = c("Female", "Male")),
               year = as.character(year) %>% as.numeric(),
               dist_i = factor(dist),
               age = factor(age,
                               levels = c(0, 1, seq(5, 100, 5)),
                            # set 95 yo as max age to reduce noise
                               labels = c(0, 1, seq(5, 95, 5), 95))) %>% 
        filter(year %in% 2015:2020) %>% 
        group_by(age, sex, year, dist) %>% 
        summarise(dth = sum(dth),
                  pop = sum(pop),
                  exp = sum(exp)) %>% 
        ungroup()


# Load 2021 pop data from statbel
# to replace previous extrapolation
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
                            right = FALSE),
               sex = factor(sex,
                            levels = c("F", "M"),
                            labels = c("Female", "Male"))) %>% 
        group_by(age.gp, sex ,dist) %>% 
        summarise(pop = sum(pop)) %>% 
        ungroup() %>% 
        # st dist name match btw sources
        mutate(dist = gsub("’", "'", dist))
# Check that dist names match
# all(unique(df.pop.2021$dist) %in% unique(df$dist))

# Extrapolation led to tot exp of 11,521,002 in 2020
sum(df$exp[df$year == 2020])

# Careful !! Following code is highly dependent 
# on order
# Compute exposure
exp.2020 = round((df$pop[df$year == 2020] + df.pop.2021$pop)/2,0)

# Replace new exp in data set
df$exp[df$year == 2020] = exp.2020

# 15,000 difference with previous extrapolation
sum(df$exp[df$year == 2020])

# Dimensions of various variables
Nage = length(unique(df$age))
Ndist = length(unique(df$dist))
Nyear = length(unique(df$year))
# District names in order of df
Namedist = unique(df$dist) 

# Belgian national mortality rates
mx_Bel = df %>% 
        group_by(age, year) %>% 
        summarise(dth = sum(dth),
                  exp = sum(exp)) %>%
        ungroup() %>% 
        mutate(mx = dth/exp) %>% 
        arrange(year, age)
# Convert to matrix
Mx = matrix(mx_Bel$mx, c(Nage, Nyear))

# Check shape
# matplot(log(Mx), type = "l")

# Convert yearly deaths by district to matrix
D = df %>% 
        group_by(year, dist) %>% 
        summarise(dth = sum(dth)) %>% 
        ungroup() %>% 
        pivot_wider(values_from = dth, names_from = year) %>% 
        select(!dist) %>% 
        as.matrix()
# Check implausible values
# matplot(y= 1:43, x =D,
#           type = "p",
#           pch = 14:20)
                

# Compute yearly expectation by district
# according to belgian ASMR and district pop
exposures = df %>% 
        group_by(year, dist, age) %>% 
        summarise(exp = sum(exp)) %>% 
        ungroup() %>% 
        arrange(year, age, dist)
# store in array
N = array(exposures$exp, c(Ndist, Nage, Nyear))
# Check implausible values
# N_check = apply(N, c(1,3), sum)
# matplot(y = 1:43, x = N_check,
#         type = "p",
#         pch = 14:20)

# matrix multiplication to obtain district-level
# expected deaths, had the district experienced 
# national mortality#
EXP = sapply(1:Nyear, function(t){ N[,,t] %*% Mx[,t]})
# Check implausible values
# matplot(y= 1:43, x =EXP,
#         type = "p",
#         pch = 14:20)




####################################################################################################################
# ----- EDA --------------------------------------------------------------------------------------------------
####################################################################################################################


# Change in SMR sd over time?
SMR_raw = D/EXP
SMR_raw %>% 
        as.data.frame() %>% 
        mutate(dist = Namedist) %>% 
        pivot_longer(!dist, names_to = "year", values_to = "smr") %>% 
        ggplot(aes(x = year, y = smr)) +
        geom_point()

# Volatility of smr over time
SMR_raw %>% 
        as.data.frame() %>% 
        mutate(dist = Namedist) %>% 
        pivot_longer(!dist, names_to = "year", values_to = "smr") %>% 
        mutate(year = as.numeric(year)) %>% 
        ggplot(aes(x = smr, y = dist, col = year)) +
        geom_point() +
        geom_vline(xintercept = 1, 
                   linetype = "dashed") +
        theme_bw() +
        theme(axis.title.y = element_blank(),
              legend.position = c(0.9, 0.9),
              legend.title = element_blank()) +
        labs(y = "SMR")




#####################################################################################################################
# ----- Modeling in STAN ---------------------------------------------------------------------------------------------------
#####################################################################################################################


stan_data <- list(A = Nage,
                  T = Nyear,
                  R = Ndist,
                  D = D,
                  EXP = EXP)
nchain = 4
niter = 4000
nwarmup = 500
nsim = (niter-nwarmup)*nchain

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)

fit <- stan("./code/stan/smr_subnational.stan",
            data = stan_data,
            chains = nchain,
            iter=niter, 
            warmup=nwarmup, 
            seed = 444)


#####################################################################################################################
# ----- Models outputs ---------------------------------------------------------------------------------------------------
#####################################################################################################################

# SMR of all districts
theta = array(as.matrix(fit,'theta'), c(nsim, Ndist, Nyear),
          dimnames=list(paste0('sample',1:nsim),
                        Namedist, 
                        2015:2020))

df.theta = apply(theta, c(2, 3), quantile, probs = c(0.025, 0.5, 0.975)) %>% 
        as.data.frame.table() %>% 
        pivot_wider(names_from = Var1, values_from = Freq) %>% 
        rename(dist = Var2,
               year = Var3,
               lower = `2.5%`,
               median = `50%`,
               upper = `97.5%`)


# Save estimates for maping
# saveRDS(df.theta,
#         "./data/estimates/smr_for_maps.rda")

# Load estimates (avoid running previous stan code)
# df.theta = readRDS("./data/estimates/smr_for_maps.rda")

# get order of districts by SMR value in 2020 
# for ploting purposes
order_theta = df.theta %>% filter(year == 2020) %>% 
        mutate(dist = as.character(dist),
               dist = case_when(substr(dist, 16, 17) == "de" ~ substr(dist, 19, nchar(dist)),
                                substr(dist, 16, 17) == "d'" ~ substr(dist, 18, nchar(dist)))) %>% 
        arrange(median) %>% 
        pull(dist)

df.theta %>% 
        # add plot details
        mutate(year_c = ifelse(year == 2020, "2020", "2015-2019"),
               dist = as.character(dist),
               dist = case_when(substr(dist, 16, 17) == "de" ~ substr(dist, 19, nchar(dist)),
                                substr(dist, 16, 17) == "d'" ~ substr(dist, 18, nchar(dist))),
               dist = factor(dist,
                             labels = order_theta,
                             levels = order_theta),
               alpha = ifelse(year == 2020, 1, 0.95)) %>% 
        ggplot(aes(x = median, y = dist)) +
        geom_pointrange(aes(x=median, xmin=lower, xmax=upper, col = year, alpha = alpha, shape = year),
                        position = position_dodge(width = 0.55)) +        
        geom_vline(xintercept = 1, linetype = "dashed") +
        ggforestplot::geom_stripes(odd = "#33333333", even = "#00000000") +
        theme(legend.title = element_blank(),
              legend.position = c(0.14, 0.8),
              axis.title.y = element_blank(),
              axis.text.y = element_text(face = "bold",
                                         size = 10,
                                         hjust = 0),
              plot.margin = grid::unit(c(0,0,0,0), "mm")) +
        scale_y_discrete(position = "right") +
        scale_alpha_continuous(guide = FALSE,
                               range = c(0.4, 1)) +
        scale_color_manual(values = c("2020" = "#460B6AFF",
                                      "2019" = "#238A8DFF",
                                      "2018" = "#20A386FF",
                                      "2017" = "#3CBC75FF",
                                      "2016" = "#74D055FF",
                                      "2015" = "#B8DE29FF")) +
        scale_shape_manual(values = c("2020" = 17,
                                      "2019" = 16,
                                      "2018" = 16,
                                      "2017" = 16,
                                      "2016" = 16,
                                      "2015" = 16)) +
        labs(x = expression(hat(SMR)))

# ggsave(filename = "./figures/smr.pdf",
#        device = "pdf")


# Obtain highest vs lowest SMR ratio for each year

# lapply(2015:2020, function(y) {
# 
#         range.t <- df.theta %>% filter(year == y) %>%
#                 pull(median) %>%
#                 range(.)
#         dist.t <- df.theta %>% filter(year == y,
#                                       median %in% range.t)
#         ratio.t <- range.t[2]/range.t[1]
#         return(list(ratio = ratio.t,
#                     dist = dist.t))
# })

# Obtain district with small temporal variation

# ratio.D <- lapply(unique(df.theta$dist), function(d) {
#         
#         ratio.d <- df.theta %>% filter(dist == d,
#                                        year != 2020) %>% 
#                 pull(median) %>% 
#                 range(.)
#         return(c(ratio = ratio.d[2]/ratio.d[1],
#                     ratio.d = ratio.d))
# })
# range.D <- tibble(ratio = do.call(rbind, ratio.D)[, 1],
#                   lowest = do.call(rbind, ratio.D)[, 2],
#                   highest = do.call(rbind, ratio.D)[, 3],
#                   dist = unique(df.theta$dist)) %>% 
#         arrange(ratio)

# SMR posterior ranking distribution

# Empty rank matrix
rank_mat <- matrix(0, ncol = 43, nrow = 43)
dist_id <- 1:43
# Algorithm to obtain posterior rank probabilities
for (k in 1:nsim){
        # cbind dist_id and theta (6 = 2020)
        temp <- cbind(dist_id, theta[k, , 6])
        # order dist_id by theta and store dist_id in rank order
        rank <- temp[order(temp[,2], decreasing = TRUE), 1]
        # Add +1 to the rank of each district 
        # in this particular k posterior draw
        for (j in 1:43){
                
                rank_mat[rank[j], j] <- rank_mat[rank[j], j] + 1      
        }
}
# Check: all(rowSums(rank_mat) == nsim)

# Convert rank_mat into rank_prob_mat
rank_prob_mat <- rank_mat/nsim
# Obtain expected rank for each district
exp_rank <- as.numeric(rank_prob_mat %*% c(1:43))
# store in a data set
df_exp_rank <- tibble(dist_i = 1:43,
                      exp_rank = exp_rank,
                      dist = Namedist) 
# Exctract order according to expected rank for 
# ploting purposes
order_dist <- df_exp_rank %>% 
        mutate(dist = case_when(substr(dist, 16, 17) == "de" ~ substr(dist, 19, nchar(dist)),
                                substr(dist, 16, 17) == "d'" ~ substr(dist, 18, nchar(dist)))) %>% 
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
               dist = Namedist[dist_i]) %>% 
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
              legend.position = c(0.9, 0.8),
              plot.margin = grid::unit(c(0,0,0,2), "mm")) +
        labs(x = "Rank",
             fill = "Probability")

# ggsave(filename = "./figures/posterior_rank_prob_smr.pdf",
#        device = "pdf")



# Overall theta distribution:
# not centered over 1 because different
# size of district.
mu_theta = array(as.matrix(fit,'mu_theta'), c(nsim, Nyear),
                 dimnames=list(paste0('sample',1:nsim),
                               2015:2020))
# Plot posterior distribution mean of the SMR
# distribution in each year
mu_theta %>% 
        as_tibble() %>% 
        mutate(sim = paste0("sim", 1:nsim)) %>% 
        pivot_longer(!sim, values_to = "mu_theta", names_to = "year") %>% 
        mutate(year = as.numeric(year)) %>% 
        ggplot(aes(x = mu_theta, col = year, fill = year, group = year)) +
        geom_density(alpha=0.2) +
        theme_bw() +
        theme(legend.position = c(0.8, 0.8),
              legend.title = element_blank()) +
        labs(x = expression(mu[theta]))




# Std error in the SMR distribution
sigma_theta = array(as.matrix(fit,'sigma_theta'), c(nsim, Nyear),
                 dimnames=list(paste0('sample',1:nsim),
                               2015:2020)) 
# quantiles of std error distributions
apply(sigma_theta, 2, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
        as.data.frame.table() %>% 
        rename("quantile" = Var1,
               "year" = Var2,
               "value" = Freq) 
# plot SMR distributions over 2015-2020
sigma_theta %>% 
        as_tibble() %>% 
        mutate(sim = paste0("sim", 1:nsim)) %>% 
        pivot_longer(!sim, values_to = "sigma_theta", names_to = "year") %>% 
        mutate(year = case_when(year == "2015" ~ "2015 -- 0.11 (0.09-0.15)",
                                year == "2016" ~ "2016 -- 0.12 (0.10-0.16)",
                                year == "2017" ~ "2017 -- 0.11 (0.09-0.14)",
                                year == "2018" ~ "2018 -- 0.12 (0.09-0.15)",
                                year == "2019" ~ "2019 -- 0.12 (0.10-0.15)",
                                year == "2020" ~ "2020 -- 0.15 (0.12-0.19)")) %>% 
        # mutate(year = as.numeric(year)) %>% 
        ggplot(aes(x = sigma_theta, col = year, fill = year, group = year)) +
        geom_density(alpha=0.2) +
        scale_color_manual(values = c("2020 -- 0.15 (0.12-0.19)" = "#460B6AFF",
                                      "2019 -- 0.12 (0.10-0.15)" = "#238A8DFF",
                                      "2018 -- 0.12 (0.09-0.15)" = "#20A386FF",
                                      "2017 -- 0.11 (0.09-0.14)" = "#3CBC75FF",
                                      "2016 -- 0.12 (0.10-0.16)" = "#74D055FF",
                                      "2015 -- 0.11 (0.09-0.15)" = "#B8DE29FF")) +
        scale_fill_manual(values = c("2020 -- 0.15 (0.12-0.19)" = "#460B6AFF",
                                     "2019 -- 0.12 (0.10-0.15)" = "#238A8DFF",
                                     "2018 -- 0.12 (0.09-0.15)" = "#20A386FF",
                                     "2017 -- 0.11 (0.09-0.14)" = "#3CBC75FF",
                                     "2016 -- 0.12 (0.10-0.16)" = "#74D055FF",
                                     "2015 -- 0.11 (0.09-0.15)" = "#B8DE29FF")) +
        theme_bw() +
        theme(legend.position = c(0.65, 0.75),
              legend.title = element_blank(),
              panel.grid = element_blank(),
              plot.margin = grid::unit(c(0,0,0,2), "mm")) +
        labs(y = expression(f(sigma[theta])),
             x = expression(sigma[theta])) +
        annotate("text", x = 0.2, y = 28.5, label = "Year -- Posterior median (95% C.I)", 
                 size = 4, hjust = .5) 

# ggsave(filename = "./figures/posterior_density_sigma_theta.pdf",
#        device = "pdf")
        
