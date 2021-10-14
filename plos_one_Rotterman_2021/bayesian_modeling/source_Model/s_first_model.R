setwd("~/Desktop/FirstModel/")

library(magrittr)
library(dplyr)
library(ggplot2)
library(rstan)
library(tidybayes)
library(emmeans)
library(broom)
library(brms)
library(modelr)
library(forcats)
library(cowplot)
library(RColorBrewer)
library(gganimate)
library(ggdist)
library(loo)

theme_set(theme_tidybayes() + panel_border())

# https://bookdown.org/marklhc/notes_bookdown/group-comparisons.html


scaled.popdata %>%
  ggplot(aes(x = Dis_2_Soma, y = Muscle)) +
  geom_point(alpha = 0.5) +
  ylab("Muscle")
## do not exclude these options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#compile model
model <- stan_model('~/Desktop/FirstModel/first_model.stan')

# run model
m = sampling(model, 
             data = compose_data(na.omit_popdata), 
             control = list(adapt_delta=0.999,max_treedepth = 15),
             chains = 4)


m %<>% recover_types(na.omit_popdata)

y <- na.omit_popdata$Dis_2_Soma
yrep <- posterior_predict(m, draws = 500)
dim(yrep_poisson)


pp_check(m, nreps = 500)

# graph
m %>%
  spread_draws(Muscle_mean[Muscle]) %>%
  ggplot(aes(x = Muscle_mean, y = Muscle)) +
  stat_halfeye()

# estimates
m %>%
  spread_draws(Muscle_mean[Muscle]) %>%
  median_qi(Muscle_mean)
