# Origination Version of August 13, 2020.
# Stephen N. Housley 
# housley.nick@gmail.com
# 

# This work is licensed under the licenses 
# Paper: Creative Commons Attribution 3.0 Unported License 
# Code: GPL-3 License 
# Depends: R (>= 3.5.0)
# Version: 0.3
# Description: code to run analytics and graphic functions associated with 
# Rotterman et al 2021 Axon initial segment geometry in relation to motoneuron recruitability
# PONE-D-21-28122 PLOS ONE

#
# This program is believed to be free of errors, but it comes with no guarantee! 
# The user bears all responsibility for interpreting the results. 
#

## version Hx
# v0.1- original

# v0.2 - add custom function for correlation graphing 

# Commit comments
## added chart.Correlation.linear function to tailor chart.corr plotting function to travis's needs
## must call chart.Correlation.linear in place of chart.Correlation

# v0.3 - add frequentist analysis for pop data

# v0.4 - added  graphing for Rheo~Conductance graphing to compare AIS MN to all other MNs from the lab


# v0.5 - added  graphing for Rheo~Conductance graphing to compare AIS MN to all other MNs from the lab


# v0.6 - final refining for deposite in Dryad (10/11/2021)

# v1.0 - remove personal paths and generalize (10/11/2021)

########################### prelims ########################### 

invisible(rm(list = ls()))
invisible(gc())

########################### timer start ########################### 
ptm <- proc.time()

########################### set dirs ########################### 
# data_Name_filled<- 'AIS_NB fillings data_TR Dario_Final_CA_2.xlsx'
data_Name_filled<- 'AIS_BioPhys_MorphParam.xlsx'   ### after the input resistance and rheo have be reconciled 
mainDir <- "~/Desktop/plos_one_Rotterman_2021" ### set this filepath to the main directory 
subDir <- "final_data" ## this is the subdirectory where data is held
figFolder <-"figFolder" ## this will be created if not already in existance
saveFolder <- "saveFolder" ## this will be created if not already in existance
invisible(  ifelse(!dir.exists(file.path(mainDir, subDir, figFolder)), dir.create(file.path(mainDir, subDir,figFolder)), FALSE))
invisible(  ifelse(!dir.exists(file.path(mainDir, subDir, saveFolder)), dir.create(file.path(mainDir, subDir,saveFolder)), FALSE))

#setwd(file.path(mainDir, subDir, figFolder))


########################### load general dependencies ########################### 
packagesS<-c("devtools",
               "dplyr",
               "parallel",
               "ggplot2",
               "readxl",
             "crayon",
             "PerformanceAnalytics",
             "Hmisc",
             "tidyr",
             "GGally",
             "tibble",
             "ggplot2",
             "leaps",
             "caret",
             "rstan",
             "rstanarm",
             "coda",
             "broom",
             "tidybayes",
             "emmeans",
             "loo",
             "bayesplot",
             "magrittr",
             "ggpubr",
             "modelr",
             "broom.mixed"
)

invisible(suppressWarnings(suppressMessages(package.check <- lapply(
  packagesS,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
))))

## do not exclude these options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#invisible(suppressWarnings(suppressMessages(lapply(packagesS, function(xxx) require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))))
#invisible(lapply(packagesS, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE))))







################################################################################# 
########################### Filled cells ########################### 
################################################################################# 

########################### Custom Functions ########################### 
# http://padamson.github.io/r/ggally/ggplot2/ggpairs/2016/02/16/multiple-regression-lines-with-ggpairs.html 
#two regression lines in each plot in the lower diagonal elements in which both X and Y data are continuous.
LoesLMPlotting <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    # geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}



# https://www.statology.org/plot-confidence-interval-r/
chart.Correlation.linear <-
  function (R, histogram = TRUE, method=c("pearson", "kendall", "spearman"), ...)
  { # @author R Development Core Team
    # @author modified by Peter Carl & Marek Lahoda
    # Visualization of a Correlation Matrix. On top the (absolute) value of the correlation plus the result 
    # of the cor.test as stars. On botttom, the bivariate scatterplots, with a linear regression fit. 
    # On diagonal, the histograms with probability, density and normal density (gaussian) distribution.
    
    x = checkData(R, method="matrix")
    
    if(missing(method)) method=method[1] #only use one
    cormeth <- method
    
    # Published at http://addictedtor.free.fr/graphiques/sources/source_137.R
    panel.cor <- function(x, y, digits=2, prefix="", use="pairwise.complete.obs", method=cormeth, cex.cor, ...)
    {
      usr <- par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r <- cor(x, y, use=use, method=method) # MG: remove abs here
      txt <- format(c(r, 0.123456789), digits=digits)[1]
      txt <- paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
      
      test <- cor.test(as.numeric(x),as.numeric(y), method=method)
      # borrowed from printCoefmat
      Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                       symbols = c("***", "**", "*", ".", " "))
      # MG: add abs here and also include a 30% buffer for small numbers
      text(0.5, 0.5, txt, cex = cex * (abs(r) + .3) / 1.3)
      text(.8, .8, Signif, cex=cex, col=2)
    }
    
    #remove method from dotargs
    dotargs <- list(...)
    dotargs$method <- NULL
    rm(method)
    
    hist.panel = function (x, ...=NULL ) {
      par(new = TRUE)
      hist(x,
           col = "light gray",
           probability = TRUE,
           axes = FALSE,
           main = "",
           breaks = "FD")
      lines(density(x, na.rm=TRUE),
            col = "red",
            lwd = 1)
      # adding line representing density of normal distribution with parameters correponding to estimates of mean and standard deviation from the data 
      ax.x = seq(min(x), max(x), 0.1)                                                  # ax.x containts points corresponding to data range on x axis
      density.est = dnorm(ax.x, mean = mean(x), sd = sd(x))   # density corresponding to points stored in vector ax.x 
      lines(ax.x, density.est, col = "blue", lwd = 1, lty = 1)                                # adding line representing density into histogram
      rug(x)
    }
    
    # Linear regression line fit over points
    reg <- function(x, y, ...) {
      
      newx <- seq(min(x), max(x), length.out=100)
      preds <- predict(lm(y~x), newdata = data.frame(x=newx), interval = 'confidence')
      #fill in area between regression line and confidence interval
      polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = '#3250a8', border = NA)
      points(x,y, ...)
      
      abline(lm(y~x), col = "red") 
      #add dashed lines for confidence bands
      lines(newx, preds[ ,3], lty = 'dashed', col = 'blue')
      lines(newx, preds[ ,2], lty = 'dashed', col = 'blue')
    }
    
    # Draw the chart
    if(histogram)
      pairs(x, gap=0, lower.panel=reg, upper.panel=panel.cor, diag.panel=hist.panel)
    else
      pairs(x, gap=0, lower.panel=reg, upper.panel=panel.cor) 
  }

########################### Load Data ########################### 
allData<-read_excel(file.path(mainDir,subDir,data_Name_filled), 
                    na = "NA")

########################### Data Wrangling ########################### 

### find large correlations with Rheo (0.7)
sigCorParamsWithRheo<-allData %>% 
  filter(dendritic_Origin == 0) %>% 
  select(Reho,	
         soma_Surf,
         soma_Vol,	
         hilock_Lgth,	
         ais_Lgth,	
         ahp_amp,
         ahp_half,
         rin ,
         cond_delay,
         twt_Amp,
         ttp,
         ) %>% 
  drop_na() %>%
  cor() %>% 
  round(digits = 2) %>% as.data.frame()%>%
    rownames_to_column('params') %>%
    # filter_at(vars(contains("Reho")), .vars_predicate =  any_vars(abs(.) > 0.7))%>%
   filter(Reho >0.7 & Reho <0.99 |Reho < -0.7 & Reho > -0.99) %>%   ### 0.99 is to remove perfect correlation with self
column_to_rownames('params') %>% rownames() 



  

### subset the relevant params
allDataScreened<-allData %>% 
  select(Reho,	
         soma_Surf,
         soma_Vol,	
         hilock_Lgth,	
         ais_Lgth,	
         ahp_amp,
         ahp_half,
         rin ,
         cond_delay,
         twt_Amp,
         ttp
  ) %>%
  drop_na() %>% 
  select(Reho,any_of(sigCorParamsWithRheo)) 

### find un-correlations with Rheo (0.4)
unCorParamsWithRheo<-allData %>% 
  filter(dendritic_Origin == 0) %>% 
  select(Reho,	
         soma_Surf,
         soma_Vol,	
         hilock_Lgth,	
         # hilock_Surf,	
         # hilocka_Vol,	
         ais_Lgth,	
         # ais_Surf,	
         # ais_Vol, 
         ap,
         ahp_amp, 
         ahp_half, 
         rin ,
         # mV,
         # thresh,	
         # sd_thresh,	
         cond_delay,	
         # twt_Dealy,	
         twt_Amp,	
         # twt_half_width ,
         ttp,	
         # relax
         ) %>% drop_na() %>%
  cor() %>% 
  round(digits = 2) %>% as.data.frame()%>%
  rownames_to_column('params') %>%
  filter_at(vars(contains("Reho")), .vars_predicate =  any_vars(abs(.) < 0.4))%>%
  column_to_rownames('params') %>% rownames() 

### subset the uncorrelated  params
allDataScreened_un<-allData %>% 
  select(Reho,	
         soma_Surf,
         soma_Vol,	
         hilock_Lgth,	
         # hilock_Surf,	
         # hilocka_Vol,	
         ais_Lgth,	
         # ais_Surf,	
         # ais_Vol, 
         ap,
         ahp_amp, 
         ahp_half, 
         rin ,
         # mV,
         # thresh,	
         # sd_thresh,	
         cond_delay,	
         # twt_Dealy,	
         twt_Amp,	
         # twt_half_width ,
         ttp,	
         # relax
  ) %>% drop_na() %>%
  select(Reho,one_of(unCorParamsWithRheo)) 

#### p values
res1<-allData %>% 
  filter(dendritic_Origin == 0) %>% 
  select(Reho,	
         soma_Surf,
         soma_Vol,	
         hilock_Lgth,	
         # hilock_Surf,	
         # hilocka_Vol,	
         ais_Lgth,	
         # ais_Surf,	
         # ais_Vol, 
         ap,
         ahp_amp, 
         ahp_half, 
         rin ,
         # mV,
         # thresh,	
         # sd_thresh,	
         cond_delay,	
         # twt_Dealy,	
         twt_Amp,	
         # twt_half_width ,
         ttp,	
         # relax
  ) %>% drop_na() %>%
  as.matrix() %>% rcorr(type="pearson")

res2<-res1$r %>% 
  as.data.frame() %>% 
  rownames_to_column('params') %>%
  filter(Reho >0.7 & Reho <0.99 |Reho < -0.7 & Reho > - 0.99) %>%
  column_to_rownames('params') %>% 
  select(Reho)%>%
  round(digits = 3)%>%
  `colnames<-`(c('R2'))

res3<-res1$P %>% 
  as.data.frame() %>% 
  rownames_to_column('params') %>%
  filter(params %in% rownames(res2)) %>%
  column_to_rownames('params') %>% 
  select(Reho)%>%
  round(digits = 5)%>%
  `colnames<-`(c('pvalue'))


### filter only structural params
StructuralParamsWithRheo<-allData %>% 
  filter(dendritic_Origin == 0) %>% 
  select(Reho,	
         soma_Surf,
         soma_Vol,	
         hilock_Lgth,	
         # hilock_Surf,	
         # hilocka_Vol,	
         ais_Lgth,	
         # ais_Surf,	
         # ais_Vol
         ) %>% 
  drop_na() 




########################### quick visualization_ print messages ########################### 
cat(red("These params are sig correlated with Rheo (p <0.01 & R2<0.7)\n"))
print(cbind('R2'=res2$R2,'pvalue'=res3),   sep="\n")  

cat(red("These params are un-correlated with Rheo (<0.4)\n"))
cat(unCorParamsWithRheo,   sep="\n")  

########################### analyses/modeling ########################### 

# do any of these variables differ by Group?
responses <- paste(colnames(allDataScreened[ , !(names(allDataScreened) %in% c("Reho"))]), collapse=",") ## dynamically grab the dependent vars
myformula <- as.formula( paste0( "cbind(", responses , ")~ Reho" ) ) ## write the formula
summary(manova( myformula, data = allDataScreened )) ## run the manova

models_stepwise <- regsubsets(Reho~., data = allDataScreened, nvmax = 4,
                     method = "exhaustive")
models_stepwise.sum <- summary(models_stepwise)
cat(red("This is how you can perform model selection\n"))
cat(red("How many params are optimal?\n"))
data.frame(
  Adj.R2 = which.max(models_stepwise.sum$adjr2),
  CP = which.min(models_stepwise.sum$cp),
  BIC = which.min(models_stepwise.sum$bic)
)
cat(red("Which ones?\n"))
# test<-as.data.frame(summary(models_stepwise)$outmat)  #### code to dynamically pull the sig params out/// currently just hard coded below
# test[test == " "] <- NA
# colnames(Filter(function(x)!all(is.na(x)), test[3,]))
print(summary(models_stepwise)$outmat)

cat(red("Whats the information gained by adding more?\n"))
print(models_stepwise.sum$rsq)


models_final <- anova(aov(Reho~ hilock_Lgth+twt_Amp, data = allDataScreened))

# data.frame(
#   Adj.R2 = which.max(res.sum$adjr2),
#   CP = which.min(res.sum$cp),
#   BIC = which.min(res.sum$bic)
# )
########################### saving data ########################### 
setwd(file.path(mainDir,subDir,saveFolder))
write.csv(cbind('R2'=res2$R2,res3),file.path(mainDir,subDir,saveFolder,"data.rheo_Sig_correlation.csv"), row.names = TRUE)

write.table(summary(models_stepwise)$outmat, file.path(mainDir,subDir,saveFolder, "data.aov_model_selection.txt"), sep = "\t",
            row.names = TRUE, col.names = NA)

write.table(models_final, file.path(mainDir,subDir,saveFolder, "data.aov_model_final.txt"), sep = "\t",
            row.names = TRUE, col.names = NA)

cat(red("after selection, this is the final model predicting Rheo \n"))
print(as.data.frame(models_final),   sep="\n")  


########################### saving figures ###########################
setwd(file.path(mainDir,subDir,figFolder))

Figure_1_cor<-ggpairs(allDataScreened, 
                      lower = list(continuous = LoesLMPlotting),
                      diag=list(continuous="barDiag"))+
theme(legend.position = "none",
      panel.grid.major = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(linetype = "dashed", colour = "black", fill = NA))

invisible(suppressMessages(suppressWarnings(ggsave(
  "ggpairs.rheo_correlation.pdf",
  plot = Figure_1_cor,
  device = "pdf",
  path = file.path(mainDir, subDir,figFolder),
  scale = 1,
  width = 12,
  height = 12,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
))))

pdf(file.path(mainDir,subDir,figFolder,"chart.Correlation.rheo_correlation.pdf"))
allDataScreened%>%
  chart.Correlation.linear(histogram=TRUE, pch=19) #, use = "complete.obs"
invisible(dev.off())

Figure_1_un<-ggpairs(allDataScreened_un, lower = list(continuous = LoesLMPlotting))+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(linetype = "dashed", colour = "black", fill = NA))

invisible(suppressMessages(suppressWarnings(ggsave(
  "ggpairs.rheo_un_correlation.pdf",
  plot = Figure_1_un,
  device = "pdf",
  path = file.path(mainDir, subDir,figFolder),
  scale = 1,
  width = 12,
  height = 12,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
))))

pdf(file.path(mainDir,subDir,figFolder,"chart.Correlation.rheo_un_correlation.pdf"))
allDataScreened_un%>%
  chart.Correlation.linear(histogram=TRUE, pch=19) #, use = "complete.obs"
invisible(dev.off())







### figures for only structural params

Figure_3_cor_16_structural<-ggpairs(StructuralParamsWithRheo, lower = list(continuous = LoesLMPlotting))+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(linetype = "dashed", colour = "black", fill = NA))

invisible(suppressMessages(suppressWarnings(ggsave(
  "ggpairs.rheo_correlation_16cells.pdf",
  plot = Figure_3_cor_16_structural,
  device = "pdf",
  path = file.path(mainDir, subDir,figFolder),
  scale = 1,
  width = 12,
  height = 12,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE
))))



pdf(file.path(mainDir,subDir,figFolder,"chart.Correlation.rheo_correlation_16cells.pdf"))
StructuralParamsWithRheo%>%
  chart.Correlation.linear(histogram=TRUE, pch=19) #, use = "complete.obs"
invisible(dev.off())



cat(green("Filled cells Figs successfully created",
          "Find them here", file.path(mainDir,subDir,figFolder),   sep="\n"))

########################### Clean up ###########################
invisible(suppressWarnings(rm(Figure_1_cor,Figure_1_un,models_deltaBYrheo, models_final,models_stepwise,package.check,res1,res2,res3)))






#################################################################################
########################### Filled cells Bayesian Modeling ###########################
#################################################################################



########################### Custom Functions ###########################
# na

########################### Load Data ###########################
# hold from previous

########################### Data Wrangling ###########################
# hold from previous

### subset the relevant params
scaled.data<-allData %>%
  select(exp,
         Reho,
         soma_Surf,
         soma_Vol,
         hilock_Lgth,
         # hilock_Surf,
         # hilocka_Vol,
         ais_Lgth,
         # ais_Surf,
         # ais_Vol,
         ap,
         ahp_amp,
         ahp_half,
         rin ,
         # mV,
         # thresh,
         # sd_thresh,
         cond_delay,
         # twt_Dealy,
         twt_Amp,
         # twt_half_width ,
         ttp,
         # relax
         ) %>%
  drop_na() %>%
  select(exp,Reho,any_of(sigCorParamsWithRheo))

scaled.data[c(2:length(scaled.data))] <- lapply(scaled.data[c(2:length(scaled.data))], function(x) c(scale(x)))


########################### quick visuakization ###########################

########################### analyses/modeling ###########################
SEED<-set.seed(315789)
sink("/dev/null")

fmla_1 <- Reho ~ hilock_Lgth
# fmla_2 <- Reho ~ hilock_Surf
# fmla_3 <- Reho ~ hilocka_Vol
fmla_4 <- Reho ~ twt_Amp

fmla_5 <- Reho ~ hilock_Lgth+twt_Amp
# fmla_6 <- Reho ~ hilock_Surf+twt_Amp
# fmla_7 <- Reho ~ hilocka_Vol+twt_Amp

fmla_8 <- Reho ~ hilock_Lgth*twt_Amp
# fmla_9 <- Reho ~ hilock_Surf*twt_Amp
# fmla_10 <- Reho ~ hilocka_Vol*twt_Amp
# 
# fmla_11 <- Reho ~ hilocka_Vol+twt_Amp+hilock_Lgth

# fmla_5 <- Reho ~ hilock_Lgth + ais_Lgth +
#   (1 + hilock_Lgth + ais_Lgth  | exp)

my_prior <- normal(location = c(0), scale = c(1), autoscale = TRUE)
my_prior_intercept <- normal(location = c(0), scale = c(5),autoscale = TRUE)
my_prior_aux <- exponential(1, autoscale=TRUE)

#####Build a model and perform model comparison on different predictors
post1 <-stan_glm(formula=fmla_1,
                 data = scaled.data,
                 algorithm = "sampling",
                 family = "gaussian",
                 adapt_delta = .99,
                 prior= normal(location = c(0), scale = c(1), autoscale = TRUE),
                 prior_aux = exponential(1, autoscale=TRUE),
                 prior_intercept = normal(location = c(0), scale = c(5),autoscale = TRUE),
                 iter = 4000,
                 warmup = 400,
                 chains = 4,
                 thin = 2,
                 seed = SEED)
# post2 <- update(post1, formula = fmla_2)
# post3 <- update(post1, formula = fmla_3)
post4 <- update(post1, formula = fmla_4)


post5 <-stan_glm(formula=fmla_5,
                 data = scaled.data,
                 algorithm = "sampling",
                 family = "gaussian",
                 adapt_delta = .99,
                 prior= normal(location = c(0, 0), scale = c(1, 1), autoscale = TRUE),
                 prior_aux = exponential(1, autoscale=TRUE),
                 prior_intercept = normal(location = c(0), scale = c(5),autoscale = TRUE),
                 iter = 4000,
                 warmup = 400,
                 chains = 4,
                 thin = 2,
                 seed = SEED)
# post6 <- update(post5, formula = fmla_6)
# post7 <- update(post5, formula = fmla_7)

post8 <-stan_glm(formula=fmla_8,
                 data = scaled.data,
                 algorithm = "sampling",
                 family = "gaussian",
                 adapt_delta = .99,
                 prior= normal(location = c(0, 0, 0), scale = c(1, 1, 1), autoscale = TRUE),
                 prior_aux = exponential(1, autoscale=TRUE),
                 prior_intercept = normal(location = c(0), scale = c(5),autoscale = TRUE),
                 iter = 4000,
                 warmup = 400,
                 chains = 4,
                 thin = 2,
                 seed = SEED)
# post9 <- update(post8, formula = fmla_9)
# post10 <- update(post8, formula = fmla_10)
# post11 <- update(post8, formula = fmla_11)
# 

########################### validation ###########################

#### Leave-one-out cross validation of models
loo1 <- loo(post1, save_psis = TRUE, cores = 4, k_threshold =0.7)
# loo2 <- loo(post2, save_psis = TRUE, cores = 4, k_threshold =0.7)
# loo3 <- loo(post3, save_psis = TRUE, cores = 4, k_threshold =0.7)
loo4 <- loo(post4, save_psis = TRUE, cores = 4, k_threshold =0.7)
loo5 <- loo(post5, save_psis = TRUE, cores = 4, k_threshold =0.7)
# loo6 <- loo(post6, save_psis = TRUE, cores = 4, k_threshold =0.7)
# loo7 <- loo(post7, save_psis = TRUE, cores = 4, k_threshold =0.7)
loo8 <- loo(post8, save_psis = TRUE, cores = 4, k_threshold =0.7)
# loo9 <- loo(post9, save_psis = TRUE, cores = 4, k_threshold =0.7)
# loo10 <- loo(post10, save_psis = TRUE, cores = 4, k_threshold =0.7)
# loo11 <- loo(post11, save_psis = TRUE, cores = 4, k_threshold =0.7)

sink()


loo_compare(loo1,
            # loo2, 
            # loo3, 
            loo4, 
            loo5, 
            # loo6, 
            # loo7, 
            loo8, 
            # loo9, 
            # loo10, 
            # loo11
            )


launch_shinystan(post5)
summary(post5)
## effect size
sink("/dev/null")


mcmc = as.matrix(post5)
newdata = expand.grid(hilock_Lgth = c(min(scaled.data$hilock_Lgth), max(scaled.data$hilock_Lgth)), twt_Amp = (-2:2) *
                        sd(scaled.data$hilock_Lgth))
Xmat = model.matrix(~hilock_Lgth + twt_Amp, newdata)
coefs = mcmc[, c("(Intercept)", "hilock_Lgth", "twt_Amp")]
fit = coefs %*% t(Xmat)
s1 = seq(1, 9, b = 2)
s2 = seq(2, 10, b = 2)
## Raw effect size
RES = tidyMCMC(as.mcmc(fit[, s2] - fit[, s1]), conf.int = TRUE, conf.method = "HPDinterval")

## Cohen's D
cohenD = (fit[, s2] - fit[, s1])/sqrt(mcmc[, "sigma"])
cohenDES = tidyMCMC(as.mcmc(cohenD), conf.int = TRUE, conf.method = "HPDinterval")
sink()

## model R^2

mcmc <- as.matrix(post5)
Xmat = model.matrix(~hilock_Lgth + twt_Amp, scaled.data)
coefs = mcmc[, c("(Intercept)", "hilock_Lgth", "twt_Amp")]
fit = coefs %*% t(Xmat)
resid = sweep(fit, 2, scaled.data$Reho, "-")
var_f = apply(fit, 1, var)
var_e = apply(resid, 1, var)
R2 = var_f/(var_f + var_e)
tidyMCMC(as.mcmc(R2), conf.int = TRUE, conf.method = "HPDinterval")


summary(lm(Reho ~ hilock_Lgth+twt_Amp, scaled.data))

########################### model behavior ###########################

# pp_check(post1, nreps = 500)
# pp_check(post1, plotfun = "stat_2d", stat = c("mean", "sd"))
# pp_check(post1, plotfun = "scatter_avg") # y vs. average yrep
# bayesplot::ppc_scatter_avg(y = scaled.data$Reho, yrep = posterior_predict(post1))


########################### saving data ###########################


write.csv(tidyMCMC(post5, conf.int = TRUE, conf.method = "HPDinterval", conf.level = 0.95),
          file.path(mainDir,subDir,saveFolder, "data.rstan_model_final.csv"), row.names = TRUE)
cat(red("Bayesian Parameter estimates of important coeffs predicting Rheo (remember if they do not cross 0 they are 'sig') \n"))
tidyMCMC(post5, conf.int = TRUE, conf.method = "HPDinterval", conf.level = 0.95)

write.csv(tidyMCMC(as.mcmc(R2), conf.int = TRUE, conf.method = "HPDinterval", conf.level = 0.95),
          file.path(mainDir,subDir,saveFolder, "data.rstan_model_final_adj_r2_filled.csv"), row.names = TRUE)
cat(red("Bayesian Parameter estimates_filled_cells R^2 \n"))
tidyMCMC(as.mcmc(R2), conf.int = TRUE, conf.method = "HPDinterval", conf.level = 0.95)



write.csv(mean(cohenDES$estimate),
          file.path(mainDir,subDir,saveFolder, "data.rstan_model_CohenD.csv"), row.names = TRUE)
cat(red(
  "The Cohen's D associated change over the observed range of Rheo is " %+%
    blue$underline$bold( mean(cohenDES$estimate)) %+%
    '\n'
))

write.csv(rbind(mean(RES$estimate),min(RES$conf.low),max(RES$conf.high)),
          file.path(mainDir,subDir,saveFolder, "data.rstan_model_raw_effectSize.csv"), row.names = TRUE)
cat(red(
  "On average, rheo increases by " %+%
    blue$underline$bold(mean(RES$estimate)) %+%
    " over the observed range.
    We are 95% confident that the increase is between " %+%
    blue$underline$bold(min(RES$conf.low)) %+%
    "and " %+%
    blue$underline$bold(max(RES$conf.high)) %+%
    '\n'
))



########################### saving figures ###########################

########################### Clean up ###########################
invisible(suppressWarnings(rm(coefs,cohenD,cohenDES, fit, loo1,loo2, loo3, loo4, loo5, loo6, loo7, loo8, loo9, loo10, mcmc, my_prior, my_prior_aux, my_prior_intercept, newdata, p, package.check, post1,post2, post3, post4,post5,post6,post8,post9,post10, Xmat)))





#################################################################################
########################### Population data ###########################
#################################################################################

########################### Custom Functions ###########################

########################### Load Data ###########################
data_Name_pop<- "Population_Data.xlsx"
popData<-read_excel(file.path(mainDir,subDir,data_Name_pop),
                    na = "NA")

########################### Data Wrangling ###########################
cols.num <- c("Axon_Hil_Area","AIS_Area","AIS_Area_corrected","Soma_Diam","AIS_Length") ## get columns with NAs
suppressWarnings(popData[cols.num] <- sapply(popData[cols.num],as.numeric)) ## convert all to numeric instead of character class
na.omit_popdata<- na.omit(popData) ## remove NAs
na.omit_popdata$Muscle<- as.factor(na.omit_popdata$Muscle)

scaled.popdata<-na.omit_popdata
scaled.popdata[c(4:9)] <- lapply(na.omit_popdata[c(4:9)], function(x) c(scale(x)))
scaled.popdata$Dis_2_Soma<-scaled.popdata$Dis_2_Soma+abs(min(scaled.popdata$Dis_2_Soma))

########################### quick visuakization ###########################
# boxplot(na.omit_popdata$Dis_2_Soma ~ Muscle, na.omit_popdata)

########################### analyses/modeling ###########################
SEED<-set.seed(315789)
sink("/dev/null")

## old model
# post1 <- stan_glm(formula=Dis_2_Soma ~ Muscle,
#                   data = na.omit_popdata,
#                   adapt_delta = .99,
#                   iter = 4000,
#                   warmup = 400,
#                   chains = 4,
#                   thin = 2,
#                   seed=SEED)

## new model
## shifted_logNorm ##
setwd("~/Desktop/plos_one_Rotterman_2021/bayesian_modeling/source_Model/")
# rstan:::rstudio_stanc("shifted_logNorm.stan")
m3 <- stan("shifted_logNorm_partial_pool_sigma.stan",
           data = list(N1 = sum(na.omit_popdata$Muscle == "MG"),
                       N2 = sum(na.omit_popdata$Muscle == "SOL"),
                       y1 = na.omit_popdata$Dis_2_Soma[which(na.omit_popdata$Muscle == "MG")],
                       y2 = na.omit_popdata$Dis_2_Soma[which(na.omit_popdata$Muscle == "SOL")]),
           control = list(adapt_delta = 0.999)
)

# m4<-   stan_glmer(Soma_Diam ~ 1|Muscle, data = na.omit_popdata, iter = 4000, warmup = 500,
#                           chains = 4, refresh = 0,   adapt_delta = 0.99, seed = 1234)
# 
# cat(red("95% and 66% PPD of Soma_Diam in each motor pool  \n"))
# tidyMCMC(m4, conf.int = TRUE, conf.method = "HPDinterval", conf.level = 0.95)
# 


########### get useful items out of the model ###########
## extract PPD reps
y1 <- na.omit_popdata$Dis_2_Soma[which(na.omit_popdata$Muscle == "MG")]
y2 <- na.omit_popdata$Dis_2_Soma[which(na.omit_popdata$Muscle == "SOL")]
# Replicated data (randomly sampling 100)
y1rep <- as.matrix(m3, pars = "y1rep")[sample.int(2000, 250), ]
y2rep <- as.matrix(m3, pars = "y2rep")[sample.int(2000, 250), ]


# Extract posterior samples of the conditional means
post_sam <- as.data.frame(m3,
                          pars = c("MG", "SOL", "sigma", "sigma_1"))





sink()

########################### validation ###########################

########################### model behavior ###########################

########################### quick analysis (frequentists) ###########################
popdata_summary_manova<- summary.aov(manova( cbind(Dis_2_Soma, 
                          Axon_Hil_Area,
                          AIS_Area,
                          AIS_Length,
                          Soma_Diam,
                          AIS_Area_corrected) ~ Muscle, data = na.omit_popdata )
            ) ## run the manova

########################### saving data ###########################

### just some summaries of raw data ##
data.na.omit_popdata_summary<- na.omit_popdata %>%
  group_by(Muscle) %>%
  summarise(
    count = n(),
    mean_Dis_2_Soma = mean(Dis_2_Soma, na.rm = TRUE),
    sd_Dis_2_Soma = sd(Dis_2_Soma, na.rm = TRUE),
    mean_Axon_Hil_Area = mean(Axon_Hil_Area, na.rm = TRUE),
    sd_Axon_Hil_Area = sd(Axon_Hil_Area, na.rm = TRUE),
    mean_AIS_Length = mean(AIS_Length, na.rm = TRUE),
    sd_AIS_Length = sd(AIS_Length, na.rm = TRUE),
    mean_AIS_Area = mean(AIS_Area, na.rm = TRUE),
    sd_AIS_Area = sd(AIS_Area, na.rm = TRUE),
    mean_Soma_Diam = mean(Soma_Diam, na.rm = TRUE),
    sd_Soma_Diam = sd(Soma_Diam, na.rm = TRUE),
    mean_AIS_Area_corrected = mean(AIS_Area_corrected, na.rm = TRUE),
    sd_AIS_Area_corrected = sd(AIS_Area_corrected, na.rm = TRUE),
    mean_AIS_Area = mean(AIS_Area, na.rm = TRUE),
    sd_AIS_Area = sd(AIS_Area, na.rm = TRUE)
  )



write.csv(data.na.omit_popdata_summary,
          file.path(mainDir,subDir,saveFolder, "data.na.omit_popdata_summary.csv"), row.names = TRUE)
cat(red("Difference in MN characteristics by muscle \n"))
print(data.na.omit_popdata_summary)


sink("/dev/null")
setwd(file.path(mainDir,subDir,saveFolder))
lapply(1:length(popdata_summary_manova), function(i) write.table(popdata_summary_manova[[i]], 
                                                                 file = paste0("popdata_summary_manova",names(popdata_summary_manova[i]), ".csv"),
                                                                 sep=',',
                                                                 row.names = FALSE))
sink()

#
#
# write.csv(tidyMCMC(post1, conf.int = TRUE, conf.method = "HPDinterval", conf.level = 0.95),
#           file.path(mainDir,subDir,saveFolder, "pop.data.rstan_model_final.csv"), row.names = TRUE)
#
# cat(red("Regression coeffs by motor pool  \n"))
# tidyMCMC(post1, conf.int = TRUE, conf.method = "HPDinterval", conf.level = 0.95)

# suppressWarnings(data.na.omit_popdata_ppd_summary<-post1 %>%
# emmeans(~Muscle) %>%
# gather_emmeans_draws() %>%
# mean_qi()
# )

# write.csv(data.na.omit_popdata_ppd_summary,
#           file.path(mainDir,subDir,saveFolder, "data.na.omit_popdata_ppd_summary.csv"), row.names = TRUE)
# cat(red("95% PPD of Dis-to-soma in each motor pool  \n"))
# print(data.na.omit_popdata_ppd_summary)
#




suppressWarnings(data.na.omit_popdata_ppd_summary<-m3 %>%
                   spread_draws(MG, SOL) %>%
                   gather(group, b, MG:SOL, factor_key=TRUE)%>%
                   group_by(group) %>%
                   median_qi(condition_mean = b, .width = c(.95, .66))
)

write.csv(data.na.omit_popdata_ppd_summary,
          file.path(mainDir,subDir,saveFolder, "data.na.omit_popdata_ppd_summary.csv"), row.names = TRUE)
cat(red("95% and 66% PPD of Soma_Diam in each motor pool  \n"))
print(data.na.omit_popdata_ppd_summary)



## alternative ppd of estimates
# broom::tidy(m3, pars = c("MG", "SOL"),
#             conf.int = TRUE, conf.method = "quantile")

## formula for calculaitng conditional means from mu_1 and mu_2 is in the generated quantities block in model and references in 'ais_modeling_testing.R'
# generated quantities{
#   real condition_A = ndt_1 + exp(mu_1 + 0.5*sigma^2);
#   real condition_B = ndt_2 + exp(mu_2 + 0.5*sigma^2);
# }

########################### saving figures ###########################


setwd(file.path(mainDir,subDir,figFolder))


### graphing options for distributions
#Color by groups without histogram

cairo_ps("observed_data_box_whisker.eps")
ggplot(na.omit_popdata, aes(x = Muscle, y=Dis_2_Soma, fill=Muscle))+
  geom_boxplot(outlier.shape = NA) +
  scale_color_manual(values = c("#868686FF", "#EFC000FF"))+
  theme_classic()
invisible(suppressMessages(suppressWarnings(dev.off())))


# yreps 'generative'
cairo_ps("ppc_dens_overlay.y1.eps")
# lims <- ggplot2::lims(x = c(0, 50), y = c(0,.25))
# ppc1<-ppc_dens_overlay(y1, yrep = y1rep)
# ppc1+lims

ppc_dens_overlay(y1, yrep = y1rep)

invisible(suppressMessages(suppressWarnings(dev.off())))

cairo_ps("ppc_dens_overlay.y2.eps")
# ppc2<-ppc_dens_overlay(y2, yrep = y2rep)
# ppc2+lims

ppc_dens_overlay(y2, yrep = y2rep)

invisible(suppressMessages(suppressWarnings(dev.off())))


# Posterior density
cairo_ps("ppd_predicted_1.eps")
bayesplot::mcmc_areas(post_sam, pars = c("MG", "SOL"),
                      prob = .95)
invisible(suppressMessages(suppressWarnings(dev.off())))

cairo_ps("ppd_predicted_sigma_1.eps")
bayesplot::mcmc_areas(post_sam, pars = c("sigma", "sigma_1"),
                      prob = .95)
invisible(suppressMessages(suppressWarnings(dev.off())))

cairo_ps("ppd_predicted_2.eps")
m3 %>%
  spread_draws(MG, SOL) %>%
  gather(group, b, MG:SOL, factor_key=TRUE)%>%
  group_by(group) %>%
  # mutate(condition_mean = `(Intercept)` + b) %>%
  ggplot(aes(y = group, x = b)) +
  stat_halfeye()
invisible(suppressMessages(suppressWarnings(dev.off())))

cairo_ps("ppd_predicted_sigma_2.eps")
m3 %>%
  spread_draws(sigma, sigma_1) %>%
  gather(group, b, sigma:sigma_1, factor_key=TRUE)%>%
  group_by(group) %>%
  # mutate(condition_mean = `(Intercept)` + b) %>%
  ggplot(aes(y = group, x = b)) +
  stat_halfeye()
invisible(suppressMessages(suppressWarnings(dev.off())))


# data density
cairo_ps("observed_data_distributions_1.eps")
## https://stackoverflow.com/questions/58837773/pivot-wider-issue-values-in-values-from-are-not-uniquely-identified-output-w
na.omit_popdata %>%
  select(Muscle, Dis_2_Soma) %>%
  group_by(Muscle) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = Muscle, values_from = Dis_2_Soma) %>%
  select(-row) %>%
  na.omit()%>%
  bayesplot::mcmc_areas(pars = c("MG", "SOL"),
                        prob = .95)
invisible(suppressMessages(suppressWarnings(dev.off())))



cairo_ps("observed_data_distributions_2.eps")
na.omit_popdata %>%
  group_by(Muscle) %>%
  # mutate(condition_mean = `(Intercept)` + b) %>%
  ggplot(aes(y = Muscle, x = Dis_2_Soma)) +
  stat_halfeye()
invisible(suppressMessages(suppressWarnings(dev.off())))

cat(green("Pop data Figs successfully created",
          "Find them here", file.path(mainDir,subDir,figFolder),   sep="\n"))


## if you just want to have the lines and no densities

# m3 %>%
#   spread_draws(MG, SOL) %>%
#   gather(group, b, MG:SOL, factor_key=TRUE)%>%
#   group_by(group) %>%
#   median_qi(condition_mean = b, .width = c(.95, .66)) %>%
#   ggplot(aes(y = group, x = condition_mean, xmin = .lower, xmax = .upper)) +
#   geom_pointinterval()
########################### Clean up ###########################



#################################################################################
########################### Rheo~Conduc all cells ###########################
#################################################################################

########################### load dependencies ###########################
data_Name_full_Rheo_Conduc<- 'Rin_vs_Rheo_MNs_Combined.xls'   ### after the input resistance and rheo have be reconciled 
mainDir <- "~/Desktop/plos_one_Rotterman_2021" ### set this filepath to the main directory 
subDir <- "final_data" ## this is the subdirectory where data is held
figFolder <-"figFolder" ## this will be created if not already in existance
saveFolder <- "saveFolder" ## this will be created if not already in existance

########################### Custom Functions ###########################
########################### Load Data ###########################

full_Rheo_Conduc<-read_excel(file.path(mainDir,subDir,data_Name_full_Rheo_Conduc), 
                             na = "NA")

########################### Data Wrangling ###########################

# Column 3 - Muscle: 1 = MG, 2 = Triceps Surae 
# Column 4 - Who: 2 = Paul, 3 = Travis, 4 = Neurobiotin, 5 = Ahmed 
full_Rheo_Conduc$Who<- as.factor(full_Rheo_Conduc$Who)
full_Rheo_Conduc$Who_collapse<- as.factor(full_Rheo_Conduc$Who_collapse)
full_Rheo_Conduc$Muscle<- as.factor(full_Rheo_Conduc$Muscle)

# collapse all previous experimenters down to one factor (== 0) and call out AIS MN as (== 1)
# full_Rheo_Conduc$Who_collapse<-ifelse(full_Rheo_Conduc$Who == 4,1,0)
# full_Rheo_Conduc$Who_collapse<- as.factor(full_Rheo_Conduc$Who_collapse)

scaled.data<-full_Rheo_Conduc%>%
  # filter(Who_collapse != 0) %>%      ### change to '0' if you want the new cells and '1' to the old cells
  select(Who,
         Who_collapse,
         Rheo,	
         Input_Resistance,
         Conductance
         
  )
########################### quick visuakization ###########################
full_Rheo_Conduc%>%
  # filter(Who != 4) %>%
  select(Rheo,	
         Input_Resistance,
         Conductance,
         Who_collapse
  ) %>% 
  ggscatter(x = "Conductance", 
            y = "Rheo", 
            color = "Who_collapse",
            palette = c("#1f75FE", 
                        "#4CB7A5", 
                        # "#926EAE",
                        "#EE204D"),
            add = "reg.line",
            conf.int = TRUE)+
  stat_cor(aes(color = Who_collapse))        # Add correlation coefficient
# stat_cor()        # Add correlation coefficient

########################### analyses/modeling ###########################

post1 = stan_glm(Rheo ~ Conductance, data = scaled.data)
post2 = stan_glm(Rheo ~ Conductance + Who_collapse, data = scaled.data)
post3 = stan_glm(Rheo ~ Conductance, data = scaled.data)


loo1 <- loo(post1, save_psis = TRUE, cores = 4) #, k_threshold =0.7)
loo2 <- loo(post2, save_psis = TRUE, cores = 4) #, k_threshold =0.7)
loo3 <- loo(post3, save_psis = TRUE, cores = 4) #, k_threshold =0.7)
loo_compare(loo1,loo2)


### post 3 is the best model to compare
mcmc <- as.matrix(post3)
Xmat = model.matrix(~Conductance, scaled.data)
coefs = mcmc[, c("(Intercept)", "Conductance")]
fit = coefs %*% t(Xmat)
resid = sweep(fit, 2, scaled.data$Rheo, "-")
var_f = apply(fit, 1, var)
var_e = apply(resid, 1, var)
R2 = var_f/(var_f + var_e)
tidyMCMC(as.mcmc(R2), conf.int = TRUE, conf.method = "HPDinterval")
sqrt(median(R2))
summary(lm(Rheo ~ Conductance , data = scaled.data))


########################### saving data ###########################


write.csv(tidyMCMC(post3, conf.int = TRUE, conf.method = "HPDinterval", conf.level = 0.95),
          file.path(mainDir,subDir,saveFolder, "data.rstan_model_final_Rheo~Conduc.csv"), row.names = TRUE)
cat(red("Bayesian Parameter estimates Rheo~Conduc ) \n"))
tidyMCMC(post3, conf.int = TRUE, conf.method = "HPDinterval", conf.level = 0.95)


write.csv(tidyMCMC(as.mcmc(R2), conf.int = TRUE, conf.method = "HPDinterval", conf.level = 0.95),
          file.path(mainDir,subDir,saveFolder, "data.rstan_model_final_adj_r2_pop.csv"), row.names = TRUE)
cat(red("Bayesian Parameter estimates R^2 \n"))
tidyMCMC(as.mcmc(R2), conf.int = TRUE, conf.method = "HPDinterval", conf.level = 0.95)



########################### saving figures ###########################
setwd(file.path(mainDir,subDir,figFolder))

plot_1<-full_Rheo_Conduc%>%
  filter(Who != 4) %>%
  select(Rheo,	
         Input_Resistance,
         Conductance,
         Who
  ) %>% 
  ggscatter(x = "Conductance", 
            y = "Rheo", 
            color = "Who",
            palette = c("#1f75FE", 
                        "#4CB7A5", 
                        # "#926EAE",
                        "#EE204D"),
            add = "reg.line",
            conf.int = TRUE)+
  stat_cor(aes(color = Who))        # Add correlation coefficient
# stat_cor()        # Add correlation coefficient

plot_2<-full_Rheo_Conduc%>%
  # filter(Who != 4) %>%
  select(Rheo,	
         Input_Resistance,
         Conductance,
         Who
  ) %>% 
  ggscatter(x = "Conductance", 
            y = "Rheo", 
            color = "Who",
            palette = c("#1f75FE", 
                        "#4CB7A5", 
                        "#926EAE",
                        "#EE204D"),
            add = "reg.line",
            conf.int = TRUE)+
  stat_cor(aes(color = Who))    

plot_3<-full_Rheo_Conduc%>%
  # filter(Who != 4) %>%
  select(Rheo,	
         Input_Resistance,
         Conductance,
         Who_collapse
  ) %>% 
  ggscatter(x = "Conductance", 
            y = "Rheo", 
            # color = "Who",
            palette = c("#1f75FE", 
                        "#4CB7A5", 
                        "#926EAE",
                        "#EE204D"),
            add = "reg.line",
            conf.int = TRUE)+
  stat_cor()  

cairo_ps("observed_rheo~conduc_3.eps")
plot_1
invisible(suppressMessages(suppressWarnings(dev.off())))


cairo_ps("observed_rheo~conduc_4.eps")
plot_2
invisible(suppressMessages(suppressWarnings(dev.off())))


cairo_ps("observed_rheo~conduc_group.eps")
plot_3
invisible(suppressMessages(suppressWarnings(dev.off())))



####


plot_4<-scaled.data %>%
  group_by(Who_collapse) %>%
  data_grid(Conductance = seq_range(Conductance, n = 51)) %>%
  add_fitted_draws(post1) %>%
  ggplot(aes(x = Conductance, y = Rheo )) +
  stat_lineribbon(aes(y = .value)) +
  geom_point(data = scaled.data) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")




plot_5<-scaled.data %>%
  group_by(Who_collapse) %>%
  data_grid(Conductance = seq_range(Conductance, n = 101)) %>%
  add_fitted_draws(post2, n = 100) %>%
  ggplot(aes(x = Conductance, y = Rheo, color = ordered(Who_collapse))) +
  geom_line(aes(y = .value, group = paste(Who_collapse, .draw)),alpha = .1) +
  geom_point(data = scaled.data) +
  scale_color_brewer(palette = "Dark2")+
  theme(legend.title = element_blank())


# The palette with black:
cbbPalette <- c("#000000", "#bebebe")
cbPalette <- c("#3366FF", "#999999")


plot_6<-scaled.data %>%
  group_by(Who_collapse) %>%
  data_grid(Conductance = seq_range(Conductance, n = 101)) %>%
  add_fitted_draws(post2, n = 100) %>%
  ggplot(aes(x = Conductance, y = Rheo, color = ordered(Who_collapse))) +
  geom_line(aes(y = .value, group = paste(Who_collapse, .draw)),alpha = .1) +
  geom_point(data = scaled.data) +
  scale_colour_manual(values=cbPalette)+
  theme_classic()+
  theme(legend.title = element_blank())

  scale_color_brewer(palette = "Dark2") 


cairo_ps("BayesModel_rheo~conduc_3.eps")
plot_4
invisible(suppressMessages(suppressWarnings(dev.off())))


cairo_ps("BayesModel_rheo~conduc+who_4.eps")
plot_5
invisible(suppressMessages(suppressWarnings(dev.off())))


cairo_ps("BayesModel_rheo~conduc*who_4.eps")
plot_6
invisible(suppressMessages(suppressWarnings(dev.off())))


########################### Clean up ###########################




#################################################################################
########################### Filled cells bayes correlations ###########################
#################################################################################

########################### load dependencies ###########################
########################### Custom Functions ###########################
########################### Load Data ###########################
########################### Data Wrangling ###########################
########################### quick visuakization ###########################
########################### analyses/modeling ###########################
########################### saving data ###########################
########################### saving figures ###########################
########################### Clean up ###########################

########################### References ###########################

# https://www.r-bloggers.com/multiple-regression-lines-in-ggpairs/

# https://stackoverflow.com/questions/17175190/trouble-with-cbind-in-manova-call-in-r

# https://www.r-project.org/nosvn/pandoc/crayon.html


########################### timer stop ###########################
proc.time() - ptm




## TEMPLATE

#################################################################################
########################### Filled cells ###########################
#################################################################################

########################### load dependencies ###########################
########################### Custom Functions ###########################
########################### Load Data ###########################
########################### Data Wrangling ###########################
########################### quick visuakization ###########################
########################### analyses/modeling ###########################
########################### saving data ###########################
########################### saving figures ###########################
########################### Clean up ###########################


## TEMPLATE



