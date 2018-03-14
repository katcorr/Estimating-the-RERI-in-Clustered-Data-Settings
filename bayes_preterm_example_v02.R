#################################################################
### Kat Correia                                              ####
### 04.24.2017                                               ####
### bayesian analysis of preterm delivery application        ####
### v02 (03.05.2018) -- specifying function for initial values ##
###             to aid sampling; and setting seed so can       ##
###             reproduce results                              ##
#################################################################

current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

library("rstan")
library("shinystan")
library("Hmisc")
library("brms")
library("ggplot2")
library("dplyr")
library("mosaic")
library("geepack")
library("mvtnorm")

set.seed(1515)

arvs <- sasxport.get("arvs.xpt")

tally(~instn1, data=arvs)


stancode.1a <- make_stancode(preterm ~ nvpconc + lowfcd4 + nvpconc:lowfcd4
                             + (1|instn1)
                             , family=bernoulli(link="logit")
                             , data=arvs
                             , prior = set_prior("cauchy(0,5)", class = "sd")
                             , save_model="bayes_logbin_cauchy_hivapp.stan")

stancode.1b <- make_stancode(preterm ~ nvpconc + lowfcd4 + nvpconc:lowfcd4
                         + (1|instn1)
                         , family=bernoulli(link="logit")
                         , data=arvs
                         , prior = set_prior("gamma(2,0.1)", class = "sd")
                         , save_model="bayes_logbin_gamma_hivapp.stan")


stancode.2a <- make_stancode(preterm ~ nvpconc + lowfcd4 + nvpconc:lowfcd4
                             + black + agegp2 + agegp3
                             + (1|instn1)
                             , family=bernoulli(link="logit")
                             , data=arvs
                             , prior = set_prior("cauchy(0,5)", class = "sd")
                             , save_model="bayes_logbin_cauchy2_hivapp.stan")

stancode.2b <- make_stancode(preterm ~ nvpconc + lowfcd4 + nvpconc:lowfcd4
                             + black + agegp2 + agegp3
                             + (1|instn1)
                             , family=bernoulli(link="logit")
                             , data=arvs
                             , prior = set_prior("gamma(2,0.1)", class = "sd")
                             , save_model="bayes_logbin_gamma2_hivapp.stan")

standat1 <- make_standata(preterm ~ nvpconc + lowfcd4 + nvpconc:lowfcd4 
                          + (1|instn1)
                          , data=arvs)

standat2 <- make_standata(preterm ~ nvpconc + lowfcd4 + nvpconc:lowfcd4
                          + black + agegp2 + agegp3
                          + (1|instn1)
                         , data=arvs)

###################################################################
####      fit random intercept BAYESIAN log binomial model    #####
###################################################################

#######################
## UNADJUSTED MODELS ##
#######################

# unadjusted with half-cauchy prior
stfit1.a <- stan(file="bayes_logbin_cauchy_hivapp.stan"
               , data=standat1
               , chains=4, iter=2000
               , seed=20151031)
RERI.results1.a <- (summary(stfit1.a))[[1]][51,c(1,4,8,10)]
SD.results1.a <- (summary(stfit1.a))[[1]][5,c(1,4,8,10)]

chains1.a <- as.data.frame(stfit1.a, pars = c("sd_1[1]", "RERI"))
colnames(chains1.a)[1] <- "SD"
acf(chains1.a[3001:4000,1])
acf(chains1.a[3001:4000,2])
gf_density(~SD,data=chains1.a)
gf_density(~RERI,data=chains1.a)

# get initial values used for each chain
get_inits(stfit1.a)
# get acceptance rates for each chain
sampler_params <- get_sampler_params(stfit1.a, inc_warmup = FALSE)
mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
mean_accept_stat_by_chain
# plots
chains1.a$draw <- rep(1:1000,4)
chains1.a$chain <- rep(1:4,each=1000)
gf_line(RERI~draw,color=~chain,data=chains1.a)
gf_line(SD~draw,color=~chain,data=chains1.a)

# unadjusted with gamma prior
stfit1.b <- stan(file="bayes_logbin_gamma_hivapp.stan"
                  , data=standat1
                  , chains=4, iter=2000
                  , seed=19850503)
RERI.results1.b <- (summary(stfit1.b))[[1]][51,c(1,4,8,10)]
SD.results1.b <- (summary(stfit1.b))[[1]][5,c(1,4,8,10)]

chains1.b <- as.data.frame(stfit1.b, pars = c("sd_1[1]", "RERI"))
colnames(chains1.b)[1] <- "SD"
acf(chains1.b[1001:2000,1])
acf(chains1.b[1:1000,2])
gf_density(~SD,color="red",fill="red",data=chains1.b) %>%
  gf_density(~SD,color="darkblue",fill="darkblue"
              , alpha=0.25,data=chains1.a)
gf_density(~RERI,color="red",fill="red",data=chains1.b) %>%  
  gf_density(~RERI,color="darkblue",fill="darkblue"
              , alpha=0.25, data=chains1.a) 

# get initial values used for each chain
get_inits(stfit1.b)
# get acceptance rates for each chain
sampler_params <- get_sampler_params(stfit1.b, inc_warmup = FALSE)
mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
mean_accept_stat_by_chain
# plots
chains1.b$draw <- rep(1:1000,4)
chains1.b$chain <- rep(1:4,each=1000)
gf_line(RERI~draw,color=~chain,data=chains1.b)
gf_line(SD~draw,color=~chain,data=chains1.b)

#####################
## ADJUSTED MODELS ##
#####################

# for adjusted models, specify initial values for betas based off 
# marginal Poisson model
# for SD, specify initial values based off unadjusted model's
# SD estimates

# need to order variables as Stan data orders them (interaction last!)
geefit <- geeglm(formula = preterm ~ nvpconc + lowfcd4 + black
                 + agegp2 + agegp3 + nvpconc:lowfcd4
                 , id=instn1, data=arvs
                 , corstr="exchangeable", family="poisson"(link="log"))
summary(geefit)

geefit.est <- coef(geefit)[2:7]
geefit.cov.mat <- geefit$geese$vbeta[2:7,2:7]


init_fun <- function() { 
    bvec <- rmvnorm(n=1,mean=geefit.est,sigma=geefit.cov.mat)
    sd <- runif(n=1,min=0.107,max=0.350) # 95% credible interval from unadjusted random int. model 
    list(b=c(bvec[1], bvec[2], bvec[3], bvec[4], bvec[5], bvec[6])
         , sd_1=array(sd,dim=length(sd))) 
} 
init_fun()
init_fun()

# adjusted with half-cauchy prior
stfit2.a <- stan(file="bayes_logbin_cauchy2_hivapp.stan"
                 , data=standat2
                 , chains=4, iter=2000
                 , init = init_fun
                 , seed=19820422)
summary(stfit2.a)
RERI.results2.a <- (summary(stfit2.a))[[1]][54,c(1,4,8,10)]
SD.results2.a <- (summary(stfit2.a))[[1]][8,c(1,4,8,10)]
RERI.results2.a
SD.results2.a

chains2.a <- as.data.frame(stfit2.a, pars = c("sd_1[1]", "RERI"))
colnames(chains2.a)[1] <- "SD"
acf(chains2.a[1001:2000,1])
acf(chains2.a[1:1000,2])
gf_density(~SD,color="red",fill="red",data=chains2.a) 
gf_density(~RERI,color="red",fill="red",data=chains2.a) 

# get initial values used for each chain
get_inits(stfit2.a)
# get acceptance rates for each chain
sampler_params <- get_sampler_params(stfit2.a, inc_warmup = FALSE)
mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
mean_accept_stat_by_chain
# plots
chains2.a$draw <- rep(1:1000,4)
chains2.a$chain <- rep(1:4,each=1000)
gf_line(RERI~draw,color=~chain,data=chains2.a)
gf_line(SD~draw,color=~chain,data=chains2.a)

#Rhat for SD is 1.06, poor chain convergence ... try increasing samples
stfit2.ab <- stan(file="bayes_logbin_cauchy2_hivapp.stan"
                 , data=standat2
                 , chains=4, iter=10000
                 , init=init_fun
                 , seed=18204002)
RERI.results2.ab <- (summary(stfit2.ab))[[1]][54,c(1,4,8,10)]
SD.results2.ab <- (summary(stfit2.ab))[[1]][8,c(1,4,8,10)]
RERI.results2.ab
SD.results2.ab

chains2.ab <- as.data.frame(stfit2.ab, pars = c("sd_1[1]", "RERI"))
colnames(chains2.ab)[1] <- "SD"
acf(chains2.ab[1:1000,1])
acf(chains2.ab[1:1000,2])
gf_density(~SD,color="red",fill="red",data=chains2.ab) 
gf_density(~RERI,color="red",fill="red",data=chains2.ab) 

# get initial values used for each chain
get_inits(stfit2.ab)
# get acceptance rates for each chain
sampler_params <- get_sampler_params(stfit2.a, inc_warmup = FALSE)
mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
mean_accept_stat_by_chain
# plots
chains2.ab$draw <- rep(1:1000,4)
chains2.ab$chain <- rep(1:4,each=1000)
gf_line(RERI~draw,color=~chain,data=chains2.ab)
gf_line(SD~draw,color=~chain,data=chains2.ab)

RERI.results2.ab
mean(chains2.ab$RERI[1:10000*100],na.rm=TRUE)
SD.results2.ab
mean(chains2.ab$SD[1:10000*100],na.rm=TRUE)

# adjusted with gamma(2,0.1) prior
stfit2.b <- stan(file="bayes_logbin_gamma2_hivapp.stan"
                 , data=standat2
                 , chains=4, iter=2000
                 , init=init_fun
                 , seed=19820422)
summary(stfit2.b)
RERI.results2.b <- (summary(stfit2.b))[[1]][54,c(1,4,8,10)]
SD.results2.b <- (summary(stfit2.b))[[1]][8,c(1,4,8,10)]
RERI.results2.b
SD.results2.b

chains2.b <- as.data.frame(stfit2.b, pars = c("sd_1[1]", "RERI"))
colnames(chains2.b)[1] <- "SD"
acf(chains2.b[3001:4000,1])
acf(chains2.b[3001:4000,2])
gf_density(~SD,color="red",fill="red",data=chains2.b) 
gf_density(~RERI,color="red",fill="red",data=chains2.b) 

# get initial values used for each chain
get_inits(stfit2.b)
# get acceptance rates for each chain
sampler_params <- get_sampler_params(stfit2.b, inc_warmup = FALSE)
mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
mean_accept_stat_by_chain
# plots
chains2.b$draw <- rep(1:1000,4)
chains2.b$chain <- rep(1:4,each=1000)
gf_line(RERI~draw,color=~chain,data=chains2.b)
gf_line(SD~draw,color=~chain,data=chains2.b)

#Rhat for SD is 1.12, poor chain convergence ... try increasing samples
stfit2.bb <- stan(file="bayes_logbin_gamma2_hivapp.stan"
                 , data=standat2
                 , chains=4, iter=10000
                 , init=init_fun
                 , seed=19820422)
RERI.results2.bb <- (summary(stfit2.bb))[[1]][54,c(1,4,8,10)]
SD.results2.bb <- (summary(stfit2.bb))[[1]][8,c(1,4,8,10)]
RERI.results2.bb
SD.results2.bb

chains2.bb <- as.data.frame(stfit2.bb, pars = c("sd_1[1]", "RERI"))
colnames(chains2.bb)[1] <- "SD"
acf(chains2.bb[1:1000,1])
acf(chains2.bb[1:1000,2])
gf_density(~SD,color="red",fill="red",data=chains2.bb) %>%
  gf_density(~SD,color="darkblue",fill="darkblue"
             , alpha=0.25,data=chains2.ab) 
gf_density(~RERI,color="red",fill="red",data=chains2.bb) %>%
  gf_density(~RERI,color="darkblue",fill="darkblue",data=chains2.ab) 

# get initial values used for each chain
get_inits(stfit2.bb)
# get acceptance rates for each chain
sampler_params <- get_sampler_params(stfit2.bb, inc_warmup = FALSE)
mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
mean_accept_stat_by_chain
# plots
chains2.bb$draw <- rep(1:1000,4)
chains2.bb$chain <- rep(1:4,each=1000)
gf_line(RERI~draw,color=~chain,data=chains2.bb)
gf_line(SD~draw,color=~chain,alpha=0.25,data=chains2.bb)
