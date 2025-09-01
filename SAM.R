#########################################################################
# SAM prior for binary endpoints
#
# https://cran.r-project.org/web/packages/SAMprior/vignettes/Example_binary.html
#
# SAM: Self-adapting mixture prior to dynamically borrow information from historical data in clinical trials
# (Yang et. al 2023)
#
# CL 28/08/2025
#########################################################################

remove(list = ls())

library(ggplot2)
library(RBesT)
library(SAMprior)
library(knitr)


theme_set(theme_bw()) # sets up plotting theme
set.seed(22)

### SAM prior derivation

map_AS <- gMAP(cbind(r, n-r) ~ 1 | study,
                   family = binomial,
                   #data = ASAS20, 
                   data = AS,
                   tau.dist = "HalfNormal", 
                   tau.prior = 1,
                   beta.prior = 2)

map_automix <- automixfit(map_ASAS20)
map_automix

plot(map_automix)$mix

### SAM weight determination

n <- 35
r = 10 

wSAM <- SAM_weight(if.prior = map_automix, 
                   delta = 0.2,
                   n = n, r = r)
cat('SAM weight: ', wSAM)

wSAM <- SAM_weight(if.prior = map_automix, 
                   delta = 0.2,
                   method.w = 'PPR',
                   prior.odds = 3/7,
                   n = n, r = r)
cat('SAM weight: ', wSAM)

### SAM prior construction

SAM.prior <- SAM_prior(if.prior = map_automix, 
                       nf.prior = mixbeta(nf.prior = c(1,1,1)),
                       weight = wSAM)
SAM.prior

### Operating Characteristics: type 1 error


set.seed(123)
TypeI <- get_OC(if.prior = map_automix,       ## MAP prior from historical data
                nf.prior = mixbeta(c(1,1,1)), ## Non-informative prior for treatment arm
                delta    = 0.2,               ## CSD for SAM prior
                ## Method to determine the mixture weight for the SAM prior
                method.w = 'LRT',             
                n        = 35, n.t = 70,      ## Sample size for control and treatment arms
                ## Decisions
                decision = decision2S(0.95, 0, lower.tail=FALSE), 
                ntrial   = 1000,              ## Number of trials simulated
                if.MAP   = TRUE,              ## Output robust MAP prior for comparison
                weight   = 0.5,               ## Weight for robust MAP prior
                ## Response rates for control and treatment arms
                theta    = c(0.36, 0.36, 0.11, 0.55),
                theta.t  = c(0.34, 0.33, 0.11, 0.55)
)

kable(TypeI)

## Operating characteristics: power


Power <- get_OC(if.prior = map_automix,       ## MAP prior from historical data
                nf.prior = mixbeta(c(1,1,1)), ## Non-informative prior for treatment arm
                delta    = 0.2,               ## CSD for SAM prior
                n        = 35, n.t = 70,      ## Sample size for control and treatment arms
                ## Decisions
                decision = decision2S(0.95, 0, lower.tail=FALSE), 
                ntrial   = 1000,              ## Number of trials simulated
                if.MAP   = TRUE,           ## Output robust MAP prior for comparison
                weight   = 0.5,               ## Weight for robust MAP prior
                ## Response rates for control and treatment arms
                theta    = c(0.37, 0.34, 0.16, 0.11),
                theta.t  = c(0.57, 0.54, 0.36, 0.31)
)
kable(Power)


### Decision making

## Sample size and number of responses for treatment arm
n_t <- 70; x_t <- 22 

## first obtain posterior distributions...
post_SAM <- postmix(priormix = SAM.prior,         ## SAM Prior
                    r = r,   n = n)
post_trt <- postmix(priormix = mixbeta(c(1,1,1)), ## Non-informative prior
                    r = x_t, n = n_t)

## Define the decision function
decision = decision2S(0.95, 0, lower.tail=FALSE)

## Decision-making
decision(post_trt, post_SAM)
