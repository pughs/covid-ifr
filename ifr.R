library(rstan)
library(splines)
library(cmdstanr)

# Read in and format data ----

#### Read in data ----

studydatifr <- data.frame(readxl::read_xlsx(
  "./IFR Data.xlsx",sheet="Seroprevalence Study Summary"))
seroprevdat <- data.frame(readxl::read_xlsx(
  "./IFR Data.xlsx",sheet="Seroprevalence Data"))
deathdat <- data.frame(readxl::read_xlsx(
  "./IFR Data.xlsx",sheet="Death Data"))
testkitdats <- data.frame(readxl::read_xlsx(
  "./IFR Data.xlsx",sheet="Test Assay Data"))
popdat <- data.frame(readxl::read_xlsx(
  "./IFR Data.xlsx",sheet="Age Distribution"))
natdat <- data.frame(readxl::read_xlsx(
  "./IFR Data.xlsx",sheet="National Age Distribution"))

seroprevdat$location_id <- factor(seroprevdat$location_id, studydatifr$location_id)
deathdat$location_id <- factor(deathdat$location_id, studydatifr$location_id)


#### Expand population using national data ----

ages <- 0:99
num_age <- length(ages)
num_loc <- length(unique(seroprevdat$location_id))

ifrlocs <- studydatifr$location_id

expandNatDat <- function(subdat1) {
  ## Must be 5 year age bins 
  
  expandedf <- numeric(length(0:84)+1)
  
  subdat1$percent_of_pop <- subdat1$percent_of_pop/sum(subdat1$percent_of_pop)
  
  ## Repeat those less than 85 uniformly
  subdat <- subdat1[subdat1$lower_age_bin < 85,]
  for (j in 1:(nrow(subdat))) {
    up <- min(subdat$upper_age_bin[j], 84)
    expandedf[subdat$lower_age_bin[j]:up+1] <- subdat$percent_of_pop[j]/
      (up-subdat$lower_age_bin[j]+1)
  }
  
  ## Make the 85+ bin
  expandedf[length(0:84)+1] <- sum(subdat1$percent_of_pop[subdat$lower_age>=85,])
  
  return(expandedf)
}

locind <- 1
countryind <- 1

expandedf <- matrix(NA,86, length(ifrlocs))
expandedNat <- matrix(NA,86, length(ifrlocs))

for (i in 1:length(ifrlocs)) {
  subdat <- popdat[popdat$location_id==ifrlocs[i],]
  subdat$percent_of_pop <- subdat$percent_of_pop/sum(subdat$percent_of_pop)
  nat <- natdat[natdat$Country==studydatifr$country[i],]
  
  ## expand national
  nat_expanded <- c(rep(as.numeric(nat[-c(1:2,ncol(nat))])/5, each=5), as.numeric(nat[ncol(nat)]))
  nat_expanded <- nat_expanded/sum(nat_expanded)
  expandedNat[,i] <- nat_expanded
  
  ## expand others wrt national
  binsizes <- subdat$upper_age_bin-subdat$lower_age_bin
  endd <- max(which(subdat$upper_age_bin<85))
  
  for (j in 1:endd) {
    if (binsizes[j]<=4) {
      expandedf[subdat$lower_age_bin[j]:subdat$upper_age_bin[j]+1,i] <- subdat$percent_of_pop[j]/
        (subdat$upper_age_bin[j]-subdat$lower_age_bin[j]+1)
    } else {
      expandedf[subdat$lower_age_bin[j]:subdat$upper_age_bin[j]+1,i] <- subdat$percent_of_pop[j]*
        nat_expanded[subdat$lower_age_bin[j]:subdat$upper_age_bin[j]+1]/
        sum(nat_expanded[subdat$lower_age_bin[j]:subdat$upper_age_bin[j]+1])
    }
  }
  
  if (any(subdat$lower_age_bin==85)) {
    expandedf[86, i] <- sum(subdat$percent_of_pop[subdat$lower_age_bin>=85])
  } else {
    ind <- min(which(subdat$upper_age_bin>85))
    natOver85 <- nat_expanded[length(nat_expanded)]/
      sum(nat_expanded[subdat$lower_age_bin[ind]:length(nat_expanded)])
    over85 <- sum(subdat$percent_of_pop[ind:nrow(subdat)])*natOver85
    under85 <- sum(subdat$percent_of_pop[ind:nrow(subdat)])*(1-natOver85)
    expandedf[86, i] <- over85
    expandedf[subdat$lower_age_bin[ind]:84 +1, i] <- under85*nat_expanded[subdat$lower_age_bin[ind]:84+1]/
      sum(nat_expanded[subdat$lower_age_bin[ind]:84+1])
  }
}

expandedfBy1yr <- expandedf

h <- 0.5 # integration width

ages <- seq(0, 100, by=h)

expandedf <- apply(expandedf, 2, function(x) {
  rep(x*h, each=1/h)
})


## Repeat for those >85 for expandedf ##
junk <- expandedf[nrow(expandedf),]
numUnif <- length(85:99)*1/h+1
expandedf[(nrow(expandedf)-1/h+1):nrow(expandedf),] <- matrix(junk/numUnif/h, 1/h, num_loc, byrow=T)
expandedf <- rbind(expandedf, matrix(junk/numUnif/h, 1/h*14+1, num_loc, byrow=T)) 
# +1 is to get the 100+ in one instead of expanded

## Repeat for >85 for exandedfBy1yr ##
junk <- expandedfBy1yr[nrow(expandedfBy1yr),]
numUnif <- length(85:99)*1/1+1
expandedfBy1yr[(nrow(expandedfBy1yr)-1/1+1):nrow(expandedfBy1yr),] <- matrix(junk/numUnif/1, 1/1, num_loc, byrow=T)
expandedfBy1yr <- rbind(expandedfBy1yr, matrix(junk/numUnif/1, 1/1*14+1, num_loc, byrow=T)) 
# +1 is to get the 100+ in one instead of expanded

#### Smooth population age distribution ----

fineages <- seq(0,100, by=.5)

getExpandedfNew <- function(i) {
  dat <- data.frame(y=expandedfBy1yr[,i]/2, x=0:100)
  dat <- dat[dat$x<85,]
  dat <- rbind(dat,c(0,100))
  preddat <- data.frame(y=1, x=ages)
  preddat2 <- data.frame(y=1, x=fineages)
  # fit <- loess(y~x,dat, span=.1) 
  fit <- loess(y~x,dat, span=.3) 
  preds <- predict(fit, preddat)
  preds2 <- predict(fit, preddat2)
  
  plot(dat$x, dat$y, xlab="age", ylab="density", main=
         paste(studydatifr$location_id[i], studydatifr$country[i],
               studydatifr$location[i]))
  lines(ages, preds)
  # lines(seq(0,100, by=.001), preds2, col="red", lty=2)
  if (any(preds<0)) {
    preds2 <- preds2-min(preds)+.000001
    preds <- preds-min(preds)+.000001
  }
  preds2 <- preds2/sum(preds)
  preds <- preds/sum(preds)
  lines(ages, preds, col="red")
  lines(fineages, preds2, col="blue", lty=2)
  sum(preds)
  
  return(list(coarse=preds, fine=preds2))
}

out <- sapply(1:nrow(studydatifr), getExpandedfNew, simplify=F)
newexpandedf <- sapply(out, function(x) x$coarse)
fineexpandedf <- sapply(out, function(x) x$fine)

expandedf <- newexpandedf


## For stan

num_bin_dstar <- as.numeric(table(deathdat$location_id)) 
length_dstar <- sum(num_bin_dstar)

#### Dstar and Rstar lengths ----

num_bin_dstar <- as.numeric(table(deathdat$location_id)) 
num_bin_rstar <- as.numeric(table(seroprevdat$location_id)) 
length_dstar <- sum(num_bin_dstar)
length_rstar <- sum(num_bin_rstar)

D_star <- as.numeric(deathdat$deaths)
N <- round(as.numeric(deathdat$population))
R_star <- round(as.numeric(seroprevdat$num_positive))
n <- round(as.numeric(seroprevdat$sample_size))

#### Study match ----

study_match_dstar <- as.numeric(deathdat$location_id)
study_match_rstar <- as.numeric(seroprevdat$location_id)

#### Country indicators ----
studydatifr$country <- factor(studydatifr$country, levels=unique(studydatifr$country))
country_inds <- as.numeric(studydatifr$country)
num_country <- length(unique(country_inds))

tab <- table(country_inds)
multLocInd <- which(country_inds %in% names(tab[tab>1]))
notMultLocInd <- which(!(country_inds %in% names(tab[tab>1])))

country_inds <- country_inds[multLocInd]
country_inds <- as.numeric(as.factor(country_inds))

#### Location indicators/total sero inds ----

loc_ind <- as.numeric(deathdat$location_id)
sero_to_studydat <- sapply(studydatifr$location_id, function(x) which(seroprevdat$location_id==x))
min_match_total_sero <- sapply(sero_to_studydat, function(x) min(x))
max_match_total_sero <- sapply(sero_to_studydat, function(x) max(x))

#### Covariates----

# h <- 0.5
ages <- seq(0,100, by=h) 

## For IFR
X <- cbind(1,ns(ages, knots=c(10, 60), Boundary.knots=c(0,80)))
X[,-1] <- scale(X[,-1])

min_match_X_to_death <- numeric(nrow(deathdat))
max_match_X_to_death <- numeric(nrow(deathdat))

for (i in 1:nrow(deathdat)) {
  min_match_X_to_death[i] <- which(ages == deathdat$lower_age[i]) 
  upage <- deathdat$upper_age[i]+1
  upage <- min(upage, 100)
  if (upage==100) {
    max_match_X_to_death[i] <- length(ages)
  } else {
    max_match_X_to_death[i] <- which(ages == upage)
  }
}

## For sero
Z <- cbind(1,ns(ages,knots=60, Boundary.knots = c(10,80)))
Z[,2:ncol(Z)] <- scale(Z[,2:ncol(Z)])

min_match_Z_to_sero <- numeric(nrow(seroprevdat))
max_match_Z_to_sero <- numeric(nrow(seroprevdat))

for (i in 1:nrow(seroprevdat)) {
  min_match_Z_to_sero[i] <- which(ages == seroprevdat$lower_age[i]) 
  
  upage <- seroprevdat$upper_age[i]+1
  upage <- min(upage, 100)
  if (upage==100) {
    max_match_Z_to_sero[i] <- length(ages)
  } else {
    max_match_Z_to_sero[i] <- which(ages == upage)
  }
}

#### Normalizing constants ----

f_normalizing_constants_dstar <- numeric(nrow(deathdat))
for (i in 1:nrow(deathdat)) {
  lwlim <- which(fineages==deathdat$lower_age[i])
  uplim <- which(fineages==min(deathdat$upper_age[i]+1,100))
  fs <- fineexpandedf[,as.numeric(deathdat$location_id)[i]] 
  # f_normalizing_constants_dstar[i] <- 1/sum(fs)
  f_normalizing_constants_dstar[i] <- (0.5*(fs[lwlim]+fs[uplim]) + 
                                         sum(fs[lwlim:(uplim-1)])) * (fineages[uplim]-fineages[lwlim])/(uplim-lwlim+1)
  f_normalizing_constants_dstar[i] <- 1/f_normalizing_constants_dstar[i] * 
    (fineages[uplim]-fineages[lwlim])/(max_match_X_to_death[i]-min_match_X_to_death[i]+1)
}

f_normalizing_constants_rstar <- numeric(sum(seroprevdat$location_id %in% ifrlocs))
for (i in which(seroprevdat$location_id %in% ifrlocs)) {
  lwlim <- which(fineages==seroprevdat$lower_age[i])
  uplim <- which(fineages==min(seroprevdat$upper_age[i]+1,100))
  fs <- fineexpandedf[, as.numeric(seroprevdat$location_id)[i]]
  f_normalizing_constants_rstar[i] <- (0.5*(fs[lwlim]+fs[uplim]) + 
                                         sum(fs[lwlim:(uplim-1)])) * (fineages[uplim]-fineages[lwlim])/(uplim-lwlim+1)
  f_normalizing_constants_rstar[i] <- 1/f_normalizing_constants_rstar[i] * 
    (fineages[uplim]-fineages[lwlim])/(max_match_Z_to_sero[i]-min_match_Z_to_sero[i]+1)
}


#### Match test kits ----

study_match_pi <- rep(99, nrow(studydatifr))
for (i in 1:nrow(studydatifr)){
  study_match_pi[i] <- which(studydatifr$test_id[i]==testkitdats$test_id)
}
study_match_pi


#### Make the list ----

stan_dat <- list(
  num_loc=nrow(studydatifr),
  num_country=num_country,
  
  length_dstar=length_dstar,
  D_star=D_star,
  N=as.numeric(N),
  
  length_rstar=length_rstar,
  R_star=R_star,
  n=n,
  
  num_test=nrow(testkitdats),
  sens_n=round(as.numeric(testkitdats$sen_tested_positives)),
  sens_x=round(as.numeric(testkitdats$sen_num_confirm_positives)),
  spec_n=round(as.numeric(testkitdats$spec_tested_negatives)),
  spec_x=round(as.numeric(testkitdats$spec_num_confirmed_negative)),
  
  n_int=nrow(X),
  h=h, 
  f_expanded=expandedf,
  f_normalizing_constants_dstar=f_normalizing_constants_dstar,
  f_normalizing_constants_rstar=f_normalizing_constants_rstar,
  
  ncolX=ncol(X),
  X=X,
  min_match_X_to_death=min_match_X_to_death,
  max_match_X_to_death=max_match_X_to_death,
  
  ncolZ=ncol(Z),
  Z=Z,
  min_match_Z_to_sero = min_match_Z_to_sero,
  max_match_Z_to_sero = max_match_Z_to_sero,
  
  study_match_pi=study_match_pi,
  study_match_dstar=study_match_dstar,
  study_match_rstar=study_match_rstar,
  
  sens=as.numeric(testkitdats$sensitivity_estimate)/100,
  spec=as.numeric(testkitdats$specificity_estimate)/100,
  
  num_mult_loc=length(multLocInd),
  num_not_mult_loc=length(notMultLocInd),
  mult_loc_ind=multLocInd,
  not_mult_loc_ind=notMultLocInd,
  num_country_mult_loc=length(unique(country_inds)),
  country_inds=country_inds,
  
  gamma_world=c(0,0),
  tau=c(0.05, 0.05)
)


# Initialization ----

#### Estimate IFR for initialization ----

## Estimate raw IFR ##
pis <- as.numeric(seroprevdat$num_positive)/as.numeric(seroprevdat$sample_size)
testkitdats$sensitivity_estimate <- testkitdats$sen_num_confirm_positives/testkitdats$sen_tested_positives
testkitdats$specificity_estimate <- testkitdats$spec_num_confirmed_negative/testkitdats$spec_tested_negatives
senss <- unlist(sapply(1:nrow(studydatifr), function(i) {
  as.numeric(testkitdats$sensitivity_estimate[rep(study_match_pi[i], each=num_bin_rstar[i])])
}))
specs <-  unlist(sapply(1:nrow(studydatifr), function(i) {
  as.numeric(testkitdats$specificity_estimate[rep(study_match_pi[i], each=num_bin_rstar[i])])
}))
rgadj <- (pis+specs-1)/(senss+specs-1)
rgadj[rgadj<.00001] <- .00001 

ifr_raw <- matrix(NA, 101, num_loc)
seros <- rep(NA,101)
deaths <- rep(NA,101)

for (ind in 1:num_loc) {
  loc <- studydatifr$location_id[ind]
  
  d <- deathdat[deathdat$location_id==loc,]
  s <- seroprevdat[seroprevdat$location_id==loc,]
  rg <- rgadj[seroprevdat$location_id==loc]
  
  for (i in 1:nrow(s)) {
    seros[s$lower_age[i]:s$upper_age[i]+1] <- rg[i]
  }
  for (i in 1:nrow(d)) {
    deaths[d$lower_age[i]:d$upper_age[i]+1] <- d$deaths[i]/d$population[i]
  }
  
  ifr_raw[,ind] <- deaths/seros
}

ifr_raw[ifr_raw<=0] <- .0000001
ifr_raw[ifr_raw>.9] <- .9 

#### IFR function ----

xs <- 0:100
matchX <- sapply(xs, function(x) which.min((ages-x)^2))

logIFR <- log(ifr_raw)

beta <- matrix(NA, ncol(X), num_loc)
for (i in 1:nrow(studydatifr)) {
  beta[,i] <- coef(lm(logIFR[,i] ~ X[matchX,]-1))
}

delta <- beta

#### Sero function ----

xsz <- apply(seroprevdat[,c("lower_age","upper_age")],1,mean)
matchZ <- sapply(xsz, function(x) which.min((ages-x)^2)) 
logitsero <- qlogis(rgadj)

gamma <- matrix(NA, ncol(Z), num_loc)
for (i in 1:nrow(studydatifr)) {
  ind <- which(seroprevdat$location_id==studydatifr$location_id[i])
  if (length(ind)==1) {
    gamma[,i] <- c(logitsero[ind],rep(0, ncol(Z)-1))
  } else {
    if (min(ages[matchZ[ind]])>18) {
      gamma[,i] <- coef(lm(c(mean(logitsero[ind]),logitsero[ind]) ~ Z[c(1,matchZ[ind]),]-1))
    } else {
      gamma[,i] <- coef(lm(logitsero[ind] ~ Z[matchZ[ind],]-1))
    }
  }
}

gamma[which(is.na(gamma),arr.ind = T)] <- 0

#### Initial values ----

init <- vector(mode="list", length=4)

init[[1]]$sens <- stan_dat$sens_x/stan_dat$sens_n-.00001
init[[1]]$spec <- stan_dat$spec_x/stan_dat$spec_n-.00001

init[[1]]$sigma <- rep(1,ncol(X)) 
init[[1]]$tau <- stan_dat$tau

init[[1]]$beta_world <- apply(beta, 1, median)
init[[1]]$gamma_world <- stan_dat$gamma_world

init[[1]]$gamma_raw <- apply(gamma, 2, function(x) (x-c(0,init[[1]]$gamma_world))/
                               c(5, init[[1]]$tau))

init[[1]]$beta_country <- rep(0, num_loc) 
init[[1]]$sigma_country <- .5

init[[1]]$beta_raw <- apply(delta, 2, function(x) (x-init[[1]]$beta_world)/
                              init[[1]]$sigma) 
init[[1]]$beta_raw[1,stan_dat$mult_loc_ind] <- (delta[1,stan_dat$mult_loc_ind]-
                                                  init[[1]]$beta_world[1]-init[[1]]$beta_country[country_inds])/init[[1]]$sigma[1]
init[[1]]$beta_raw[1,stan_dat$not_mult_loc_ind] <- (delta[1,stan_dat$not_mult_loc_ind]-
                                                      init[[1]]$beta_world[1])/
  sqrt(init[[1]]$sigma[1]^2 + init[[1]]$sigma_country^2)

init <- list(init[[1]])

## Run stan ----

m0 <- cmdstan_model("ifr.stan")

res <- m0$sample(data=stan_dat, chains=3, 
                 adapt_delta=0.8, iter_warmup=2500, iter_sampling=3000, refresh=100, max_treedepth=15,
                 init=list(init[[1]], init[[1]], init[[1]]))

