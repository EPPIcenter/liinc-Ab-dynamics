## ----packages and data---------------------------------------------------

library(dplyr)
library(tidyverse)
library(readxl)
library(lme4)
library(broom)
library(rstan)
library(tibble)

## Read in the raw data
dat_ind_URL <- "https://www.medrxiv.org/content/medrxiv/early/2021/03/05/2021.03.03.21251639/DC12/embed/media-12.xlsx"
dat_visit_URL <- "https://www.medrxiv.org/content/medrxiv/early/2021/03/05/2021.03.03.21251639/DC13/embed/media-13.xlsx"

tmp_ind <- tempfile(fileext=".xlsx")
tmp_visit <- tempfile(fileext=".xlsx")

download.file(dat_ind_URL, destfile=tmp_ind, mode="wb")
download.file(dat_visit_URL, destfile=tmp_visit, mode="wb")

dat_ind <- read_excel(tmp_ind, sheet=1)
dat_visit <- read_excel(tmp_visit, sheet=1)

## ----transform-----------------------------------------------------------

## Get the analytical scale from Supplementary Table 1
trans_natural <- c("N-Abbott", "S-Ortho Ig", "S-Ortho IgG", "S-DiaSorin")
trans_log <- c("N-Roche", "Neut-Monogram", "RBD-Split Luc", "N-Split Luc", "RBD-LIPS", "S-LIPS", "N-LIPS", "S-Lum", "RBD-Lum", "N(full)-Lum", "N(frag)-Lum")

## Transform the data
dat_visit %>%
	mutate(`N-Abbott` = ifelse(`N-Abbott`==0, 0.01, `N-Abbott`)) %>%
	mutate_at(vars(contains(trans_log)), log) %>%
	mutate(days_since_seroconv = days_since_onset-21) -> dat_visit

## Join with individual-level data
dat_visit %>%
	left_join(dat_ind, by="participant_ID") -> dat_visit

## ----mixed effects-------------------------------------------------------

dat_visit %>%
	pivot_longer(`N-Abbott`:`Neut-Monogram`, names_to="assay", values_to="response") %>%
	mutate(assay = factor(assay)) -> dat_visit_long

## Fit the lmer model to each assay
dat_visit_long %>%
	group_by(assay) %>%
	nest() %>%
	mutate(
		fit = map(data, ~ lmer(response ~ factor(hosp_status) + days_since_seroconv + (1|participant_ID), REML=TRUE, data = .x)),
		tidied = map(fit, tidy)) %>% 
	unnest(tidied) -> fit_lmer

## ----time to seroreversion-----------------------------------------------

data.frame(assay=levels(dat_visit_long$assay)) %>%
	mutate(trans = ifelse(assay %in% trans_natural, "natural", "log")) -> dat_assays

## Get the cutoff for positivity from Supplementary Table 1
dat_assays$cutoff <- c(1.4, 125000, 1, 83.1, 0.02684, 0.02473, 40, 52000, 0.0396, 45.9, 15, 45000, 0.0426, 1, 1)
dat_assays %>% mutate(cutoff = ifelse(trans=="log", log(cutoff), cutoff)) -> dat_assays

## Add meta-data
dat_visit_long %>%
	mutate(hosp = ifelse(hosp_status=="Yes", 1, 0)) %>%
	mutate(hosp = factor(hosp)) -> dat_visit_long_with_meta

## To save output
fit_days_to_reversion <- vector(mode="list", length=nrow(dat_assays))
fit_days_to_reversion_CI <- vector(mode="list", length=nrow(dat_assays))

for(i in 1:length(unique(dat_assays$assay))) {
	
	## Get the data for a single assay
	dat_visit_long_with_meta %>%
		filter(assay==as.character(unique(dat_assays$assay)[i])) %>%
		filter(!is.na(response)) -> dat_single
	
	## Fit the lmer model
	fit_lmer_single <- lmer(response ~ -1 + factor(hosp_status) + days_since_seroconv + (1|participant_ID), REML=TRUE, data=dat_single)
	
	## Function to estimate time to sero-reversion
	func_days_to_reversion <- function(fitted_lmer_model, cut_point=as.numeric(dat_assays$cutoff[i])) {
		
		c(Days_to_reversion = (cut_point - fixef(fitted_lmer_model)[1]) / fixef(fitted_lmer_model)[3],
		Days_to_reversion = (cut_point - fixef(fitted_lmer_model)[2]) / fixef(fitted_lmer_model)[3])
		
	}
	
	## Get bootstrap CIs
	fit_days_to_reversion[[i]] <- bootstrapped_days_to_reversion <- bootMer(x=fit_lmer_single, FUN=func_days_to_reversion, nsim=10000, type="semiparametric", use.u=TRUE)
	
	fit_days_to_reversion_CI[[i]] <- confint(bootstrapped_days_to_reversion)
	
}

## ----sensitivity---------------------------------------------------------

## Get the specificity from Supplementary Table 1
dat_assays$Sp <- c(1, 1, 0.998, 0.99, 1, 1, 0.988, 1, 1, 1, 0.993, 1, 1, 1, 1)

## To save output
fit_sensitivity <- vector(mode="list", length=nrow(dat_assays))
fit_NPV <- vector(mode="list", length=nrow(dat_assays))

for(i in 1:length(unique(dat_assays$assay))) {
	
	## Get the data for a single assay
	dat_visit_long_with_meta %>%
		filter(assay==as.character(unique(dat_assays$assay)[i])) %>%
		filter(!is.na(response)) %>%
		droplevels() %>%
		group_by(participant_ID) %>%
		mutate(participant_ID_index = group_indices()) %>%
		ungroup() %>%
		arrange(participant_ID_index, days_since_seroconv) %>%
		select(participant_ID, hosp, participant_ID_index, days_since_seroconv, response) %>%
		as.data.frame() -> dat_for_Stan_single
	
	## Fit the model in Stan
	fit_Stan_single <- stan(file="Estimate_Se.stan",
		data=list(
			N = nrow(dat_for_Stan_single),
			J = length(unique(dat_for_Stan_single$participant_ID_index)),
			participant_ID_index = as.integer(dat_for_Stan_single$participant_ID_index),
			days_since_seroconv = dat_for_Stan_single$days_since_seroconv,
			response = dat_for_Stan_single$response,
			## Severities
			hosp_by_ind = as.integer(dat_for_Stan_single %>% group_by(participant_ID) %>% slice(1) %>% ungroup() %>% select(hosp) %>% pull()),
			## How many times to sample hierarchical params?
			N_rep = as.integer(500),
			## Which dates to estimate sensitivity at?
			K = 4,
			days_for_Se = c(0,2,4,6)*30,
			## What is the cutoff?
			cutoff = dat_assays$cutoff[i]),
		chains=4, iter=2000, thin=1, init=0, seed=1234,
		control=list(stepsize=0.1, adapt_delta=0.9, max_treedepth=15))

	## Save output
	summary(fit_Stan_single)$summary %>%
		as.data.frame() %>%
		rownames_to_column(var="varname") %>%
		filter(str_detect(varname,"sensitivity") & !str_detect(varname,"ratio")) %>%
		mutate(assay = as.character(unique(dat_assays$assay)[i])) %>%
		mutate(hosp = c(rep(0,4),rep(1,4))) %>%
		mutate(month = rep(c(0,2,4,6), times=2)) -> fit_sensitivity[[i]]
	
	## Look at NPV, for commercial assays
	if(as.character(unique(dat_assays$assay)[i]) %in% c("N-Abbott", "N-Roche", "S-Ortho Ig", "S-Ortho IgG", "S-DiaSorin")) {
		
		calc_NPV <-'
		data {
			int len_prev;
			real prev[len_prev];
			real Sp;
		}
		parameters {
			matrix<lower=0.0, upper=1.0>[2,4] sensitivity_by_severity_time;
		}
		generated quantities {
			real<lower=0.0, upper=1.0> NPV[2,4,len_prev];
			
			for(i in 1:2) {
				for(j in 1:4) {
					for(k in 1:len_prev) {
						NPV[i,j,k] = (Sp*(1.0-prev[k])) / ((Sp*(1.0-prev[k])) + ((1.0-sensitivity_by_severity_time[i,j])*prev[k]));
					}
				}
			}
		}
		'
	
		## Read in model to generate NPVs
		model_NPV <- stan_model(model_code=calc_NPV)
		
		## Sweep over prevalences
		PREV <- seq(0.05, 0.5, by=0.05)
		
		## Get the Sp
		SPEC <- dat_assays %>% filter(assay==as.character(unique(dat_assays$assay)[i])) %>% select(Sp) %>% pull()
		
		## Ordered as NPV[hosp, month, prevalence]
		NPV <- gqs(model_NPV, data=list(len_prev=length(PREV), prev=PREV, Sp=SPEC), draws=as.matrix(fit_Stan_single))
		
		as.data.frame(summary(NPV)$summary) %>%
			rownames_to_column() %>%
			mutate(Assay=as.character(unique(dat_assays$assay)[i])) -> tmp
	
		tmp$hosp <- rep(c(0,1), each=nrow(tmp)/2)
		tmp$month <- rep(rep(c(0,2,4,6), each=length(PREV)), times=2)
		tmp$prevalence <- rep(PREV, times=nrow(tmp)/length(PREV))
		
		fit_NPV[[i]] <- tmp
		
	}
	
}
