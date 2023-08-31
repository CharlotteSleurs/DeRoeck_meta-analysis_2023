###############################################
#### ANALYSES FOR META-ANALYSIS ON GLIOMA #####
###############################################

rm(list = ls())

################
### PACKAGES ###
################

library(readxl)
library(metafor)

#################
### FUNCTIONS ###
#################

### Function for computing Hedges' g and its sampling variance in case of 
# pre and posttest design. Formulas are from sections 11.2.2.2 of 
# Cooper, Hedges, and Valentine (2019). The Handbook of Research Synthesis and 
# Meta-Analysis (3rd ed.)
es_long <- function(change, sd_change, m1, sd1, m2, sd2, n1, n2, zval1, zval2, r)
{
  
  n <- 1/mean(1/c(n1, n2)) # Compute harmonic mean of sample sizes
  
  if (is.na(change) == FALSE & is.na(sd_change) == FALSE)
  { # If change scores are reported 
    swithin <- sd_change/sqrt(2*(1-r))
    d <- change/swithin
  } else if (is.na(m1) == FALSE & is.na(sd1) == FALSE & is.na(m2) == FALSE & 
             is.na(sd2) == FALSE)
  {
    swithin <- sqrt((sd1^2+sd2^2)/2)
    d <- (m2-m1)/swithin 
    
  } else if (is.na(zval1) == FALSE & is.na(zval2) == FALSE & is.na(sd1) == FALSE & 
             is.na(sd2) == FALSE)
  {
    yi <- zval2-zval1
    vi <- (sd1^2 + sd2^2 - 2*r*sd1*sd2)/n
  } 
  
  if ((is.na(change) == FALSE & is.na(sd_change) == FALSE) | 
      (is.na(m1) == FALSE & is.na(sd1) == FALSE & is.na(m2) == FALSE & 
       is.na(sd2) == FALSE))
  {
    vd <- (1/n+d^2/(2*n))*2*(1-r)
    m <- n-1
    ### Hedges' g correction factor (exact rather than approximation)
    J <- exp(lgamma(m/2) - log(sqrt(m/2)) - lgamma((m-1)/2)) 
    yi <- J*d
    vi <- J^2*vd
  } 
  
  return(data.frame(yi = yi, vi = vi))
}

### Function for computing Hedges' g and its sampling variance in case of 
# a cross-sectional design. Formulas are from section 11.2.2.1 of Cooper, Hedges,
# and Valentine (2019). The Handbook of Research Synthesis and Meta-Analysis (3rd ed.)
es_cross <- function(m1, sd1, m2, sd2, n1, n2, zval2)
{
  
  if (is.na(m1) == FALSE & is.na(sd1) == FALSE & is.na(m2) == FALSE & 
      is.na(sd2) == FALSE)
  {
    spool <- sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
    d <- (m2-m1)/spool
    
  } else if (is.na(zval2) == FALSE & is.na(sd2) == FALSE)
  {
    yi <- zval2
    vi <- sd2^2/n2
  } else 
  {
    yi <- vi <- NA
  }
  
  if (is.na(m1) == FALSE & is.na(sd1) == FALSE & is.na(m2) == FALSE & 
      is.na(sd2) == FALSE)
  {
    vd <- (n1+n2)/(n1*n2) + d^2/(2*(n1+n2))
    m <- n1+n2-2
    ### Hedges' g correction factor (exact rather than approximation)
    J <- exp(lgamma(m/2) - log(sqrt(m/2)) - lgamma((m-1)/2)) 
    yi <- J*d
    vi <- J^2*vd
  } 
  
  return(data.frame(yi = yi, vi = vi))
}

################################################################################

####################
### LONGITUDINAL ###
####################

### Load data
dat_long <- read_excel("Meta-analyse_300722_final_2.xlsx", sheet = "Longitudinaal")

yi <- vi <- es_type <- ID <- RiskOfBias2Risks <- n_baseline <- n_post <- 
  type_test <- Testnumber <- Domain_number <- Majority_HGG <- NA

b <- 1 # Counter for storing results

for (i in 1:nrow(dat_long))
{
  
  m1 <- m2 <- sd1 <- sd2 <- change <- sd_change <- zval1 <- zval2 <- NA
  
  ### Store sample sizes of baseline and post measurement
  n1 <- dat_long$n_baseline[i]
  n2 <- dat_long$n_post[i]
  
  if (is.na(dat_long$mean_baseline[i]) == FALSE & 
      is.na(dat_long$mean_post[i]) == FALSE & 
      is.na(dat_long$sd_baseline[i]) == FALSE & 
      is.na(dat_long$sd_post[i]) == FALSE)
  { # If group means and standard deviations are reported
    m1 <- dat_long$mean_baseline[i]
    m2 <- dat_long$mean_post[i]
    sd1 <- dat_long$sd_baseline[i]
    sd2 <- dat_long$sd_post[i]
    
    ### Compute Hedges' g and sampling variance given correlation is 0.5
    es <- es_long(change = change, sd_change = sd_change, m1 = m1, sd1 = sd1, 
                  m2 = m2, sd2 = sd2, n1 = n1, n2 = n2, zval1 = zval1, 
                  zval2 = zval2, r = 0.5)
    
    ### Multiply effect size with -1 if a lower score implies better performance
    yi[b] <- ifelse(dat_long$`Interpretation_scores_raw-scores`[i] == 0, -1*es$yi, es$yi)
    vi[b] <- es$vi
    es_type[b] <- "means"

    ### Store data
    ID[b] <- dat_long$ID[i]
    RiskOfBias2Risks[b] <- dat_long$"RiskOfBias>2Risks"[i]
    n_baseline[b] <- dat_long$n_baseline[i]
    n_post[b] <- dat_long$n_post[i]
    type_test[b] <- dat_long$type_test[i]
    Testnumber[b] <- dat_long$"Test number"[i]
    Domain_number[b] <- dat_long$Domain_number[i]
    Majority_HGG[b] <- dat_long$Majority_HGG[i]

    b <- b+1 # Add one to counter for storing results
  } 
  
  if (is.na(dat_long$mean_diff[i]) == FALSE &
             is.na(dat_long$sd_diff[i]) == FALSE)
  { # If change scores are reported
    change <- dat_long$mean_diff[i]
    sd_change <- dat_long$sd_diff[i]
    
    ### Compute Hedges' g and sampling variance given correlation is 0.5
    es <- es_long(change = change, sd_change = sd_change, m1 = m1, sd1 = sd1, 
                  m2 = m2, sd2 = sd2, n1 = n1, n2 = n2, zval1 = zval1, 
                  zval2 = zval2, r = 0.5)
    
    ### Multiply effect size with -1 if a lower score implies better performance
    yi[b] <- ifelse(dat_long$`Interpretation_scores_raw-scores`[i] == 0, -1*es$yi, es$yi)
    vi[b] <- es$vi
    es_type[b] <- "change"
    
    ### Store data
    ID[b] <- dat_long$ID[i]
    RiskOfBias2Risks[b] <- dat_long$"RiskOfBias>2Risks"[i]
    n_baseline[b] <- dat_long$n_baseline[i]
    n_post[b] <- dat_long$n_post[i]
    type_test[b] <- dat_long$type_test[i]
    Testnumber[b] <- dat_long$"Test number"[i]
    Domain_number[b] <- dat_long$Domain_number[i]
    Majority_HGG[b] <- dat_long$Majority_HGG[i]

    b <- b+1 # Add one to counter for storing results
  } 
    
  if (is.na(dat_long$"Z-score_mean_baseline"[i]) == FALSE &
             is.na(dat_long$"Z-score_mean_post"[i]) == FALSE)
  { # If z-values are reported
    zval1 <- dat_long$"Z-score_mean_baseline"[i]
    zval2 <- dat_long$"Z-score_mean_post"[i]
    sd1 <- dat_long$"Z-score_SD_baseline"[i]
    sd2 <- dat_long$"Z-score_SD_post"[i]
    
    ### Compute mean difference in z-scores and sampling variance given 
    # correlation is 0.5
    es <- es_long(change = change, sd_change = sd_change, m1 = m1, sd1 = sd1, 
                  m2 = m2, sd2 = sd2, n1 = n1, n2 = n2, zval1 = zval1, 
                  zval2 = zval2, r = 0.5)
    
    ### Multiply effect size with -1 if a lower score implies better performance
    yi[b] <- ifelse(dat_long$`Interpretation_scores_z-scores`[i] == 0, -1*es$yi, es$yi)
    vi[b] <- es$vi
    es_type[b] <- "z-score"
    
    ### Store data
    ID[b] <- dat_long$ID[i]
    RiskOfBias2Risks[b] <- dat_long$"RiskOfBias>2Risks"[i]
    n_baseline[b] <- dat_long$n_baseline[i]
    n_post[b] <- dat_long$n_post[i]
    type_test[b] <- dat_long$type_test[i]
    Testnumber[b] <- dat_long$"Test number"[i]
    Domain_number[b] <- dat_long$Domain_number[i]
    Majority_HGG[b] <- dat_long$Majority_HGG[i]

    b <- b+1 # Add one to counter for storing results
  }
}

### Create dataframe
dat_long <- data.frame(ID = ID, RiskOfBias2Risks = RiskOfBias2Risks, 
                       n_baseline = n_baseline, n_post = n_post, 
                       type_test = type_test, Testnumber = Testnumber, 
                       Domain_number = Domain_number, 
                       Majority_HGG = Majority_HGG, yi = yi,
                       vi = vi, es_type = es_type)

### Save data
save(dat_long, file = "dat_long.RData")

################################################################################

#######################
### CROSS-SECTIONAL ###
#######################

### Load data
dat_cross <- read_excel("Meta-analyse_300722_final_2.xlsx", sheet = "Cross-sectioneel")

yi <- vi <- es_type <- ID <- RiskOfBias2Risks <- n_contr <- n_interv <- 
  type_test <- Domain_number <- Majority_HGG <- NA

b <- 1 # Counter for storing results

for (i in 1:nrow(dat_cross))
{
  
  m1 <- m2 <- sd1 <- sd2 <- zval1 <- zval2 <- NA
  
  ### Store sample sizes of baseline and post measurement
  n1 <- dat_cross$n_contr[i]
  n2 <- dat_cross$n_interv[i]
  
  if (is.na(dat_cross$mean_contr[i]) == FALSE & 
      is.na(dat_cross$mean_interv[i]) == FALSE & 
      is.na(dat_cross$sd_contr[i]) == FALSE & 
      is.na(dat_cross$sd_int[i]) == FALSE)
  { # If group means and standard deviations are reported
    m1 <- dat_cross$mean_contr[i]
    m2 <- dat_cross$mean_interv[i]
    sd1 <- dat_cross$sd_contr[i]
    sd2 <- dat_cross$sd_int[i]
    
    ### Compute Hedges' g and sampling variance
    es <- es_cross(m1 = m1, sd1 = sd1, m2 = m2, sd2 = sd2, n1 = n1, n2 = n2, 
                   zval2 = zval2)
    
    ### Multiply effect size with -1 if a lower score implies better performance
    yi[b] <- ifelse(dat_cross$`Interpretation_scores_raw-scores`[i] == 0, -1*es$yi, es$yi)
    vi[b] <- es$vi
    es_type[b] <- "means"
    
    ### Store data
    ID[b] <- dat_cross$ID[i]
    RiskOfBias2Risks[b] <- dat_cross$"RiskOfBias>2Risks"[i]
    n_interv[b] <- dat_cross$n_interv[i]
    n_contr[b] <- dat_cross$n_contr[i]
    type_test[b] <- dat_cross$type_test[i]
    Domain_number[b] <- dat_cross$Domain_number[i]
    Majority_HGG[b] <- dat_cross$Majority_HGG[i]

    b <- b+1 # Add one to counter for storing results
  }
    
  if (is.na(dat_cross$"Z-score_mean_interv"[i]) == FALSE)
  { # If z-values are reported
    zval2 <- dat_cross$"Z-score_mean_interv"[i]
    sd2 <- dat_cross$"Z-score_SD_interv"[i]
    
    ### Compute Hedges' g and sampling variance
    es <- es_cross(m1 = m1, sd1 = sd1, m2 = m2, sd2 = sd2, n1 = n1, n2 = n2, 
                   zval2 = zval2)
    
    ### Multiply effect size with -1 if a lower score implies better performance
    yi[b] <- ifelse(dat_cross$`Interpretation_scores_z-scores`[i] == 0, -1*es$yi, es$yi)
    vi[b] <- es$vi
    es_type[b] <- "z-score"
    
    ### Store data
    ID[b] <- dat_cross$ID[i]
    RiskOfBias2Risks[b] <- dat_cross$"RiskOfBias>2Risks"[i]
    n_interv[b] <- dat_cross$n_interv[i]
    n_contr[b] <- dat_cross$n_contr[i]
    type_test[b] <- dat_cross$type_test[i]
    Domain_number[b] <- dat_cross$Domain_number[i]
    Majority_HGG[b] <- dat_cross$Majority_HGG[i]

    b <- b+1 # Add one to counter for storing results
  }
}

### Create dataframe
dat_cross <- data.frame(ID = ID, RiskOfBias2Risks = RiskOfBias2Risks, 
                       n_interv = n_interv, n_contr = n_contr, 
                       type_test = type_test, Domain_number = Domain_number, 
                       Majority_HGG = Majority_HGG, yi = yi, vi = vi, 
                       es_type = es_type)

### Save data
save(dat_cross, file = "dat_cross.RData")

################################################################################
################################################################################
################################################################################

##### SUMMARIZE LONGITUDINAL STUDIES PER TYPE OF TEST #####

### Fix some spelling mistakes in the data
dat_long$type_test <- ifelse(dat_long$type_test == "MMSE ", 
                             "MMSE", dat_long$type_test)

### Create dummy variable to denote whether the posttest was the first one 
# (score = 0) or a later one (score = 1). The score is NA if it is unknown.

### Empty object that is going to be filled
dat_long$multiple_posts <- NA

for (b in 1:nrow(dat_long))
{
  if (dat_long$Testnumber[b] == "T1")
  {
    dat_long$multiple_posts[b] <- 0
  } else if (dat_long$Testnumber[b] == "Unkown")
  {
    dat_long$multiple_posts[b] <- NA
  } else
  {
    dat_long$multiple_posts[b] <- 1
  }
}

### Split data per type of test
splits <- split(dat_long, dat_long$type_test)

### Empty objects for storing results
long_means <- data.frame(type_test = character(length(splits)), 
                         domain_number = rep(NA, length(splits)),
                         k = rep(NA, length(splits)), 
                         sum_n = rep(NA, length(splits)), 
                         est = rep(NA, length(splits)), 
                         se = rep(NA, length(splits)),
                         lb = rep(NA, length(splits)), 
                         ub = rep(NA, length(splits)),
                         zval = rep(NA, length(splits)), 
                         pval = rep(NA, length(splits)),
                         tau2 = rep(NA, length(splits)),
                         se_tau2 = rep(NA, length(splits)),
                         lb_tau2 = rep(NA, length(splits)), 
                         ub_tau2 = rep(NA, length(splits)), 
                         Qval = rep(NA, length(splits)), 
                         pval_Q = rep(NA, length(splits)), 
                         I2 = rep(NA, length(splits)),
                         est_ee = rep(NA, length(splits)), 
                         se_ee = rep(NA, length(splits)),
                         lb_ee = rep(NA, length(splits)), 
                         ub_ee = rep(NA, length(splits)),
                         zval_ee = rep(NA, length(splits)), 
                         pval_ee = rep(NA, length(splits)),
                         k_risk = rep(NA, length(splits)), 
                         sum_n_risk = rep(NA, length(splits)), 
                         est_risk = rep(NA, length(splits)), 
                         se_risk = rep(NA, length(splits)),
                         lb_risk = rep(NA, length(splits)), 
                         ub_risk = rep(NA, length(splits)),
                         zval_risk = rep(NA, length(splits)), 
                         pval_risk = rep(NA, length(splits)),
                         tau2_risk = rep(NA, length(splits)),
                         se_tau2_risk = rep(NA, length(splits)),
                         lb_tau2_risk = rep(NA, length(splits)), 
                         ub_tau2_risk = rep(NA, length(splits)), 
                         Qval_risk = rep(NA, length(splits)), 
                         pval_Q_risk = rep(NA, length(splits)), 
                         I2_risk = rep(NA, length(splits)),
                         k1 = rep(NA, length(splits)), 
                         b0 = rep(NA, length(splits)), 
                         b0_se = rep(NA, length(splits)),
                         b0_lb = rep(NA, length(splits)), 
                         b0_ub = rep(NA, length(splits)),
                         b0_zval = rep(NA, length(splits)), 
                         b0_pval = rep(NA, length(splits)),
                         b1 = rep(NA, length(splits)), 
                         b1_se = rep(NA, length(splits)),
                         b1_lb = rep(NA, length(splits)), 
                         b1_ub = rep(NA, length(splits)),
                         b1_zval = rep(NA, length(splits)), 
                         b1_pval = rep(NA, length(splits)),
                         b0_ee = rep(NA, length(splits)), 
                         b0_se_ee = rep(NA, length(splits)),
                         b0_lb_ee = rep(NA, length(splits)), 
                         b0_ub_ee = rep(NA, length(splits)),
                         b0_zval_ee = rep(NA, length(splits)), 
                         b0_pval_ee = rep(NA, length(splits)),
                         b1_ee = rep(NA, length(splits)), 
                         b1_se_ee = rep(NA, length(splits)),
                         b1_lb_ee = rep(NA, length(splits)), 
                         b1_ub_ee = rep(NA, length(splits)),
                         b1_zval_ee = rep(NA, length(splits)), 
                         b1_pval_ee = rep(NA, length(splits)),
                         k1_rb = rep(NA, length(splits)), 
                         b0_rb = rep(NA, length(splits)), 
                         b0_se_rb = rep(NA, length(splits)),
                         b0_lb_rb = rep(NA, length(splits)), 
                         b0_ub_rb = rep(NA, length(splits)),
                         b0_zval_rb = rep(NA, length(splits)), 
                         b0_pval_rb = rep(NA, length(splits)),
                         b1_rb = rep(NA, length(splits)), 
                         b1_se_rb = rep(NA, length(splits)),
                         b1_lb_rb = rep(NA, length(splits)), 
                         b1_ub_rb = rep(NA, length(splits)),
                         b1_zval_rb = rep(NA, length(splits)), 
                         b1_pval_rb = rep(NA, length(splits)))

long_zscore <- data.frame(type_test = character(length(splits)), 
                          domain_number = rep(NA, length(splits)),
                          k = rep(NA, length(splits)), 
                          sum_n = rep(NA, length(splits)), 
                          est = rep(NA, length(splits)), 
                          se = rep(NA, length(splits)),
                          lb = rep(NA, length(splits)), 
                          ub = rep(NA, length(splits)),
                          zval = rep(NA, length(splits)), 
                          pval = rep(NA, length(splits)),
                          tau2 = rep(NA, length(splits)),
                          se_tau2 = rep(NA, length(splits)),
                          lb_tau2 = rep(NA, length(splits)), 
                          ub_tau2 = rep(NA, length(splits)), 
                          Qval = rep(NA, length(splits)), 
                          pval_Q = rep(NA, length(splits)), 
                          I2 = rep(NA, length(splits)),
                          est_ee = rep(NA, length(splits)), 
                          se_ee = rep(NA, length(splits)),
                          lb_ee = rep(NA, length(splits)), 
                          ub_ee = rep(NA, length(splits)),
                          zval_ee = rep(NA, length(splits)), 
                          pval_ee = rep(NA, length(splits)),
                          k_risk = rep(NA, length(splits)), 
                          sum_n_risk = rep(NA, length(splits)), 
                          est_risk = rep(NA, length(splits)), 
                          se_risk = rep(NA, length(splits)),
                          lb_risk = rep(NA, length(splits)), 
                          ub_risk = rep(NA, length(splits)),
                          zval_risk = rep(NA, length(splits)), 
                          pval_risk = rep(NA, length(splits)),
                          tau2_risk = rep(NA, length(splits)),
                          se_tau2_risk = rep(NA, length(splits)),
                          lb_tau2_risk = rep(NA, length(splits)), 
                          ub_tau2_risk = rep(NA, length(splits)), 
                          Qval_risk = rep(NA, length(splits)), 
                          pval_Q_risk = rep(NA, length(splits)), 
                          I2_risk = rep(NA, length(splits)),
                          k1_rb = rep(NA, length(splits)), 
                          b0_rb = rep(NA, length(splits)), 
                          b0_se_rb = rep(NA, length(splits)),
                          b0_lb_rb = rep(NA, length(splits)), 
                          b0_ub_rb = rep(NA, length(splits)),
                          b0_zval_rb = rep(NA, length(splits)), 
                          b0_pval_rb = rep(NA, length(splits)),
                          b1_rb = rep(NA, length(splits)), 
                          b1_se_rb = rep(NA, length(splits)),
                          b1_lb_rb = rep(NA, length(splits)), 
                          b1_ub_rb = rep(NA, length(splits)),
                          b1_zval_rb = rep(NA, length(splits)), 
                          b1_pval_rb = rep(NA, length(splits)))

for (m in 1:length(splits))
{
  tmp <- splits[[m]]
  
  ### Select data based on means or change scores
  sub_means <- subset(tmp, es_type == "means" | es_type == "change")
  
  ### Store type of test and number of effect sizes
  long_means$type_test[m] <- names(splits)[[m]]
  long_means$k[m] <- nrow(sub_means)
  
  if (nrow(sub_means) > 0)
  {
    ### Conduct a random-effects meta-analysis
    out <- rma(yi = yi, vi = vi, dat = sub_means)
    
    ### Store results
    long_means$est[m] <- out$b[1]
    long_means$se[m] <- out$se
    long_means$lb[m] <- out$ci.lb
    long_means$ub[m] <- out$ci.ub
    long_means$zval[m] <- out$zval
    long_means$pval[m] <- out$pval
    long_means$tau2[m] <- out$tau2
    long_means$se_tau2[m] <- out$se.tau2
    long_means$Qval[m] <- out$QE
    long_means$pval_Q[m] <- out$QEp
    long_means$I2[m] <- out$I2
    
    ### Compute sum of harmonic mean of sample sizes
    long_means$sum_n[m] <- round(sum(apply(rbind(sub_means$n_baseline, sub_means$n_post), 
                                           2, function(x) 1/mean(1/x))))
    
    ### Store domain number for ordering rows in table
    long_means$domain_number[m] <- unique(sub_means$Domain_number)
    
    if (out$k > 1)
    {
      ci <- confint(out) # 95% CI tau^2
      long_means$lb_tau2[m] <- paste0(ci$lb.sign, ci$random["tau^2","ci.lb"])
      long_means$ub_tau2[m] <- paste0(ci$ub.sign, ci$random["tau^2","ci.ub"])
    }
    
    ### Conduct an equal-effect meta-analysis
    out <- rma(yi = yi, vi = vi, dat = sub_means, method = "EE")
    
    ### Store results
    long_means$est_ee[m] <- out$b[1]
    long_means$se_ee[m] <- out$se
    long_means$lb_ee[m] <- out$ci.lb
    long_means$ub_ee[m] <- out$ci.ub
    long_means$zval_ee[m] <- out$zval
    long_means$pval_ee[m] <- out$pval
    
    ### Meta-analysis without studies at high risk of bias
    sub_means_risk <- subset(sub_means, sub_means$RiskOfBias2Risks == 0)
    
    ### Store number of studies in meta-analysis without studies at high risk of bias
    long_means$k_risk[m] <- nrow(sub_means_risk)
    
    if (nrow(sub_means_risk) > 0)
    {
      
      ### Conduct a random-effects meta-analysis
      out <- rma(yi = yi, vi = vi, dat = sub_means_risk)
      
      ### Store results
      long_means$est_risk[m] <- out$b[1]
      long_means$se_risk[m] <- out$se
      long_means$lb_risk[m] <- out$ci.lb
      long_means$ub_risk[m] <- out$ci.ub
      long_means$zval_risk[m] <- out$zval
      long_means$pval_risk[m] <- out$pval
      long_means$tau2_risk[m] <- out$tau2
      long_means$se_tau2_risk[m] <- out$se.tau2
      long_means$Qval_risk[m] <- out$QE
      long_means$pval_Q_risk[m] <- out$QEp
      long_means$I2_risk[m] <- out$I2
      
      ### Compute sum of harmonic mean of sample sizes
      long_means$sum_n_risk[m] <- round(sum(apply(rbind(sub_means_risk$n_baseline, 
                                                        sub_means_risk$n_post), 
                                                  2, function(x) 1/mean(1/x))))
      
      if (out$k > 1)
      {
        ci <- confint(out) # 95% CI tau^2
        long_means$lb_tau2_risk[m] <- paste0(ci$lb.sign, ci$random["tau^2","ci.lb"])
        long_means$ub_tau2_risk[m] <- paste0(ci$ub.sign, ci$random["tau^2","ci.ub"])
      }
    }
    
    ### Number of studies that score 1 on "multiple_posts"
    long_means$k1[m] <- sum(sub_means$multiple_posts == 1, na.rm = TRUE)
    
    ### If there are studies with scores 0 and 1 on the moderator "multiple_posts"
    # and there are more than two studies, conduct the moderator analysis
    if (all(c(0,1) %in% unique(sub_means$multiple_posts)) & 
        (sum(sub_means$multiple_posts == 0, na.rm = TRUE)+
         sum(sub_means$multiple_posts == 1, na.rm = TRUE)) > 2)
    {
      ### Conduct a random-effects meta-analysis with moderator "multiple_posts"
      out <- rma(yi = yi, vi = vi, mods = ~ multiple_posts, dat = sub_means)
      
      ### Store results
      long_means$b0[m] <- out$b[1]
      long_means$b0_se[m] <- out$se[1]
      long_means$b0_lb[m] <- out$ci.lb[1]
      long_means$b0_ub[m] <- out$ci.ub[1]
      long_means$b0_zval[m] <- out$zval[1]
      long_means$b0_pval[m] <- out$pval[1]
      long_means$b1[m] <- out$b[2]
      long_means$b1_se[m] <- out$se[2]
      long_means$b1_lb[m] <- out$ci.lb[2]
      long_means$b1_ub[m] <- out$ci.ub[2]
      long_means$b1_zval[m] <- out$zval[2]
      long_means$b1_pval[m] <- out$pval[2]
      
      ### Conduct an equal-effect meta-analysis with moderator "multiple_posts"
      out <- rma(yi = yi, vi = vi, mods = ~ multiple_posts, dat = sub_means, 
                 method = "EE")
      
      ### Store results
      long_means$b0_ee[m] <- out$b[1]
      long_means$b0_se_ee[m] <- out$se[1]
      long_means$b0_lb_ee[m] <- out$ci.lb[1]
      long_means$b0_ub_ee[m] <- out$ci.ub[1]
      long_means$b0_zval_ee[m] <- out$zval[1]
      long_means$b0_pval_ee[m] <- out$pval[1]
      long_means$b1_ee[m] <- out$b[2]
      long_means$b1_se_ee[m] <- out$se[2]
      long_means$b1_lb_ee[m] <- out$ci.lb[2]
      long_means$b1_ub_ee[m] <- out$ci.ub[2]
      long_means$b1_zval_ee[m] <- out$zval[2]
      long_means$b1_pval_ee[m] <- out$pval[2]
    }
    
    ### Number of studies that score 1 on "Majority_HGG"
    long_means$k1_rb[m] <- sum(sub_means$Majority_HGG == 1, na.rm = TRUE)
    
    ### If there are studies with scores 0 and 1 on the moderator "Majority_HGG"
    # and there are more than two studies, conduct the moderator analysis
    if (all(c(0,1) %in% unique(sub_means$Majority_HGG)) & 
        (sum(sub_means$Majority_HGG == 0, na.rm = TRUE)+
         sum(sub_means$Majority_HGG == 1, na.rm = TRUE)) > 2)
    {
      ### Conduct a random-effects meta-analysis with moderator "Majority_HGG"
      out <- rma(yi = yi, vi = vi, mods = ~ Majority_HGG, dat = sub_means)
      
      ### Store results
      long_means$b0_rb[m] <- out$b[1]
      long_means$b0_se_rb[m] <- out$se[1]
      long_means$b0_lb_rb[m] <- out$ci.lb[1]
      long_means$b0_ub_rb[m] <- out$ci.ub[1]
      long_means$b0_zval_rb[m] <- out$zval[1]
      long_means$b0_pval_rb[m] <- out$pval[1]
      long_means$b1_rb[m] <- out$b[2]
      long_means$b1_se_rb[m] <- out$se[2]
      long_means$b1_lb_rb[m] <- out$ci.lb[2]
      long_means$b1_ub_rb[m] <- out$ci.ub[2]
      long_means$b1_zval_rb[m] <- out$zval[2]
      long_means$b1_pval_rb[m] <- out$pval[2]
    }
    
  }
  
  ##############################################################################
  
  ### Select data based on z-scores
  sub_zscore <- subset(tmp, es_type == "z-score")
  
  ### Store type of test and number of effect sizes
  long_zscore$type_test[m] <- names(splits)[[m]]
  long_zscore$k[m] <- nrow(sub_zscore)
  
  if (nrow(sub_zscore) > 0)
  {
    ### Conduct a random-effects meta-analysis
    out <- rma(yi = yi, vi = vi, dat = sub_zscore)
    
    ### Store results
    long_zscore$est[m] <- out$b[1]
    long_zscore$se[m] <- out$se
    long_zscore$lb[m] <- out$ci.lb
    long_zscore$ub[m] <- out$ci.ub
    long_zscore$zval[m] <- out$zval
    long_zscore$pval[m] <- out$pval
    long_zscore$tau2[m] <- out$tau2
    long_zscore$se_tau2[m] <- out$se.tau2
    long_zscore$Qval[m] <- out$QE
    long_zscore$pval_Q[m] <- out$QEp
    long_zscore$I2[m] <- out$I2
    
    ### Compute sum of harmonic mean of sample sizes
    long_zscore$sum_n[m] <- round(sum(apply(rbind(sub_zscore$n_baseline, sub_zscore$n_post), 
                                            2, function(x) 1/mean(1/x))))
    
    ### Store domain number for ordering rows in table
    long_zscore$domain_number[m] <- unique(sub_zscore$Domain_number)
    
    if (out$k > 1)
    {
      ci <- confint(out) # 95% CI tau^2
      long_zscore$lb_tau2[m] <- paste0(ci$lb.sign, ci$random["tau^2","ci.lb"])
      long_zscore$ub_tau2[m] <- paste0(ci$ub.sign, ci$random["tau^2","ci.ub"])
    }
    
    ### Conduct an equal-effect meta-analysis
    out <- rma(yi = yi, vi = vi, dat = sub_zscore, method = "EE")
    
    ### Store results
    long_zscore$est_ee[m] <- out$b[1]
    long_zscore$se_ee[m] <- out$se
    long_zscore$lb_ee[m] <- out$ci.lb
    long_zscore$ub_ee[m] <- out$ci.ub
    long_zscore$zval_ee[m] <- out$zval
    long_zscore$pval_ee[m] <- out$pval
  }
  
  ### Meta-analysis without studies at high risk of bias
  sub_zscore_risk <- subset(sub_zscore, sub_zscore$RiskOfBias2Risks == 0)
  
  ### Store number of studies in meta-analysis without studies at high risk of bias
  long_zscore$k_risk[m] <- nrow(sub_zscore_risk)
  
  if (nrow(sub_zscore_risk) > 0)
  {
    
    ### Conduct a random-effects meta-analysis
    out <- rma(yi = yi, vi = vi, dat = sub_zscore_risk)
    
    ### Store results
    long_zscore$est_risk[m] <- out$b[1]
    long_zscore$se_risk[m] <- out$se
    long_zscore$lb_risk[m] <- out$ci.lb
    long_zscore$ub_risk[m] <- out$ci.ub
    long_zscore$zval_risk[m] <- out$zval
    long_zscore$pval_risk[m] <- out$pval
    long_zscore$tau2_risk[m] <- out$tau2
    long_zscore$se_tau2_risk[m] <- out$se.tau2
    long_zscore$Qval_risk[m] <- out$QE
    long_zscore$pval_Q_risk[m] <- out$QEp
    long_zscore$I2_risk[m] <- out$I2
    
    ### Compute sum of harmonic mean of sample sizes
    long_zscore$sum_n_risk[m] <- round(sum(apply(rbind(sub_zscore_risk$n_baseline, 
                                                      sub_zscore_risk$n_post), 
                                                2, function(x) 1/mean(1/x))))
    
    if (out$k > 1)
    {
      ci <- confint(out) # 95% CI tau^2
      long_zscore$lb_tau2_risk[m] <- paste0(ci$lb.sign, ci$random["tau^2","ci.lb"])
      long_zscore$ub_tau2_risk[m] <- paste0(ci$ub.sign, ci$random["tau^2","ci.ub"])
    }
  }
  
  ### Number of studies that score 1 on "Majority_HGG"
  long_zscore$k1_rb[m] <- sum(sub_zscore$Majority_HGG == 1, na.rm = TRUE)
  
  ### If there are studies with scores 0 and 1 on the moderator "Majority_HGG"
  # and there are more than two studies, conduct the moderator analysis
  if (all(c(0,1) %in% unique(sub_zscore$Majority_HGG)) & 
      (sum(sub_zscore$Majority_HGG == 0, na.rm = TRUE)+
       sum(sub_zscore$Majority_HGG == 1, na.rm = TRUE)) > 2)
  {
    ### Conduct a random-effects meta-analysis with moderator "Majority_HGG"
    out <- rma(yi = yi, vi = vi, mods = ~ Majority_HGG, dat = sub_zscore)
    
    ### Store results
    long_zscore$b0_rb[m] <- out$b[1]
    long_zscore$b0_se_rb[m] <- out$se[1]
    long_zscore$b0_lb_rb[m] <- out$ci.lb[1]
    long_zscore$b0_ub_rb[m] <- out$ci.ub[1]
    long_zscore$b0_zval_rb[m] <- out$zval[1]
    long_zscore$b0_pval_rb[m] <- out$pval[1]
    long_zscore$b1_rb[m] <- out$b[2]
    long_zscore$b1_se_rb[m] <- out$se[2]
    long_zscore$b1_lb_rb[m] <- out$ci.lb[2]
    long_zscore$b1_ub_rb[m] <- out$ci.ub[2]
    long_zscore$b1_zval_rb[m] <- out$zval[2]
    long_zscore$b1_pval_rb[m] <- out$pval[2]
  }
  
}

### Save results
save(long_means, file = "long_means.RData")
save(long_zscore, file = "long_zscore.RData")

################################################################################

##### SUMMARIZE CROSS-SECTIONAL STUDIES PER TYPE OF TEST #####

### Fix some spelling mistakes in the data
dat_cross$type_test <- ifelse(dat_cross$type_test == "digit span forward", 
                              "Digit span forward", dat_cross$type_test)

### Split data per type of test
splits <- split(dat_cross, dat_cross$type_test)

### Empty objects for storing results
cross_means <- data.frame(type_test = character(length(splits)), 
                          domain_number = rep(NA, length(splits)),
                          k = numeric(length(splits)), 
                          sum_n = numeric(length(splits)), 
                          est = numeric(length(splits)), 
                          se = numeric(length(splits)),
                          lb = numeric(length(splits)), 
                          ub = numeric(length(splits)),
                          zval = numeric(length(splits)), 
                          pval = numeric(length(splits)),
                          tau2 = numeric(length(splits)), 
                          se_tau2 = numeric(length(splits)), 
                          lb_tau2 = rep(NA, length(splits)), 
                          ub_tau2 = rep(NA, length(splits)), 
                          Qval = numeric(length(splits)), 
                          pval_Q = numeric(length(splits)), 
                          I2 = numeric(length(splits)),
                          est_ee = rep(NA, length(splits)), 
                          se_ee = rep(NA, length(splits)),
                          lb_ee = rep(NA, length(splits)), 
                          ub_ee = rep(NA, length(splits)),
                          zval_ee = rep(NA, length(splits)), 
                          pval_ee = rep(NA, length(splits)),
                          k_risk = numeric(length(splits)), 
                          sum_n_risk = numeric(length(splits)), 
                          est_risk = numeric(length(splits)), 
                          se_risk = numeric(length(splits)),
                          lb_risk = numeric(length(splits)), 
                          ub_risk = numeric(length(splits)),
                          zval_risk = numeric(length(splits)), 
                          pval_risk = numeric(length(splits)),
                          tau2_risk = numeric(length(splits)), 
                          se_tau2_risk = numeric(length(splits)), 
                          lb_tau2_risk = rep(NA, length(splits)), 
                          ub_tau2_risk = rep(NA, length(splits)), 
                          Qval_risk = numeric(length(splits)), 
                          pval_Q_risk = numeric(length(splits)), 
                          I2_risk = numeric(length(splits)),
                          k1_rb = rep(NA, length(splits)), 
                          b0_rb = rep(NA, length(splits)), 
                          b0_se_rb = rep(NA, length(splits)),
                          b0_lb_rb = rep(NA, length(splits)), 
                          b0_ub_rb = rep(NA, length(splits)),
                          b0_zval_rb = rep(NA, length(splits)), 
                          b0_pval_rb = rep(NA, length(splits)),
                          b1_rb = rep(NA, length(splits)), 
                          b1_se_rb = rep(NA, length(splits)),
                          b1_lb_rb = rep(NA, length(splits)), 
                          b1_ub_rb = rep(NA, length(splits)),
                          b1_zval_rb = rep(NA, length(splits)), 
                          b1_pval_rb = rep(NA, length(splits)))

cross_zscore <- data.frame(type_test = character(length(splits)), 
                           domain_number = rep(NA, length(splits)),
                           k = rep(NA, length(splits)), 
                           sum_n = numeric(length(splits)), 
                           est = rep(NA, length(splits)), 
                           se = rep(NA, length(splits)),
                           lb = rep(NA, length(splits)), 
                           ub = rep(NA, length(splits)),
                           zval = rep(NA, length(splits)), 
                           pval = rep(NA, length(splits)),
                           tau2 = rep(NA, length(splits)), 
                           se_tau2 = numeric(length(splits)), 
                           lb_tau2 = rep(NA, length(splits)), 
                           ub_tau2 = rep(NA, length(splits)), 
                           Qval = rep(NA, length(splits)), 
                           pval_Q = rep(NA, length(splits)), 
                           I2 = rep(NA, length(splits)),
                           est_ee = rep(NA, length(splits)), 
                           se_ee = rep(NA, length(splits)),
                           lb_ee = rep(NA, length(splits)), 
                           ub_ee = rep(NA, length(splits)),
                           zval_ee = rep(NA, length(splits)), 
                           pval_ee = rep(NA, length(splits)),
                           k_risk = numeric(length(splits)), 
                           sum_n_risk = numeric(length(splits)), 
                           est_risk = numeric(length(splits)), 
                           se_risk = numeric(length(splits)),
                           lb_risk = numeric(length(splits)), 
                           ub_risk = numeric(length(splits)),
                           zval_risk = numeric(length(splits)), 
                           pval_risk = numeric(length(splits)),
                           tau2_risk = numeric(length(splits)), 
                           se_tau2_risk = numeric(length(splits)), 
                           lb_tau2_risk = rep(NA, length(splits)), 
                           ub_tau2_risk = rep(NA, length(splits)), 
                           Qval_risk = numeric(length(splits)), 
                           pval_Q_risk = numeric(length(splits)), 
                           I2_risk = numeric(length(splits)),
                           k1_rb = rep(NA, length(splits)), 
                           b0_rb = rep(NA, length(splits)), 
                           b0_se_rb = rep(NA, length(splits)),
                           b0_lb_rb = rep(NA, length(splits)), 
                           b0_ub_rb = rep(NA, length(splits)),
                           b0_zval_rb = rep(NA, length(splits)), 
                           b0_pval_rb = rep(NA, length(splits)),
                           b1_rb = rep(NA, length(splits)), 
                           b1_se_rb = rep(NA, length(splits)),
                           b1_lb_rb = rep(NA, length(splits)), 
                           b1_ub_rb = rep(NA, length(splits)),
                           b1_zval_rb = rep(NA, length(splits)), 
                           b1_pval_rb = rep(NA, length(splits)))

for (m in 1:length(splits))
{
  tmp <- splits[[m]]
  
  ### Select data based on means
  sub_means <- subset(tmp, es_type == "means")
  
  ### Store type of test and number of effect sizes
  cross_means$type_test[m] <- names(splits)[[m]]
  cross_means$k[m] <- nrow(sub_means)
  
  if (nrow(sub_means) > 0)
  {
    ### Conduct a random-effects meta-analysis
    out <- rma(yi = yi, vi = vi, dat = sub_means)
    
    ### Store results
    cross_means$est[m] <- out$b[1]
    cross_means$se[m] <- out$se
    cross_means$lb[m] <- out$ci.lb
    cross_means$ub[m] <- out$ci.ub
    cross_means$zval[m] <- out$zval
    cross_means$pval[m] <- out$pval
    cross_means$tau2[m] <- out$tau2
    cross_means$se_tau2[m] <- out$se.tau2
    cross_means$Qval[m] <- out$QE
    cross_means$pval_Q[m] <- out$QEp
    cross_means$I2[m] <- out$I2
    
    ### Sum of the sample sizes in both groups and across studies
    cross_means$sum_n[m] <- sum(sub_means$n_contr+sub_means$n_interv)
    
    ### Store domain number for ordering rows in table
    cross_means$domain_number[m] <- unique(sub_means$Domain_number)
    
    if (out$k > 1)
    {
      ci <- confint(out) # 95% CI tau^2
      cross_means$lb_tau2[m] <- paste0(ci$lb.sign, ci$random["tau^2","ci.lb"])
      cross_means$ub_tau2[m] <- paste0(ci$ub.sign, ci$random["tau^2","ci.ub"])
    }
    
    ### Conduct an equal-effect meta-analysis
    out <- rma(yi = yi, vi = vi, dat = sub_means, method = "EE")
    
    ### Store results
    cross_means$est_ee[m] <- out$b[1]
    cross_means$se_ee[m] <- out$se
    cross_means$lb_ee[m] <- out$ci.lb
    cross_means$ub_ee[m] <- out$ci.ub
    cross_means$zval_ee[m] <- out$zval
    cross_means$pval_ee[m] <- out$pval
    
    ### Meta-analysis without studies at high risk of bias
    sub_means_risk <- subset(sub_means, sub_means$"RiskOfBias2Risks" == 0)
    
    ### Store number of studies in meta-analysis without studies at high risk of bias
    cross_means$k_risk[m] <- nrow(sub_means_risk)
    
    if (nrow(sub_means_risk) > 0)
    {
      
      ### Conduct a random-effects meta-analysis
      out <- rma(yi = yi, vi = vi, dat = sub_means_risk)
      
      ### Store results
      cross_means$est_risk[m] <- out$b[1]
      cross_means$se_risk[m] <- out$se
      cross_means$lb_risk[m] <- out$ci.lb
      cross_means$ub_risk[m] <- out$ci.ub
      cross_means$zval_risk[m] <- out$zval
      cross_means$pval_risk[m] <- out$pval
      cross_means$tau2_risk[m] <- out$tau2
      cross_means$se_tau2_risk[m] <- out$se.tau2
      cross_means$Qval_risk[m] <- out$QE
      cross_means$pval_Q_risk[m] <- out$QEp
      cross_means$I2_risk[m] <- out$I2
      
      ### Sum of the sample sizes in both groups and across studies
      cross_means$sum_n_risk[m] <- sum(sub_means_risk$n_contr+sub_means_risk$n_interv)
      
      if (out$k > 1)
      {
        ci <- confint(out) # 95% CI tau^2
        cross_means$lb_tau2_risk[m] <- paste0(ci$lb.sign, ci$random["tau^2","ci.lb"])
        cross_means$ub_tau2_risk[m] <- paste0(ci$ub.sign, ci$random["tau^2","ci.ub"])
      }
    }
    
    ### Number of studies that score 1 on "Majority_HGG"
    cross_means$k1_rb[m] <- sum(sub_means$Majority_HGG == 1, na.rm = TRUE)
    
    ### If there are studies with scores 0 and 1 on the moderator "Majority_HGG"
    # and there are more than two studies, conduct the moderator analysis
    if (all(c(0,1) %in% unique(sub_means$Majority_HGG)) & 
        (sum(sub_means$Majority_HGG == 0, na.rm = TRUE)+
         sum(sub_means$Majority_HGG == 1, na.rm = TRUE)) > 2)
    {
      ### Conduct a random-effects meta-analysis with moderator "Majority_HGG"
      out <- rma(yi = yi, vi = vi, mods = ~ Majority_HGG, dat = sub_means)
      
      ### Store results
      cross_means$b0_rb[m] <- out$b[1]
      cross_means$b0_se_rb[m] <- out$se[1]
      cross_means$b0_lb_rb[m] <- out$ci.lb[1]
      cross_means$b0_ub_rb[m] <- out$ci.ub[1]
      cross_means$b0_zval_rb[m] <- out$zval[1]
      cross_means$b0_pval_rb[m] <- out$pval[1]
      cross_means$b1_rb[m] <- out$b[2]
      cross_means$b1_se_rb[m] <- out$se[2]
      cross_means$b1_lb_rb[m] <- out$ci.lb[2]
      cross_means$b1_ub_rb[m] <- out$ci.ub[2]
      cross_means$b1_zval_rb[m] <- out$zval[2]
      cross_means$b1_pval_rb[m] <- out$pval[2]
    }
    
  }
  
  ##############################################################################
  
  ### Select data based on z-scores
  sub_zscore <- subset(tmp, es_type == "z-score")
  
  ### Store type of test and number of effect sizes
  cross_zscore$type_test[m] <- names(splits)[[m]]
  cross_zscore$k[m] <- nrow(sub_zscore)
  
  if (nrow(sub_zscore) > 0)
  {
    ### Conduct a random-effects meta-analysis
    out <- rma(yi = yi, vi = vi, dat = sub_zscore)
    
    ### Store results
    cross_zscore$est[m] <- out$b[1]
    cross_zscore$se[m] <- out$se
    cross_zscore$lb[m] <- out$ci.lb
    cross_zscore$ub[m] <- out$ci.ub
    cross_zscore$zval[m] <- out$zval
    cross_zscore$pval[m] <- out$pval
    cross_zscore$tau2[m] <- out$tau2
    cross_zscore$se_tau2[m] <- out$se.tau2
    cross_zscore$Qval[m] <- out$QE
    cross_zscore$pval_Q[m] <- out$QEp
    cross_zscore$I2[m] <- out$I2
    
    ### Sum of the sample sizes of the intervention group across studies
    cross_zscore$sum_n[m] <- sum(sub_zscore$n_interv, na.rm = TRUE)
    
    ### Store domain number for ordering rows in table
    cross_zscore$domain_number[m] <- unique(sub_zscore$Domain_number)
    
    if (out$k > 1)
    {
      ci <- confint(out) # 95% CI tau^2
      cross_zscore$lb_tau2[m] <- paste0(ci$lb.sign, ci$random["tau^2","ci.lb"])
      cross_zscore$ub_tau2[m] <- paste0(ci$ub.sign, ci$random["tau^2","ci.ub"])
    }
    
    ### Conduct an equal-effect meta-analysis
    out <- rma(yi = yi, vi = vi, dat = sub_zscore, method = "EE")
    
    ### Store results
    cross_zscore$est_ee[m] <- out$b[1]
    cross_zscore$se_ee[m] <- out$se
    cross_zscore$lb_ee[m] <- out$ci.lb
    cross_zscore$ub_ee[m] <- out$ci.ub
    cross_zscore$zval_ee[m] <- out$zval
    cross_zscore$pval_ee[m] <- out$pval
    
    ### Meta-analysis without studies at high risk of bias
    sub_zscore_risk <- subset(sub_zscore, sub_zscore$"RiskOfBias2Risks" == 0)
    
    ### Store number of studies in meta-analysis without studies at high risk of bias
    cross_zscore$k_risk[m] <- nrow(sub_zscore_risk)
    
    if (nrow(sub_zscore_risk) > 0)
    {
      
      ### Conduct a random-effects meta-analysis
      out <- rma(yi = yi, vi = vi, dat = sub_zscore_risk)
      
      ### Store results
      cross_zscore$est_risk[m] <- out$b[1]
      cross_zscore$se_risk[m] <- out$se
      cross_zscore$lb_risk[m] <- out$ci.lb
      cross_zscore$ub_risk[m] <- out$ci.ub
      cross_zscore$zval_risk[m] <- out$zval
      cross_zscore$pval_risk[m] <- out$pval
      cross_zscore$tau2_risk[m] <- out$tau2
      cross_zscore$se_tau2_risk[m] <- out$se.tau2
      cross_zscore$Qval_risk[m] <- out$QE
      cross_zscore$pval_Q_risk[m] <- out$QEp
      cross_zscore$I2_risk[m] <- out$I2
      
      ### Sum of the sample sizes in both groups and across studies
      cross_zscore$sum_n_risk[m] <- sum(sub_zscore_risk$n_interv, na.rm = TRUE)
      
      if (out$k > 1)
      {
        ci <- confint(out) # 95% CI tau^2
        cross_zscore$lb_tau2_risk[m] <- paste0(ci$lb.sign, ci$random["tau^2","ci.lb"])
        cross_zscore$ub_tau2_risk[m] <- paste0(ci$ub.sign, ci$random["tau^2","ci.ub"])
      }
    }
    
    ### Number of studies that score 1 on "Majority_HGG"
    cross_zscore$k1_rb[m] <- sum(sub_zscore$Majority_HGG == 1, na.rm = TRUE)
    
    ### If there are studies with scores 0 and 1 on the moderator "Majority_HGG"
    # and there are more than two studies, conduct the moderator analysis
    if (all(c(0,1) %in% unique(sub_zscore$Majority_HGG)) & 
        (sum(sub_zscore$Majority_HGG == 0, na.rm = TRUE)+
         sum(sub_zscore$Majority_HGG == 1, na.rm = TRUE)) > 2)
    {
      ### Conduct a random-effects meta-analysis with moderator "Majority_HGG"
      out <- rma(yi = yi, vi = vi, mods = ~ Majority_HGG, dat = sub_zscore)
      
      ### Store results
      cross_zscore$b0_rb[m] <- out$b[1]
      cross_zscore$b0_se_rb[m] <- out$se[1]
      cross_zscore$b0_lb_rb[m] <- out$ci.lb[1]
      cross_zscore$b0_ub_rb[m] <- out$ci.ub[1]
      cross_zscore$b0_zval_rb[m] <- out$zval[1]
      cross_zscore$b0_pval_rb[m] <- out$pval[1]
      cross_zscore$b1_rb[m] <- out$b[2]
      cross_zscore$b1_se_rb[m] <- out$se[2]
      cross_zscore$b1_lb_rb[m] <- out$ci.lb[2]
      cross_zscore$b1_ub_rb[m] <- out$ci.ub[2]
      cross_zscore$b1_zval_rb[m] <- out$zval[2]
      cross_zscore$b1_pval_rb[m] <- out$pval[2]
    }
    
  }
}

### Save results
save(cross_means, file = "cross_means.RData")
save(cross_zscore, file = "cross_zscore.RData")