
lapply(c("data.table", "tidyverse", "magrittr", "survey", "stargazer", "sandwich", "lmtest"), library, ch = TRUE)
options(scipen = 999)

full.DT <- readxl::read_xlsx("PNAS_toShare.xlsx", sheet = "violent_crime") %>% as.data.table()

#----
# Correlation Matrix
cor.mat.vars <- full.DT[, .(log.vcr.08, log.hr.08, mob_adjusted, poverty.rate.08, gini.09, leo_per_100k_08, unemp.rate.08)]

cor.mat <- Hmisc::rcorr(as.matrix(cor.mat.vars))$r
cor.mat

xtable::xtable(cor.mat) %>% print(type = "html")

#----
# Create categorical variables for mobility and poverty

full.DT[, .(quantile(mob_adjusted, c(0.33, 0.67), na.rm = TRUE))]
full.DT[, .(quantile(poverty.rate.08, c(0.33, 0.67), na.rm = TRUE))]
full.DT[, .(quantile(log.vcr.08, c(0.33, 0.67), na.rm = TRUE))]


full.DT[, ':='(mob_cat = ifelse(mob_adjusted < 43.544, "Low", ifelse(mob_adjusted < 47.595, "MEDIUM", ifelse(mob_adjusted >= 47.595, "High", "NA"))),
               viol_cat = ifelse(log.vcr.08 < 4.880, "Low", ifelse(log.vcr.08 < 5.710, "Medium", ifelse(log.vcr.08 >= 5.710, "High", "NA"))))]

full.DT[viol_cat == "High" & mob_cat == "Low"][, .(uniqueN(FIPS))]
full.DT[viol_cat == "Low" & mob_cat == "High"][, .(uniqueN(FIPS))]
full.DT[is.na(mob_adjusted)]

#-----
full.DT[, ':='(mob_cat = ifelse(mob_adjusted < 43.544, "Low", ifelse(mob_adjusted < 47.595, "MEDIUM", ifelse(mob_adjusted >= 47.595, "High", "NA"))),
               pov_cat = ifelse(poverty.rate.08 < 11.8, "Low", ifelse(poverty.rate.08 < 16.8, "Medium", ifelse(poverty.rate.08 >= 16.8, "High", "NA"))))]

full.DT[pov_cat == "High" & mob_cat == "Low"][, .(uniqueN(FIPS))]
full.DT[pov_cat == "Low" & mob_cat == "High"][, .(uniqueN(FIPS))]
full.DT[is.na(mob_adjusted)]

#--------------------------------------------------
# summary statistics
summ.stats <- full.DT[, .(mean = lapply(.SD, mean, na.rm = TRUE), sd = lapply(.SD, sd, na.rm = TRUE)), .SDcols = setdiff(names(full.DT), c("FIPS", "violent_crime_rate", "homicide_rate", "COVIND", "state"))]
summ.stats <- cbind(setdiff(names(full.DT), c("FIPS", "violent_crime_rate", "homicide_rate", "COVIND", "state")), summ.stats)
summ.stats


#---- Set up OLS models ----
# robust standard errors function
se_robust <- function(x){
  coeftest(x, vcov = vcovCL, type = "HC1")[, "Std. Error"]
}

#---- state fixed effects ----
se_robust_fe <- function(x){
  coeftest(x, vcov = vcovCL, type = "HC1", cluster = ~stateFIPS, df = uniqueN(x$model$`factor(stateFIPS)`) - 1)[, "Std. Error"]
}

#----
choose_robust_func <- function(model){
  if (str_detect(as.character(model$call[2]), "factor")) {
    se_robust_fe(model)
  }
  else {
    se_robust(model)
  }
}

# Hommel adjustments for p-values ----
p.adjust.hommel <- function(model){
  orig.pval <- summary(model)$coefficients[,4] %>% as.data.table(keep.rownames = TRUE)
  setnames(orig.pval, new = c("rn", "orig.p"))
  adjusted.p <- stats::p.adjust(orig.pval$orig.p, method = "hommel")
  names(adjusted.p) <- orig.pval$rn
  return(adjusted.p)
}



#---- violent crime rate ----
model1.vcr <- lm(log.vcr.08 ~ mob_adjusted, data = full.DT)
model2.vcr <- lm(log.vcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + urban06 + gini.09 + pct.singleHoH.09, data = full.DT)
model3.vcr <- lm(log.vcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + urban06 + gini.09 + pct.singleHoH.09 + pct.pop.diff.0708 + pct.samehouse.08.09 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08, data = full.DT)
model4.vcr <- lm(log.vcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + urban06 + gini.09 + pct.singleHoH.09 + pct.pop.diff.0708 + pct.samehouse.08.09 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08 + log_vcr_rollmean_1997, data = full.DT)
model.4.vcr.fe <- lm(log.vcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + urban06 + gini.09 + pct.singleHoH.09 + pct.pop.diff.0708 + pct.samehouse.08.09 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08 + log_vcr_rollmean_1997 + factor(stateFIPS), data = full.DT)

models.all.vcr <- list(model2.vcr, model3.vcr, model4.vcr, model.4.vcr.fe)

list.tot.se.vcr <- rlist::list.append(lapply(models.all.vcr, choose_robust_func))

list.tot.p.vcr <- lapply(models.all.vcr, p.adjust.hommel)

stargazer(models.all.vcr, type = "text", align = TRUE,
          # uncomment to write table to file
          # out = "model_tables/vcr_full_122123.html",
          omit = "state",
          se = c(list.tot.se.vcr),
          p = c(list.tot.p.vcr),
          dep.var.caption = "Log(Violent Crime Rate + 1)",
          dep.var.labels.include = FALSE,
          order = c(1,2,9),
          covariate.labels = c("Mobility",  "Poverty",  "GINI", "LEOs Per 100,000", "Unemployment Rate",  "Population Per 100,000",  "Percent Non-Hisp. Black", "Percent Male 15-24 Pop.", "Urban", "Percent Single Parent Head of Household",
                               "Percent Difference in Population from the Previous Year", "Percent in Same County from the Previous Year",  "Median Number of Years Lived in Current House",  "Percent Owner Households", "Percent Occupied Households",  "Percent on Housing Vouchers",  "3-Year Average of Past Violence (1996 - 1998)"),
          add.lines = list(c("State Fixed-Effects", rep("No", 3), "Yes")),
          title = "Estimation Results for OLS Models Predicting Violent Crime in 2008 (per Uniform Crime Reports)",
          notes = "Reports robust SEs and uses Hommel adjustments for p-values")



#---- homicide rate ----
model1.hr <- lm(log.hr.08 ~ mob_adjusted, data = full.DT)
model2.hr <- lm(log.hr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + urban06 + gini.09 + pct.singleHoH.09, data = full.DT)
model3.hr <- lm(log.hr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + urban06 + gini.09 + pct.singleHoH.09 + pct.pop.diff.0708 + pct.samehouse.08.09 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08, data = full.DT)
model4.hr <- lm(log.hr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + urban06 + gini.09 + pct.singleHoH.09 + pct.pop.diff.0708 + pct.samehouse.08.09 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08 + log_hr_rollmean_1997, data = full.DT)
model.4.hr.fe <- lm(log.hr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + urban06 + gini.09 + pct.singleHoH.09 + pct.pop.diff.0708 + pct.samehouse.08.09 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08 + log_hr_rollmean_1997 + factor(stateFIPS), data = full.DT)

models.all.hr <- list(model2.hr, model3.hr, model4.hr, model.4.hr.fe)
list.tot.se.hr <- rlist::list.append(lapply(models.all.hr, choose_robust_func))
list.tot.p.hr <- lapply(models.all.hr, p.adjust.hommel)

stargazer(models.all.hr, type = "text", align = TRUE,
          # uncomment to write table to file
          # out = "model_tables/hr_full_122123.html",
          omit = "state",
          se = c(list.tot.se.hr),
          p = c(list.tot.p.hr),
          dep.var.caption = "Log(Homicide Rate + 1)",
          dep.var.labels.include = FALSE,
          order = c(1,2,9),
          covariate.labels = c("Mobility",  "Poverty",  "GINI", "LEOs Per 100,000", "Unemployment Rate",  "Population Per 100,000",  "Percent Non-Hisp. Black", "Percent Male 15-24 Pop.", "Urban", "Percent Single Parent Head of Household",
                               "Percent Difference in Population from the Previous Year", "Percent in Same County from the Previous Year",  "Median Number of Years Lived in Current House",  "Percent Owner Households", "Percent Occupied Households",  "Percent on Housing Vouchers",  "3-Year Average of Past Violence (1996 - 1998)"),
          add.lines = list(c("State Fixed-Effects", rep("No", 3), "Yes")),
          title = "Estimation Results for OLS Models Predicting Homicide Rates in 2008 (per Uniform Crime Reports)",
          notes = "Reports robust SEs and uses Hommel adjustments for p-values")

rm(model1.vcr, model2.vcr, model3.vcr, model4.vcr, model.4.vcr.fe,
   model1.hr, model2.hr, model3.hr, model4.hr, model.4.hr.fe)

#---- Bivariate Models ----
# Violent crime rate
model1.vcr <- lm(log.vcr.08 ~ mob_adjusted, data = full.DT)
model2.vcr <- lm(log.vcr.08 ~ poverty.rate.08, data = full.DT)
model3.vcr <- lm(log.vcr.08 ~ gini.09, data = full.DT)
model4.vcr <- lm(log.vcr.08 ~ leo_per_100k_08, data = full.DT)
model5.vcr <- lm(log.vcr.08 ~ unemp.rate.08, data = full.DT)
model6.vcr <- lm(log.vcr.08 ~ pct.nh.black.08, data = full.DT)

# Homicide rate
model1.hr <- lm(log.hr.08 ~ mob_adjusted, data = full.DT)
model2.hr <- lm(log.hr.08 ~ poverty.rate.08, data = full.DT)
model3.hr <- lm(log.hr.08 ~ gini.09, data = full.DT)
model4.hr <- lm(log.hr.08 ~ leo_per_100k_08, data = full.DT)
model5.hr <- lm(log.hr.08 ~ unemp.rate.08, data = full.DT)
model6.hr <- lm(log.hr.08 ~ pct.nh.black.08, data = full.DT)

models.all <- list(model1.vcr, model2.vcr, model3.vcr, model4.vcr, model5.vcr, model6.vcr,
                   model1.hr, model2.hr, model3.hr, model4.hr, model5.hr, model6.hr)

# list.tot.se.vcr <- lapply(models.all.vcr, choose_robust_func, clust_level = "state")
list.tot.se <- rlist::list.append(lapply(models.all, choose_robust_func))

list.tot.p <- lapply(models.all, p.adjust.hommel)

stargazer(models.all, type = "html", align = TRUE,
          # out = "model_tables/all_bivar_122123.html",
          omit = "state",
          se = c(list.tot.se),
          p = c(list.tot.p),
          # dep.var.caption = "Log(Violent Crime Rate + 1)",
          dep.var.labels.include = TRUE,
          # order = c(1,2,9),
          covariate.labels = c("Mobility",  "Poverty",  "GINI", "LEOs Per 100,000", "Unemployment Rate", "Proportion Non-Hispanic Black"),
          # add.lines = list(c("State Fixed-Effects", rep("No", 6), rep("Yes", 6))),
          title = "Estimation Results for OLS Models Predicting Violent Crime in 2008 (per Uniform Crime Reports)",
          notes = "Reports robust SEs and uses Hommel adjustments for p-values")



#---- Urban Counties ----
#---- violent crime rate ----
model1.vcr <- lm(log.vcr.08 ~ mob_adjusted, data = full.DT[urban06 == 1])
model2.vcr <- lm(log.vcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + gini.09 + pct.singleHoH.09, data = full.DT[urban06 == 1])
model3.vcr <- lm(log.vcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + gini.09 + pct.singleHoH.09 + pct.samehouse.08.09 + pct.pop.diff.0708 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08, data = full.DT[urban06 == 1])
model4.vcr <- lm(log.vcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + gini.09 + pct.singleHoH.09 + pct.samehouse.08.09 + pct.pop.diff.0708 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08 + log_vcr_rollmean_1997, data = full.DT[urban06 == 1])
model.4.vcr.fe <- lm(log.vcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + gini.09 + pct.singleHoH.09 + pct.samehouse.08.09 + pct.pop.diff.0708 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08 + log_vcr_rollmean_1997 + factor(stateFIPS), data = full.DT[urban06 == 1])

models.all.vcr <- list(model2.vcr, model3.vcr, model4.vcr, model.4.vcr.fe)

list.tot.se.vcr <- rlist::list.append(lapply(models.all.vcr, choose_robust_func))

list.tot.p.vcr <- lapply(models.all.vcr, p.adjust.hommel)

stargazer(models.all.vcr, type = "text", align = TRUE,
          # uncomment to write table to file
          # out = "model_tables/vcr_popscaled_urb_121823.html",
          keep = "mob_adjusted",
          se = c(list.tot.se.vcr),
          p = c(list.tot.p.vcr),
          dep.var.caption = "Log(Violent Crime Rate + 1)",
          dep.var.labels.include = FALSE,
          order = c(1,2,9),
          covariate.labels = c("Mobility",  "Poverty",  "GINI", "LEOs Per 100,000", "Unemployment Rate",  "Population Per 100,000",  "Percent Non-Hisp. Black", "Percent Male 15-24 Pop.", "Percent Single Parent Head of Household",
                               "Percent in Same House from 2008 to 2009",  "Percent Difference in Population from 2007 to 2008",  "Median Number of Years Lived in Current House",  "Percent Owner Households", "Percent Occupied Households",  "Percent on Housing Vouchers",  "3-Year Average of Past Violence (1996 - 1998)"),
          add.lines = list(c("State Fixed-Effects", rep("No", 3), "Yes")),
          title = "Estimation Results for OLS Models Predicting Violent Crime in 2008 for Urban Counties (per Uniform Crime Reports)",
          notes = "Reports robust SEs and uses Hommel adjustments for p-values")



#---- homicide rate ----
model1.hr <- lm(log.hr.08 ~ mob_adjusted, data = full.DT[urban06 == 1])
model2.hr <- lm(log.hr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + gini.09 + pct.singleHoH.09, data = full.DT[urban06 == 1])
model3.hr <- lm(log.hr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + gini.09 + pct.singleHoH.09 + pct.samehouse.08.09 + pct.pop.diff.0708 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08, data = full.DT[urban06 == 1])
model4.hr <- lm(log.hr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + gini.09 + pct.singleHoH.09 + pct.samehouse.08.09 + pct.pop.diff.0708 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08 + log_hr_rollmean_1997, data = full.DT[urban06 == 1])
model.4.hr.fe <- lm(log.hr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + gini.09 + pct.singleHoH.09 + pct.samehouse.08.09 + pct.pop.diff.0708 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08 + log_hr_rollmean_1997 + factor(stateFIPS), data = full.DT[urban06 == 1])

models.all.hr <- list(model2.hr, model3.hr, model4.hr, model.4.hr.fe)
list.tot.se.hr <- rlist::list.append(lapply(models.all.hr, choose_robust_func))
list.tot.p.hr <- lapply(models.all.hr, p.adjust.hommel)

stargazer(models.all.hr, type = "text", align = TRUE,
          # uncomment to write table to file
          # out = "model_tables/hr_popscaled_urb_121823.html", omit = "state",
          keep = "mob_adjusted",
          se = c(list.tot.se.hr),
          p = c(list.tot.p.hr),
          dep.var.caption = "Log(Homicide Rate + 1)",
          dep.var.labels.include = FALSE,
          order = c(1,2,9),
          covariate.labels = c("Mobility",  "Poverty",  "GINI", "LEOs Per 100,000", "Unemployment Rate",  "Population Per 100,000",  "Percent Non-Hisp. Black", "Percent Male 15-24 Pop.", "Percent Single Parent Head of Household",
                               "Percent in Same House from 2008 to 2009",  "Percent Difference in Population from 2007 to 2008",  "Median Number of Years Lived in Current House",  "Percent Owner Households", "Percent Occupied Households",  "Percent on Housing Vouchers",  "3-Year Average of Past Homicide (1996 - 1998)"),
          add.lines = list(c("State Fixed-Effects", rep("No", 3), "Yes")),
          title = "Estimation Results for OLS Models Predicting Homicide Rates in 2008 For Urban Counties (per Uniform Crime Reports)",
          notes = "Reports robust SEs and uses Hommel adjustments for p-values")

rm(model1.vcr, model2.vcr, model3.vcr, model4.vcr, model.4.vcr.fe,
   model1.hr, model2.hr, model3.hr, model4.hr, model.4.hr.fe)

#---- Rural Counties ----
#--- violent crime rate ----
model1.vcr <- lm(log.vcr.08 ~ mob_adjusted, data = full.DT[urban06 == 0])
model2.vcr <- lm(log.vcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + gini.09 + pct.singleHoH.09, data = full.DT[urban06 == 0])
model3.vcr <- lm(log.vcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + gini.09 + pct.singleHoH.09 + pct.samehouse.08.09 + pct.pop.diff.0708 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08, data = full.DT[urban06 == 0])
model4.vcr <- lm(log.vcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + gini.09 + pct.singleHoH.09 + pct.samehouse.08.09 + pct.pop.diff.0708 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08 + log_vcr_rollmean_1997, data = full.DT[urban06 == 0])
model.4.vcr.fe <- lm(log.vcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + gini.09 + pct.singleHoH.09 + pct.samehouse.08.09 + pct.pop.diff.0708 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08 + log_vcr_rollmean_1997 + factor(stateFIPS), data = full.DT[urban06 == 0])

models.all.vcr <- list(model2.vcr, model3.vcr, model4.vcr, model.4.vcr.fe)

list.tot.se.vcr <- rlist::list.append(lapply(models.all.vcr, choose_robust_func))

list.tot.p.vcr <- lapply(models.all.vcr, p.adjust.hommel)

stargazer(models.all.vcr, type = "text", align = TRUE,
          # uncomment to write table to file
          # out = "model_tables/vcr_popscaled_rur_121823.html", omit = "state",
          keep = "mob_adjusted",
          se = c(list.tot.se.vcr),
          p = c(list.tot.p.vcr),
          dep.var.caption = "Log(Violent Crime Rate + 1)",
          dep.var.labels.include = FALSE,
          order = c(1,2,9),
          covariate.labels = c("Mobility",  "Poverty",  "GINI", "LEOs Per 100,000", "Unemployment Rate",  "Population Per 100,000",  "Percent Non-Hisp. Black", "Percent Male 15-24 Pop.", "Percent Single Parent Head of Household",
                               "Percent in Same House from 2008 to 2009",  "Percent Difference in Population from 2007 to 2008",  "Median Number of Years Lived in Current House",  "Percent Owner Households", "Percent Occupied Households",  "Percent on Housing Vouchers",  "3-Year Average of Past Violence (1996 - 1998)"),
          add.lines = list(c("State Fixed-Effects", rep("No", 3), "Yes")),
          title = "Estimation Results for OLS Models Predicting Violent Crime in 2008 for Rural Counties (per Uniform Crime Reports)",
          notes = "Reports robust SEs and uses Hommel adjustments for p-values")



#---- homicide rate ----
model1.hr <- lm(log.hr.08 ~ mob_adjusted, data = full.DT[urban06 == 0])
model2.hr <- lm(log.hr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + gini.09 + pct.singleHoH.09, data = full.DT[urban06 == 0])
model3.hr <- lm(log.hr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + gini.09 + pct.singleHoH.09 + pct.samehouse.08.09 + pct.pop.diff.0708 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08, data = full.DT[urban06 == 0])
model4.hr <- lm(log.hr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + gini.09 + pct.singleHoH.09 + pct.samehouse.08.09 + pct.pop.diff.0708 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08 + log_hr_rollmean_1997, data = full.DT[urban06 == 0])
model.4.hr.fe <- lm(log.hr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + gini.09 + pct.singleHoH.09 + pct.samehouse.08.09 + pct.pop.diff.0708 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08 + log_hr_rollmean_1997 + factor(stateFIPS), data = full.DT[urban06 == 0])

models.all.hr <- list(model2.hr, model3.hr, model4.hr, model.4.hr.fe)
list.tot.se.hr <- rlist::list.append(lapply(models.all.hr, choose_robust_func))
list.tot.p.hr <- lapply(models.all.hr, p.adjust.hommel)

stargazer(models.all.hr, type = "text", align = TRUE,
          # uncomment to write table to file
          # out = "model_tables/hr_popscaled_rur_121823.html", omit = "state",
          keep = "mob_adjusted",
          se = c(list.tot.se.hr),
          p = c(list.tot.p.hr),
          dep.var.caption = "Log(Homicide Rate + 1)",
          dep.var.labels.include = FALSE,
          order = c(1,2,9),
          covariate.labels = c("Mobility",  "Poverty",  "GINI", "LEOs Per 100,000", "Unemployment Rate",  "Population Per 100,000",  "Percent Non-Hisp. Black", "Percent Male 15-24 Pop.", "Percent Single Parent Head of Household",
                               "Percent in Same House from 2008 to 2009",  "Percent Difference in Population from 2007 to 2008",  "Median Number of Years Lived in Current House",  "Percent Owner Households", "Percent Occupied Households",  "Percent on Housing Vouchers",  "3-Year Average of Past Homicide (1996 - 1998)"),
          add.lines = list(c("State Fixed-Effects", rep("No", 3), "Yes")),
          title = "Estimation Results for OLS Models Predicting Homicide Rates in 2008 for Rural Counties (per Uniform Crime Reports)",
          notes = "Reports robust SEs and uses Hommel adjustments for p-values")

rm(model1.vcr, model2.vcr, model3.vcr, model4.vcr, model.4.vcr.fe,
   model1.hr, model2.hr, model3.hr, model4.hr, model.4.hr.fe)

#---- Black Counties ----
blk.DT <- full.DT[pct.nh.black.08 > 16] # counties with a black population percentage of at least 16%

#---- violent crime rate ----
model1.vcr <- lm(log.vcr.08 ~ mob_adjusted, data = blk.DT)
model2.vcr <- lm(log.vcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.male.15.24.08 + urban06 + gini.09, data = blk.DT)
model3.vcr <- lm(log.vcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.male.15.24.08 + urban06 + gini.09 + pct.samehouse.08.09 + pct.pop.diff.0708 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08, data = blk.DT)
model4.vcr <- lm(log.vcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.male.15.24.08 + urban06 + gini.09 + pct.samehouse.08.09 + pct.pop.diff.0708 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08 + log_vcr_rollmean_1997, data = blk.DT)
model.4.vcr.fe <- lm(log.vcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.male.15.24.08 + urban06 + gini.09 + pct.samehouse.08.09 + pct.pop.diff.0708 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08 + log_vcr_rollmean_1997 + factor(stateFIPS), data = blk.DT)

models.all.vcr <- list(model2.vcr, model3.vcr, model4.vcr, model.4.vcr.fe)

list.tot.se.vcr <- rlist::list.append(lapply(models.all.vcr, choose_robust_func))

list.tot.p.vcr <- lapply(models.all.vcr, p.adjust.hommel)

stargazer(models.all.vcr, type = "text", align = TRUE,
          # uncomment to write table to file
          # out = "model_tables/vcr_popscaled_blk_121923.html",
          keep = "mob_adjusted",
          se = c(list.tot.se.vcr),
          p = c(list.tot.p.vcr),
          dep.var.caption = "Log(Violent Crime Rate + 1)",
          dep.var.labels.include = FALSE,
          order = c(1,2,9),
          covariate.labels = c("Mobility",  "Poverty",  "GINI", "LEOs Per 100,000", "Unemployment Rate",  "Population Per 100,000",  "Percent Male 15-24 Pop.", "Urban",
                               "Percent in Same House from 2008 to 2009",  "Percent Difference in Population from 2007 to 2008",  "Median Number of Years Lived in Current House",  "Percent Owner Households", "Percent Occupied Households",  "Percent on Housing Vouchers",  "3-Year Average of Past Violence (1996 - 1998)"),
          add.lines = list(c("State Fixed-Effects", rep("No", 3), "Yes")),
          title = "Estimation Results for OLS Models Predicting Violent Crime in 2008 for Counties with a Significant Black Population (per Uniform Crime Reports)",
          notes = "Reports robust SEs and uses Hommel adjustments for p-values")



#---- homicide rate ----
model1.hr <- lm(log.hr.08 ~ mob_adjusted, data = blk.DT)
model2.hr <- lm(log.hr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.male.15.24.08 + urban06 + gini.09, data = blk.DT)
model3.hr <- lm(log.hr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.male.15.24.08 + urban06 + gini.09 + pct.samehouse.08.09 + pct.pop.diff.0708 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08, data = blk.DT)
model4.hr <- lm(log.hr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.male.15.24.08 + urban06 + gini.09 + pct.samehouse.08.09 + pct.pop.diff.0708 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08 + log_hr_rollmean_1997, data = blk.DT)
model.4.hr.fe <- lm(log.hr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.male.15.24.08 + urban06 + gini.09 + pct.samehouse.08.09 + pct.pop.diff.0708 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08 + log_hr_rollmean_1997 + factor(stateFIPS), data = blk.DT)

models.all.hr <- list(model2.hr, model3.hr, model4.hr, model.4.hr.fe)
list.tot.se.hr <- rlist::list.append(lapply(models.all.hr, choose_robust_func))
list.tot.p.hr <- lapply(models.all.hr, p.adjust.hommel)

stargazer(models.all.hr, type = "text", align = TRUE,
          # uncomment to write table to file
          # out = "model_tables/hr_popscaled_blk_121923.html",
          keep = "mob_adjusted",
          se = c(list.tot.se.hr),
          p = c(list.tot.p.hr),
          dep.var.caption = "Log(Homicide Rate + 1)",
          dep.var.labels.include = FALSE,
          order = c(1,2,9),
          covariate.labels = c("Mobility",  "Poverty",  "GINI", "LEOs Per 100,000", "Unemployment Rate",  "Population Per 100,000",  "Percent Male 15-24 Pop.", "Urban",
                               "Percent in Same House from 2008 to 2009",  "Percent Difference in Population from 2007 to 2008",  "Median Number of Years Lived in Current House",  "Percent Owner Households", "Percent Occupied Households",  "Percent on Housing Vouchers",  "3-Year Average of Past Homicide (1996 - 1998)"),
          add.lines = list(c("State Fixed-Effects", rep("No", 3), "Yes")),
          title = "Estimation Results for OLS Models Predicting Homicide Rates in 2008 for Counties with a Significant Black Population (per Uniform Crime Reports)",
          notes = "Reports robust SEs and uses Hommel adjustments for p-values")

rm(model1.vcr, model2.vcr, model3.vcr, model4.vcr, model.4.vcr.fe,
   model1.hr, model2.hr, model3.hr, model4.hr, model.4.hr.fe)


#---- PROPERTY CRIME ----
property.DT <- readxl::read_xlsx("PNAS_toShare.xlsx", sheet = "property_crime") %>% as.data.table()

# MODELS
model1.pcr <- lm(log.pcr.08 ~ mob_adjusted, data = property.DT)
model2.pcr <- lm(log.pcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + urban06 + gini.09 + pct.singleHoH.09, data = property.DT)
model3.pcr <- lm(log.pcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + urban06 + gini.09 + pct.singleHoH.09 + pct.pop.diff.0708 + pct.samehouse.08.09 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08, data = property.DT)
model4.pcr <- lm(log.pcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + urban06 + gini.09 + pct.singleHoH.09 + pct.pop.diff.0708 + pct.samehouse.08.09 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08 + log_pcr_rollmean_1997, data = property.DT)
model.4.pcr.fe <- lm(log.pcr.08 ~ mob_adjusted + poverty.rate.08 + leo_per_100k_08 + unemp.rate.08 + TOT_POP_08_scaled + pct.nh.black.08 + pct.male.15.24.08 + urban06 + gini.09 + pct.singleHoH.09 + pct.pop.diff.0708 + pct.samehouse.08.09 + median_yrs_in_current_house_09 + pct.own.09 + pct.occupied.09 + pct.pop.HUD.08 + log_pcr_rollmean_1997 + factor(stateFIPS), data = property.DT)

models.all.pcr <- list(model2.pcr, model3.pcr, model4.pcr, model.4.pcr.fe)

list.tot.se.pcr <- rlist::list.append(lapply(models.all.pcr, choose_robust_func))

list.tot.p.pcr <- lapply(models.all.pcr, p.adjust.hommel)

stargazer(models.all.pcr, type = "html", align = TRUE,
          # uncomment to write table to file
          # out = "model_tables/pcr_08_122123.html",
          omit = "state",
          se = c(list.tot.se.pcr),
          p = c(list.tot.p.pcr),
          # out = "model_tables/pcr04_noFE.html",
          dep.var.caption = "Log(Property Crime Rate + 1)",
          dep.var.labels.include = FALSE,
          order = c(1,2,9),
          covariate.labels = c("Mobility",  "Poverty",  "GINI", "LEOs Per 100,000", "Unemployment Rate",  "Population Per 100,000",  "Percent Non-Hisp. Black", "Percent Male 15-24 Pop.", "Urban", "Percent Single Parent Head of Household",
                               "Percent Difference in Population from the Previous Year", "Percent in Same County from the Previous Year",  "Median Number of Years Lived in Current House",  "Percent Owner Households", "Percent Occupied Households",  "Percent on Housing Vouchers",  "3-Year Average of Past Property Crime (1996 - 1998)"),
          add.lines = list(c("State Fixed-Effects", rep("No", 3), "Yes")),
          title = "Estimation Results for OLS Models Predicting Property Crime (per Uniform Crime Reports)",
          notes = "Reports robust SEs and uses Hommel adjustments for p-values")





