#' Author: TJ Sipin
#' Date: 06-11-2025

#Set working directory
setwd("/home/tjsipin/network-storage/meta_analysis/")
#Read in packages
library(orchaRd)
library(patchwork)
library(tidyverse)
require(metafor)
require(meta)
library(stringr)
require(ggplot2)
require(metaviz)
library(RColorBrewer)

#Reorder model files
model_files = list.files(
    "Data/1_DataProcessing/meta_analyses/models",
    pattern = ".rds",
    full.names = T
)


## Read in files and gather model fit statistics
modelFitStats = function(fi){
    #Read in model
    model = readRDS(fi)
    
    aic_i = AIC.rma(model)
    r2_i = r2_ml(model)
    QE_i = model$QE
    QM_i = model$QM
    tau2 <- sum(model$sigma2)   # total heterogeneity variance
    mean_se2 <- mean(model$vi)  # mean sampling variance
    
    I2_total <- (tau2 / (tau2 + mean_se2)) * 100
    
    #do the same for the random effects
    tau2_components <- model$sigma2
    names(tau2_components) <- c("RecordNo", "groupID", "OR_method")
    
    total_variance <- sum(tau2_components) + mean_se2
    
    I2_each <- (tau2_components / total_variance) * 100
    
    #Get moderator(s) from formula
    mod = as.character(model$formula.mods)[2] %>% 
        str_split_i(pattern = " \\+ ", i=2)
    
    if(is.na(mod)){
        mod="NULL"
    }
    
    tibble(
        mod = mod,
        aic = aic_i,
        I2_total = I2_total,
        I2_RecordNo = I2_each[[1]],
        I2_groupID = I2_each[[2]],
        I2_OR_method = I2_each[[3]],
        R2_marginal = r2_i[[1]],
        R2_conditional = r2_i[[2]],
        QE_i = QE_i,
        QM_i = QM_i,
        mean_se2 = mean_se2, 
        tau2_total = tau2,
        tau2_RecordNo = tau2_components[[1]],
        tau2_groupID = tau2_components[[2]],
        tau2_OR_method = tau2_components[[3]],
        total_variance = total_variance
    )
}

#Gather model fit statistics ordered by AIC
model_fit_stats = map(model_files, modelFitStats) %>% 
    bind_rows() %>% 
    mutate(model_name = NA)

#Get intercept only model
model_fit_stats_intercept_only = model_fit_stats %>% 
    filter(mod=="NULL") %>% 
    mutate(model_name = "Intercept only")

#Get all other models
model_fit_stats_nonintercept = model_fit_stats %>% 
    filter(mod!="NULL") %>% 
    #sort by AIC
    arrange(aic) %>% 
    mutate(model_name = str_replace(mod, "\\...", replacement = " : "))

#Overall effect as first row, all double else sorted by AIC
model_fit_stats_all = rbind(model_fit_stats_intercept_only, model_fit_stats_nonintercept) %>% 
    dplyr::select(model_name, everything())

dir.create("Data/1_DataProcessing/meta_analyses/stats", recursive=T)
write_csv(model_fit_stats_all, "Data/1_DataProcessing/meta_analyses/stats/model_fit_stats.csv")