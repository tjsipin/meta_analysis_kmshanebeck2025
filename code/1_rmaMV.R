#' Author: TJ Sipin
#' Date: 06-11-2025

################# Install the following packages if necessary #################
# install.packages("devtools")
# install.packages("tidyverse")
# install.packages("metafor")
# install.packages("patchwork")
# install.packages("R.rsp")
# install.packages("remotes")
# remotes::install_github("daniel1noble/metaAidR")
# devtools::install_github("itchyshin/orchard_plot", subdir = "orchaRd", force = TRUE, build_vignettes = TRUE)
# devtools::install_github("daniel1noble/orchaRd", ref = "main", force = TRUE)
# devtools::install_github("daniel1noble/orchaRd", force = TRUE)
# install.packages('meta')
# install.packages('metaviz')

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

#Create data directory
dir.create("results/models", recursive = T)

#Read in data
madata <- read_csv("data/cleaned_data_6_11_25_data_cleaned_collapse_risk_groups.csv")

#Declare factors
##Factors to feed into single factor models
single_factors = c(
    "NULL", "Year", "marg_group", "group_proxy", "disease_type", "pathogen",
    "transmission", "risk_group", "risk_proxy", "country", 
    "scale", "data_source", "broad_region"
)

#Factors to feed into double factor models
double_factors_grid = combn(
    c("marg_group", 'group_proxy', "disease_type", 
      "transmission", "risk_group", "risk_proxy", 
      "scale", "data_source", "broad_region"),
    2, FUN = NULL, simplify = T
) %>% 
    matrix(ncol = 2, byrow = T) %>% 
    as.data.frame() 
names(double_factors_grid) = c("mod1", "mod2")

#Create function to make meta-analysis regression models
rmaMvFunc = function(mod1, mod2=NULL){
    #Print moderators
    print(paste(mod1, mod2))
    
    #Make model based on input moderators
    if(mod1=="NULL"){ ##If mod1 is NULL, then model only uses intercept
        mod2 = "NULL"
        mod_formula = "~ 1"
        output_filename = paste0("results/models/intercept_only.rds")
    } else if(is.null(mod2)){ ##If mod1 is NULL and mod2 is NULL, then model only mod1
        mod2 = "NULL"
        mod_formula = paste0("~ 1 + ", mod1)
        output_filename = paste0("results/models/", mod1, ":", mod2, ".rds")
        
    } else if(mod1 == mod2){ ##In the rare case that mod1 == mod2 then return NULL (nonsensical)
        return(NULL)
    } else { ##This is the case of two interacting moderators
        #Create interaction variable combining mod1 and mod2
        interterm = paste0(mod1, "...", mod2)
        madata[, interterm] = paste0(madata[, mod1][[1]], ":", madata[, mod2][[1]])
        
        #Create mod formula for interaction terms
        mod_formula = paste0("~ 1 + ", interterm)
        output_filename = paste0("results/models/", mod1, ":", mod2, ".rds")
    }
    
    #Cancel iteration if output file exists
    # if(file.exists(output_filename)){
    #     return(NULL)
    # }
    
    #Create model
    meta.analysis = rma.mv(
        yi = LogOR, #yi
        V = LogOR_v, #V
        mods = formula(mod_formula),
        #Random effects
        random = list(~ 1 | RecordNo, ~ 1 | groupID, ~ 1 | OR_method),
        #Model fitted with restricted maximum likelihood estimation
        method = "REML",
        data = madata
    )
    
    #Save model
    saveRDS(meta.analysis, output_filename)
    
    return(meta.analysis)
    
}

#Run on all single factors
single_factor_run = map(
    .x = single_factors,
    ~ rmaMvFunc(mod1 = .x)
)

#Run on all possible combinations of factors
two_factor_run = map2(
    .x = double_factors_grid$mod1,
    .y = double_factors_grid$mod2,
    ~ rmaMvFunc(.x, .y)
)