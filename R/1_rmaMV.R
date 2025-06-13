#' Author: TJ Sipin
#' Date: 03-20-2025

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

#Create data directory
dir.create("Data/1_DataProcessing/meta_analyses/models", recursive = T)

#Read in data
madata <- read_csv("Data/1_DataProcessing/cleaned_data_6_11_25_data_cleaned_collapse_risk_groups.csv")

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
        output_filename = paste0("Data/1_DataProcessing/meta_analyses/models/intercept_only.rds")
    } else if(is.null(mod2)){ ##If mod1 is NULL and mod2 is NULL, then model only mod1
        mod2 = "NULL"
        mod_formula = paste0("~ 1 + ", mod1)
        output_filename = paste0("Data/1_DataProcessing/meta_analyses/models/", mod1, ":", mod2, ".rds")
        
    } else if(mod1 == mod2){ ##In the rare case that mod1 == mod2 then return NULL (nonsensical)
        return(NULL)
    } else { ##This is the case of two interacting moderators
        #Create interaction variable combining mod1 and mod2
        interterm = paste0(mod1, "...", mod2)
        madata[, interterm] = paste0(madata[, mod1][[1]], ":", madata[, mod2][[1]])
        
        #Create mod formula for interaction terms
        mod_formula = paste0("~ 1 + ", interterm)
        output_filename = paste0("Data/1_DataProcessing/meta_analyses/models/", mod1, ":", mod2, ".rds")
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


funnelPlotFunc = function(fi){
    #Read in model
    model = readRDS(fi)
    
    #Get moderator(s)
    mod = as.character(model$formula.mods)[2] %>% 
        str_split_i(pattern = " \\+ ", i=2) %>% 
        str_replace("\\...", ":")
    #Set moderator to "intercept" if no moderator detected
    if(is.na(mod)){
        mod="intercept"
    }
    
    #Set plot path
    png(paste0("Data/1_DataProcessing/meta_analyses/png/funnel/se/", mod, "_se.png"))
    funnel_se = funnel(
        model, main="Standard Error", refline=0, level=c(90, 95, 99), 
        shade=c("white", "gray", "darkgray")
    )
    dev.off()
    
    png(paste0("Data/1_DataProcessing/meta_analyses/png/funnel/vi/", mod, "_vi.png"))
    funnel_vi = funnel(
        model, yaxis="vi", main="Sampling Variance", refline=0, level=c(90, 95, 99), 
        shade=c("white", "gray", "darkgray")
    )
    dev.off()
    
    png(paste0("Data/1_DataProcessing/meta_analyses/png/funnel/seinv/", mod, "_seinv.png"))
    funnel_seinv = funnel(
        model, yaxis="seinv", main="Inverse Standard Error", refline=0, level=c(90, 95, 99), 
        shade=c("white", "gray", "darkgray")
    )
    dev.off()
    
    png(paste0("Data/1_DataProcessing/meta_analyses/png/funnel/vinv/", mod, "_vinv.png"))
    funnel_vinv = funnel(
        model, yaxis="vinv", main="Inverse Sampling Variance", refline=0, level=c(90, 95, 99), 
        shade=c("white", "gray", "darkgray")
    )
    dev.off()
}

orchardPlotFunc = function(fi){
    #Read in model
    model = readRDS(fi)
    
    #Get moderator(s)
    mod = as.character(model$formula.mods)[2] %>% 
        str_split_i(pattern = " \\+ ", i=2)
    mod_name = mod %>% 
        str_replace("\\...", ":")
    
    
    #Color palette
    coul <- brewer.pal(n=9, name="Pastel1") 
    coul <- colorRampPalette(coul)(80)
    #Intercept-only plot
    if(is.na(mod)){
        output_filename = paste0("Data/1_DataProcessing/meta_analyses/png/orchard/intercept_orchard.png")
        # if(file.exists(output_filename)) return(NULL)
        orchaRd::orchard_plot(
            model, mod='1', group='groupID',
            xlab = "Log Odds Ratio", angle = 45, g = FALSE, legend.pos = "top.left",
            k = T, transfm = 'none', cb = F
        ) + 
            scale_color_manual(values = coul) +
            theme(legend.direction = "vertical")
        
        ggsave(output_filename, scale=3, width=1080, height=1080, units='px')
        return(NULL)
    }
    
    #Declare output filename
    output_filename = paste0("Data/1_DataProcessing/meta_analyses/png/orchard/", mod, "_orchard.png")
    # if(file.exists(output_filename)) return(NULL)
    
    #If there is an error, return NULL and a message
    tryCatch({
        HetModel <- orchaRd::mod_results(
            model, group = "groupID", 
            mod = mod, 
            weights = "prop"
        )
        
        
        oplot = orchaRd::orchard_plot(
            HetModel,
            xlab = "Log Odds Ratio", angle = 0, g = FALSE, legend.pos = "top.left",
            condition.lab = small_mod, k = T, transfm = 'none', cb = F
        ) + 
            scale_color_manual(values = coul) +
            theme(legend.direction = "vertical")
        oplot
        
        ggsave(output_filename, scale=3, width=1080, height=1080, units='px')
    }, error = function(e){
        message(
            cat(paste0(mod, ": \n", e))
        )
        return(message(
            cat(paste0(mod, ": \n", e))
        ))
    })
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
dir.create("Data/1_DataProcessing/meta_analyses/png/funnel/se", recursive=T)
dir.create("Data/1_DataProcessing/meta_analyses/png/funnel/seinv", recursive=T)
dir.create("Data/1_DataProcessing/meta_analyses/png/funnel/vi", recursive=T)
dir.create("Data/1_DataProcessing/meta_analyses/png/funnel/vinv", recursive=T)
dir.create("Data/1_DataProcessing/meta_analyses/png/orchard", recursive=T)

write_csv(model_fit_stats_all, "Data/1_DataProcessing/meta_analyses/stats/model_fit_stats.csv")

#Map funnel and orchard plot functions on models
map(model_files, funnelPlotFunc)
map(model_files, orchardPlotFunc)