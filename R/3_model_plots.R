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

dir.create("Data/1_DataProcessing/meta_analyses/png/funnel/se", recursive=T)
dir.create("Data/1_DataProcessing/meta_analyses/png/funnel/seinv", recursive=T)
dir.create("Data/1_DataProcessing/meta_analyses/png/funnel/vi", recursive=T)
dir.create("Data/1_DataProcessing/meta_analyses/png/funnel/vinv", recursive=T)
dir.create("Data/1_DataProcessing/meta_analyses/png/orchard", recursive=T)

#Map funnel and orchard plot functions on models
map(model_files, funnelPlotFunc)
map(model_files, orchardPlotFunc)