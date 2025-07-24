#' Author: TJ Sipin
#' Date: 06-11-2025

#Read in packages
library(orchaRd)
require(metaviz)
require(metafor)
require(meta)
library(RColorBrewer)
library(dplyr)


dir.create("results/png/orchard/main/", recursive=T)

orchardPlotFunc = function(fi, coul){
    model = readRDS(fi)
    
    mod = as.character(model$formula.mods)[2] %>% 
        str_split_i(pattern = " \\+ ", i=2)
    mod_name = mod %>% 
        str_replace("\\...", ":")
    
    output_filename = paste0("results/png/orchard/main/", mod, "_orchard.png")
    
    tryCatch({
        
        mod_res <- orchaRd::mod_results(
            model, group = "groupID", 
            mod = mod, 
            weights = "prop"
        )
        
        #For each moderator combination, we want specific tree orders
        if(mod_name == "risk_group:risk_proxy"){
            oplot = orchaRd::orchard_plot(
                mod_res,
                xlab = "Log Odds Ratio", angle = 85, g = FALSE, legend.pos = "top.left",
                condition.lab = small_mod, k = T, transfm = 'none', cb = F, flip = T,
                k.pos = 'right',
                tree.order = c(
                    "Disease_risk:diagnosis",
                    "Disease_risk:health_outcome",
                    "Disease_risk:hospitalizations",
                    
                    "Infection_risk:seroprevalance",
                    "Infection_risk:reported_cases",
                    
                    "Exposure_risk:behavior",
                    "Exposure_risk:knowledge",
                    "Exposure_risk:vector_proximity"
                )
            ) 
        } else if(mod_name == "marg_group:risk_group"){
            oplot = orchaRd::orchard_plot(
                mod_res,
                xlab = "Log Odds Ratio", angle = 85, g = FALSE, legend.pos = "top.left",
                condition.lab = small_mod, k = T, transfm = 'none', cb = F, flip = T,
                k.pos = 'right',
                tree.order = c(
                    "Race:disease_risk",
                    "Race:infection_risk",
                    "Race:exposure_risk",
                    
                    "Socioeconomic:disease_risk",
                    "Socioeconomic:infection_risk",
                    "Socioeconomic:exposure_risk"
                )
            ) 
        } else{
            oplot = orchaRd::orchard_plot(
                mod_res,
                xlab = "Log Odds Ratio", angle = 85, g = FALSE, legend.pos = "top.left",
                condition.lab = small_mod, k = T, transfm = 'none', cb = F, flip = T,
                k.pos = 'right'
            ) 
        }
        
        oplot_out = oplot + 
            coord_cartesian(ylim=c(-2, 4), clip = "on", expand = F) +
            coord_flip(ylim=c(-2, 4), clip = "on", expand = F) +
            scale_fill_manual(values = coul) +
            ggplot2::theme(
                legend.direction = "vertical",
                text = element_text(size = 12),
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 6),
                axis.text.r = element_text(size = 12)
            )
        
        oplot_out
        
        #Save plots individually
        ggsave(output_filename, scale=2, width=1080, height=1080, units='px')
        return(oplot_out)
    }, error = function(e){
        message(
            cat(paste0(mod, ": \n", e))
        )
        return(NULL)
    })
}

model_files = c(
    "results/models/marg_group:group_proxy.rds",
    "results/models/risk_group:risk_proxy.rds",
    "results/models/marg_group:risk_group.rds",
    # extra standalone plots
    "results/models/transmission:broad_region.rds"
)

#Create color palettes for each orchard plot
coul1 <- brewer.pal(n=8, name="Reds") 
coul1 <- colorRampPalette(coul1)(readRDS(model_files[1])$coef %>% length())

coul2 <- brewer.pal(n=8, name="Blues") 
coul2 <- colorRampPalette(coul2)(readRDS(model_files[2])$coef %>% length())

coul3 <- brewer.pal(n=8, name="Oranges") 
coul3 <- colorRampPalette(coul3)(readRDS(model_files[3])$coef %>% length())

o1 = orchardPlotFunc(model_files[1], coul=coul1)
o2 = orchardPlotFunc(model_files[2], coul=coul2)
o3 = orchardPlotFunc(model_files[3], coul=coul3)

#Arrange into a 1x3 faceted plot
ggpubr::ggarrange(
    o1 + 
        labs(y = "marg_group:group_proxy", caption = " ") +
        theme(
            axis.text.y = element_text(margin = margin(20, 20, 20, 20)),
            axis.title.x = element_text(size = 8)
        ), 
    o2 + 
        labs(y = "risk_group:risk_proxy") +
        theme(
            axis.text.y = element_text(margin = margin(20, 20, 20, 20)),
            axis.title.x = element_text(size = 8)
        ), 
    o3 + 
        labs(y = "marg_group:risk_group", caption = " ") +
        theme(
            axis.text.y = element_text(margin = margin(20, 20, 20, 20)),
            axis.title.x = element_text(size = 8)
        ), 
    nrow = 1
)
ggsave(
    "results/png/orchard/main/main_fig_facet.png", scale=2, width=2040, height=1080, units='px'
)


#Same as the first orchard plot function but making sure the K is visible
orchardPlotFuncWithK = function(fi, coul){
    #Read in model
    model = readRDS(fi)
    
    #Pull in moderator names
    mod = as.character(model$formula.mods)[2] %>% 
        str_split_i(pattern = " \\+ ", i=2)
    mod_name = mod %>% 
        str_replace("\\...", ":")
    
    #Declare output file
    output_filename = paste0("results/png/orchard/main/", mod, "_orchard_withK.png")
    
    tryCatch({
        #Create 
        mod_res <- orchaRd::mod_results(
            model, group = "groupID", 
            mod = mod, 
            weights = "prop"
        )
        #For each moderator combination, we want specific tree orders
        if(mod_name == "risk_group:risk_proxy"){
            oplot = orchaRd::orchard_plot(
                mod_res,
                xlab = "Log Odds Ratio", angle = 85, g = FALSE, legend.pos = "top.left",
                condition.lab = small_mod, k = T, transfm = 'none', cb = F, flip = T,
                k.pos = 'right',
                tree.order = c(
                    "Disease_risk:diagnosis",
                    "Disease_risk:health_outcome",
                    "Disease_risk:hospitalizations",
                    
                    "Infection_risk:seroprevalance",
                    "Infection_risk:reported_cases",
                    
                    "Exposure_risk:behavior",
                    "Exposure_risk:knowledge",
                    "Exposure_risk:vector_proximity"
                )
            ) 
        } else if(mod_name == "marg_group:risk_group"){
            oplot = orchaRd::orchard_plot(
                mod_res,
                xlab = "Log Odds Ratio", angle = 85, g = FALSE, legend.pos = "top.left",
                condition.lab = small_mod, k = T, transfm = 'none', cb = F, flip = T,
                k.pos = 'right',
                tree.order = c(
                    "Race:disease_risk",
                    "Race:infection_risk",
                    "Race:exposure_risk",
                    
                    "Socioeconomic:disease_risk",
                    "Socioeconomic:infection_risk",
                    "Socioeconomic:exposure_risk"
                )
            ) 
        } else if(mod_name == "marg_group:group_proxy"){
            oplot = orchaRd::orchard_plot(
                mod_res,
                xlab = "Log Odds Ratio", angle = 85, g = FALSE, legend.pos = "top.left",
                condition.lab = small_mod, k = T, transfm = 'none', cb = F, flip = T,
                k.pos = 'right'
            ) 
        } else if(mod_name == "transmission:broad_region"){
            ##This is a different one, so no color palette here
            oplot = orchaRd::orchard_plot(
                mod_res,
                xlab = "Log Odds Ratio", angle = 85, g = FALSE, legend.pos = "top.left",
                condition.lab = small_mod, k = T, transfm = 'none', cb = F, flip = T,
                k.pos = 'right',
                tree.order = c(
                    "Lice_flea:Non_US",
                    "Lice_flea:West",
                    "Lice_flea:South",
                    "Lice_flea:National",
                    
                    "Mosquito:Non_US",
                    "Mosquito:Territories",
                    "Mosquito:West",
                    "Mosquito:South",
                    "Mosquito:Midwest",
                    "Mosquito:Northeast",
                    "Mosquito:National",
                    
                    "Tick:Non_US",
                    "Tick:South",
                    "Tick:Midwest",
                    "Tick:Northeast",
                    "Tick:National",
                    
                    "Zoonotic:Non_US",
                    "Zoonotic:West",
                    "Zoonotic:South",
                    "Zoonotic:Midwest",
                    "Zoonotic:Northeast",
                    "Zoonotic:National"
                )
            ) 
            
            oplot_out = oplot + 
                coord_cartesian(ylim=c(-2, 5), clip = "on", expand = T) +
                coord_flip(ylim=c(-2, 5), clip = "on", expand = T) +
                ggplot2::theme(
                    legend.direction = "vertical",
                    text = element_text(size = 12),
                    axis.text.x = element_text(size = 12),
                    axis.text.y = element_text(size = 3),
                    axis.text.r = element_text(size = 12)
                )
            oplot_out
            
            ggsave(output_filename, scale=2, width=1080, height=1080, units='px')
            ##Jump out of function when done, we don't want further adjustments to our plot
            return(oplot_out)
        }
        
        #
        oplot_out = oplot + 
            coord_cartesian(ylim=c(-2, 5), clip = "on", expand = T) +
            coord_flip(ylim=c(-2, 5), clip = "on", expand = T) +
            scale_fill_manual(values = coul) +
            ggplot2::theme(
                legend.direction = "vertical",
                text = element_text(size = 12),
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 6),
                axis.text.r = element_text(size = 12)
            )
        oplot_out
        
        #Save plots individually
        ggsave(output_filename, scale=2, width=1080, height=1080, units='px')
        return(oplot_out)
    }, error = function(e){
        message(
            cat(paste0(mod, ": \n", e))
        )
        return(NULL)
    })
}

o1_withK = orchardPlotFuncWithK(model_files[1], coul=coul1)
o2_withK = orchardPlotFuncWithK(model_files[2], coul=coul2)
o3_withK = orchardPlotFuncWithK(model_files[3], coul=coul3)
o4_withK = orchardPlotFuncWithK(model_files[4])

#Arrange plots for MS
ggpubr::ggarrange(
    o1_withK + 
        labs(y = "marg_group:group_proxy", caption = " ") +
        theme(
            axis.text.y = element_text(margin = margin(20, 20, 20, 20)),
            axis.title.x = element_text(size = 8)
        ), 
    o2_withK + 
        labs(y = "risk_group:risk_proxy") +
        theme(
            axis.text.y = element_text(margin = margin(20, 20, 20, 20)),
            axis.title.x = element_text(size = 8)
        ), 
    o3_withK + 
        labs(y = "marg_group:risk_group", caption = " ") +
        theme(
            axis.text.y = element_text(margin = margin(20, 20, 20, 20)),
            axis.title.x = element_text(size = 8)
        ), 
    nrow = 1
)
ggsave(
    "results/png/orchard/main/main_fig_facet_withK.png", scale=2, width=2040, height=1080, units='px'
)