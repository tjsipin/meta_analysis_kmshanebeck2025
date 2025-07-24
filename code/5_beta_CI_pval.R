#' Author: TJ Sipin
#' Date: 06-11-2025

#Read in packages
library(ggplot2)
library(metafor)
library(tidymodels)
tidymodels_prefer()

#Create directory for log odds ratio plots
dir.create("results/png/log_odds_ratio", recursive=T)

#Plot order
# A) overall effect (intercept only) + marg group, B) Risk Group, C) Transmission, D) Country, E) Data source, F) scale

#Read in model files of interest
ma.intercept = readRDS("results/models/intercept_only.rds")
ma.marg_group = readRDS("results/models/marg_group:NULL.rds")
ma.risk_group = readRDS("results/models/risk_group:NULL.rds")
ma.transmission = readRDS("results/models/transmission:NULL.rds")
ma.country = readRDS("results/models/country:NULL.rds")
ma.data_source = readRDS("results/models/data_source:NULL.rds")
ma.scale = readRDS("results/models/scale:NULL.rds")
ma.group_proxy = readRDS("results/models/group_proxy:NULL.rds")

#Calculate p value from estimate, upper CI, and lower CI
pvalFunc = function(est, uci, lci){
    se = (uci - lci)/(2*qnorm(.975, lower.tail=T))
    z = est/se
    p = 2*(1-pnorm(abs(z), lower.tail = T))
    return(p)
}

#Function for obtaining point estimates, confidence intervals, and p values
estConfIntPval <- function(model, mod='1'){
    #Create data frame with point estimates, conf ints, and p vals
    mod_res = orchaRd::mod_results(model, mod = mod, group='groupID')$mod_table %>% 
        #Get p value
        mutate(p.value = pvalFunc(est=estimate, uci=upperCL, lci=lowerCL) %>% round(4))
    
    #Get k for each country, maintaining the capitalization
    if(mod=="country"){
        #Get k using table()
        mod_data_k = model$data[,mod][[1]] %>% 
            table() %>% 
            data.frame() %>% 
            rename(estimate = ".", k = Freq) 
        
        output = mod_res %>% 
            select(
                estimate=name, mean=estimate, lower=lowerCL, upper=upperCL, p.value
            ) %>% 
            full_join(mod_data_k)
    #Get k for each moderator, capitalizing the first letter of each name
    } else if(mod!="1"){
        #Get k using table()
        mod_data_k = model$data[,mod][[1]] %>% 
            table() %>% 
            data.frame() %>% 
            rename(estimate = ".", k = Freq) %>% 
            mutate(estimate = str_to_sentence(estimate))
        
        output = mod_res %>% 
            select(
                estimate=name, mean=estimate, lower=lowerCL, upper=upperCL, p.value
            ) %>% 
            full_join(mod_data_k) %>% 
            mutate(estimate = str_to_sentence(estimate))
    #For mod='1'
    } else{
        output = mod_res %>% 
            select(
                estimate=name, mean=estimate, lower=lowerCL, upper=upperCL, p.value
            ) %>% 
            mutate(k = nrow(model$data))
    }
    
    return(output)
}

#Apply function to models of interest
ma.intercept.ecip <- estConfIntPval(ma.intercept) %>% 
    #call analysis marginalized group so it can be joined with the marg group output
    mutate(
        analysis = "Marginalized group"
    ) %>%
    mutate(estimate = factor(
        estimate, 
        levels = c(
            "Race",
            "Socioeconomic",
            "Intrcpt"
        )
    )) %>% 
    mutate(estimate = fct_recode(
        estimate, 
        #Rename intercept for publication
        "Overall\nEffect" = "Intrcpt"
    ))
ma.marg_group.ecip <- estConfIntPval(ma.marg_group, "marg_group") %>% 
    mutate(
        analysis = "Marginalized group"
    ) %>% 
    mutate(estimate = factor(
        estimate, 
        levels = c(
            "Race",
            "Socioeconomic",
            "Intrcpt"
        )
    )) %>% 
    full_join(ma.intercept.ecip) %>%
    mutate(
        #rename and order estimates for publication
        estimate = fct_recode(
            estimate, 
            "Race/\nEthnicity" = "Race",
            "Socioeconomic\nStatus" = "Socioeconomic"
        ) %>% 
            ordered(c('Race/\nEthnicity', 'Socioeconomic\nStatus', 'Overall\nEffect'))
    )
ma.risk_group.ecip <- estConfIntPval(ma.risk_group, "risk_group") %>% 
    mutate(
        analysis = "Risk group"
    ) %>% 
    #Reorder according to author preference for publication
    mutate(
        estimate = fct_recode(
            estimate,
            "Exposure Risk" = "Exposure_risk",
            "Infection Risk" = "Infection_risk",
            "Disease Risk" = "Disease_risk"
        ) %>% 
            ordered(c("Exposure Risk", "Infection Risk", "Disease Risk"))
    )
ma.transmission.ecip <- estConfIntPval(ma.transmission, "transmission") %>% 
    mutate(
        analysis = "Transmission"
    ) %>% 
    #Rename and order for publication
    mutate(
        estimate = fct_recode(
            estimate,
            "Lice/flea" = "Lice_flea"
        ) %>% 
            ordered(c('Lice/flea', 'Mosquito', 'Tick', 'Zoonotic'))
    ) 
ma.country.ecip <- estConfIntPval(ma.country, "country") %>% 
    mutate(
        analysis = "Country"
    ) %>% 
    mutate(
        #Order for publication
        estimate = ordered(
            estimate,
            c("USA", "Europe", "Canada")
        )
    )
ma.scale.ecip <- estConfIntPval(ma.scale, "scale") %>% 
    mutate(
        analysis = "Scale"
    ) %>% 
    #Order for publication
    mutate(
        estimate = ordered(estimate, c("Local", "Regional", "National"))
    )
ma.data_source.ecip <- estConfIntPval(ma.data_source, "data_source") %>% 
    mutate(
        analysis = "Data source"
    ) %>% 
    #Rename and order for publication
    mutate(
        estimate = fct_recode(
            estimate,
            "Secondary-\nSource" = "Secondary",
            "Directly\nSampled" = "Direct"
        ) %>% ordered(c("Survey", "Secondary-\nSource", "Directly\nSampled"))
    )

#overall + marg group, data source, scale, risk group, country, and transmission
full.database <- rbind(
    ma.marg_group.ecip,
    ma.risk_group.ecip,
    ma.transmission.ecip,
    ma.country.ecip,
    ma.scale.ecip,
    ma.data_source.ecip
) %>% 
    distinct() %>% 
    #Change order of analysis
    mutate(
        analysis = ordered(
            analysis,
            c("Marginalized group", "Risk group", "Transmission",
              "Country", "Scale", "Data source")
        )
    ) %>% 
    mutate(
        estimate = ordered(
            estimate, 
            c(
                'Directly\nSampled', 'Secondary-\nSource', 'Survey',
                'National', 'Regional', 'Local',
                'Canada', 'Europe', 'USA',
                'Lice/flea', 'Mosquito', 'Tick', 'Zoonotic',
                'Disease Risk', 'Infection Risk', 'Exposure Risk', 
                'Race/\nEthnicity', 'Socioeconomic\nStatus', 'Overall\nEffect'
            )
        )
    )

facet_fig = full.database %>% 
    ggplot() +
    #Construct lines for confidence intervals
    geom_errorbarh(aes(xmin = lower,
                       xmax = upper, y = estimate), 
                   height = 0, show.legend = F, position = position_dodge(width=0.9),
                   color=rgb(0,0,0, 1)) +
    #Add point estimates
    geom_point(aes(x = mean, y = estimate, size=2),
               shape=20, fill = rgb(0,0,0, 0.9), position = position_dodge(width=0.9), show.legend = F) + 
    #Add k
    geom_text(aes(x = mean, y = estimate,label = k),vjust = 0,nudge_y = 0.12) +
    #Add dashed line at x = 0
    geom_vline(xintercept = 0, linetype = 2, colour = "black") + 
    ylab("")+
    xlab("Log Odds Ratio")+
    theme_classic() +
    theme(panel.spacing = unit(0.1, "lines"),
          text = element_text(size=12),
          panel.border= element_blank(),
          axis.line=element_line(), 
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(), 
          legend.text = element_text(size=16), 
          legend.title=element_text(size=16, face = "bold"),
          axis.title.x = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.y = element_text(size = 8, hjust = 0.35, margin = margin(r=5)),
          axis.text.y = element_text(angle = 0, color="black",hjust=0.95),
          axis.text.x = element_text(color="black"),
          plot.title = element_text(size = 16),
          plot.margin = unit(c(1.5,0,0,-0.5), "lines")) +
    #Clip plot 
    coord_cartesian(clip = "on", xlim = c(-0.4, 1.3))

#Save plot as both a 3x2 or 2x3 for publication
facet_fig_3x2 = facet_fig +
    facet_wrap(~analysis, nrow = 3, scales = "free_y") 
ggsave("results/png/log_odds_ratio/facet_fig_3x2.png", scale = 4, units = 'px', width=2040, height=1080)

facet_fig_2x3 = facet_fig +
    facet_wrap(~analysis, nrow = 2, scales = "free_y") 
ggsave("results/png/log_odds_ratio/facet_fig_2x3.png", scale = 4, units = 'px',width=1080, height=2040)

#Individual plots
figure1.1 = full.database %>% 
    filter(analysis=="Marginalized group") %>%
    ggplot() +
    geom_errorbarh(aes(xmin = lower,
                       xmax = upper, y = estimate), 
                   height = 0, show.legend = F, position = position_dodge(width=0.9),
                   color=rgb(0,0,0, 1)) +
    geom_point(aes(x = mean, y = estimate, size=2),
               shape=20, fill = rgb(0,0,0, 0.9), position = position_dodge(width=0.9), show.legend = F) + 
    geom_text(aes(x = mean, y = estimate,label = k),vjust = 0,nudge_y = 0.12) +
    geom_vline(xintercept = 0, linetype = 2, colour = "black") + 
    ylab("")+
    xlab("Log Odds Ratio")+
    theme_classic() +
    theme(panel.spacing = unit(0.1, "lines"),
          text = element_text(size=12),
          panel.border= element_blank(),
          axis.line=element_line(), 
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(), 
          legend.text = element_text(size=16), 
          legend.title=element_text(size=16, face = "bold"),
          axis.title.x = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.y = element_text(size = 8, hjust = 0.35, margin = margin(r=5)),
          axis.text.y = element_text(angle = 0, color="black",hjust=0.95),
          axis.text.x = element_text(color="black"),
          plot.title = element_text(size = 16),
          plot.margin = unit(c(1.5,0,0,-0.5), "lines")) +
    coord_cartesian(clip = "on", xlim = c(-0.4, 1.3)) +
    scale_x_continuous(breaks = seq(-0.5, 1.3, 0.2))

figure1.2 = full.database %>% 
    filter(analysis=="Risk group") %>%
    ggplot() +
    geom_errorbarh(aes(xmin = lower,
                       xmax = upper, y = estimate), 
                   height = 0, show.legend = F, position = position_dodge(width=0.9),
                   color=rgb(0,0,0, 1))+
    geom_point(aes(x = mean, y = estimate, size=2),
               shape=20, fill = rgb(0,0,0, 0.9), position = position_dodge(width=0.9), show.legend = F) + 
    geom_text(aes(x = mean, y = estimate,label = k),vjust = 0,nudge_y = 0.12) +
    geom_vline(xintercept = 0, linetype = 2, colour = "black") + 
    ylab("")+
    xlab("Log Odds Ratio")+
    theme_classic() +
    theme(panel.spacing = unit(0.1, "lines"),
          text = element_text(size=12),
          panel.border= element_blank(),
          axis.line=element_line(), 
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(), 
          legend.text = element_text(size=16), 
          legend.title=element_text(size=16, face = "bold"),
          axis.title.x = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.y = element_text(size = 8, hjust = 0.35, margin = margin(r=5)),
          axis.text.y = element_text(angle = 0, color="black",hjust=0.95),
          axis.text.x = element_text(color="black"),
          plot.title = element_text(size = 16),
          plot.margin = unit(c(1.5,0,0,-0.5), "lines")) +
    coord_cartesian(clip = "on", xlim = c(-0.5, 1.3)) +
    scale_x_continuous(breaks = seq(-0.5, 1.3, 0.2))

figure1.3 = full.database %>% 
    filter(analysis=="Transmission") %>%
    ggplot() +
    geom_errorbarh(aes(xmin = lower,
                       xmax = upper, y = estimate), 
                   height = 0, show.legend = F, position = position_dodge(width=0.9),
                   color=rgb(0,0,0, 1))+
    geom_point(aes(x = mean, y = estimate, size=2),
               shape=20, fill = rgb(0,0,0, 0.9), position = position_dodge(width=0.9), show.legend = F) + 
    geom_text(aes(x = mean, y = estimate,label = k),vjust = 0,nudge_y = 0.12) +
    geom_vline(xintercept = 0, linetype = 2, colour = "black") + 
    ylab("")+
    xlab("Log Odds Ratio")+
    theme_classic() +
    theme(panel.spacing = unit(0.1, "lines"),
          text = element_text(size=12),
          panel.border= element_blank(),
          axis.line=element_line(), 
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(), 
          legend.text = element_text(size=16), 
          legend.title=element_text(size=16, face = "bold"),
          axis.title.x = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.y = element_text(size = 8, hjust = 0.35, margin = margin(r=5)),
          axis.text.y = element_text(angle = 0, color="black",hjust=0.95),
          axis.text.x = element_text(color="black"),
          plot.title = element_text(size = 16),
          plot.margin = unit(c(1.5,0,0,-0.5), "lines")) +
    coord_cartesian(clip = "on", xlim = c(-0.5, 1.3)) +
    scale_x_continuous(breaks = seq(-0.5, 1.3, 0.2))

figure2.1 = full.database %>% 
    filter(analysis=="Country") %>%
    ggplot() +
    geom_errorbarh(aes(xmin = lower,
                       xmax = upper, y = estimate), 
                   height = 0, show.legend = F, position = position_dodge(width=0.9),
                   color=rgb(0,0,0, 1))+
    geom_point(aes(x = mean, y = estimate, size=2),
               shape=20, fill = rgb(0,0,0, 0.9), position = position_dodge(width=0.9), show.legend = F) + 
    geom_text(aes(x = mean, y = estimate,label = k),vjust = 0,nudge_y = 0.12) +
    geom_vline(xintercept = 0, linetype = 2, colour = "black") + 
    ylab("")+
    xlab("Log Odds Ratio")+
    theme_classic() +
    theme(panel.spacing = unit(0.1, "lines"),
          text = element_text(size=12),
          panel.border= element_blank(),
          axis.line=element_line(), 
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(), 
          legend.text = element_text(size=16), 
          legend.title=element_text(size=16, face = "bold"),
          axis.title.x = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.y = element_text(size = 8, hjust = 0.35, margin = margin(r=5)),
          axis.text.y = element_text(angle = 0, color="black",hjust=0.95),
          axis.text.x = element_text(color="black"),
          plot.title = element_text(size = 16),
          plot.margin = unit(c(1.5,0,0,-0.5), "lines")) +
    coord_cartesian(clip = "on", xlim = c(-0.5, 1.3)) +
    scale_x_continuous(breaks = seq(-0.5, 1.3, 0.2))


figure2.2 = full.database %>% 
    filter(analysis=="Data source") %>%
    ggplot() +
    geom_errorbarh(aes(xmin = lower,
                       xmax = upper, y = estimate), 
                   height = 0, show.legend = F, position = position_dodge(width=0.9),
                   color=rgb(0,0,0, 1))+
    geom_point(aes(x = mean, y = estimate, size=2),
               shape=20, fill = rgb(0,0,0, 0.9), position = position_dodge(width=0.9), show.legend = F) + 
    geom_text(aes(x = mean, y = estimate,label = k),vjust = 0,nudge_y = 0.12) +
    geom_vline(xintercept = 0, linetype = 2, colour = "black") + 
    ylab("")+
    xlab("Log Odds Ratio")+
    theme_classic() +
    theme(panel.spacing = unit(0.1, "lines"),
          text = element_text(size=12),
          panel.border= element_blank(),
          axis.line=element_line(), 
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(), 
          legend.text = element_text(size=16), 
          legend.title=element_text(size=16, face = "bold"),
          axis.title.x = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.y = element_text(size = 8, hjust = 0.35, margin = margin(r=5)),
          axis.text.y = element_text(angle = 0, color="black",hjust=0.95),
          axis.text.x = element_text(color="black"),
          plot.title = element_text(size = 16),
          plot.margin = unit(c(1.5,0,0,-0.5), "lines")) +
    coord_cartesian(clip = "on", xlim = c(-0.5, 1.3)) +
    scale_x_continuous(breaks = seq(-0.5, 1.3, 0.2))

figure2.3 = full.database %>% 
    filter(analysis=="Scale") %>% 
    ggplot() +
    geom_errorbarh(aes(xmin = lower,
                       xmax = upper, y = estimate), 
                   height = 0, show.legend = F, position = position_dodge(width=0.9),
                   color=rgb(0,0,0, 1))+
    geom_point(aes(x = mean, y = estimate, size=2),
               shape=20, fill = rgb(0,0,0, 0.9), position = position_dodge(width=0.9), show.legend = F) + 
    geom_text(aes(x = mean, y = estimate,label = k),vjust = 0,nudge_y = 0.12) +
    geom_vline(xintercept = 0, linetype = 2, colour = "black") + 
    ylab("")+
    xlab("Log Odds Ratio")+
    theme_classic() +
    theme(panel.spacing = unit(0.1, "lines"),
          text = element_text(size=12),
          panel.border= element_blank(),
          axis.line=element_line(), 
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(), 
          legend.text = element_text(size=16), 
          legend.title=element_text(size=16, face = "bold"),
          axis.title.x = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.y = element_text(size = 8, hjust = 0.35, margin = margin(r=5)),
          axis.text.y = element_text(angle = 0, color="black",hjust=0.95),
          axis.text.x = element_text(color="black"),
          plot.title = element_text(size = 16),
          plot.margin = unit(c(1.5,0,0,-0.5), "lines")) +
    coord_cartesian(clip = "on", xlim = c(-0.5, 1.3)) +
    scale_x_continuous(breaks = seq(-0.5, 1.3, 0.2))

#Display and save each figure
figure1.1
ggsave("results/png/log_odds_ratio/intercept_marg_group.png", scale = 2.22)
figure1.2
ggsave("results/png/log_odds_ratio/risk_group.png", scale = 2.22)
figure1.3
ggsave("results/png/log_odds_ratio/transmission.png", scale = 2.22)
figure2.1
ggsave("results/png/log_odds_ratio/country.png", scale = 2.22)
figure2.2
ggsave("results/png/log_odds_ratio/data_source.png", scale = 2.22)
figure2.3
ggsave("results/png/log_odds_ratio/scale.png", scale = 2.22)


# Table for each of the main effects --------------------------------------

#KMS needed this for manuscript
ma.group_proxy.transmission = readRDS("results/models/group_proxy:transmission.rds")
ma.group_proxy.transmission.ecip <- estConfIntPval(ma.group_proxy.transmission, "group_proxy...transmission") %>% 
    mutate(analysis = "Group proxy : Transmission")
ma.group_proxy.ecip = estConfIntPval(ma.group_proxy, "group_proxy") %>% 
    mutate(analysis = "Group proxy") 

table_061225 = rbind(
    ma.intercept.ecip %>% 
        mutate(analysis = "Intercept only"),
    ma.marg_group.ecip,
    ma.group_proxy.ecip,
    ma.risk_group.ecip,
    ma.transmission.ecip,
    ma.country.ecip,
    ma.scale.ecip,
    ma.data_source.ecip,
    ma.group_proxy.transmission.ecip
) %>% 
    #remove extra row that was specified for the plot
    filter(!(estimate == "Overall\nEffect" & (analysis == "Marginalized group")))

write_csv(table_061225, "results/stats/main_effects_table.csv") 
