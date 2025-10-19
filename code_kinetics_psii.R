###...................................###
# Bison NN, Michaletz ST. The kinetic basis of photosynthetic heat tolerance
# Code S1: R code for reproducing analyses and figures. 
# Prepared by Nicole Bison (nicole.bison@ubc.ca), August 2025.
###...................................###

##### Load packages ####
library(ggplot2)
library(tidyverse)
library("devtools")
#devtools::install_github("jinyizju/V.PhyloMaker")
library(V.PhyloMaker)
library(ggpubr)
library(gghalves)
library(phytools)
library(ggtree)
library(smatr)
library(scales)
library(caper)
library(phylolm)
library(gghalves)
library(ggpubr)
library(dplyr)
library(broom)
library(multcomp)
library(lubridate)

# Load psii heat tolerance data 
leaf_data_all <- read.csv("./clean_data/leaf_data_clean.csv")

#Get species,family counts
length(unique(leaf_data_all$species))
length(unique(leaf_data_all$family))

#### Fig 2a  ####
sp_list <- leaf_data_all[,c(4,3,2)]
tree <- phylo.maker(sp_list, tree = GBOTB.extended, output.tree = T, r = 1)
tree3 <- tree$scenario.3

ferns <- subset(leaf_data_all[,4], leaf_data_all$clade == "pteridophyte")
findMRCA(tree3, tips = ferns)

gymnosperms <- subset(leaf_data_all[,4], leaf_data_all$clade == "gymnosperm")
findMRCA(tree3, tips = gymnosperms)

angiosperms <- subset(leaf_data_all[,4], leaf_data_all$clade == "angiosperm")
findMRCA(tree3, tips = angiosperms)
# Get the root node of the entire tree
root_node <- getMRCA(tree3, tree3$tip.label)

fig2a <- ggtree(tree3, layout = "circular") +
  # Grey highlight behind the entire phylogeny
  geom_hilight(
    node = root_node,
    fill = "#4D4D4D",
    type = "rect",
    alpha = 1,            # Fully opaque
    extendto = 410
  ) +
  # Foreground highlights for the three clades
  geom_hilight(node = 326, fill = '#003f5c', type = "rect", alpha = 1) +
  geom_hilight(node = 340, fill = '#bc5090', type = "rect", alpha = 1) +
  geom_hilight(node = 180, fill = '#ffa600', type = "rect", alpha = 1) +
  # Draw the phylogeny *on top* with white branches
  geom_tree(color = "black", size = 0.3)

fig2a

#### Fig 2b  ####
leaf_data_all_combined <- leaf_data_all %>%
  mutate(clade = "All")

clade_colors <- c(
  "#ffa600",  # clade 1
  "#003f5c",  # clade 2
  "#bc5090"   # clade 3
)

label_e_sup_1digit <- function(x) {
  sapply(x, function(val) {
    # Handle NA or NaN safely
    if (is.na(val) || !is.finite(val)) return(NA)
    # Handle zero separately
    if (val == 0) return("0")
    # Format with 1 significant digit before exponent
    sci <- formatC(val, format = "e", digits = 0)
    parts <- strsplit(sci, "e", fixed = TRUE)[[1]]
    coef <- parts[1]
    exp <- gsub("^(-?)0+", "\\1", parts[2]) # remove leading zeros
    parse(text = paste0(coef, "*e^", exp))
  })
}
fig2b <- ggplot(leaf_data_all, aes(x = clade, y = ea, color = clade)) +
  scale_color_manual(values = clade_colors) +
  geom_half_violin(position = position_nudge(x = 0.3, y = 0),
                   side = 'R', adjust = 1, trim = FALSE, alpha = 0.5, width = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.5, size = 1, fill = NA) + 
  geom_jitter(size = 1.5, alpha = 0.5, width = 0.2) +
  geom_hline(yintercept = 298.45, col = "black", linetype = "dashed", size = 1) +
  geom_hline(yintercept = 61, col = "black", linetype = "dashed", size = 1)

fig2b <- fig2b +
  geom_half_violin(
    data = leaf_data_all_combined,
    aes(x = "All", y = ea),
    side = 'R',
    adjust = 1,
    trim = FALSE,
    alpha = 0.5,
    width = 0.5,
    fill = NA,            
    color = "#4D4D4D",
    position = position_nudge(x = 0.3),
    inherit.aes = FALSE
  )+
  geom_boxplot(
    data = leaf_data_all_combined,
    aes(x = "All", y = ea),
    fill = NA,          
    color = "#4D4D4D",    
    width = 0.4,
    alpha = 0.6,
    size = 1,
    outlier.shape = NA,
    inherit.aes = FALSE
  ) +
  geom_jitter(
    data = leaf_data_all_combined,
    aes(x = "All", y = ea),
    size = 1.5,
    alpha = 0.5,
    width = 0.2,
    color = "#4D4D4D",
    inherit.aes = FALSE
  )

# Adjust x-axis so "All" is on the right
fig2b <- fig2b +
  scale_x_discrete(limits = c(levels(factor(leaf_data_all$clade)), "All")) +
  xlab("Clade") +
  ylab(expression(italic(E)[a]~"(kJ mol"^{-1}*")")) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_blank(),
     legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

fig2b

#### Fig 2c ####
protein_only <- read.csv("clean_data/protein_ea.csv")
protein_only <- protein_only[,c(2,4,11)]

colnames(protein_only) <- c("species","ea","a")
protein_only$clade = "protein"


leaf_data_proteins <- leaf_data_all[,c(1,4:6)]
leaf_data_proteins$a <- exp(leaf_data_proteins$ln_a)
leaf_data_proteins <- leaf_data_proteins[,c(1,2,3,5)]
colnames(leaf_data_proteins)
colnames(protein_only)
#leaf_data_proteins$type <- "psii"
leaf_data_proteins <- rbind(leaf_data_proteins,protein_only)
leaf_data_proteins$clade <- factor(leaf_data_proteins$clade, levels = c("angiosperm", "protein","gymnosperm","pteridophyte"))


fig2c<- ggplot(leaf_data_proteins, aes(x = ea, y = a, color = clade)) +
  #scale_color_manual(values = c("#ffa600", "black", "#003f5c", "#bc5090")) +
  geom_point(
    data = subset(leaf_data_proteins, clade == "angiosperm"),
    color = "#ffa600", size = 3, alpha = 0.5
  ) +
  geom_point(
    data = subset(leaf_data_proteins, clade == "gymnosperm"),
    color = "#003f5c", size = 3, alpha = 0.5
  ) +
  geom_point(
    data = subset(leaf_data_proteins, clade == "pteridophyte"),
    color = "#bc5090", size = 3, alpha = 0.5
  )+
  geom_point(
    data = subset(leaf_data_proteins, clade == "protein"),
    color = "black", size = 3, alpha = 0.5
  )+
  # Regression line for proteins (black)
  geom_smooth(
    data = subset(leaf_data_proteins, clade == "protein"),
    method = "lm", se = FALSE,
    linewidth = 0.5, alpha = 0.5, color = "black"
  ) +
  
  # Regression line for plants (grey)
  geom_smooth(
    data = subset(leaf_data_proteins, clade %in% c("angiosperm", "gymnosperm", "pteridophyte")),
    method = "lm", se = FALSE,
    linewidth = 0.5, alpha = 0.5, linetype = "solid", color = "black"
  ) +
  
  annotate("text", label = "R2 = 0.9998") +
  
  ylab(expression(ln(italic(A))~"(s"^{-1}*")"))+
  xlab(expression(italic(E)[a]~"(kJ mol"^{-1}*")"))+
  scale_y_continuous(
    trans = "log",
    labels = label_e_sup_1digit
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )
fig2c

### Fig 2d ####
#create a data frame of temperatures (40°C to 70°C, in Kelvin)
arrhenius_plot <- leaf_data_all[,c(1,4:6)]
arrhenius_plot$ln_a_pred <- 0.3716727*arrhenius_plot$ea + 1

temperature_data <- data.frame(
  temp_C = seq(30, 70, by = 1)  
)
temperature_data$temp_K <- temperature_data$temp_C + 273.15  #convert to Kelvin

R= 0.0083145 #universal gas constant
plot_data <- do.call(rbind, lapply(1:nrow(arrhenius_plot), function(i) {
  data.frame(
    species = arrhenius_plot$species[i],
    clade = arrhenius_plot$clade[i],
    temp_C = temperature_data$temp_C,
    RT = R * temperature_data$temp_K,
    ln_k = arrhenius_plot$ln_a_pred[i] - (arrhenius_plot$ea[i] / (R * temperature_data$temp_K))
    
  )
}))


leaf_subset <- plot_data

leaf_subset$temp_K <- leaf_subset$temp_C +273.15
leaf_subset$temp_axis <- -1/(leaf_subset$temp_K*R)
temp_axis=-1/(leaf_subset$temp_K*R)
leaf_subset <- subset(leaf_subset, leaf_subset$temp_C <70)
leaf_subset <- subset(leaf_subset, leaf_subset$temp_C >30)
leaf_subset$k <- exp(leaf_subset$ln_k)

#filter to non-zero k values (since log(0) = -Inf) — adjust filter as needed
k_vals <- leaf_subset$k
k_vals <- k_vals[k_vals > 0 & !is.na(k_vals)]

#get y-axis range for log scale
ymin_val <- min(k_vals)
ymax_val <- max(k_vals)

label_e_sup <- function(x) {
  s <- scales::scientific_format()(x)        
  s <- gsub("e-0*", "e-", s)                 
  s <- gsub("e", "*e^", s)                   
  parse(text = s)                            
}



fig2d <- ggplot() +
  geom_line(data = subset(leaf_subset, clade == "angiosperm"), 
            aes(x = temp_axis, y = k, group = species),
            linewidth = 0.5, alpha = 0.2, color = "#ffa600") +
  geom_line(data = subset(leaf_subset, clade == "gymnosperm"), 
            aes(x = temp_axis, y = k, group = species), 
            linewidth = 0.5, alpha = 0.2, color = "#003f5c") +
  geom_line(data = subset(leaf_subset, clade == "pteridophyte"), 
            aes(x = temp_axis, y = k, group = species), 
            linewidth = 0.5, alpha = 0.2, color = "#bc5090") +
  
  scale_x_continuous(
    name = expression("1 /"~italic(RT)~"(mol"~J^{-1}*")"),
    breaks = c(-0.39, -0.38, -0.37, -0.36),
    labels = scales::number_format(accuracy = 0.01),  
    sec.axis = sec_axis(
      ~ -(1 / (. * R)) - 273.15,
      name = "Temperature (°C)",
      breaks = seq(30, 70, by = 5)
    )
  ) +
  scale_y_continuous(
    trans = "log",
    labels = label_e_sup_1digit,
    expand = c(0, 0)
  )+
  annotate("rect", xmin = -0.3696931, xmax = -0.3734591, 
           ymin = ymin_val, ymax = ymax_val, 
           alpha = 0.2, fill = "black") +
  geom_vline(xintercept = -0.3716727, col = "black", linetype = "dashed", size = 1) +
  ylab(expression(italic(k)~(s^{-1}))) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.line = element_line(color = "black"),
    legend.position = "none"
  )
fig2d

##### Fig 2e ####
lit_review<-read.csv("./clean_data/lit_review_cleaned_oct2025_tcrit_only.csv")
length(unique(lit_review$species))
length(unique(lit_review$unique_id))
length(lit_review$heat_tolerance)

# Get world map data
min(lit_review$latitude,na.rm = T)
max(lit_review$latitude,na.rm = T)

df_summary <- lit_review %>%
  group_by(longitude, latitude, species) %>%
  summarise(count = n(), .groups = "drop")

data_complete_map <- df_summary[complete.cases(df_summary),]

world_map <- map_data("world")
world_map

# Plot map with points
fig2e<- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "transparent", color = "black") + 
  geom_point(data = df_summary, aes(x = longitude, y = latitude, size = count), 
             colour = "#539987" ,alpha = 0.9,position = position_jitter(width = 0.5, height = 0.5)) + 
  #geom_point(data = subset(df_summary, df_summary$clade == "gymnosperm"), aes(x = longitude, y = latitude, size = count), 
  #        colour = "#003f5c" ,alpha = 0.7, position = position_jitter(width = 0.5, height = 0.5)) + 
  # geom_point(data = subset(df_summary, df_summary$clade == "pteridophyte"), aes(x = longitude, y = latitude, size = count), 
  #           colour = "#bc5090" ,alpha = 0.7, position = position_jitter(width = 0.5, height = 0.5)) + 
  theme_minimal() +
  scale_size_continuous(range = c(2, 8)) +  # Adjust point sizes
  
  theme(panel.background = element_blank(),  # Removes the background color
        panel.grid.major = element_blank(),  # Removes major grid lines
        panel.grid.minor = element_blank(),  # Removes minor grid lines
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 14),  # Change axis labels text size
        axis.text = element_text(size = 12),
        # axis.text.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))+
  labs(title = "", x = "Longitude", y = "Latitude")
fig2e

#### Fig 2f ####
fig2f<-ggplot(lit_review, aes(x=heat_tolerance)) +
  #geom_histogram( binwidth=2, fill="#93B7BE", color="black", alpha=0.9) +
  geom_histogram( binwidth=2, fill="#539987", color="black", alpha=0.9) +
  geom_vline(xintercept = 50.47, col="black", linetype="dashed", linewidth=1.2)+
  ylab("Count")+
  xlim(30,70)+
  # ylim(-35,5)+
  xlab("Temperature (C)")+
  theme(panel.background = element_blank(),  # Removes the background color
        panel.grid.major = element_blank(),  # Removes major grid lines
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 14),  # Change axis labels text size
        axis.text = element_text(size = 12), # Removes minor grid lines
        axis.line = element_line(color = "black"),
        #axis.text.x = element_blank(),
        legend.position = "none")
fig2f

#### Fig 2 ####
ggarrange(fig2a,fig2b, fig2c, fig2d, fig2e, fig2f, 
          labels = c("a", "b", "c","d","e","f"),
          ncol = 2, nrow = 3)
#ggsave("./figures/final/fig2_v19.pdf", last_plot(), height = 25, width = 25, units = "cm")

####### Figure 3 #######
#Load kinetics, heat tolerance data from seasonal dataset 
df_seasons<-read.csv("./clean_data/seasonal_data.csv")
df_seasons$season <- factor(df_seasons$season, levels = c("spring","summer","fall"))

#Load climate data
df_seasons_climate <- read.csv("./clean_data/seasonal_with_climate.csv")
  
#### Fig 3a ####
fig3a<- ggplot(df_seasons, aes(x= season, y = tcrit))+
  scale_color_manual(values = c("#A8E6CF","#FFD54F","#FF8C42"))+
  geom_jitter(aes(color = season, shape = Species),size =2,alpha =0.5,width =0.2)+ 
  geom_boxplot(aes(group = season),outlier.shape = NA, width =0.5,alpha=0.5,size=1, fill = NA)+
  xlab("Season")+
  ylab(expression(italic(T)[crit]~"(°C)"))+
  theme(panel.background = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 14),  
        axis.text = element_text(size = 12),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
fig3a


df_seasons_climate$season <- factor(df_seasons_climate$season, levels = c("spring","summer","fall"))

#### Fig 3b ####
fig3b <- ggplot(df_seasons_climate, aes(x= Heat_Stress_Days_Cumulative, y = tcrit))+
  scale_color_manual(values = c("#A8E6CF","#FFD54F","#FF8C42"))+
  geom_jitter(aes(color = season),size =2,alpha =0.5,width =0.2)+ 
  geom_smooth(aes(fill= season), method = "lm", color = "black", fill = "black")+
  scale_fill_manual(values = c("#A8E6CF","#FFD54F","#FF8C42"))+
  xlab("Cumulative days > 21°C")+
  ylab(expression(italic(T)[crit]~"(°C)"))+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
fig3b

#### Fig 3c ####
fig3c <- ggplot(df_seasons, aes(x= season, y = ea))+
  scale_color_manual(values = c("#A8E6CF","#FFD54F","#FF8C42"))+
  geom_jitter(aes(color = season, shape = Species),size =2,alpha =0.5,width =0.2)+ 
  geom_boxplot(aes(group = season),outlier.shape = NA, width =0.5,alpha=0.5,size=1, fill = NA)+
  scale_fill_manual(values = c("#A8E6CF","#FFD54F","#FF8C42"))+
  xlab("Season")+
  ylab(expression(italic(E)[a]~"(kJ mol"^{-1}*")"))+
  theme(panel.background = element_blank(),  
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 14),  
        axis.text = element_text(size = 12),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
fig3c

#### Fig 3d ####
fig3d <-ggplot(df_seasons_climate, aes(x= Heat_Stress_Days_Cumulative, y = ea))+
  scale_color_manual(values = c("#A8E6CF","#FFD54F","#FF8C42"))+
  geom_jitter(aes(color = season),size =2,alpha =0.5,width =0.2)+ 
  geom_smooth(aes(fill= season), method = "lm", color = "black", fill = "black")+
  scale_fill_manual(values = c("#A8E6CF","#FFD54F","#FF8C42"))+
  xlab("Cumulative days > 21°C")+
  ylab(expression(italic(E)[a]~"(kJ mol"^{-1}*")"))+
  theme(panel.background = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))

#### Fig 3e ####
fig3e <- ggplot(df_seasons, aes(x= season, y = ln_a))+
  scale_color_manual(values = c("#A8E6CF","#FFD54F","#FF8C42"))+
  geom_jitter(aes(color = season, shape = Species),size =2,alpha =0.5,width =0.2)+ 
  geom_boxplot(aes(group = season),outlier.shape = NA, width =0.5,alpha=0.5,size=1, fill = NA)+
  scale_fill_manual(values = c("#A8E6CF","#FFD54F","#FF8C42"))+
  xlab("Season")+
  ylab(expression(ln(italic(A))~"(s"^{-1}*")"))+
  theme(panel.background = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 14),  
        axis.text = element_text(size = 12),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
fig3e

#### Fig 3f ####
fig3f <-ggplot(df_seasons_climate, aes(x= Heat_Stress_Days_Cumulative, y = exp(ln_a)))+
  scale_color_manual(values = c("#A8E6CF","#FFD54F","#FF8C42"))+
  scale_y_continuous(
    trans = "log",
    labels = function(x) parse(text = paste0("e^", round(log(x), 1)))
  )+
  geom_jitter(aes(color = season),size =2,alpha =0.5,width =0.2)+ 
  scale_fill_manual(values = c("#A8E6CF","#FFD54F","#FF8C42"))+
  xlab("Cumulative days > 21°C")+
  ylab(expression(ln(italic(A))~"(s"^{-1}*")"))+
  theme(panel.background = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
fig3f

##### Fig 3g ####
df_seasons <- df_seasons %>%
  mutate(season_lumped = case_when(
    season == "spring" ~ "spring",
    season %in% c("summer", "fall") ~ "summer_fall",
    TRUE ~ NA_character_
  ))

fig3g <- ggplot(df_seasons, aes(x = ea, y = exp(ln_a), color = season)) +
  scale_color_manual(values = c("#A8E6CF","#FFD54F","#FF8C42")) +
  
  geom_point(size = 2.5, alpha = 0.5, aes(color = season)) +
  
  geom_smooth(aes(group = season_lumped), color = "black", method = "lm", se = FALSE,
              linewidth = 0.5, alpha = 0.5, linetype = "solid")+
  
  annotate("text", label = "R2 = 0.9998") +
  
  ylab(expression(ln(italic(A))~"(s"^{-1}*")"))+
  xlab(expression(italic(E)[a]~"(kJ mol"^{-1}*")"))+
  scale_y_continuous(
    trans = "log",
    labels = function(x) parse(text = paste0("e^", round(log(x), 1)))
  )+
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )
fig3g

#Filter to spring, and summer+fall as slopes were not signficantly different between summer and fall
spring      <- df_seasons %>% filter(season == "spring")
summer_fall <- df_seasons %>% filter(season %in% c("summer", "fall"))

seasonal_data <- list(
  spring = spring,
  summer_fall = summer_fall
)

#bootstrap function with SMA regression
bootstrap_sma <- function(data, n_boot = 1000) {
  replicate(n_boot, {
    sample_data <- data[sample(nrow(data), replace = TRUE), ]
    
    model <- try(sma(ln_a ~ ea, data = sample_data), silent = TRUE)
    if (inherits(model, "try-error")) return(c(slope = NA, Tiso = NA))
    
    slope_vec <- coef(model)
    if (!("slope" %in% names(slope_vec))) return(c(slope = NA, Tiso = NA))
    
    slope <- slope_vec["slope"]
    Tiso <- (1 / (0.008314 * slope)) - 273.15
    
    c(slope = slope, Tiso = Tiso)
  }, simplify = TRUE)
}

#apply bootstrap
boot_results_list <- lapply(seasonal_data, bootstrap_sma, n_boot = 1000)

#convert to dataframe with season labels
boot_df_spring <- as.data.frame(t(boot_results_list$spring))
boot_df_spring$season <- "spring"
boot_df_spring$season_new <- "spring"

boot_df_summer_fall <- as.data.frame(t(boot_results_list$summer_fall))
boot_df_summer_fall$season <- "summer_fall"

boot_df_summer <- boot_df_summer_fall
boot_df_summer$season_new <-"summer"

boot_df_fall <- boot_df_summer_fall
boot_df_fall$season_new <-"fall"

#combine
seasons_Tiso <- rbind(boot_df_spring, boot_df_summer,boot_df_fall)
seasons_Tiso <- seasons_Tiso[,c(1,2,4)]
colnames(seasons_Tiso) <- c("slope", "Tiso", "season")

seasons_Tiso$season <- factor(seasons_Tiso$season, levels = c("spring","summer","fall"))

##### Fig 3h ####
fig3h <-ggplot(seasons_Tiso, aes(x= season, y = Tiso))+
  geom_jitter(aes(color = season),size =2,alpha =0.1,width =0.2)+
  geom_boxplot(aes(group = season),outlier.shape = NA, width =0.5,alpha=0.5,size=1, fill = NA)+
  scale_color_manual(values = c("#A8E6CF","#FFD54F","#FF8C42"))+
  geom_hline(yintercept = 55.1, col="black", linetype="dashed", size=1)+
  geom_hline(yintercept = 58.7, col="black", linetype="dashed", size=1)+
  geom_hline(yintercept = 43.7, col="black", linetype="dashed", size=1)+
  geom_hline(yintercept = 47.4, col="black", linetype="dashed", size=1)+
  ylab(expression(italic(T)[iso]~"(°C)"))+
  theme(panel.background = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 14),  
        axis.text = element_text(size = 12),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
fig3h

##### Fig 3 ####
fig3<-ggarrange(fig3a,fig3b, fig3c,fig3d, fig3e, fig3f,fig3g,fig3h,
                labels = c("a", "b", "c","d","e","f","g","h"),
                ncol = 2, nrow = 4)
fig3
#ggsave("./figures/final/fig3_v19.pdf", last_plot(), height = 26, width = 15, units = "cm")

####### Extended data #####
protein_only <- read.csv("clean_data/protein_only.csv")
protein_only$ln_a <- log(protein_only$a)
protein_only <- protein_only[,c(2:4,6)]
protein_only$clade <- "NA"

lipid_only <- read.csv("clean_data/lipid_only.csv")
lipid_only <- lipid_only[,c(2:5)]
lipid_only$clade <- "NA"

leaf_data <- leaf_data_all[,c(1,5,6)]
leaf_data$type <- "psii"
leaf_data$material <- "NA"

leaf_proteins_lipids <- rbind(leaf_data,lipid_only,protein_only)

#### Table 1 ####
#summarize across clades and type
df_summary_clade <- leaf_data %>%
  group_by(clade) %>%
  summarise(
    n = n(),
    min = min(ea),
    max = max(ea),
    median = median(ea),
    mean = mean(ea),
    sd = sd(ea),
    se = sd / sqrt(n),
    ci95 = qt(0.975, df = n - 1) * se,
    lower = mean - ci95,
    upper = mean + ci95
  )

#Now for summary for proteins, lipids and for PSII across ALL clades
df_summary<- leaf_proteins_lipids %>%
  group_by(type) %>%
  summarise(
    n = n(),
    min = min(ea),
    max = max(ea),
    median = median(ea),
    mean = mean(ea),
    sd = sd(ea),
    se = sd / sqrt(n),
    ci95 = qt(0.975, df = n - 1) * se,
    lower = mean - ci95,
    upper = mean + ci95
  )

#### Table 2 ####
#Get general data and protein data
leaf_data_proteins <- rbind(leaf_data,protein_only)
leaf_data_proteins <- leaf_data_proteins[,c(4,2,3)]

#Get seasonal data by season
spring <- df_seasons %>%
  filter(season == "spring") %>%
  dplyr::select(ea, ln_a)
colnames(spring) <- c("ea","ln_a")
spring$type <- "spring"

summer <- df_seasons %>%
  filter(season == "summer") %>%
  dplyr::select(ea, ln_a)
colnames(summer) <- c("ea","ln_a")
summer$type <- "summer"

fall <- df_seasons %>%
  filter(season == "fall") %>%
  dplyr::select(ea, ln_a)
colnames(fall) <- c("ea","ln_a")
fall$type <- "fall"

leaf_proteins_seasons <- rbind(leaf_data_proteins,spring,summer,fall)

#apply bootstrap to all data
datasets_all <- list(
  spring = spring,
  summer = summer,
  fall=fall,
  summer_fall = summer_fall,
  psii = subset(leaf_data_proteins, leaf_data_proteins$type == "psii"),
  protein = subset(leaf_data_proteins, leaf_data_proteins$type == "Protein")
)

set.seed(123)  #for reproducibility
boot_results_list <- lapply(datasets_all, bootstrap_sma, n_boot = 1000)

boot_results_list <- lapply(boot_results_list, function(x) {
  # If the element is a data frame or matrix, rename its columns
  if (is.data.frame(x) || is.matrix(x)) {
    colnames(x) <- sub("\\.slope$", "", colnames(x))
    return(x)
  }
  # If it's a list, recurse deeper
  if (is.list(x)) {
    return(lapply(x, function(y) {
      if (is.data.frame(y) || is.matrix(y)) {
        colnames(y) <- sub("\\.slope$", "", colnames(y))
        return(y)
      } else {
        return(y)
      }
    }))
  }
  # Otherwise, leave untouched
  return(x)
})

#convert to dataframe with dataset label
boot_df_spring <- as.data.frame(t(boot_results_list$spring))
colnames(boot_df_spring) <- gsub("\\.slope$", "", colnames(boot_df_spring))
boot_df_spring$season <- "spring"
boot_df_spring$season_new <- "spring"

boot_df_summer_fall <- as.data.frame(t(boot_results_list$summer_fall))
colnames(boot_df_summer_fall) <- gsub("\\.slope$", "", colnames(boot_df_summer_fall))
boot_df_summer_fall$season <- "summer_fall"

boot_df_summer <- as.data.frame(t(boot_results_list$summer))
boot_df_summer$season_new <-"summer"

boot_df_fall <- as.data.frame(t(boot_results_list$fall))
boot_df_fall$season_new <-"fall"

boot_df_psii<- as.data.frame(t(boot_results_list$psii))
colnames(boot_df_psii) <- gsub("\\.slope$", "", colnames(boot_df_psii))
boot_df_psii$type <- "psii"

boot_df_protein<- as.data.frame(t(boot_results_list$protein))
colnames(boot_df_protein) <- gsub("\\.slope$", "", colnames(boot_df_protein))
boot_df_protein$type <- "protein"

#Now get Tiso for each df
#Spring
quantile(boot_df_spring$Tiso, probs = c(0.025, 0.975), na.rm = TRUE)
mean(boot_df_spring$Tiso)

#Summer
quantile(boot_df_summer$Tiso, probs = c(0.025, 0.975), na.rm = TRUE)
mean(boot_df_summer$Tiso)

#Fall
quantile(boot_df_fall$Tiso, probs = c(0.025, 0.975), na.rm = TRUE)
mean(boot_df_fall$Tiso)

# Summer - Fall
quantile(boot_df_summer_fall$Tiso, probs = c(0.025, 0.975), na.rm = TRUE)
mean(boot_df_summer_fall$Tiso)

#PSII
quantile(boot_df_psii$Tiso, probs = c(0.025, 0.975), na.rm = TRUE)
mean(boot_df_psii$Tiso)

#Protein
quantile(boot_df_protein$Tiso, probs = c(0.025, 0.975), na.rm = TRUE)
mean(boot_df_protein$Tiso)

#Now get slope for each df
#Spring
quantile(boot_df_spring$slope, probs = c(0.025, 0.975), na.rm = TRUE)
mean(boot_df_spring$slope)

#Summer
quantile(boot_df_summer$slope, probs = c(0.025, 0.975), na.rm = TRUE)
mean(boot_df_summer$slope)

#Fall
quantile(boot_df_fall$slope, probs = c(0.025, 0.975), na.rm = TRUE)
mean(boot_df_fall$slope)

# Summer - Fall
quantile(boot_df_summer_fall$slope, probs = c(0.025, 0.975), na.rm = TRUE)
mean(boot_df_summer_fall$slope)

#PSII
quantile(boot_df_psii$slope, probs = c(0.025, 0.975), na.rm = TRUE)
mean(boot_df_psii$slope)

#Protein
quantile(boot_df_protein$slope, probs = c(0.025, 0.975), na.rm = TRUE)
mean(boot_df_protein$slope)

#### Table 3 ####
# ******Note need to change R version to 4.0.2 in this step to see output with pairwise p-values
sma_model <- sma(ln_a ~ ea * type, data = leaf_proteins_seasons, multcomp = TRUE)
sma_model


##### Table 4 #######
climate_raw <- read.csv("./clean_data/climate_data_2023.csv")
climate_raw$date <- gsub("-", "", climate_raw$Date.Time)
climate_raw <- climate_raw %>%
  mutate(date = mdy(date),
         date = format(date, "%Y%m%d")) 

climate <- climate_raw[, c(8, 9, 10, 12)]
colnames(climate)[1:4] <- c("max_temp", "min_temp", "mean_temp", "Date")
# Define thresholds to test
thresholds <- 10:26
slope_results <- data.frame()

climate <- climate %>%
  mutate(Date = as.integer(Date))  # make Date column integer globally

# Loop through each threshold
for (thresh in thresholds) {
  
  # Recalculate cumulative heat stress days for current threshold
  climate_thresh <- climate %>%
    arrange(Date) %>%
    mutate(
      Max_Temp_Lag1 = lag(max_temp, 1),
      Heat_Stress_Days_Cumulative = cumsum(
        if_else(Max_Temp_Lag1 > thresh, 1L, 0L, missing = 0L)
      )
    ) %>%
    dplyr::select(Date, Heat_Stress_Days_Cumulative)
  
  # Merge with trait data
  model_data <- df_seasons_climate %>%
    dplyr::select(Date, ea, Species) %>%
    left_join(climate_thresh, by = "Date") %>%
    filter(!is.na(ea), !is.na(Heat_Stress_Days_Cumulative), !is.na(Species))
  
  # Fit linear model with interaction
  if (nrow(model_data) >= 5) {
    lm_model <- lm(ea ~ Heat_Stress_Days_Cumulative, data = model_data)
    
    # Extract only the main effect of heat stress (not species-specific terms)
    tidy_terms <- tidy(lm_model) %>%
      filter(term == "Heat_Stress_Days_Cumulative") %>%
      mutate(threshold = thresh) %>%
      dplyr::select(threshold, slope = estimate, p.value)
    
    slope_results <- bind_rows(slope_results, tidy_terms)
  }
}

# View final summary table
print(slope_results)

##### Table 5 #######
prothermdb_a_thaliana_hsp <- read.csv("clean_data/prothermdb_a_thaliana_hsp.csv")
prothermdb_a_thaliana_hsp <- prothermdb_a_thaliana_hsp[,c(2,29)]
head(prothermdb_a_thaliana_hsp)

prothermdb_a_thaliana_hsp <- prothermdb_a_thaliana_hsp %>%
  mutate(HSP_group = case_when(
    grepl("20", PROTEIN, ignore.case = TRUE) ~ "HSP20",
    grepl("60", PROTEIN, ignore.case = TRUE) ~ "HSP60",
    grepl("70", PROTEIN, ignore.case = TRUE) ~ "HSP70",
    grepl("90|91", PROTEIN, ignore.case = TRUE) ~ "HSP90",
    TRUE ~ NA_character_
  ))
prothermdb_a_thaliana_hsp <- prothermdb_a_thaliana_hsp %>%
  mutate(Tm_numeric = as.numeric(sub(" .*", "", Tm_.C.)))

prothermdb_a_thaliana_hsp <- prothermdb_a_thaliana_hsp %>%
  filter(!is.na(HSP_group))

summary_hsp <- prothermdb_a_thaliana_hsp %>%
  group_by(HSP_group) %>%
  summarise(
    mean_Tm = mean(Tm_numeric, na.rm = TRUE),
    n = n()
  )

#### Figure 2 ####
ggplot(lit_review, aes(x=heat_tolerance)) +
  #geom_histogram( binwidth=2, fill="#93B7BE", color="black", alpha=0.9) +
  geom_histogram( binwidth=2, fill="#539987", color="black", alpha=0.9) +
  geom_vline(xintercept = 47.4, col="black", linetype="dashed", linewidth=1.2)+ #HSP90
  geom_vline(xintercept = 43.7, col="black", linetype="dashed", linewidth=1.2)+ #HSP70
  geom_vline(xintercept = 58.7, col="black", linetype="dashed", linewidth=1.2)+ #HSP20
  geom_vline(xintercept = 55.1, col="black", linetype="dashed", linewidth=1.2)+ #HSP60
  
  ylab("Count")+
  xlim(30,70)+
  # ylim(-35,5)+
  xlab("Temperature (C)")+
  theme(panel.background = element_blank(),  # Removes the background color
        panel.grid.major = element_blank(),  # Removes major grid lines
        panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 14),  # Change axis labels text size
        axis.text = element_text(size = 12), # Removes minor grid lines
        axis.line = element_line(color = "black"),
        #axis.text.x = element_blank(),
        legend.position = "none")
ggsave("ed_Figure2_hsp.pdf",last_plot(), height = 12, width =15, units = "cm")

#### Figure 3 ####
protein_lipids <- subset(leaf_proteins_lipids, leaf_proteins_lipids$type != "psii")
protein_lipids_summary <- protein_lipids %>%
  group_by(type) %>%
  summarise(
    mean = mean(ea),
    median = median(ea),
    n = n(),
    sd = sd(ea),
    se = sd / sqrt(n),
    ci95 = qt(0.975, df = n - 1) * se,
    lower = mean - ci95,
    upper = mean + ci95
  )
ggplot(protein_lipids_summary, aes(x = type, y = mean,fill = type,color = type))+
  #scale_fill_manual(values = c("red","blue"))+
  #scale_colour_manual(values = c("red","blue"))+
  scale_color_manual(values = c("black","black"))+
  geom_errorbar(data = protein_lipids_summary, 
                aes(x = type, ymin = lower, ymax = upper, color = type), width = 0.1,linewidth = 1.5) +
  geom_point(size = 4) +
  xlab("")+
  # xlim(0,900)+
  ylab(expression(italic(E)[a]~"(kJ mol"^{-1}*")"))+
  theme(panel.background = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12),
        #axis.text.y = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))

#### Figure 4 ####
d <- df_seasons[, 2:10]
colnames(d)[3] <- "Date"

results <- data.frame()

#loop over averaging windows
for (j in 1:120) {
  n <- j
  x <- lag(climate$mean_temp, 1)
  len <- length(x)
  y <- rep(NA, len)
  
  for (i in 1:len) {
    if (i < n) next
    y[i] <- mean(x[(i - (n - 1)):(i)], na.rm = TRUE)
  }
  
  climate$Rolling_Mean_C <- y
  clim_cut <- climate %>%
    dplyr::select(Date, Rolling_Mean_C) %>%
    mutate(Rolling_Mean_Days = j)
  
  d$Date <- as.integer(d$Date)
  clim_cut$Date <- as.integer(clim_cut$Date)
  traits_with_weather <- merge(d, clim_cut, by = "Date")
  
  results <- bind_rows(results, traits_with_weather)
}

#pivot to wide format
results_wider <- results %>%
  pivot_wider(names_from = Rolling_Mean_Days, values_from = Rolling_Mean_C)

summary(results_wider)

#fit linear models for each trait and window size
fit_temp_models <- function(data, temp_cols = as.character(1:120)) {
  responses <- c("tcrit", "ea", "ln_a")
  
  map_dfr(responses, function(resp) {
    map_dfr(temp_cols, function(temp_col) {
      formula <- as.formula(paste(resp, "~", paste0("`", temp_col, "`")))
      fit <- tryCatch(lm(formula, data = data), error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      
      coef_info <- tidy(fit) %>% filter(term != "(Intercept)")
      model_info <- glance(fit)
      
      tibble(
        response = resp,
        temp_col = temp_col,
        slope = coef_info$estimate,
        std_error = coef_info$std.error,
        p_value = coef_info$p.value,
        r_squared = model_info$r.squared,
        adj_r_squared = model_info$adj.r.squared
      )
    })
  })
}

model_results <- fit_temp_models(results_wider)

#convert temp_col to numeric for plotting
model_results <- model_results %>%
  mutate(temp_col = factor(temp_col, levels = as.character(1:120)),
         averaging_window = as.numeric(as.character(temp_col)))

# Filter significant models
significant_models <- model_results %>%
  filter(p_value < 0.05) %>%
  mutate(temp_col = factor(temp_col, levels = as.character(1:120)))

# Plot 2: All models by response variable
ggplot(model_results, aes(x = averaging_window, y = adj_r_squared, group = response, color = response)) +
  geom_point(size = 1, alpha = 0.75) +
  geom_line(linewidth = 2, alpha = 0.25) +
  xlab("Temperature averaging window (days)") +
  ylab(expression(R^2)) +
  theme(
    text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.background = element_rect(fill = 'white'),
    plot.background = element_rect(fill = 'white'),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent')
  )

##### Statistical tests ####
##### Figure 2 a
#Is Ea phylogenetically structured?
trait <- setNames(leaf_data_all$ea, c(leaf_data_all$species))
phylosig(tree3, trait, method="lambda", test=TRUE, nsim=999)

##### Figure 2 b

#Variation in Ea across clades
anova_model <- aov(ea ~ clade, data = leaf_data)
summary(anova_model)

#Is there a difference between Ea for PSII inactivation and protein denaturation?
leaf_protein <- subset(leaf_proteins_lipids,leaf_proteins_lipids$type!= "Lipid bilayer")
anova_model <- aov(ea ~ type, data = leaf_protein)
summary(anova_model)

#What if we compore each clade to protein as a baseline?
#post-hoc with protein as baseline
leaf_protein$clade[leaf_protein$type == "Protein"] <- "Protein"
leaf_protein$clade <- relevel(factor(leaf_protein$clade), ref = "Protein")

#fit ANOVA
anova_model <- aov(ea ~ clade, data = leaf_protein)

#pairwise comparisons vs protein
posthoc <- emmeans(anova_model, pairwise ~ clade)

summary(posthoc$contrasts, infer = TRUE, adjust = "none") %>%
  as.data.frame() %>%
  filter(grepl("Protein", contrast)) %>%
  mutate(p.value = format(p.value, digits = 10, scientific = TRUE))

#Is there a difference between Ea for PSII inactivation and lipid bilayer breakdown?
leaf_lipid <- subset(leaf_proteins_lipids,leaf_proteins_lipids$type!= "Protein")
anova_model <- aov(ea ~ type, data = leaf_lipid)
summary(anova_model)

#What if we compore each clade to lipid bilayer as a baseline?
#post-hoc with lipid bilayer  as baseline
leaf_lipid$clade[leaf_lipid$type == "Lipid bilayer"] <- "Lipid bilayer"
leaf_lipid$clade <- relevel(factor(leaf_lipid$clade), ref = "Lipid bilayer")

#fit ANOVA
anova_model <- aov(ea ~ clade, data = leaf_lipid)

#pairwise comparisons vs protein
posthoc <- emmeans(anova_model, pairwise ~ clade)

summary(posthoc$contrasts, infer = TRUE, adjust = "none") %>%
  as.data.frame() %>%
  filter(grepl("Lipid bilayer", contrast)) %>%
  mutate(p.value = format(p.value, digits = 10, scientific = TRUE))

##### Figure 2c
#Is there a diff in intercepts between proteins and psii?
protein_leaf <- subset(leaf_proteins_lipids, leaf_proteins_lipids$type != "Lipid bilayer")
sma_model <- elev.com(y=protein_leaf$ln_a, x=protein_leaf$ea, groups=protein_leaf$type)
sma_model


##### Figure 2e
## Proportion of Tcrit observations above estimated Tiso of 50C
# Total observations
length(subset(lit_review$species, lit_review$heat_tolerance > 50.422))/length(lit_review$species)

#Grouped by species
species <- lit_review %>%
  group_by(species) %>%
  summarise(heat_tolerance = mean(heat_tolerance, na.rm = T))
length(subset(species$species, species$heat_tolerance > 50.422))/length(species$species)

## Proportion of Tcrit observations above estimated Tm from HSP60 59C
# Total observations
length(subset(lit_review$species, lit_review$heat_tolerance > 59))/length(lit_review$species)
#Grouped by species
length(subset(species$species, species$heat_tolerance > 59))/length(species$species)

#####Figure 3
#What is average Tcrit by season?
df_seasons %>%
  group_by(season) %>%
  summarise(mean_Tcrit = mean(tcrit, na.rm = TRUE))

#Do parameters vary seasonally?
#Tcrit
anova_model <- aov(tcrit ~ season, data = df_seasons)
summary(anova_model)

#Now posthoc comparison
tukey_glht <- glht(anova_model, linfct = mcp(season = "Tukey"))
summary(tukey_glht)

#Ea
res_aov <- aov(ea ~ season,
               data = df_seasons
)
summary(res_aov)

#ln(A)
res_aov <- aov(ln_a ~ season,
               data = df_seasons
)
summary(res_aov)

#Do parameters correlate with cumulative heat exposure?
summary(lm(tcrit~Heat_Stress_Days_Cumulative, data=df_seasons_climate))
summary(lm(ea~Heat_Stress_Days_Cumulative, data=df_seasons_climate))
summary(lm(ln_a~Heat_Stress_Days_Cumulative, data=df_seasons_climate))

#Is there a difference between bootsrapped Tiso for proteins and PSII?

#testing for difference between Tiso for proteins and PSII
#extract Tiso values
T1 <- boot_df_protein$Tiso
T2 <- boot_df_psii$Tiso

#drop NAs
T1 <- T1[!is.na(T1)]
T2 <- T2[!is.na(T2)]

#observed difference in means
obs_diff <- mean(T1) - mean(T2)

# Generate bootstrap distribution of differences
set.seed(123)  # For reproducibility
n_boot <- 10000
boot_diffs <- replicate(n_boot, {
  m1 <- mean(sample(T1, replace = TRUE))
  m2 <- mean(sample(T2, replace = TRUE))
  m1 - m2
})

#two-tailed p-value
mean(abs(boot_diffs) >= abs(obs_diff))

