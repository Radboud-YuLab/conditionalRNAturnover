# compiles progeny scores with QC figures

library(tidyverse)
library(tidymodels)
library(here)
library(pheatmap)

# loading the progeny scores
bazzini <- readRDS(file = here("4_processed_data/wholeset/progeny/bazzini/20240822_Bazzini_progeny_log10_cpm.RDS"))
rownames(bazzini) <- c("Bazzini_4sU_K562_1_A","Bazzini_4sU_K562_1_B","Bazzini_4sU_K562_1_C")
Bazzini_4sU_K562 <- colMeans(bazzini)

cramer <- readRDS(file = here("4_processed_data/wholeset/progeny/cramer1/20240822_cramer1_progeny_log10_tpm.RDS"))
Cramer_4sU_K562 <- colMeans(cramer)

dieterich <- readRDS(file = here("4_processed_data/wholeset/progeny/dieterich/20240822_dieterich_progeny_log10_tpm.RDS"))
row.names(dieterich) <- c("Dieterich_4sU_HEK293","Dieterich_4sU_MCF7")
dieterich

mortavazi <- readRDS(file = here("4_processed_data/wholeset/progeny/mortavazi/20240822_mortavazi_progeny_log10_tpm.RDS"))
Mortazavi_4sU_GM12878 <- mortavazi[1,]

rissland <- readRDS(file = here("4_processed_data/wholeset/progeny/rissland/20240822_rissland_progeny_log10_tpm.RDS"))
Rissland_4sU_HEK293 <- colMeans(rissland)

mulder <- readRDS(file = here("4_processed_data/wholeset/progeny/mulder/20240822_mulder_progeny_log10_tpm.RDS"))
mulder[c(6,7,8),]

Mulder_4sU_KC.AG1478 <- colMeans(mulder[c(1,2),])
Mulder_4sU_KC.BMP27 <- colMeans(mulder[c(3,4,5),])
Mulder_4sU_KC.Vehicle <- colMeans(mulder[c(6,7,8),])


# these scores used for modelling 
progeny_results <- rbind(dieterich,
                         Bazzini_4sU_K562,
                         Cramer_4sU_K562,
                         Mortazavi_4sU_GM12878,
                         Rissland_4sU_HEK293,
                         Mulder_4sU_KC.AG1478,
                         Mulder_4sU_KC.BMP27,
                         Mulder_4sU_KC.Vehicle
)

progeny_results_unscaled <- progeny_results
# saveRDS(progeny_results_unscaled, file = "4_processed_data/wholeset/progeny/unscaled/20240910_progeny_results_all_mean_per_condition_unscaled.RDS")


# with replicates
progeny_results_with_replicates <-
  rbind(Bazzini_4sU_K562, # only one half life value
        Cramer_4sU_K562, 
        dieterich,
        Mortazavi_4sU_GM12878,
        rissland,
        mulder)
rownames(progeny_results_with_replicates)[seq(6,15)]
rownames(progeny_results_with_replicates)[seq(6,15)] <- c("Rissland_4sU_HEK293_2","Rissland_4sU_HEK293_3",
                                                          "Mulder_4sU_KC.AG1478_1","Mulder_4sU_KC.AG1478_2",
                                                          "Mulder_4sU_KC.BMP27_1","Mulder_4sU_KC.BMP27_2",
                                                          "Mulder_4sU_KC.BMP27_3","Mulder_4sU_KC.Vehicle_1",
                                                          "Mulder_4sU_KC.Vehicle_2","Mulder_4sU_KC.Vehicle_3")

# saveRDS(progeny_results_with_replicates[seq(6,15),], file = "4_processed_data/wholeset/progeny/20240902_progeny_output_only_replicates.RDS")


# create figure to check the sd of replicates
progeny_reps <- scale(progeny_results_with_replicates[seq(6,15),])
data <- as.data.frame(progeny_reps)
data <- rownames_to_column(data, var = "condition")
data_long <- pivot_longer(data, cols = all_of(colnames(data)[-1]))
data_long <- data_long %>% mutate(group = stringr::str_sub(condition, end = -3L))

sd_data_long <- data_long %>% group_by(group,name) %>% 
  summarise(sd = sd(value, na.rm=T))


ggplot(data_long, aes(x = name, y = value)) +
  geom_point() +
  facet_wrap(~group) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
    element_text(size = 2)
  ) + 
  labs(x = "", y = "Pathway activity scores")
ggsave(filename = "6_results/QC/20240902_progeny_replicates_check.png",dpi = 600, height = 5, width = 7)


# create heatmap

progeny_results_scaled <- scale(progeny_results_unscaled)
rownames(progeny_results_scaled)
rownames(progeny_results_scaled)[5] <- "Mortazavi_4sU_LCL"

paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

png(filename = "6_results/QC/20240902_heatmap_progeny_output_mean_per_condition.png",
    width = 13,
    height = 15,
    units = "cm",
    res = 500)

progeny_hmap <- pheatmap(t(progeny_results_scaled),fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor,, angle_col = 45,
                         treeheight_col = 0,  border_color = NA)
dev.off()



