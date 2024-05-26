library(tidyverse)
library(data.table)
source("utils.R")
library(RColorBrewer)

# where the CBSS results are saved (after running the `stats/sum_data.R` file)
df_dir <- "/home/yutong/Dropbox/Projects/CBSS/results/" 
# save the analysis results
save_dir <- "/home/yutong/Dropbox/Projects/CBSS/reports/figures/"
# demographic features
demo_path <- "/home/yutong/Data/CBSS/demographics/ZhejiangAD_metrics.csv"

# cognitive features to be studied and feature order
cogn <- c("Age", "Sex", "MMSE", "MOCA", "STM", "LTM", "TMTA")
feature_order <- c("Age", "MMSE", "MOCA", "STM", "LTM", "TMTA")
# the following 2 regions may not be captured in IVIM and should be removed
remove_reg <- c("premedullary cistern", "cerebellomedullary angle cistern")

# ====================cognition correlation====================
# please adjust the codes according to the demographic attributes in your data
demo_df <- fread(demo_path) %>%
    mutate(WMH_norm = log(total_WMH/ICV + 1)) %>%
    rename(STM = ShortMem) %>%
    rename(LTM = LateMem) %>%
    mutate(TBV = GM + WM) %>%
    mutate(TBV_frac = (GM + WM)/(GM + WM + vent + CSF))

diff_frame <- fread(paste(df_dir, "all_diff.csv", sep="/"))
# average left and right-sided pseudodiffusivity weighted by sulcal volume
diff_frame <- get_sided(diff_frame, vol_weight=T) %>%
    mutate(labID = ifelse(labID == "vent", "sulcus", labID)) %>%
    filter(!grepl("inferior", Region))
all_df <- diff_frame %>% merge(demo_df, by="ID")

all_cor <- find_cog_corr(all_df, cogn, "median_val",
                         covars=c("Age", "Sex", "edu"))

# plotting
all_cor %>%
    filter(labID == "sulcus") %>%  # only analyse the sulcus
    mutate(Region = factor(Region, levels=rev(get_sulcus_order()))) %>%
    mutate(feature = factor(feature, levels=cogn)) %>%
    ggplot(aes(x=feature, y=Region, fill=corr)) +
    geom_tile() +
    geom_text(aes(label=gtools::stars.pval(p.adj))) +
    scale_fill_distiller(palette="RdYlBu") +
    theme_publication() +
    xlab("") + ylab("") + labs(fill="Slope") +
    extra_theme()
ggsave(paste(save_dir, "/cogn_corr.jpg", sep="/"), width=7, height=9, dpi=500)

# ----------Average level in each CSF region----------
diff_frame %>% 
    filter(labID == "sulcus") %>%
    group_by(Region) %>%
    summarize(mean_val = mean(median_val)) %>%
    fwrite(paste(save_dir, "/sum_stats.csv", sep="/"))
# Now you can run `stats/show_map.py` to display average values in each sulcus

diff_frame %>%
    filter(labID == "sulcus") %>%
    filter(!Region %in% remove_reg) %>%
    mutate(Region = factor(Region, levels=rev(get_sulcus_order()))) %>%
    ggplot(aes(x=median_val, y=Region)) +
    geom_boxplot() +
    xlab("Median pseudodiffusivity") + ylab("") +
    theme_publication()
ggsave(paste(save_dir, "/median_val_box.jpg", sep="/"), width=9, height=7)

# ----------correlation between sulci and cisterns----------
cor_df <- diff_frame %>% 
    filter(labID == "sulcus") %>%
    ungroup() %>%
    select(ID, Region, median_val)

cor_plot_df <- cor_heat(cor_df)
cor_plot_df %>% 
    mutate(pval = gtools::stars.pval(padj)) %>%
    ggplot(aes(x=regionA, y=regionB, fill=corr)) +
    geom_tile() +
    geom_text(aes(label=pval)) +
    scale_fill_distiller(palette="RdYlBu") +
    theme_publication() +
    xlab("") + ylab("") + labs(fill="Correlation") +
    extra_theme()
ggsave(paste(save_dir, "/sulcus_corr.jpg", sep="/"), width=11, height=9, dpi=500)

# ----------Ventricular zones fine mapping----------
plot_df <- all_cor %>% 
    filter(labID == "vent_zone") %>%
    filter(feature != "Sex") %>%
    filter(!grepl("fourth", Region)) %>%
    mutate(psym = gtools::stars.pval(p.adj)) %>%
    assign_dist() %>%
    mutate(feature = factor(feature, levels=feature_order)) %>%
    mutate(label_pos = corr + 0.03*sign(corr))

max_y <- max(plot_df$label_pos)
ggplot(plot_df, aes(x=Region, y=corr, color=feature)) +
    geom_line() +
    geom_text(aes(label=psym, y=label_pos), size=4, show.legend=F) +
    labs(color="") +
    xlab("Distance from the start of the third ventricle (mm)") +
    ylab("Linear regression slope") +
    geom_vline(xintercept=0, linetype="dashed") +
    annotate("text", x=10, y=max_y, label="Third ventricle") +
    annotate("text", x=-10, y=max_y, label="Lateral ventricle") +
    theme_publication()
ggsave(paste(save_dir, "/Monro.jpg", sep="/"), width=7, height=7, dpi=500)

# ====================association with perfusion====================
perf_frame <- fread(paste(df_dir, "all_perf.csv", sep="/"))
perf_frame <- get_sided(perf_frame, vol_weight=T) %>%
    filter(labID == "CBF") %>%
    mutate(Region = gsub("fasciculus", "fissure", Region)) %>%
    filter(!grepl("inferior", Region))

perf_val <- perf_frame %>%
    mutate(regionID = paste(ID, Region)) %>%
    rename(CBF = median_val) %>%
    select(regionID, CBF)

diff_val <- all_df %>%
    mutate(regionID = paste(ID, Region)) %>%
    filter(labID == "sulcus") %>%
    rename(Diff = median_val) %>%
    select(regionID, Region, ID, Diff, total_vol, Age, Sex)

all_val <- diff_val %>% merge(perf_val, by="regionID") %>%
    filter(!Region.x %in% remove_reg)
all_regions <- unique(all_val$Region.x)
all_res <- as.list(all_regions) %>% lapply(function(xx) {
    one_val <- all_val %>% filter(Region.x == xx) %>%
        mutate_if(is.numeric, scale2)
    mod <- summary(lm(Diff ~ CBF + total_vol + Age + Sex, data=one_val))
    data.frame(mod$coefficients) %>%
        rownames_to_column("Variable") %>%
        rename(p = Pr...t..) %>%
        select(Variable, Estimate, p) %>%
        filter(Variable == "CBF") %>%
        mutate(Region = xx) %>%
        return()
})

all_res <- do.call(rbind, all_res) %>%
    mutate(p.adj = p.adjust(p, method="fdr")) %>%
    mutate(feature = gsub("CBF", "Regional CBF", Variable)) %>%
    rename(corr = Estimate)

imaging <- c("CBF", "ALPS", "PVS_CSO", "PVS_BG")
imaging_order <- c("DTI-ALPS", "Cortical CBF", "Regional CBF", "PVS_CSO", "PVS_BG")
all_cor <- find_cog_corr(all_df, imaging, "median_val",
                         covars=c("Age", "Sex"))

# plotting
plyr::rbind.fill(all_cor %>% filter(labID == "sulcus"), all_res) %>%
    mutate(feature = gsub("^CBF$", "Cortical CBF", feature)) %>%
    mutate(feature = gsub("^ALPS$", "DTI-ALPS", feature)) %>%
    mutate(feature = factor(feature, levels=imaging_order)) %>%
    mutate(Region = factor(Region, levels=rev(get_sulcus_order()))) %>%
    ggplot(aes(x=feature, y=Region, fill=corr)) +
    geom_tile() +
    geom_text(aes(label=gtools::stars.pval(p.adj))) +
    scale_fill_distiller(palette="RdYlBu") +
    theme_publication() +
    xlab("") + ylab("") + labs(fill="Slope") +
    extra_theme()
ggsave(paste(save_dir, "/imaging_corr.jpg", sep="/"), width=6, height=9, dpi=500)
