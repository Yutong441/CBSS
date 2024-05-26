# summarize image processing results
library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
# the ${data_dir}/processed, i.e., where is the processed data stored
root <- paste(args[1], "/processed", sep="/")
df_dir <- args[2] # where to save the output files

# ====================key function====================
collect_csf <- function(root, filename="sum_stat_CSF.csv", labID="CSF",
                        convert_root="../data/", folder="stats") {
    all_files <- list.files(root)
    csf_list <- as.list(all_files) %>% lapply(function(x) {
        file_path <- paste(root, "/", x, "/CBSS/", folder, "/", filename, sep="")
        if (file.exists(file_path)) {
            if (dim(one_df)[2] > 1) {
                one_df %>% filter(vol > 10) %>%
                    mutate(ID = x) %>%
                    rename(Region = V1) %>% 
                    filter(Region != 3000) %>%
                    return()
            }
        }
    })
    csf_df <- do.call(rbind, csf_list)

    if (labID == "sulcus") {
        conversion <- fread(paste(convert_root, "/CBSS_atlas.csv", sep="/")) %>%
            mutate(label = gsub("C", "R", label))
    } else {
        conversion <- fread(paste(convert_root, "/synthseg_atlas.csv", sep="/")) %>%
            mutate(label = paste("R", label, sep=""))
    }
    if (labID != "vent_zone") {
        csf_df$Region <- paste("R", csf_df$Region, sep="")
        csf_df$Region <- conversion$atlas[match(csf_df$Region, conversion$label)]
    }
    csf_df$labID <- labID
    return(csf_df)
}


# ====================CBSS data====================
all_files <- c("sum_stat_CSF.csv", "sum_stat_vent.csv",
               "sum_stat_sulcus.csv", "sum_vent_region.csv")
labIDs <- c("CSF", "vent", "sulcus", "vent_zone")
result_list <- as.list(1:length(all_files)) %>% lapply(function(xx) {
    collect_csf(root, filename=all_files[xx], labID=labIDs[xx])   
})

do.call(rbind, result_list) %>%
    select(!one_of(c("max_val", "min_val", "sd"))) %>%
    fwrite(paste(df_dir, "/all_diff.csv", sep="/"))

# ====================perfusion statistics====================
# ----------regional perfusion
all_files <- c("stats_sulcus_CBF.csv", "stats_sulcus_CVR.csv",
               "stats_sulcus_CVRlag.csv",
               "stats_sulcus_Fivim.csv", "stats_sulcus_delay.csv",
               "stats_sulcus_NVC.csv", "stats_sulcus_ALFF.csv"
               )
labIDs <- c("CBF", "CVR", "CVRlag", "F", "Delay", "NVC", "ALFF")
result_list <- as.list(1:length(all_files)) %>% lapply(function(xx) {
    one_df <- collect_csf(root, filename=all_files[xx], labID="sulcus",
                         folder="blood_stats")
    one_df$labID <- labIDs[xx]
    return(one_df)
})

all_res <- do.call(rbind, result_list) %>%
    select(!one_of(c("max_val", "min_val", "sd"))) %>%
    fwrite(paste(df_dir, "/all_perf.csv", sep="/"))

# ----------global cortical perfusion
all_files <- c("stats_gm_CBF.csv", "stats_gm_CVR.csv",
               "stats_gm_CVRlag.csv",
               "stats_gm_Fivim.csv", "stats_gm_delay.csv",
               "stats_gm_ALFF.csv"
)
labIDs <- c("CBF", "CVR", "CVRlag", "F", "Delay", "ALFF")
result_list <- as.list(1:length(all_files)) %>% lapply(function(xx) {
    all_IDs <- list.files(root)
    one_type <- as.list(all_IDs) %>% lapply(function(yy) {
        filename <- paste(root, "/", yy, "/CBSS/blood_stats/", all_files[xx], sep="")
        if (file.exists(filename)) {
            one_df <- fread(filename)
            if (dim(one_df)[2] > 1) {
                val <- one_df$median_val
                return(data.frame(ID=yy, median_val=val))
            }
        }
    })
    one_type <- do.call(rbind, one_type)
    one_type$labID <- labIDs[xx]
    return(one_type)
})

do.call(rbind, result_list) %>%
    spread(labID, median_val) %>%
    rename(Fivim = `F`) %>%
    fwrite(paste(df_dir, "/global_perf.csv", sep="/"))
