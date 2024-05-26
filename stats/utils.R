scale2 <- function(xx) {(xx - mean(xx, na.rm=T))/sd(xx, na.rm=T) %>% return()}

get_sided <- function(one_df, vol_weight=FALSE) {
    if (vol_weight) {
        one_df %>%
            dplyr::mutate(Side = stringr::str_extract(Region, "(left|right)")) %>%
            dplyr::mutate(Region = gsub("(left|right) ", "", Region)) %>%
            dplyr::group_by(Region, labID, ID) %>%
            dplyr::summarize(mean_val = mean(mean_val*vol),
                      median_val = mean(median_val*vol),
                      total_vol = sum(vol)) %>%
            dplyr::mutate(mean_val = mean_val/total_vol) %>%
            dplyr::mutate(median_val = median_val/total_vol) %>%
            return()
    } else {
        one_df %>%
            dplyr::mutate(Side = stringr::str_extract(Region, "(left|right)")) %>%
            dplyr::mutate(Region = gsub("(left|right) ", "", Region)) %>%
            dplyr::group_by(Region, labID, ID) %>%
            dplyr::summarize(mean_val = mean(mean_val), median_val = mean(median_val)) %>%
            return()
    }
}


#' Correlation between CSF diffusivity and continuous clinical variables
#' 
#' @description use linear regression
#' @param all_df dataframe obtained from `stats/sum_data.R`. This dataframe has
#' the following columns: "Region", "ID", "labID", "mean_val", "median_val".
#' @param features any imaging or cognitive features to correlate with
#' @param analyze_thres Only analyze a particular sulcus if there are at least a
#' certain number of patients
#' @param dif_col can be "median_val" or "mean_val"
find_cog_corr <- function(all_df, features, dif_col="median_val",
                          analyze_thres=50, covars=c("Age", "Sex"),
                          method="lm") {
    all_reg <- unique(all_df$Region)
    all_fea <- as.list(features) %>% lapply(function(y) {
        all_cor <- as.list(all_reg) %>% lapply(function(x) {
            one_df <- all_df %>% filter(Region == x)
            N <- sum(!is.na(one_df %>% pull(one_of(dif_col))))
            if (N > analyze_thres) {
                if (y != "Sex") {
                    formu <- paste(y, "~", dif_col)
                    if (!y %in% c("Age", "Sex")) {
                        covars_formu <- paste(covars, collapse="+")
                        formu <- paste(formu, "+", covars_formu)
                    } else {
                        sel_cov <- covars[!covars %in% c("Age", "Sex")]
                        if (length(sel_cov) > 0) {
                            covars_formu <- paste(sel_cov, collapse="+")
                            formu <- paste(formu, "+", covars_formu)
                        }
                    }
                    one_df <- one_df %>% mutate_if(is.numeric, scale2)
                    if (method == "lm") {
                        mod <- summary(lm(as.formula(formu), data=one_df))
                        stat_val <- mod$coefficients[dif_col, "Estimate"]
                        p_val <- mod$coefficients[dif_col, "Pr(>|t|)"]
                        pearson <- cor(one_df[, dif_col], one_df[, y], use="na.or.complete")
                    } else {

                    }
                } else {
                    formu <- as.formula(paste(dif_col, "~Sex"))
                    mod <- wilcox.test(formu, data=one_df)
                    stat_val <- mod$statistic
                    p_val <- mod$p.value
                    pearson <- NA
                }
                return(data.frame(Region=x, p=p_val, corr=stat_val,
                                  labID=one_df$labID[1], pearson=pearson))
                } 
        })
        all_cor <- do.call(rbind, all_cor)
        # `match` is essential in combining dataframes later
        all_cor[match(all_reg, all_cor$Region), ] %>%
            mutate(Region = all_reg) %>%
            mutate(feature = y) %>%
            group_by(labID) %>%
            mutate(p.adj = p.adjust(p, method="fdr")) %>%
            ungroup() %>% return()
    })
    all_fea <- do.call(rbind, all_fea)
        # mutate(Region = factor(Region, levels=get_sulcus_order()))
    return(all_fea)
}

#' Correlation between CSF diffusivity and nominal clinical variables
#' 
#' @description use linear regression to remove the influence of covariates,
#' then use Wilcoxon to do comparison
#' @param all_df dataframe obtained from `stats/sum_data.R`. This dataframe has
#' the following columns: "Region", "ID", "labID", "mean_val", "median_val".
#' @param features any imaging or cognitive features to correlate with
#' @param analyze_thres Only analyze a particular sulcus if there are at least a
#' certain number of patients
#' @param dif_col can be "median_val" or "mean_val"
find_cog_corr_nominal <- function(all_df, features, dif_col="median_val",
                                  analyze_thres=50, covars=c("Age", "Sex"),
                                  method="lm") {
    all_reg <- unique(all_df$Region)
    all_fea <- as.list(features) %>% lapply(function(y) {
        all_cor <- as.list(all_reg) %>% lapply(function(x) {
            one_df <- all_df %>% filter(Region == x)
            N <- sum(!is.na(one_df %>% pull(one_of(dif_col))))
            if (N > analyze_thres) {
                formu <- paste(dif_col, "~")
                covars_formu <- paste(covars, collapse="+")
                formu <- paste(formu, "+", covars_formu)
                one_df <- one_df %>% mutate_if(is.numeric, scale2)
                # linear regression to regress out covariates
                model <- lm(as.formula(formu), data=one_df)
                pred <- predict(model, one_df %>% dplyr::select(dplyr::all_of(covars)))
                one_df$residuals <- one_df %>% pull(dplyr::one_of(dif_col)) - pred

                formu <- as.formula(paste("residuals ~", y))
                mod <- wilcox.test(formu, data=one_df)
                stat_val <- mod$statistic
                p_val <- mod$p.value
                return(data.frame(Region=x, p=p_val, corr=stat_val,
                                  labID=one_df$labID[1]))
                } 
        })
        all_cor <- do.call(rbind, all_cor)
        # `match` is essential in combining dataframes later
        all_cor[match(all_reg, all_cor$Region), ] %>%
            mutate(Region = all_reg) %>%
            mutate(feature = y) %>%
            group_by(labID) %>%
            mutate(p.adj = p.adjust(p, method="fdr")) %>%
            ungroup() %>% return()
    })
    all_fea <- do.call(rbind, all_fea)
        # mutate(Region = factor(Region, levels=get_sulcus_order()))
    return(all_fea)
}

change_zone_name <- function(all_df) {
    all_df %>%
        mutate(Region = gsub("Monro_zone1", "lateral_zone4", Region)) %>%
        mutate(Region = gsub("Monro_zone2", "lateral_zone3", Region)) %>%
        mutate(Region = gsub("Monro_zone3", "lateral_zone2", Region)) %>%
        mutate(Region = gsub("Monro_zone4", "lateral_zone1", Region)) %>%
        mutate(Region = gsub("_", " ventricle ", Region)) %>%
        return()
}


assign_dist <- function(all_df) {
    dist_vec <- c(Monro_zone1=-2.5, Monro_zone2=-7.5,
                  Monro_zone3=-12.5, Monro_zone4=-17.5,
                  third_zone1=2.5, third_zone2=7.5,
                  third_zone3=12.5, third_zone4=17.5)
    all_df$Region <- dist_vec[match(all_df$Region, names(dist_vec))]
    return(all_df)
}

get_sulcus_order <- function() {
    return(c("premedullary cistern", "cerebellomedullary angle cistern",
             "prepontine cistern", "cerebellopontine angle cistern",
             "suprasellar cistern", "crural cistern",
             "ambient cistern", "quadrigeminal cistern",
             "Sylvian fissure", "superior temporal sulcus",
             "precentral sulcus", "central sulcus", 
             "postcentral sulcus", "intraparietal sulcus",
             "superior longitudinal fissure (anterior part)",
             "superior longitudinal fissure (posterior part)",
             "lateral ventricle", "third ventricle", "fourth ventricle"
    )
    )
}

theme_publication <- function(base_size=15, base_family="arial") {
      (ggthemes::theme_foundation(base_size=base_size, base_family=base_family)
       + ggplot2::theme(
               plot.title = ggplot2::element_text(
                    face = "bold", size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = ggplot2::element_rect(colour = NA),
               plot.background = ggplot2::element_rect(colour = NA),
               panel.border = ggplot2::element_rect(colour = NA),
               axis.title = ggplot2::element_text(face = "bold", size = rel(1)),
               axis.title.y = ggplot2::element_text(angle = 90, vjust = 2),
               axis.title.x = ggplot2::element_text(vjust = -0.2),
               axis.text = ggplot2::element_text(),
               axis.line = ggplot2::element_line(colour = "black"),
               axis.ticks = ggplot2::element_line(),
               panel.grid.major = ggplot2::element_line(colour = "#f0f0f0"),
               panel.grid.minor = ggplot2::element_blank(),
               legend.key = ggplot2::element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size = unit(0.2, "cm"),
               legend.margin = unit(0, "cm"),
               legend.title = ggplot2::element_text(face = "italic"),
               plot.margin = unit(c(10, 5, 5, 5), "mm"),
               strip.background = ggplot2::element_rect(colour="#f0f0f0",
                                                        fill="#f0f0f0"),
               strip.text = ggplot2::element_text(face = "bold")
          ))
}


extra_theme <- function() {
    theme(axis.text.x = element_text(angle=45, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right",
          legend.key.height = unit(1, 'cm'),
          legend.direction = "vertical"
    ) %>% return()
}


funct_regress <- function(all_df, demo_df, cogn,
                          x_name=NULL, x_range=c(-20, 10)) {
    sel_ID <- demo_df %>% 
        dplyr::filter(!is.na(!!as.symbol(cogn))) %>%
        dplyr::pull(ID)
    x_beta <- all_df %>% 
        dplyr::filter(ID %in% sel_ID) %>%
        dplyr::select(ID, Region, median_val) %>%
        tidyr::spread(Region, median_val) %>%
        tibble::column_to_rownames("ID") %>% 
        dplyr::select(!one_of("labID")) %>%
        t()

    y_vec <- demo_df %>% dplyr::select(dplyr::one_of("ID", cogn)) %>% 
        tibble::deframe() %>% as.numeric()
    names(y_vec) <- demo_df %>% dplyr::pull(ID)
    y_vec <- y_vec[colnames(x_beta)]

    if (is.null(x_name)) {
        x_name <- c(Monro_zone4=-0.5, Monro_zone3=-1.5, Monro_zone2=-2.5,
                    Monro_zone1=-3.5, third_zone1=0.5, third_zone2=1.5)
        x_range <- c(-3.5, 1.5)
    }

    smallbasis  <- fda::create.fourier.basis(x_range, length(x_name) - 1)
    xfd <- fda::smooth.basis(x_name, x_beta, smallbasis)$fd
    # y_xfd <- fda::fRegress(y_vec ~ xfd)
    # return(y_xfd)
    return(list(xfd, y_vec))
}


#' Correlation heatmap matrix
#'
cor_heat <- function(all_df, padj_method="fdr") {
    # remove columns with lots of NA
    sum_df <- all_df %>%
        mutate(NA_dif = ifelse(is.na(median_val), 0, 1)) %>%
        group_by(Region) %>%
        summarize(NA_region = sum(NA_dif))

    all_IDs <- all_df %>% pull(ID) %>% unique()

    keep_region <- sum_df %>% 
        mutate(NA_region = NA_region/length(all_IDs)) %>%
        filter(NA_region > 0.5) %>% pull(Region)

    sel_df <- all_df %>% 
        mutate(Region = factor(Region, levels=get_sulcus_order())) %>%
        filter(Region %in% keep_region)

    N <- length(keep_region)
    keep_region <- sel_df %>% filter(!duplicated(Region)) %>% 
        arrange(Region) %>% pull(Region) %>% as.character()

    all_xx <- as.list(1:N) %>% lapply(function(x) {
        all_yy <- as.list(1:N) %>% lapply(function(y) {
            if (x > y) {
                xx <- keep_region[x]
                yy <- keep_region[y]

                val_x <- sel_df %>% filter(Region == xx) %>%
                    select(ID, median_val) %>% deframe()
                val_y <- sel_df %>% filter(Region == yy) %>%
                    select(ID, median_val) %>% deframe()
                val_x <- val_x[all_IDs]
                val_y <- val_y[all_IDs]
                mod <- cor.test(val_x, val_y, use="na.or.complete")
                data.frame(corr=mod$estimate, p=mod$p.value,
                           regionA=xx, regionB=yy) %>% return()
            }
        })
        return(do.call(rbind, all_yy))
    })
    do.call(rbind, all_xx) %>%
        mutate(regionA = factor(regionA, levels=get_sulcus_order())) %>%
        mutate(regionB = factor(regionB, levels=get_sulcus_order())) %>%
        mutate(padj = p.adjust(p, method=padj_method)) %>%
        return()
}


#' Cross correlation between sulcal statistics
cross_cor <- function(df1, df2) {
    all_regions1 <- unique(df1$Region)
    all_regions2 <- unique(df2$Region)
    all_IDs <- intersect(unique(df1$ID), unique(df2$ID))

    all_res1 <- as.list(all_regions1) %>% lapply(function(xx) {
        all_res2 <- as.list(all_regions2) %>% lapply(function(yy) {
            one_val1 <- data.frame(df1) %>% filter(Region == xx) %>%
                select(ID, median_val) %>% deframe()
            one_val2 <- data.frame(df2) %>% filter(Region == yy) %>%
                select(ID, median_val) %>% deframe()
            mod <- cor.test(one_val1[all_IDs], one_val2[all_IDs],
                            use="na.or.complete")
            data.frame(corr=mod$estimate, p=mod$p.value,
                       regionA=xx, regionB=yy) %>% return()
        })
        return(do.call(rbind, all_res2))
    })
    do.call(rbind, all_res1) %>%
        mutate(regionA = factor(regionA, levels=get_sulcus_order())) %>%
        mutate(regionB = factor(regionB, levels=get_sulcus_order())) %>%
        mutate(padj = p.adjust(p, method="fdr")) %>%
        return()
}
