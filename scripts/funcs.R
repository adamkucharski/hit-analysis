
# 1. LOAD LIBRARIES AND DATA
## 1.1 load libraries
packs <- c("tidyverse", "readxl", "patchwork", "metR", "shades", "ggdist", "ggrepel", "countrycode", "lubridate")
lapply(packs, library, character.only = TRUE)

## 1.2 load mean and ci data for R0 (measles, flu, sars-cov-2) and vac eff from excel
data_xl <- read_excel("data/data_full.xlsx") %>%
    mutate(eff_mean = eff_mean / 100, eff_lb = eff_lb / 100, eff_ub = eff_ub / 100)

## 1.3 samples list R0 for bootstrapping
samples_bs_r0 <- list(
    mumps = c(4.5, 4.3, 3.6, 4, 4.2, 4.4),
    rubella = c(3.7, 6.4, 3.4, 4.2, 7.8, 4.2, 3.7),
    varicella = c(6.47, 3.83, 4.85, 5.46, 5.22, 7.71, 3.31, 8.26, 16.91, 5.72, 3.91)
)

## 1.4 Calculate the HIT values
hit_df <- data.frame(r0 = seq(1, 21, 0.1), threshold = 1 - 1 / seq(1, 21, 0.1))

# 2. HELPLER FITTING FUNCTIONS
# 2.1 Function to fit to a distributions using mean and interval
fitdist_ci <- function(pars, data, dist) {
    qs <- c(0.025, 0.5, 0.975); 1:3 %>% map_dbl(~ (dist(qs[.x],  pars[1], pars[2]) - data[.x])^2) %>% sum
}
fitdist_iqr <- function(pars, data, dist) {
    qs <- c(0.25, 0.5, 0.75); 1:3 %>% map_dbl(~ (dist(qs[.x],  pars[1], pars[2]) - data[.x])^2) %>% sum
}

# 3. FUNCTION TO GET SAMPLES AND DATA
## 3.1 Get samples by fitting mean and ci to a distribution.
get_samples_fit <- function(data, dist, param, init = c(1, 1)) {
    dist_pdf <- str2lang(paste0("q", dist))
    dist_rng <- str2lang(paste0("r", dist))
    samples_list <- vector(mode = "list", length = nrow(data))
    for (i in seq_len(nrow(data))) {
        sample <- data[i, ]
        mean_ci <- c(sample[[paste0(param, "_lb")]], sample[[paste0(param, "_mean")]], sample[[paste0(param, "_ub")]])
        par_fit <- optim(init, fitdist_ci, data = mean_ci, dist = eval(dist_pdf))$par
        samples_list[[i]] <- data.frame(pathogen = sample$pathogen,
            metric = param,
            samples = eval(dist_rng)(1000, par_fit[1], par_fit[2]))
    }
    samples_list
}

## 3.2 Get samples by bootstrapping some samples
get_samples_bs <- function(data, samples, param) {
    samples_list <- vector(mode = "list", length = nrow(data))
    for (pathogen in names(samples)) {
        samples_list[[pathogen]] <- data.frame(pathogen = pathogen,
            metric = param,
            samples = map_dbl(1:1000, ~mean(sample(samples[[pathogen]], size = 20, replace = TRUE))))
    }
    samples_list
}

## 3.3 Get samples of HIT given some background immunity, vaccine effectiveness and r0
get_samples_nat_imm <- function(data, pathogen_str, vac_eff) {
    # trim values for plotting
    trim_val <- function(val, oob_l, oob_u) {
        if (val < 0) { return(oob_l)}
        else if(val > 1) { return(oob_u)}
        else { return(val)}
    }

    r0_samples <- data %>% filter(pathogen == pathogen_str) %>% pull(r0)
    xx <- seq(0.01, 0.99, 0.01); j <- 1;
    hit_list <- vector(mode = "list", length = length(xx))
    for (i in seq_len(length(xx))) {
        for (cov in vac_eff) {
            samples <- ((1 - 1 / r0_samples) - xx[i]) / ((1 - xx[i]) * cov)
            hit_list[[j]] <- data.frame(
                x = i / 100,
                y = median(samples) %>% trim_val(NA, NA),
                ymin = quantile(samples, c(0.025)) %>% trim_val(0, 1),
                ymax = quantile(samples, c(0.975)) %>% trim_val(NA, 1),
                coverage = cov
            )
            j <- j + 1
        }
    }
    hit_df <- bind_rows(hit_list)
    hit_df$coverage <- factor(hit_df$coverage, levels = vac_eff)
    attr(hit_df, "variant") <- pathogen_str
    hit_df
}

## 3.4 Get the proportion of persons aged 15 years and older in each country
get_prop_15p_country <- function() {
    # http://www.healthdata.org/sites/default/files/files/Projects/GBD/GBDRegions_countries.pdf
    # https://data.worldbank.org/indicator/SP.POP.TOTL?name_desc=false
    # https://data.worldbank.org/indicator/SP.POP.0014.TO.ZS?name_desc=false
    country_info <- read.csv("data/country_info.csv")
    prop_child <- read.csv("data/prop_child_country.csv")
    pop <- read.csv("data/pop.csv")

    prop_child_trim <- prop_child %>% select(country = Country.Name, iso3c = Country.Code, prop = X2019) %>% na.omit
    country_info_trim <- country_info %>% select(country = TableName, iso3c = Country.Code, who_reg = Region, gbd_reigon = GBD_super, income = IncomeGroup) 
    pop_trim <- pop %>% select(iso3c = Country.Code, pop = X2019)

    combine_df <- left_join(prop_child_trim, country_info_trim)
    combine_df[combine_df == ""] = NA
    combine_df_trim <- combine_df[complete.cases(combine_df), ]
    df_full <- left_join(combine_df_trim, pop_trim)
}


# 3.5 Function to get seroprevalence data from serotracker
get_sero_prev <- function() {
    sero_prev_info <- read.csv("data/serotracker_data.csv")
    # rename and keep relevant columns
    sero_prev_info_trim <- sero_prev_info %>% select(dates = Study.Dates, scope = Grade.of.Estimate.Scope, country = Country, 
        risk = Overall.Risk.of.Bias..JBI., 
        study_pop = Sample.Frame..groups.of.interest., pop = Denominator.Value, prev = Seroprevalence..95..CI.) %>%
        mutate(iso3c = countrycode(country, origin = 'country.name', destination = 'iso3c'))
    # sort out dates
    sero_prev_info_dates <- sero_prev_info_trim %>%
        separate(dates, c("start_date", "end_date"), sep = " - ") %>%
        mutate(start_date = case_when((as.numeric(substr(start_date, 1, 2)) < 3) ~ paste0(start_date, "/21"),
        (as.numeric(substr(start_date, 1, 2)) >= 3) ~ paste0(start_date, "/20")), 
        end_date = case_when((as.numeric(substr(end_date, 1, 2)) < 3) ~ paste0(end_date, "/21"),
        (as.numeric(substr(end_date, 1, 2)) >= 3) ~ paste0(end_date, "/20")), 
        ) %>% mutate(start_date = mdy(start_date), end_date = mdy(end_date) )
    # sep prev estimates
    sero_prev_info_dates <- sero_prev_info_dates %>% separate(prev, c("prev", "uncert"), sep = "% ") %>%
        mutate(prev = as.numeric(prev))

    sero_prev_global <- sero_prev_info_dates %>% 
        filter(study_pop %in% c("Household and community samples", "Blood donors", "Tissue donor", "Multiple populations"),
        scope %in% c("National", "Regional", "Local"), 
        end_date > ymd("2020-03-01")) 

    sero_prev_global$study_pop <- factor(sero_prev_global$study_pop, 
        levels = c("Household and community samples", "Blood donors",  "Tissue donor", "Multiple populations"))
    sero_prev_global$scope <- factor(sero_prev_global$scope, 
        levels = c("National", "Regional", "Local"))

    sero_prev_global[complete.cases(sero_prev_global), ]
}

# 3.6. Function to combine the seroprevalence data and the proportion over 15
get_sero_prev_prop <- function() { 
    latest_survey <- function(prev, ss, date) {
        latest_date <-  max(date)
        prev[which.max(date)]
        tibble(prev = prev[which.max(date)], ss = ss[which.max(date)], date = latest_date)
    }

    sero_prev_global <- get_sero_prev() %>% select(end_date, iso3c, pop, scope, prev)
    sero_prev_summ <- sero_prev_global %>% group_by(iso3c, scope) %>% summarise(latest_survey(prev, pop, end_date))
    prop_15p <- get_prop_15p_country()
    combined_country_plot <- left_join(prop_15p, sero_prev_summ, by = "iso3c") %>% na.omit %>%
        mutate(prop = (100 - prop)/100, prev = prev/100, 
            income =  factor(income,  levels = c("High income", "Upper middle income", "Lower middle income", "Low income")))
}

# 4. PLOTTING FUNCTION
## 4.1 Plot Figure 1a
plot_fig1a <- function(samples_full, samples_median) {
    sc2_cols <- c("#046927", "#4CBB17")
    flu_cols <- c("blue") %>% brightness(0.8) %>% saturation(seq(0.2, 1, 0.4))
    pathogen_colours <- c("red", "purple", "orange", "gray", sc2_cols, flu_cols)

    samples_full %>%
        ggplot(aes(x = r0, y = eff)) +
            geom_bin2d(aes(fill = pathogen), alpha = 0.7, bins = 250) +
            geom_point(data = samples_median, aes(fill = pathogen),
                shape = 21, size = 3, colour = "black") +
            scale_fill_manual(values = pathogen_colours) +
            theme_bw() + theme(aspect.ratio = 0.7, legend.position = "top") +
            xlim(0, 20) +
            geom_line(data = hit_df, aes(x = r0, y = threshold)) +
            geom_label_repel(data = samples_median, aes(x = r0, y = eff, label = pathogen), size = 1.4, alpha = 0.85) +
            annotate("label", x = 11, y = 0.65,  alpha = 0.8,
                label = "Herd immunity threshold", size = 2.5, angle = 0) +
            annotate("segment", x = 8.5, xend = 7.0, y = 0.7, yend = 0.85, colour = "black", size =  0.5, alpha = 0.8,
                arrow = arrow(length = unit(4, "mm"))) +
            labs(x = expression(R["0"]), y = "Vaccine effectiveness", fill = "Pathogen", title = "A)") +
            guides(fill = guide_legend(override.aes = list(size = 4))) +
                theme(legend.title = element_text(size = 8),
                    legend.text  = element_text(size = 6),
                    legend.key.size = unit(1, "lines"),
                    aspect.ratio = 0.5)
}

## 4.2 Plot Figure 1b/c
plot_fig1bc <- function(pltdata, prop_15p, let) {
    # stuff to truncate the values between the 0 to 1 range
    pltdata %>%
        ggplot(aes(x = x, y = y)) +
            geom_hline(data = prop_15p, aes(yintercept = prop), size = 0.3) + 
            geom_text_repel(data = prop_15p, aes(x = 0.9, y = prop, label = income), size = 2, alpha = 1, 
                min.segment.length = 0) + 
            geom_line(aes(color = coverage), size = 2) +
            geom_ribbon(aes(ymin = ymin, ymax = ymax, color = coverage), alpha = 0.05, size = 0.2) +
            scale_fill_continuous(type = "viridis") +
            theme_bw() + theme(aspect.ratio = 0.5, legend.position = "top",
                legend.title = element_text(size = 12),
                legend.text  = element_text(size = 10),
                legend.key.size = unit(1, "lines")) +
            labs(x = "Proportional reduction in transmission as a result of\n
                    naturally acquired immunity to SARS-CoV-2",
                y = "Vaccination coverage", colour = "Coverage") +
            coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
            ggtitle(paste0(let,") " ,attributes(pltdata)$variant))
}

## 4.3 Plot Figure 2
plot_fig2 <- function(data, thresholds, title) {
    add_sat <- function(data) {
        stdscale <- function(x){(x - (min(x - 80)))/(max(x) - (min(x - 80)))}
        data %>% mutate(date_num  = stdscale(as.numeric(date))) %>% 
            mutate(
            plt_col = case_when(income == "High income"~saturation(brightness("red", 0.9), date_num),
                income == "Upper middle income"~saturation(brightness("orange", 0.9), date_num),
                income == "Lower middle income"~saturation(brightness("green", 0.9), date_num),
                income == "Low income"~saturation(brightness("blue", 0.9), date_num)
                ) ) 
    }

    add_lab <- function(data) {
        data %>% 
            mutate(label = countrycode(data$iso3c, origin = "iso3c", destination = "country.name")) %>%
            mutate(label_inc = case_when(
                ((scope == "National") & (date > as.Date("2020-08-20")))~TRUE,
                ((scope == "National") & (income == "Upper middle income"))~TRUE,
                (prev > 0.4)~TRUE,
                (prop < 0.65)~TRUE,
                TRUE~FALSE
                ))
    }
    
    data <- data %>% add_sat %>% add_lab
    data_labs <- data %>% filter(label_inc == TRUE)

    data %>% 
        ggplot() + 
        geom_line(data = sc2_wild, aes(x = x, y = y, linetype = coverage), size = 1) + 
        geom_point(aes(x = prev, y = prop, shape = scope, size = ss, fill = income), alpha = 0.7) + 
        geom_point(aes(x = prev, y = prop, shape = scope, size = ss), fill = data$plt_col, alpha = 0.7) + 
        geom_text_repel(data = data_labs, aes(x = prev, y = prop, label = country), color = "black", alpha = 0.8, size = 3, 
            min.segment.length = 0, label.padding = 0.15, box.padding = 0.2, , max.overlaps = Inf,
            segment.curvature = -0.1,
            segment.ncp = 3,
            segment.angle = 20,
            xlim  = c(0.15, NA),
            bg.color = "white", bg.r = 0.15) +
        scale_shape_manual(values = c(21, 22, 23, 24)) + 
        scale_linetype_manual(values = c("solid", "dashed", "dotted")) + 
        theme_bw() + 
        theme(aspect.ratio = 0.8) + 
        labs(x = "Estimated seroprevalence", y = "Proportion of population 15 years and over", fill = "Income", 
            shape = "Study type", size = "Study size", linetype = "HIT at different vaccine\n effectiveness values", 
            title = title) +
        coord_cartesian(xlim = c(0, 0.72), ylim = c(0.55, 0.9)) + 
        scale_size(guide = "none") + 
        guides(linetype = guide_legend(order = 1),
            fill = guide_legend(order = 2, override.aes = 
                list(color = c("red", "orange", "green", "blue"), size = 4)),
            shape = guide_legend(order = 3, override.aes = list(size = 4)))
}
