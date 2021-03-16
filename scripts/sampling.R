
# 1. LOAD LIBRARIES AND DATA
# load libraries
ren
packs <- c("tidyverse", "readxl", "patchwork", "metR", "shades", "ggdist", "ggrepel")
lapply(packs, library, character.only = TRUE)
# load mean and ci from excel
data_xl <- read_excel("data/data_full.xlsx")

# 2. USEFUL FITTING FUNCTIONS
# function to fit to a distributions using mean and ci
fitdist <- function(pars, data, dist) {
    qs <- c(0.025, 0.5, 0.975)
    1:3 %>% map_dbl(~ (dist(qs[.x],  pars[1], pars[2]) - data[.x])^2) %>% sum
}
fitbeta <- function(pars, data) {
    qs <- c(0.025, 0.5, 0.975)
    1:3 %>% map_dbl(~ (qbeta(qs[.x],  pars[1], pars[2]) - data[.x])^2) %>% sum
}
# function to fit to a normal distributions using mean and ci
fitlnorm <- function(pars, data) {
    qs <- c(0.025, 0.5, 0.975)
    1:3 %>% map_dbl(~ (qlnorm(qs[.x],  pars[1], pars[2]) - data[.x])^2) %>% sum
}
# function to fit to a log-normal distributions using mean and ci
fitnorm <- function(pars, data) {
    qs <- c(0.025, 0.5, 0.975)
    1:3 %>% map_dbl(~ (qnorm(qs[.x],  pars[1], pars[2]) - data[.x])^2) %>% sum
}

# 3. GET SAMPLES FOR EFFICACY
# all fitted to a beta distribution
vac_eff_samples_list <- vector(mode = "list", length = nrow(data_xl))
for (i in seq_len(nrow(data_xl))) {
    sample <- data_xl[i, ]
    mean_ci <- c(sample$eff_lb / 100, sample$eff_mean / 100, sample$eff_ub / 100)
    par_fit <- optim(c(10, 10), fitdist, data = mean_ci, dist = qbeta)$par
    vac_eff_samples_list[[i]] <- data.frame(disease = sample$pathogen,
        metric = "eff",
        samples = rbeta(1000, par_fit[1], par_fit[2]))
}

# 4. GET SAMPLES FOR R0
# influenza, measles are all log-normal
data_xl_flu <- data_xl %>% filter(str_detect(pathogen, "^flu|measles"))
r0_flu_samples_list <- vector(mode = "list", length = nrow(data_xl_flu))
for (i in seq_len(nrow(data_xl_flu))) {
    sample <- data_xl_flu[i, ]
    mean_ci <- c(sample$r0_lb, sample$r0_mean, sample$r0_ub)
    par_fit <- optim(c(1, 1), fitdist, data = mean_ci, dist = qlnorm)$par
    r0_flu_samples_list[[i]] <- data.frame(disease = sample$pathogen,
        metric = "r0",
        samples = rlnorm(1000, par_fit[1], par_fit[2]))
}
# sars are all normal
data_xl_sars <- data_xl %>% filter(str_detect(pathogen, "^sars"))
r0_sars_samples_list <- vector(mode = "list", length = nrow(data_xl_sars))
for (i in seq_len(nrow(data_xl_sars))) {
    sample <- data_xl_sars[i, ]
    mean_ci <- c(sample$r0_lb, sample$r0_mean, sample$r0_ub)
    par_fit <- optim(c(1, 1), fitdist, data = mean_ci, dist = qnorm)$par
    r0_sars_samples_list[[i]] <- data.frame(disease = sample$pathogen,
        metric = "r0",
        samples = rnorm(1000, par_fit[1], par_fit[2]))
}
# other pathogens via bootstrapping
samplelist <- list(
    mumps = c(4.5, 4.3, 3.6, 4, 4.2, 4.4),
    rubella = c(3.7, 6.4, 3.4, 4.2, 7.8, 4.2, 3.7),
    varicella = c(6.47, 3.83, 4.85, 5.46, 5.22, 7.71, 3.31, 8.26, 16.91, 5.72, 3.91)
)
data_xl_other <- data_xl %>% filter(!(pathogen %in% data_xl_flu$pathogen))
r0_other_samples_list <- vector(mode = "list", length = nrow(data_xl_other))
for (pathogen in names(samplelist)) {
    r0_other_samples_list[[pathogen]] <- data.frame(disease = pathogen,
        metric = "r0",
        samples = map_dbl(1:1000, ~mean(sample(samplelist[[pathogen]], size = 20, replace = TRUE))))
}
# combine all r0 values together
vac_eff_samples_df <- bind_rows(vac_eff_samples_list, r0_other_samples_list, r0_flu_samples_list,
        r0_sars_samples_list) %>%
    arrange(disease) %>%
    group_by(metric) %>%
    mutate(row = row_number()) %>%
    pivot_wider(names_from = metric, values_from = samples) %>%
    as.data.frame

# 5. PLOT FIGURE 1A
# get HIT
hit_df <- data.frame(r0 = seq(1, 21, 0.1),
    threshold = 1 - 1 / seq(1, 21, 0.1))
# plotting names
pathogen_levels <- c("measles", "mumps", "rubella", "varicella", "sars_cov_2_wt", "sars_cov_2_b117",
"flu_h1n1pmd09", "flu_h3n2", "flu_b")
pathogen_labels <- c("Measles", "Mumps", "Rubella", "Varicella-Zoster", "SARS-CoV-2 (pre-B.1.1.7 variants)",
    "SARS-CoV-2 (B.1.1.7)", "Flu A(H1N1)",  "Flu A(H3N2)", "Flu B")
vac_eff_samples_df$disease <- factor(vac_eff_samples_df$disease, levels = pathogen_levels, labels = pathogen_labels)
# plotting colours
sc2_cols <- c("#046927", "#4CBB17")
flu_cols <- c("blue") %>% brightness(0.8) %>% saturation(seq(0.2, 1, 0.4))
pathogen_colours <- c("red", "purple", "orange", "gray", sc2_cols, flu_cols)
# find median points
vac_eff_samples_mean <- vac_eff_samples_df %>% group_by(disease) %>% summarise(eff = median(eff), r0 = median(r0))

p1 <- vac_eff_samples_df %>%
    ggplot(aes(x = r0, y = eff)) +
    geom_bin2d(aes(fill = disease), alpha = 0.7, bins = 250) +
    geom_point(data = vac_eff_samples_mean, aes(fill = disease),
        shape = 21, size = 3, colour = "black") +
    scale_fill_manual(values = pathogen_colours) +
    theme_bw() + theme(aspect.ratio = 0.7, legend.position = "top") +
    xlim(0, 20) +
    geom_line(data = hit_df, aes(x = r0, y = threshold)) +
    geom_label_repel(data = vac_eff_samples_mean, aes(x = r0, y = eff, label = disease), size = 1.4, alpha = 0.85,
        check_overlap = T) +
    annotate("label", x = 11, y = 0.65,  alpha = 0.8,
        label = "Herd immunity threshold", size = 2.5, angle = 0) +
    annotate("segment", x = 8.5, xend = 7.0, y = 0.7, yend = 0.85, colour = "black", size =  0.5, alpha = 0.8,
        arrow = arrow(length = unit(4, "mm"))) +
    labs(x = expression(R["0"]), y = "Vaccine effectiveness", fill = "Pathogen") +
    guides(fill = guide_legend(override.aes = list(size = 4))) +
        theme(legend.title = element_text(size = 8),
              legend.text  = element_text(size = 6),
              legend.key.size = unit(1, "lines"),
              aspect.ratio = 0.5)
# find proportion above HIT threshold
prop_path_hi <- vac_eff_samples_df %>% group_by(disease) %>% summarise(prop = mean(eff > (1 - 1 / r0)))

# 6. PLOT FIGURE 1b/c
clean_data_fig2 <- function(pathogen_str) {
    sars_cov_2_r0 <- vac_eff_samples_df %>% filter(disease == pathogen_str) %>% pull(r0)
    xx <- seq(0.01, 0.99, 0.01)
    sars_cov_2_hit <- vector(mode = "list", length = length(xx))
    j <- 1
    for (i in seq_len(length(xx))) {
        for (cov in c(0.5, 0.7, 0.9)) {
            samples <- ((1 - 1 / sars_cov_2_r0) - xx[i]) / ((1 - xx[i]) * cov)
            sars_cov_2_hit[[j]] <- data.frame(
                x = i / 100,
                y = median(samples),
                ymin = quantile(samples, c(0.025)),
                ymax = quantile(samples, c(0.975)),
                coverage = cov
            )
            j <- j + 1
        }
    }
    sars_cov_2_hit_bind <- bind_rows(sars_cov_2_hit)
    sars_cov_2_hit_bind$coverage <- factor(sars_cov_2_hit_bind$coverage, levels = c(0.5, 0.7, 0.9))
    attr(sars_cov_2_hit_bind, "variant") <- pathogen_str
    sars_cov_2_hit_bind
}

sc2_wild <- clean_data_fig2("SARS-CoV-2 (pre-B.1.1.7 variants)")
sc2_b117 <- clean_data_fig2("SARS-CoV-2 (B.1.1.7)")

plot_fig2 <- function(pltdata) {
    # stuff to truncate the values between the 0 to 1 range
    pltdata$ymax[pltdata$ymax > 0.999] <- 1
    pltdata$y[pltdata$y > 0.999] <- NA
    pltdata$y[pltdata$y < 0.001] <- NA
    pltdata$ymin[pltdata$ymin < 0.001] <- 0
    pltdata <- pltdata[pltdata$ymax > 0.001, ]
    pltdata <- pltdata[pltdata$ymin < 0.999, ]
    pltdata %>%
        ggplot(aes(x = x, y = y)) +
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
            ggtitle(attributes(pltdata)$variant)
}

# 7. MAKE FIGURE AND SAVE
p2 <- plot_fig2(sc2_wild)
p3 <- plot_fig2(sc2_b117)
pb <- guide_area() / p2 / p3 + plot_layout(guides = "collect", heights = c(0.1, 3, 3), widths = 2)
p1 / pb + plot_layout(heights = c(1, 2.5))
ggsave(file = "outputs/fig1.png", width = 8, height = 12, dpi = 300)
ggsave(file = "outputs/fig1.pdf", width = 8, height = 12, dpi = 300)