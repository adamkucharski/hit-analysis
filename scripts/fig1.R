# R script to make Figure 1
source(here::here("scripts", "funcs.R"))
# 1. RUN SAMPLING AND PLOT
# 1.1 Get the samples for fitted distributions or bootstraps
sample_eff_all <- get_samples_fit(data_xl, dist = "beta", param = "eff", init = c(10, 10)) ## i) efficacy
sample_r0_flu_mea <- get_samples_fit(data_xl %>% filter(str_detect(pathogen, "^flu|mea")), dist = "lnorm", param = "r0") ## ii) R0 for flu and measles
sample_r0_sar <- get_samples_fit(data_xl %>% filter(str_detect(pathogen, "^sar")), dist = "norm", param = "r0") ## iii) R0 for sars-cov-2
sample_r0_other <- get_samples_bs(data_xl %>% filter(str_detect(pathogen, "^mum|var|rub")), samples_bs_r0, param = "r0") ## iv) R0 for mumps, measles, rubella

# 1.2 Combine all samples together, rename pathogens and save
pathogen_levels <- c("measles", "mumps", "rubella", "varicella", "sars_cov_2_wt", "sars_cov_2_b117",
"flu_h1n1pmd09", "flu_h3n2", "flu_b")
pathogen_labels <- c("Measles", "Mumps", "Rubella", "Varicella-Zoster", "SARS-CoV-2 (pre-B.1.1.7 variants)",
    "SARS-CoV-2 (B.1.1.7)", "Flu A(H1N1)",  "Flu A(H3N2)", "Flu B")
samples_full <- 
    bind_rows(sample_eff_all, sample_r0_flu_mea, sample_r0_sar, sample_r0_other) %>%
    arrange(pathogen) %>%
    group_by(metric) %>%
    mutate(row = row_number()) %>%
    pivot_wider(names_from = metric, values_from = samples) %>%
    as.data.frame %>%
    mutate(pathogen = factor(pathogen, levels = pathogen_levels, labels = pathogen_labels))
samples_median <- samples_full %>% group_by(pathogen) %>% summarise(eff = median(eff), r0 = median(r0)) # Get medians
save(samples_full, file = here::here("outputs", "fig_data", "samples_fig1a.RDS"))

# 1.3 Get data for Figures 1b/c
sc2_wild <- get_samples_nat_imm(samples_full, "SARS-CoV-2 (pre-B.1.1.7 variants)", vac_eff = c(50, 70, 90))
sc2_b117 <- get_samples_nat_imm(samples_full, "SARS-CoV-2 (B.1.1.7)", vac_eff =  c(50, 70, 90))
prop_15p_in <- get_prop_15p_country() %>%
    group_by(income) %>%
    summarise(prop = 1 - weighted.mean(prop, pop)/100) %>%
    mutate(income = factor(income, levels = c("High income", "Upper middle income", "Lower middle income", "Low income")))
save(sc2_wild, file = here::here("outputs", "fig_data", "samples_fig1b.RDS"))
save(sc2_b117, file = here::here("outputs", "fig_data", "samples_fig1c.RDS"))

# 1.4. Make Figure 1 and save
p1a <- plot_fig1a(samples_full, samples_median)
p1b <- plot_fig1bc(sc2_wild, prop_15p_in, "B")
p1c <- plot_fig1bc(sc2_b117, prop_15p_in, "C")
p1bc <-  p1b / p1c / guide_area() + plot_layout(guides = "collect", heights = c(3, 3, 0.8))
p1a / p1bc + plot_layout(heights = c(1, 2.5))
ggsave(file = "outputs/fig1.png", width = 8, height = 12, dpi = 300)
ggsave(file = "outputs/fig1.pdf", width = 8, height = 12, dpi = 300)
#p1a + p1bc + plot_layout(widths = c(1, 1))
#ggsave(file = "outputs/fig1.png", width = 8, height = 8, dpi = 300)
#ggsave(file = "outputs/fig1.pdf", width = 8, height = 8, dpi = 300)


# 2. USEFUL METRICS
# 2.1. Proportion of samples above HIT threshold for each pathogen
prop_path_hi <- samples_full %>%
     group_by(pathogen) %>%
     summarise(prop = mean(eff > (1 - 1 / r0))) %>%
     mutate(lower = prop - 2*sqrt(prop*(1-prop)/1000), upper = prop + 2*sqrt(prop*(1-prop)/1000))