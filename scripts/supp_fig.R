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

sc2_wild_cov75 <- get_samples_nat_imm(samples_full, "SARS-CoV-2 (pre-B.1.1.7 variants)", vac_eff = c(50, 70, 90)*0.75)
sc2_b117_cov75 <- get_samples_nat_imm(samples_full, "SARS-CoV-2 (B.1.1.7)", vac_eff =  c(50, 70, 90)*0.75)

sc2_wild_cov50 <- get_samples_nat_imm(samples_full, "SARS-CoV-2 (pre-B.1.1.7 variants)", vac_eff = c(50, 70, 90)*0.5)
sc2_b117_cov50 <- get_samples_nat_imm(samples_full, "SARS-CoV-2 (B.1.1.7)", vac_eff =  c(50, 70, 90)*0.5)

prop_15p_in <- get_prop_15p_country() %>%
    group_by(income) %>%
    summarise(prop = 1 - weighted.mean(prop, pop)/100) %>%
    mutate(income = factor(income, levels = c("High income", "Upper middle income", "Lower middle income", "Low income")))


p1a <- plot_fig1bc(sc2_wild_cov75, prop_15p_in, "A")
p1b <- plot_fig1bc(sc2_b117_cov75, prop_15p_in, "B")

p2a <- plot_fig1bc(sc2_wild_cov50, prop_15p_in, "C")
p2b <- plot_fig1bc(sc2_b117_cov50, prop_15p_in, "D")

pt1 <- (p1a + p1b) & theme(legend.position = "top")
pt2 <- (p2a + p2b) & theme(legend.position = "top")

p1 <- pt1 / guide_area() + plot_layout(heights = c(1, 0.02), nrow = 2, guides = "collect")
p2 <- pt2 / guide_area() + plot_layout(heights = c(1, 0.02), nrow = 2, guides = "collect") 

(p1 / p2) + plot_layout(heights = c(1.2, 0.02, 1.2, 0.2), nrow = 4)
ggsave(file = "outputs/figS1.png", width = 10, height = 8, dpi = 300)