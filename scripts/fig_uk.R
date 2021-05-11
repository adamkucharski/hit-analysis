source(here::here("scripts", "funcs.R"))
sample_r0_sar <- get_samples_fit(data_xl %>% filter(str_detect(pathogen, "^sar")), dist = "norm", param = "r0") ## iii) R0 for sars-cov-2

r0_sample_b117 <- sample_r0_sar[[2]] %>% pull(samples)

values <- (1 - 1/r0_sample_b117) %>% quantile(c(0.05, 0.25, 0.5, 0.75, 0.95))
threshold_df <- data.frame(
    quantile = c(0.05, 0.25, 0.5, 0.75, 0.95),
    threshold = values
)

prev <- c(12.8, 15.5, 16.9, 17.8, 20.2, 18.7, 23.2, 28.3, 30.8, 34.8, 36.9, 43.0, 47.6, 51.5, 52.2, 53.1, 61.5, 68.3, 90)
dates <- as.Date("2020-12-13") + (1:length(prev) - 1)*7

df_plot <- data.frame(date = dates,
`All children have antibodies` = (prev / 100),
`No children have antibodies` = (prev / 100)*(1 - 0.213)) %>%
    pivot_longer(!date, names_to = "population", values_to = "prev")

df_plot %>% 
    ggplot(aes(x = date, y = prev)) + 
    geom_line(aes(color = population)) + 
    theme_bw() + theme(aspect.ratio = 0.8) + coord_cartesian(ylim = c(0, 1)) + 
    geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.6417625, ymax = 0.8334208,
                   fill = "pink", alpha = 0.03) + 
    geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.73859025, ymax = 0.8047842,
                   fill = "pink", alpha = 0.03) + 
    geom_hline(data = threshold_df, aes(yintercept = threshold), linetype = c("dotted", "dashed", "solid", "dashed", "dotted")) + 
    geom_text(data=data.frame(x=as.Date("2021-01-01"), y=values, labels = c(0.05, 0.25, 0.5, 0.75, 0.95)), aes(x, y, label=labels),  vjust=0) + 
    labs(x = "date", y = "Proportion with antibodies", title = "Seroprevalence in UK and relationship to HIT")


    