# R script to make Figure 2
source(here::here("scripts", "funcs.R"))
# 1. GET DATA AND SAVE
data_sero_prev <- get_sero_prev_prop()
save(data_sero_prev, file = here::here("outputs", "fig_data", "samples_fig2.RDS"))

# 2. GET DATA AND SAVE
p2a <- plot_fig2(data_sero_prev, sc2_wild, "a) Pre-B.1.1.7 variants")
p2b <- plot_fig2(data_sero_prev, sc2_b117, "b) B.1.1.7 variants")
p2a / p2b + plot_layout(guides = "collect", heights = c(3, 3))
ggsave(file = "outputs/fig2.png", width = 8, height = 10, dpi = 300)
ggsave(file = "outputs/fig2.pdf", width = 8, height = 10, dpi = 300)