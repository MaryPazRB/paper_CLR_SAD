# ============================================================
# Time-to-assess (Old vs New) — compact workflow
# Data: Google Sheet (3 tabs = fields A, B, C)
# Outcome: time (seconds) per branch assessment
# Model: Gamma GLMM (log link), random: field/evaluator2/branch
# Outputs: plots, effect size (ratio + absolute), total time (3 fields)
# ============================================================

# ----------------------------
# 0) Packages
# ----------------------------
library(gsheet)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(viridis)
library(ggridges)
library(cowplot)
library(patchwork)
library(lme4)
library(emmeans)
library(DHARMa)

# ----------------------------
# 1) Import + reshape
# ----------------------------
sheet_url <- "https://docs.google.com/spreadsheets/d/1_fO1nLXZzxPsKPBDb4JY3d76Ddthaf3bbTo7Xw76sVY/edit"

read_field <- function(gid, field_id, evaluator_map) {
  gsheet2tbl(paste0(sheet_url, "?gid=", gid, "#gid=", gid)) |>
    select(evaluator, plant, branch, leaf, S1, TIME_OLD, TIME_NEW) |>
    rename(Old = TIME_OLD, New = TIME_NEW) |>
    group_by(evaluator, plant, branch) |>
    slice(1) |>
    ungroup() |>
    pivot_longer(c(Old, New), names_to = "sad", values_to = "time") |>
    mutate(
      field = field_id,
      evaluator2 = recode(evaluator, !!!evaluator_map, .default = evaluator)
    )
}

field_all <- bind_rows(
  read_field(0,          "A", c("DEBORA" = "A", "EMME" = "B")),
  read_field(1502996524, "B", c("HUGO"   = "A", "PAZ"  = "B")),
  read_field(1841716165, "C", c("HUGO"   = "A", "PAZ"  = "B"))
) |>
  mutate(
    sad = factor(sad, levels = c("Old", "New")),
    field = factor(field),
    evaluator2 = factor(evaluator2)
  )

write_csv(field_all, "field_all.csv")

# ----------------------------
# 2) QC summary (optional)
# ----------------------------
field_all |>
  group_by(field) |>
  summarise(
    meanS1 = mean(S1, na.rm = TRUE),
    sdS1   = sd(S1, na.rm = TRUE),
    .groups = "drop"
  )

# ----------------------------
# 3) Ridgeline: time per branch (seconds)
# ----------------------------
df_plot <- field_all |>
  filter(is.finite(time))

means_sec <- df_plot |>
  group_by(sad) |>
  summarise(mean_sec = mean(time, na.rm = TRUE), .groups = "drop")

p_ridges <- ggplot(df_plot, aes(x = time, y = sad, fill = sad)) +
  geom_density_ridges(
    alpha = 0.7,
    scale = 1.4,
    rel_min_height = 0.01,
    color = "white",
    linewidth = 0.25
  ) +
  geom_vline(
    data = means_sec,
    aes(xintercept = mean_sec, color = sad),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  scale_fill_viridis_d(begin = 0.5, end = 1) +
  scale_color_viridis_d(begin = 0.5, end = 1, guide = "none") +
  theme_half_open() +
  theme(legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),      # <- ESSENCIAL
    axis.line.x  = element_line(),      # mantém a base
    panel.grid   = element_blank()
  )+
  labs(x = "Assessment time/branch (sec)", y = NULL)

p_ridges

make_box <- function(f, ev, subtitle, y_lim = c(0, 100)) {
  field_all |>
    filter(field == f, evaluator2 == ev, is.finite(time)) |>
    ggplot(aes(sad, time, color = sad)) +
    geom_boxplot() +
    scale_color_viridis_d(begin = 0.5, end = 1) +
    theme_half_open() +
    theme(legend.position = "none") +
    coord_cartesian(ylim = y_lim) +
    labs(subtitle = subtitle, x = NULL, y = "Time (s)")
}

row1 <- make_box("A","A","Field 1") + make_box("A","B","Field 1")
row2 <- make_box("B","A","Field 2") + make_box("B","B","Field 2")
row3 <- make_box("C","A","Field 3") + make_box("C","B","Field 3")

p_boxes <- (row1 / row2 / row3) + plot_annotation(tag_levels = "A")

# ----------------------------
# 4) Gamma GLMM + diagnostics + mean plot
# ----------------------------
m_time <- glmer(
  time ~ sad + (1 | field/evaluator2/branch),
  data = field_all,
  family = Gamma(link = "log")
)

plot(simulateResiduals(m_time))

emm <- emmeans(m_time, ~ sad, type = "response")
emm_df <- as.data.frame(emm)

p_means <- ggplot(emm_df, aes(sad, response, color = sad)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = .05) +
  scale_color_viridis_d(begin = 0.5, end = 1) +
  theme_half_open() +
  theme(legend.position = "none") +
  ylab("") +
  xlab(NULL)

# Combined view
(p_boxes | (p_ridges / p_means)) + plot_annotation(tag_levels = "A")

ggsave("figs/fig_time_branch.png", width = 10, height = 8)

# ----------------------------
# 5) Effect sizes + total time (Old vs New)
# ----------------------------
# Relative effect: ratio New/Old (primary for Gamma log-link)
ratio_eff <- contrast(
  emm,
  method = list("New/Old" = c(Old = -1, New = 1)),
  ratio  = TRUE
)
ratio_ci <- as.data.frame(summary(ratio_eff, infer = TRUE)) |>
  mutate(pct_change = (ratio - 1) * 100)

ratio_ci  # ratio < 1 => New faster

# Absolute effect: difference New - Old (seconds) on response scale
diff_abs <- contrast(
  regrid(emm),
  method = list("New - Old" = c(Old = -1, New = 1))
)
diff_ci <- as.data.frame(summary(diff_abs, infer = TRUE))

diff_ci  # negative => New faster (seconds saved)





# garantir nomes padronizados
diff_ci2 <- diff_ci |>
  rename(diff_sec = estimate)

ratio_ci2 <- ratio_ci  # já tem a coluna 'ratio'



# ----------------------------
# Forest plots (ratio + diff) — stable for patchwork tags
# ----------------------------
p_ratio <- ggplot(ratio_ci2, aes(x = ratio, y = 1)) +
  geom_errorbarh(aes(xmin = asymp.LCL, xmax = asymp.UCL), height = 0, linewidth = 1.2) +
  geom_point(size = 4) +
  geom_vline(xintercept = 1, linetype = 2) +
  labs(x = "Time ratio (New / Old)", y = NULL) +
  theme_half_open() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),      # <- ESSENCIAL
    axis.line.x  = element_line(),      # mantém a base
    panel.grid   = element_blank()
  )

p_diff <- ggplot(diff_ci2, aes(x = diff_sec, y = 1)) +
  geom_errorbarh(aes(xmin = asymp.LCL, xmax = asymp.UCL), height = 0, linewidth = 1.2) +
  geom_point(size = 4) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "Time diff (New − Old, sec)", y = NULL) +
  theme_half_open() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),      # <- ESSENCIAL
    axis.line.x  = element_line(),      # mantém a base
    panel.grid   = element_blank()
  )


# Top row
p_top <- (p_ridges | p_means) + plot_layout(widths = c(1.07, 0.93))


# Forest row
p_forest <- (p_ratio | p_diff) +
  plot_layout(widths = c(1, 1))


# Final figure (tags ONLY here)
fig <- (p_top / p_forest) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag.position = c(0, 1)
  )

fig

ggsave("figs/figura_plots.png", width = 7, height= 5)




