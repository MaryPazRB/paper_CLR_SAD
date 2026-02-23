
# ============================================================
# Field validation – ratings vs gold standard (plain R script)
# Short, organized, minimal redundancy
# ============================================================

# ----------------------------
# 1) Packages (only what is used)
# ----------------------------
library(gsheet)
library(dplyr)
library(tidyr)
library(ggplot2)
library(epiR)
library(cowplot)
library(patchwork)
library(ggthemes)  # theme_minimal_grid(), scale_color_few()

# ----------------------------
# 2) Load data (wide) from Google Sheets
# ----------------------------
field_wide <- gsheet2tbl(
  "https://docs.google.com/spreadsheets/d/1rVFraYtTIUxoIfk5w7F4LdoniNPRXu6R/edit?gid=1394248356#gid=1394248356"
)

# ----------------------------
# 3) Tidy to long + define method + errors
# ----------------------------
df_long <- field_wide %>%
  pivot_longer(
    cols = -c(image, GS),
    names_to = "rater",
    values_to = "estimate"
  ) %>%
  rename(actual = GS) %>%
  mutate(
    method = case_when(
      rater %in% c("R1OLD", "R2OLD") ~ "aidedold",
      rater %in% c("R1NEW", "R2NEW") ~ "aidednew",
      TRUE                          ~ "Algorithm"
    ),
    error = estimate - actual,
    abs_error = abs(error)
  )

# Keep only human raters (OLD/NEW) and set order
df_ratings <- df_long %>%
  filter(rater %in% c("R1OLD", "R2OLD", "R1NEW", "R2NEW")) %>%
  mutate(method = factor(method, levels = c("aidedold", "aidednew")))

# Remove problematic images (keep as a vector for easy maintenance)
drop_images <- c(
  "set_16_13_obj1.png",
  "set_15_20_obj3.png",
  "set_16_64_obj3.png",
  "set_16_60_obj3.png"
)

df_ratings2 <- df_ratings %>%
  filter(!image %in% drop_images)

# ----------------------------
# 4) Quick diagnostic plot: error vs actual by method
# ----------------------------
p_err <- df_ratings2 %>%
  ggplot(aes(actual, estimate - actual, color = method)) +
  geom_point(alpha = 0.8) +
  geom_smooth(se = FALSE, color = "black") +
  geom_hline(yintercept = 0, linetype = 2) +
  cowplot::theme_minimal_grid() +
  scale_color_few() +
  coord_cartesian(ylim = c(-50, 50)) +
  facet_wrap(~ method) +
  labs(x = "Actual severity (%)", y = "Error (pp)") +
  theme(legend.position = "bottom")

# print(p_err)

# ----------------------------
# 5) CCC components by method x rater (Agreement, Bias, Precision)
# ----------------------------
ccc_components <- function(df) {
  df %>%
    group_by(method, rater) %>%
    summarise(
      {
        d <- dplyr::pick(actual, estimate) %>% tidyr::drop_na()
        if (nrow(d) < 2) {
          tibble(
            Agreement = NA_real_,
            `Bias coefficient` = NA_real_,
            Precision = NA_real_
          )
        } else {
          fit <- epi.ccc(d$estimate, d$actual)
          Agreement <- as.numeric(fit$rho.c[1])
          Bias      <- as.numeric(fit$C.b)
          tibble(
            Agreement = Agreement,
            `Bias coefficient` = Bias,
            Precision = Agreement / Bias
          )
        }
      },
      .groups = "drop"
    )
}

df_raters <- ccc_components(df_ratings2)

# Create annotation labels
df_labels <- df_raters %>%
  mutate(
    label = sprintf(
      "Agreement = %.2f\nBias = %.2f\nPrecision = %.2f",
      Agreement, `Bias coefficient`, Precision
    )
  ) %>%
  select(method, rater, label)

# ----------------------------
# 6) Scatter plot function (method x rater) + annotation (top-left)
# ----------------------------
cols_method <- c(
  aidedold = "#21908C",
  aidednew = "#FDE725"
)

plot_mr <- function(dat, labels_df, m, r, lim = 60) {
  
  lab <- labels_df %>%
    filter(method == m, rater == r) %>%
    pull(label)
  
  dat %>%
    filter(method == m, rater == r) %>%
    ggplot(aes(actual, estimate)) +
    geom_point(color = cols_method[[m]], alpha = 0.85) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
    annotate(
      "text",
      x = lim * 0.03, y = lim * 0.97,
      label = lab, hjust = 0, vjust = 1, size = 3.5
    ) +
    coord_cartesian(xlim = c(0, lim), ylim = c(0, lim)) +
    theme_half_open() +
    labs(title = "", x = "Actual severity (%)", y = "Estimate (%)")
}

# ----------------------------
# 7) Build 2x2 panel (one plot per method x rater) + save
# ----------------------------
comb <- tibble::tribble(
  ~method,    ~rater,
  "aidedold", "R1OLD",
  "aidedold", "R2OLD",
  "aidednew", "R1NEW",
  "aidednew", "R2NEW"
)

p_list <- purrr::pmap(comb, ~plot_mr(df_ratings2, df_labels, ..1, ..2, lim = 60))



glimpse(df_ratings2)



# ----------------------------
# 5) Coverage probability (HYBRID rule) + mixed-effects logistic model
# ----------------------------
# Prepare rated/true in percent scale (0–100)
df_cp <- df_ratings2 |>
  transmute(
    image, rater, method,
    rated = estimate,
    true  = actual
  ) |>
  mutate(
    err  = rated - true,
    aerr = abs(err)
  )

# --- Hybrid rule parameters (tune here) ---
delta_rel_prop <- 0.20  # relative tolerance (20% of true) for moderate/high severities
low_cut_pp     <- 10    # threshold below which we use absolute tolerance
delta_abs_low  <- 2     # absolute tolerance (percentage points) when true <= low_cut_pp

# Hybrid indicator
df_cp <- df_cp |>
  mutate(
    within_hyb = as.integer(if_else(true <= low_cut_pp,
                                    aerr <= delta_abs_low,
                                    aerr <= delta_rel_prop * pmax(true, 1e-8)))
  )

# Mixed-effects logistic model (accounts for rater + leaf)
cp_hyb_glmm <- glmer(
  within_hyb ~ method + (1 | rater) + (1 | image),
  data = df_cp,
  family = binomial
)

# Marginal CP by method (probabilities + 95% CI)
hyb_emm <- emmeans(cp_hyb_glmm, ~ method, type = "response")
cp_hyb_tbl <- as.data.frame(hyb_emm) |>
  transmute(
    method,
    cp = prob,
    lo = asymp.LCL,
    hi = asymp.UCL
  )

print(cp_hyb_tbl)

# ----------------------------
# 6) Hybrid CP by severity band (descriptive plot)
# ----------------------------
breaks <- c(-Inf, 5, 10, 20, Inf)
labels <- c("≤5", "(5,10]", "(10,20]", ">20")

df_cp <- df_cp |>
  mutate(sev_band = cut(true, breaks = breaks, labels = labels, right = TRUE))

band_cp_num <- function(data, flag_var) {
  data |>
    group_by(method, sev_band) |>
    summarise(
      n  = n(),
      k  = sum(.data[[flag_var]], na.rm = TRUE),
      cp = k / pmax(n, 1),
      se = sqrt(pmax(cp * (1 - cp) / pmax(n, 1), 0)),
      lo = pmax(0, cp - 1.96 * se),
      hi = pmin(1, cp + 1.96 * se),
      .groups = "drop"
    )
}

hyb_band <- band_cp_num(df_cp, "within_hyb")

p_cp <- ggplot(hyb_band, aes(x = sev_band, y = cp, color = method, group = method)) +
  geom_pointrange(aes(ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.45), linewidth = 0.9) +
  geom_point(position = position_dodge(width = 0.45), size = 2.2) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    x = "True severity band (%)",
    y = "Coverage Prob. (%)",
    color = "Method"
  ) +
  scale_color_viridis_d(begin = 0.5, end =1)+
  cowplot::theme_half_open(font_size = 12) +
  theme(
    legend.position = c(0.98, 0.98),     # top-right inside
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    legend.key.size = unit(0.6, "lines")
  )


panel_2x2 <- ((p_list[[1]] | p_list[[2]]) /
  (p_list[[3]] | p_list[[4]]))/ p_cp  +
  plot_annotation(tag_levels = "A")



ggsave("plots_field.png", panel_2x2, bg = "white", width = 6, height = 10, dpi = 300)


