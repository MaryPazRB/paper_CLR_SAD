# ============================================================
# Pre-validation analysis (plain R script)
# Focus: Hybrid coverage probability + CCC-derived metrics
# ============================================================

# ----------------------------
# 0) Packages
# ----------------------------
library(gsheet)
library(dplyr)
library(tidyr)
library(ggplot2)
library(epiR)
library(lme4)
library(emmeans)
library(DHARMa)
library(patchwork)
library(scales)
library(cowplot)
library(pbkrtest)
library(lmerTest)
# New libraries for ICC and Excel export
library(irr)
library(writexl)


# ----------------------------
# 1) Download & reshape data
# ----------------------------
gs_id <- "141vQED6DLbzC0tQfvz16r_qR-UeQ0Hd7VOU46fJ_Z-Q"

read_tab_long <- function(gid, method_name) {
  url <- paste0("https://docs.google.com/spreadsheets/d/", gs_id, "/edit?gid=", gid, "#gid=", gid)
  gsheet2tbl(url) |>
    pivot_longer(cols = 2:22, names_to = "rater", values_to = "estimate") |>
    mutate(method = method_name)
}

unaided_h <- read_tab_long(gid = 0, method_name = "unaided")
aidedold_h <- read_tab_long(gid = 206725966, method_name = "aidedold")
aidednew_h <- read_tab_long(gid = 492462463, method_name = "aidednew")

df_all <- bind_rows(aidednew_h, aidedold_h, unaided_h)

# ----------------------------
# 2) Build "actual" (truth) and analysis dataset
# ----------------------------
# Truth is "unaided" tab, rater == "actual"
truth <- df_all |>
  filter(method == "unaided", rater == "actual") |>
  transmute(leaf, actual = estimate)

# Remove excluded raters + remove the 'actual' pseudo-rater
excluded_raters <- c("Carlos", "Beatriz", "Isabella", "actual")

df_all2 <- df_all |>
  left_join(truth, by = "leaf") |>
  filter(!rater %in% excluded_raters) |>
  mutate(
    method = factor(method, levels = c("unaided", "aidedold", "aidednew"))
  )

# Optional quick visuals (comment out if not needed)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)

# Garantir ordem dos métodos
df_all2 <- df_all2 |>
  mutate(method = factor(method, levels = c("unaided", "aidedold", "aidednew")))

# Função para gerar o gráfico de um método
plot_method <- function(m) {
  df_all2 |>
    filter(method == m) |>
    ggplot(aes(actual, estimate)) +
    geom_point(alpha = 0.7, shape = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    coord_fixed() +
    xlim(0, 100) +
    labs(
      title = "",
      x = "Actual (%)",
      y = "Estimate (%)"
    ) +
    theme_half_open()
}

# Criar os três gráficos
p_unaided <- plot_method("unaided") + geom_point(alpha = 0.5, color = "#440154")
p_aidedold <- plot_method("aidedold") + geom_point(alpha = 0.5, color = "#21908C")
p_aidednew <- plot_method("aidednew") + geom_point(alpha = 0.5, color = "#FDE725")

# Combinar com patchwork e tags A/B/C
combo1 <- (p_unaided | p_aidedold | p_aidednew) +
  plot_annotation(tag_levels = "A")

ggsave("panel_est_actual.png", bg = "white", width = 9, height = 8)
# ggplot(df_all2, aes(actual, estimate - actual)) + geom_point() + geom_hline(yintercept=0) + facet_grid(method ~ rater)

# ----------------------------
# 3) CCC by method x rater (long format)
# ----------------------------
ccc_by_method_rater_long <- function(df) {
  df |>
    group_by(method, rater) |>
    summarise(
      {
        d <- dplyr::pick(actual, estimate) |> tidyr::drop_na()

        if (nrow(d) < 2) {
          tibble(
            Agreement = NA_real_,
            `Bias coefficient` = NA_real_,
            Precision = NA_real_,
            scale_shift = NA_real_,
            location_shift = NA_real_,
            n = nrow(d)
          )
        } else {
          fit <- epi.ccc(d$estimate, d$actual)
          Agreement <- as.numeric(fit$rho.c[1])
          Bias <- as.numeric(fit$C.b)

          tibble(
            Agreement = Agreement,
            `Bias coefficient` = as.numeric(Bias),
            Precision = Agreement / Bias,
            scale_shift = as.numeric(fit$s.shift),
            location_shift = as.numeric(fit$l.shift),
            n = nrow(d)
          )
        }
      },
      .groups = "drop"
    )
}

df_raters <- ccc_by_method_rater_long(df_all2)

# ----------------------------
# 4) Mixed model plots for Agreement / Precision / Bias coefficient
#    (Model on logit scale; back-transform manually)
# ----------------------------
logit01 <- function(x, eps = 1e-6) qlogis(pmin(pmax(x, eps), 1 - eps))

fit_emm_plot <- function(data, response_var, ylab,
                         method_levels = c("unaided", "aidedold", "aidednew"),
                         eps = 1e-6,
                         ylim = c(0.7, 1.0),
                         do_dharma = FALSE) {
  dat <- data |>
    mutate(
      .eta = logit01(.data[[response_var]], eps = eps),
      method = factor(method, levels = method_levels)
    )

  fit <- lmer(.eta ~ method + (1 | rater), data = dat)

  if (do_dharma) {
    plot(simulateResiduals(fit))
  }

  em_link <- emmeans(fit, ~method) # link (logit) scale
  em_df <- as.data.frame(em_link) |>
    mutate(
      estimate = plogis(emmean),
      lower = plogis(lower.CL),
      upper = plogis(upper.CL),
      letter = case_when(
        method %in% c("aidedold", "aidednew") ~ "a",
        method == "unaided" ~ "b",
        TRUE ~ ""
      )
    )

  p <- ggplot(em_df, aes(x = method, color = method, y = estimate)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
    geom_text(aes(label = letter, y = upper), vjust = -0.6, size = 4) +
    labs(x = NULL, y = ylab) +
    coord_cartesian(ylim = ylim) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none") +
    scale_color_viridis_d(begin = 0, end = 1)

  list(fit = fit, em_df = em_df, plot = p)
}

# Independent plots
p_agree <- fit_emm_plot(df_raters, "Agreement", "Agreement",
  ylim = c(0.7, 1.0)
)$plot
p_prec <- fit_emm_plot(df_raters, "Precision", "Precision",
  ylim = c(0.7, 1.0)
)$plot
p_bias <- fit_emm_plot(df_raters, "Bias coefficient", "Bias coefficient",
  ylim = c(0.7, 1.0)
)$plot


# ----------------------------
# 5) Coverage probability (HYBRID rule) + mixed-effects logistic model
# ----------------------------
# Prepare rated/true in percent scale (0–100)
df_cp <- df_all2 |>
  transmute(
    leaf, rater, method,
    rated = estimate,
    true = actual
  ) |>
  mutate(
    err  = rated - true,
    aerr = abs(err)
  )

# --- Hybrid rule parameters (tune here) ---
delta_rel_prop <- 0.20 # relative tolerance (20% of true) for moderate/high severities
low_cut_pp <- 10 # threshold below which we use absolute tolerance
delta_abs_low <- 2 # absolute tolerance (percentage points) when true <= low_cut_pp

# Hybrid indicator
df_cp <- df_cp |>
  mutate(
    within_hyb = as.integer(if_else(true <= low_cut_pp,
      aerr <= delta_abs_low,
      aerr <= delta_rel_prop * pmax(true, 1e-8)
    ))
  )

# Mixed-effects logistic model (accounts for rater + leaf)
cp_hyb_glmm <- glmer(
  within_hyb ~ method + (1 | rater) + (1 | leaf),
  data = df_cp,
  family = binomial
)

# Marginal CP by method (probabilities + 95% CI)
hyb_emm <- emmeans(cp_hyb_glmm, ~method, type = "response")
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
      n = n(),
      k = sum(.data[[flag_var]], na.rm = TRUE),
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
    position = position_dodge(width = 0.45), linewidth = 0.9
  ) +
  geom_point(position = position_dodge(width = 0.45), size = 2.2) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    x = "True severity band (%)",
    y = "Coverage probability",
    color = "Method"
  ) +
  scale_color_viridis_d() +
  cowplot::theme_half_open(font_size = 12) +
  theme(
    legend.position = c(0.98, 0.98), # top-right inside
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.6, "lines")
  )

# ----------------------------
# 7) 2x2 panel (3 CCC plots + Hybrid CP; coverage plot is the 4th)
# ----------------------------
panel_2x2 <- ((p_agree + p_prec) / (p_bias + p_cp)) + plot_annotation(tag_levels = "A")
print(panel_2x2)

# ----------------------------
# 8) (Optional) Save outputs
# ----------------------------
ggsave("figs/panel_2x2.png", panel_2x2, width = 6, height = 5.5, dpi = 300)
# write.csv(df_raters, "df_raters_ccc.csv", row.names = FALSE)
# write.csv(cp_hyb_tbl, "cp_hybrid_by_method.csv", row.names = FALSE)

# ----------------------------
# 9) Interrater Reliability (ICC)
# ----------------------------


df_icc <- df_all2 |>
  group_by(method) |>
      select(leaf, rater, estimate) |>
      pivot_wider(names_from = rater, values_from = estimate) |>
      select(-leaf)

new <- df_icc |> 
  filter(method == "aidednew") |> 
  ungroup() |> 
  select(-method)

old <- df_icc |> 
  filter(method == "aidedold") |> 
  ungroup() |> 
  select(-method)

unaided <- df_icc |> 
  filter(method == "unaided") |> 
  ungroup() |> 
  select(-method)


library(epiR)
epi.occc(new, na.rm = FALSE, pairs = TRUE)
epi.occc(old, na.rm = FALSE, pairs = TRUE)
epi.occc(unaided, na.rm = FALSE, pairs = TRUE)


#ICC 2

library(psych)
icc_unaided <- ICC(unaided)
knitr::kable(icc_unaided$results[1:2]) # only selected columns

icc_new <- ICC(new)
knitr::kable(icc_new$results[1:2]) # only selected columns

icc_old <- ICC(old)
knitr::kable(icc_old$results[1:2]) # only selected columns




