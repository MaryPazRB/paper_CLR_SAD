# ============================================================
# Precision (min plants) + Field time
# Coffee leaf rust SAD: Old vs New
# Thesis-ready script (plain R)
#
# Inputs:
#   1) field_all.csv   : field, plant, branch, leaf, S1 (%), sad, time
#   2) df_ratings2.csv : method, actual (GS), estimate, error (= estimate-actual)
#
# Output objects:
#   - sig2_meas_by_method
#   - curves_all, minima_all
#   - recommendations
#   - times_lookup, sim_all, summary_sim
#   - p_curves, p_min, p_time, p_pct_saved, p_hours_saved
# ============================================================

# ----------------------------
# 0) Packages
# ----------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(lme4)
  library(ggplot2)
  library(cowplot)
  library(purrr)
  library(forcats)
})

# ----------------------------
# 1) User inputs (edit here)
# ----------------------------
# Precision target for field mean severity (percentage points)
h_target <- 2

# Candidate range for number of plants (p)
p_grid <- 2:80

# Critical value for 95% CI (Wald); you can swap for qt(.975, df)
tcrit <- 1.96

# Protocol definitions (branches × leaves per branch)
protocol_levels <- c(
  "4 B × 4 L",
  "4 B × 2 L",
  "2 B × 4 L",
  "2 B × 2 L"
)

protocols <- tibble::tribble(
  ~protocol,              ~r_fix, ~l_fix,
  "4 B × 4 L",     4,     4,
  "4 B × 2 L",     4,     2,
  "2 B × 4 L",     2,     4,
  "2 B × 2 L",     2,     2
) %>%
  mutate(protocol = factor(protocol, levels = protocol_levels))

# --- Time model (seconds) ---
# Movement times (assumptions)
t_move_plant  <- 15
t_move_branch <- 5

# Predicted mean time per branch when evaluating 4 leaves/branch (from GLMM)
t_branch4_old <- 28.4
t_branch4_new <- 23.6

# Helper: scale branch time from 4 leaves to l leaves
t_branch_fun <- function(t4, l) t4 * (l / 4)

# ----------------------------
# 2) Read data + harmonize method names
# ----------------------------
df_ratings2 <- read_csv("df_ratings2.csv", show_col_types = FALSE) %>%
  mutate(
    method = factor(method),
    method = fct_recode(method,
                        Old = "aidedold",
                        New = "aidednew")
  )

field_all <- read_csv("field_all.csv", show_col_types = FALSE) 
 

# Basic checks
stopifnot(all(c("Old", "New") %in% levels(df_ratings2$method)))
stopifnot(all(c("Old", "New") %in% levels(field_all$sad)))

# ----------------------------
# 3) Measurement error variance (calibration dataset)
# ----------------------------
sig2_meas_by_method <- df_ratings2 %>%
  group_by(method) %>%
  summarise(
    sig2_meas = var(error, na.rm = TRUE),
    n_calib   = sum(!is.na(error)),
    .groups = "drop"
  )

print(sig2_meas_by_method)

# ----------------------------
# 4) Variance components from field data (per field)
#    Model: S1 ~ 1 + (1|plant) + (1|plant:branch)
# ----------------------------
fit_components_field <- function(dat_field) {
  
  dat_field <- dat_field %>%
    mutate(
      plant = factor(plant),
      branch = factor(branch),
      plant_branch = interaction(plant, branch, drop = TRUE)
    )
  
  fit <- lmer(S1 ~ 1 + (1 | plant) + (1 | plant_branch),
              data = dat_field, REML = TRUE)
  
  vc <- as.data.frame(VarCorr(fit))
  
  # Extract safely
  sig2_plant  <- vc$vcov[vc$grp == "plant"][1]
  sig2_branch <- vc$vcov[vc$grp == "plant_branch"][1]
  sig2_leaf   <- attr(VarCorr(fit), "sc")^2
  
  tibble(sig2_plant = sig2_plant,
         sig2_branch = sig2_branch,
         sig2_leaf = sig2_leaf)
}

# ----------------------------
# 5) Precision model: half-width of 95% CI for field mean
#    Var(mean) = sig2_plant/p + sig2_branch/(p*r) + (sig2_leaf+sig2_meas)/(p*r*l)
# ----------------------------
half_width_fun <- function(p, r, l,
                           sig2_plant, sig2_branch, sig2_leaf, sig2_meas,
                           tcrit = 1.96) {
  
  var_mean <- sig2_plant / p +
    sig2_branch / (p * r) +
    (sig2_leaf + sig2_meas) / (p * r * l)
  
  tcrit * sqrt(var_mean)
}

# ----------------------------
# 6) Time model for a single field (minutes)
# ----------------------------
time_total_minutes <- function(p, r, l, method,
                               t_move_plant  = 15,
                               t_move_branch = 5,
                               t_branch4_old = 28.4,
                               t_branch4_new = 23.6) {
  
  t4 <- ifelse(method == "Old", t_branch4_old, t_branch4_new)
  t_branch <- t_branch_fun(t4, l)
  
  # Per plant:
  #   - assess r branches: r * t_branch
  #   - move between branches: (r-1) * t_move_branch
  #   - move between plants: t_move_plant
  T_sec <- p * (r * t_branch + (r - 1) * t_move_branch + t_move_plant)
  
  T_sec / 60
}

# ----------------------------
# 7) Analyze one field: precision curves + minima
# ----------------------------
analyze_field <- function(field_id,
                          field_all,
                          protocols,
                          p_grid,
                          sig2_meas_by_method,
                          h_target = 2,
                          tcrit = 1.96) {
  
  dat_field <- field_all %>% filter(field == field_id)
  comps <- fit_components_field(dat_field)
  
  grid <- protocols %>%
    crossing(tibble(p = p_grid)) %>%
    mutate(
      r = r_fix,
      l = l_fix,
      n_leaves = p * r * l
    )
  
  res <- sig2_meas_by_method %>%
    filter(method %in% c("Old", "New")) %>%
    crossing(grid) %>%
    mutate(
      half_width = half_width_fun(
        p = p, r = r, l = l,
        sig2_plant  = comps$sig2_plant,
        sig2_branch = comps$sig2_branch,
        sig2_leaf   = comps$sig2_leaf,
        sig2_meas   = sig2_meas,
        tcrit = tcrit
      ),
      pass  = half_width <= h_target,
      field = field_id
    )
  
  minima <- res %>%
    filter(pass) %>%
    group_by(field, method, protocol) %>%
    arrange(p) %>%
    slice(1) %>%
    ungroup()
  
  list(curves = res, minima = minima, comps = comps)
}

# ----------------------------
# 8) Run all fields
# ----------------------------
fields <- sort(unique(field_all$field))

field_labels <- field_all %>%
  group_by(field) %>%
  summarise(
    meanS1 = mean(S1, na.rm = TRUE),
    sdS1   = sd(S1,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(field_lab = paste0(field, " (mean=", sprintf("%.2f", meanS1), "%)"))

print(field_labels)

outs <- lapply(
  fields,
  analyze_field,
  field_all = field_all,
  protocols = protocols,
  p_grid = p_grid,
  sig2_meas_by_method = sig2_meas_by_method,
  h_target = h_target,
  tcrit = tcrit
)

curves_all <- bind_rows(lapply(outs, `[[`, "curves")) %>%
  left_join(field_labels %>% select(field, field_lab), by = "field")

minima_all <- bind_rows(lapply(outs, `[[`, "minima")) %>%
  mutate(protocol = factor(protocol, levels = protocol_levels)) %>%
  left_join(field_labels %>% select(field, field_lab), by = "field")

# If empty, diagnose quickly:
if (nrow(minima_all) == 0) {
  message("minima_all is empty: no designs met the target half-width. ",
          "Try increasing p_grid upper bound or relaxing h_target.")
}

print(minima_all %>% arrange(field, protocol, method))

# ----------------------------
# 9) Add time predictions for minima
# ----------------------------
minima_all <- minima_all %>%
  mutate(
    time_min = time_total_minutes(
      p = p, r = r, l = l, method = method,
      t_move_plant  = t_move_plant,
      t_move_branch = t_move_branch,
      t_branch4_old = t_branch4_old,
      t_branch4_new = t_branch4_new
    )
  )

print(minima_all %>% select(field_lab, method, protocol, p, half_width, time_min))

write_csv(minima_all, "minimaall.csv")

# ----------------------------
# 10) Plots (thesis-ready)
# ----------------------------


p_curves <- ggplot(curves_all, aes(x = p, y = half_width, color = method)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = h_target, linetype = "dashed") +
  facet_grid(field_lab ~ protocol) +
  scale_color_viridis_d(begin = 0.5, end = 1,direction = -1)+
  cowplot::theme_half_open() +
  labs(x = "Number of plants sampled (p)",
       y = "95% CI half-width (percentage points)",
       color = "SAD")

p_min <- minima_all %>%
  ggplot(aes(x = protocol, y = p, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ field_lab, nrow = 1) +
  cowplot::theme_half_open() +
  coord_flip()+
  scale_fill_viridis_d(begin = 0.5, end = 1, direction = -1)+
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
  labs(x = NULL,
       y = paste0("Minimum plants (95% CI ± ", h_target, " p.p.)"),
       fill = "SAD")

p_time <- minima_all %>%
  ggplot(aes(x = protocol, y = time_min, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ field_lab, nrow = 1) +
  coord_flip()+
  cowplot::theme_half_open() +
  scale_fill_viridis_d(begin = 0.5, end = 1, direction = -1)+
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
  labs(x = NULL,
       y = "Predicted total time per field (minutes)",
       fill = "SAD")

library(patchwork)

fig <- ((p_min  + theme_minimal_grid(font_size = 12)) /
          (p_time + theme_minimal_grid(font_size = 12))) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

ggsave("figs/fig_min_plants_time.png", fig, width = 7, height = 5, dpi = 300)

# Optional save:
# ggsave("fig_precision_curves.png", p_curves, width = 10, height = 7, dpi = 300)
# ggsave("fig_min_plants.png",       p_min,    width = 10, height = 3, dpi = 300)
# ggsave("fig_time_minutes.png",     p_time,   width = 10, height = 3, dpi = 300)

# ----------------------------
# 11) Recommended design per field & method (min time among passing)
# ----------------------------
recommendations <- minima_all %>%
  group_by(field, method) %>%
  arrange(time_min) %>%
  slice(1) %>%
  ungroup() %>%
  select(field_lab, field, method, protocol, p, r, l, half_width, time_min)

print(recommendations)

# ----------------------------
# 12) Seasonal simulation: expected time saved across many fields
#      Field class A/B/C ~ severity classes (based on your 3 fields)
# ----------------------------
times_lookup <- minima_all %>%
  select(method, protocol, field, time_min) %>%
  pivot_wider(names_from = method, values_from = time_min) %>%
  rename(time_new = New, time_old = Old)

simulate_season <- function(times_lookup,
                            protocol,
                            n_fields = 500,
                            probs = c(A = 0.3, B = 0.4, C = 0.3),
                            n_sims = 5000,
                            seed = 1) {
  
  probs <- probs[c("A","B","C")]
  if (abs(sum(probs) - 1) > 1e-8) stop("probs must sum to 1.")
  if (!protocol %in% unique(times_lookup$protocol)) stop("protocol not found in times_lookup.")
  
  lut <- times_lookup %>%
    filter(protocol == !!protocol) %>%
    mutate(field = factor(field, levels = c("A","B","C"))) %>%
    arrange(field)
  
  # Map A/B/C -> time
  time_old_vec <- lut$time_old[match(c("A","B","C"), as.character(lut$field))]
  time_new_vec <- lut$time_new[match(c("A","B","C"), as.character(lut$field))]
  
  set.seed(seed)
  draws <- replicate(
    n_sims,
    sample.int(3, size = n_fields, replace = TRUE, prob = probs),
    simplify = "matrix"
  )
  draws <- matrix(draws, nrow = n_fields, ncol = n_sims)
  
  tot_old <- colSums(matrix(time_old_vec[draws], nrow = n_fields, ncol = n_sims))
  tot_new <- colSums(matrix(time_new_vec[draws], nrow = n_fields, ncol = n_sims))
  
  tibble(
    protocol     = protocol,
    n_fields     = n_fields,
    pA           = probs["A"],
    pB           = probs["B"],
    pC           = probs["C"],
    sim          = seq_len(n_sims),
    hours_saved  = (tot_old - tot_new) / 60,
    pct_saved    = 100 * (tot_old - tot_new) / tot_old
  )
}

run_simulation_grid <- function(times_lookup,
                                n_fields = 500,
                                n_sims = 5000,
                                seed = 1) {
  
  scenarios <- tibble::tribble(
    ~scenario, ~pA,  ~pB,  ~pC,
    "Mild",     0.6,  0.3,  0.1,
    "Typical",  0.3,  0.4,  0.3,
    "Severe",   0.1,  0.3,  0.6
  )
  
  grid <- tidyr::expand_grid(protocol = unique(times_lookup$protocol), scenarios)
  
  purrr::pmap_dfr(
    .l = grid,
    .f = function(protocol, scenario, pA, pB, pC) {
      simulate_season(
        times_lookup = times_lookup,
        protocol = protocol,
        n_fields = n_fields,
        probs = c(A = pA, B = pB, C = pC),
        n_sims = n_sims,
        seed = seed
      ) %>% mutate(scenario = scenario)
    }
  )
}

summarize_sims <- function(sim_df) {
  sim_df %>%
    group_by(protocol, scenario, n_fields, pA, pB, pC) %>%
    summarise(
      pct_saved_mean   = mean(pct_saved),
      pct_saved_q025   = quantile(pct_saved, 0.025),
      pct_saved_q975   = quantile(pct_saved, 0.975),
      hours_saved_mean = mean(hours_saved),
      hours_saved_q025 = quantile(hours_saved, 0.025),
      hours_saved_q975 = quantile(hours_saved, 0.975),
      .groups = "drop"
    ) %>%
    mutate(
      protocol = factor(protocol, levels = protocol_levels),
      scenario = factor(scenario, levels = c("Mild","Typical","Severe"))
    ) %>%
    arrange(scenario, protocol)
}

# Run simulations (can be time-consuming)
sim_all     <- run_simulation_grid(times_lookup, n_fields = 500, n_sims = 10000, seed = 123)
summary_sim <- summarize_sims(sim_all)
print(summary_sim)

p_pct_saved <- summary_sim %>%
  ggplot(aes(x = scenario, y = pct_saved_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = pct_saved_q025, ymax = pct_saved_q975), width = 0.15) +
  facet_wrap(~ protocol, nrow = 1) +
  theme_bw() +
  labs(x = NULL, y = "% time saved (New vs Old)")

p_hours_saved <- summary_sim %>%
  ggplot(aes(x = scenario, y = hours_saved_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = hours_saved_q025, ymax = hours_saved_q975), width = 0.15) +
  facet_wrap(~ protocol, nrow = 1) +
  theme_bw() +
  labs(x = NULL, y = "Hours saved per season (New vs Old)")

print(p_pct_saved)
print(p_hours_saved)

# Optional save:
# write_csv(minima_all, "minima_all.csv")
# write_csv(summary_sim, "summary_sim.csv")
# ggsave("fig_pct_saved.png", p_pct_saved, width = 10, height = 3, dpi = 300)
# ggsave("fig_hours_saved.png", p_hours_saved, width = 10, height = 3, dpi = 300)



