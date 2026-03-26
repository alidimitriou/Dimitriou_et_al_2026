# R Code for:
# Dimitriou, A., Benson-Amram, S., Gaynor, K. M., Percy, M., Burton, A. C. (2026). 
# Sharing the Trail: the effect of recreation on bear behaviour in a partially closed 
# Rocky Mountain Park. Journal of Wildlife Management.

# ===============================================================
# ANALYSIS 1 — Weekly habitat use
# ===============================================================
getwd()
#setwd("~/Desktop/Dimitriou_et_al_2026")

# ---------------------------------------------------------------
# ANALYSIS 1.0 — Libraries
# ---------------------------------------------------------------
library(dplyr)
library(brms)
library(posterior)
library(tidyr)
library(ggplot2)
library(ggcorrplot)
library(car)
library(performance)

# ---------------------------------------------------------------
# ANALYSIS 1.1 — Load prepared data + output folder
# ---------------------------------------------------------------
bear_model_data <- read.csv("Data/combined_bear_model_data.csv") %>%
  filter(is.finite(Effort), Effort > 0) %>%
  mutate(
    log_effort = log(Effort),
    Species_Label = case_when(
      Species == "Ursus americanus" ~ "Black bear",
      Species == "Ursus arctos"     ~ "Grizzly bear",
      TRUE ~ Species
    )
  ) %>%
  drop_na(log_effort)

output_dir <- "Output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# ---------------------------------------------------------------
# ANALYSIS 1.2 — Build covariates (proximity + recreation)
# ---------------------------------------------------------------
#lambda_m <- 500

#bear_model_data <- bear_model_data %>%
#  mutate(
    # distance
#    prox_trail        = exp(-dist_trail / lambda_m),
#    prox_trail_scaled = as.numeric(scale(prox_trail)),
    
    # recreation
#    rec_week_log      = log1p(rec_week),
#    rec_week_scaled   = as.numeric(scale(rec_week_log))
# )


#bear_model_data <- bear_model_data %>%
#  mutate(
#    prox_trail_scaled = as.numeric(scale(prox_trail)),
#    rec_week_log      = log1p(rec_week),
#    rec_week_scaled   = as.numeric(scale(rec_week_log))
#  )

#add scaled variables to csv
#write.csv(bear_model_data, "Data/combined_bear_model_data.csv", row.names = FALSE)

# ---------------------------------------------------------------
# ANALYSIS 1.3 — Multicollinearity checks
# ---------------------------------------------------------------
lm_vif_interact <- lm(
  Bear_Detections ~ prox_trail_scaled + rec_week_scaled +
    mean_ndvi_scaled + elevation_scaled,
  data = bear_model_data %>%
    drop_na(Bear_Detections, prox_trail_scaled, rec_week_scaled,
            mean_ndvi_scaled, elevation_scaled)
)

print(car::vif(lm_vif_interact))
print(performance::check_collinearity(lm_vif_interact))

X <- bear_model_data %>%
  select(prox_trail_scaled, rec_week_scaled, mean_ndvi_scaled, elevation_scaled) %>%
  drop_na()

cor_mat <- cor(X, method = "pearson")
print(round(cor_mat, 2))
ggcorrplot(cor_mat, hc.order = TRUE, type = "lower", lab = TRUE)

# ---------------------------------------------------------------
# ANALYSIS 1.4 — Independence threshold sensitivity check
# ---------------------------------------------------------------
detections <- read.csv("Data/BLT_detection_data_urs_hom.csv")

thresholds <- c(5, 10, 15, 20, 30, 45, 60)

results <- lapply(thresholds, function(t) {
  detections %>%
    filter(Species == "Ursus americanus") %>%
    arrange(Deployment.Location.ID, Date_Time.Captured) %>%
    group_by(Deployment.Location.ID) %>%
    mutate(
      Time_Diff = as.numeric(difftime(Date_Time.Captured, lag(Date_Time.Captured), units = "mins")),
      Event_ID = cumsum(is.na(Time_Diff) | Time_Diff > t)
    ) %>%
    group_by(Deployment.Location.ID, Event_ID) %>%
    slice_min(Date_Time.Captured, n = 1) %>%
    ungroup() %>%
    summarise(
      Threshold = t,
      Total_Events = n()
    )
}) %>%
  bind_rows()

print(results)

results <- results %>%
  arrange(Threshold) %>%
  mutate(
    Prev = lag(Total_Events),
    Delta = Total_Events - Prev,
    Pct_Change = 100 * (Delta / Prev)
  )

ggplot(results, aes(x = Threshold, y = Total_Events)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 30, linetype = "dashed") +
  labs(x = "Independence threshold (minutes)", y = "Total independent events",
       title = "Event count vs. independence threshold") +
  theme_classic(base_size = 14)

# ---------------------------------------------------------------
# ANALYSIS 1.5 — Fit weekly habitat-use model
# ---------------------------------------------------------------
run_bear_model <- function(df, species_name) {
  cat("\nRunning model for:", species_name, "\n")
  
  data_fit <- df %>%
    mutate(week_start = as.Date(Week_Start)) %>%   # correct column
    drop_na(
      Bear_Detections,
      prox_trail_scaled,
      rec_week_scaled,
      elevation_scaled,
      mean_ndvi_scaled,
      log_effort,
      Deployment.Location.ID,
      Week_Start
    ) %>%
    arrange(Deployment.Location.ID, week_start) %>%
    group_by(Deployment.Location.ID) %>%
    mutate(week_id = row_number()) %>%
    ungroup()
  
  
  map_df <- data_fit %>%
    mutate(.row_fit = row_number()) %>%
    select(.row_fit, Deployment.Location.ID, week_start, Week_Label, week_id)
  
  model <- brm(
    formula = Bear_Detections ~
      prox_trail_scaled * rec_week_scaled +
      elevation_scaled + mean_ndvi_scaled +
      offset(log_effort) + (1 | Deployment.Location.ID),
    autocor = cor_ar(~ week_id | Deployment.Location.ID, p = 1),
    data = data_fit,
    family = negbinomial(),
    prior = c(
      prior(normal(0, 2), class = "b"),
      prior(normal(0, 10), class = "Intercept"),
      prior(cauchy(0, 2), class = "sd"),
      prior(exponential(1), class = "shape")
    ),
    chains = 4, iter = 10000, warmup = 4000, seed = 123,
    control = list(adapt_delta = 0.98, max_treedepth = 12)
  )
  
  post <- as_draws_df(model)
  
  base_labels <- c(
    "b_prox_trail_scaled"                 = "Proximity to trail",
    "b_rec_week_scaled"                   = "Recreation Intensity",
    "b_prox_trail_scaled:rec_week_scaled" = "Trail × Recreation",
    "b_mean_ndvi_scaled"                  = "NDVI",
    "b_elevation_scaled"                  = "Elevation"
  )
  
  effect_df <- post %>%
    select(any_of(names(base_labels))) %>%
    pivot_longer(everything(), names_to = "Variable", values_to = "Estimate") %>%
    group_by(Variable) %>%
    summarise(
      median  = median(Estimate),
      lower95 = quantile(Estimate, 0.025),
      upper95 = quantile(Estimate, 0.975),
      lower80 = quantile(Estimate, 0.10),
      upper80 = quantile(Estimate, 0.90),
      .groups = "drop"
    ) %>%
    mutate(
      Species  = species_name,
      Term_raw = Variable,
      Term     = dplyr::recode(Variable, !!!base_labels, .default = Variable)
    )
  
  
  fname_stub <- gsub(" ", "_", tolower(species_name))
  saveRDS(model,  file.path(output_dir, paste0(fname_stub, "_model.rds")))
  saveRDS(map_df, file.path(output_dir, paste0(fname_stub, "_model_data_map.rds")))
  capture.output(summary(model), file = file.path(output_dir, paste0(fname_stub, "_summary.txt")))
  
  list(model = model, effects = effect_df, var_labels = base_labels)
}

# ---------------------------------------------------------------
# ANALYSIS 1.6 — Fit models (Black bear, Grizzly bear)
# ---------------------------------------------------------------

brms_black <- run_bear_model(
  bear_model_data %>% filter(Species == "Ursus americanus"),
  "Black bear"
)

brms_grizzly <- run_bear_model(
  bear_model_data %>% filter(Species == "Ursus arctos"),
  "Grizzly bear"
)



# ---------------------------------------------------------------
# ANALYSIS 1.7 — Combine posterior effects + export
# ---------------------------------------------------------------
effects_all <- bind_rows(
  brms_black$effects,
  brms_grizzly$effects
) %>%
  select(Species, Term, Term_raw, median, lower95, upper95, lower80, upper80)

write.csv(effects_all, file.path(output_dir, "posterior_effects_both_species.csv"), row.names = FALSE)

# ---------------------------------------------------------------
# ANALYSIS 1.8 — Reload posterior effects (from CSV)
# ---------------------------------------------------------------
effects_all <- read.csv(
  file.path(output_dir, "posterior_effects_both_species.csv"),
  stringsAsFactors = FALSE
)

effects_all <- effects_all %>%
  mutate(
    Species = factor(Species, levels = c("Black bear", "Grizzly bear"))
  )

# ---------------------------------------------------------------
# ANALYSIS 1.9 — Plot posterior effects (combined species)
# ---------------------------------------------------------------
effects_all <- effects_all %>%
  mutate(
    Term = dplyr::recode(
      Term,
      "Recreation intensity" = "Recreation Intensity"
    )
  )

level_order <- effects_all %>%
  pull(Term) %>%
  unique() %>%
  rev()

plot_df <- effects_all %>%
  filter(!is.na(median), !is.na(lower80), !is.na(upper80)) %>%
  mutate(
    Term   = factor(Term, levels = level_order),
    y_base = as.numeric(Term),
    y      = case_when(
      Species == "Black bear"   ~ y_base + 0.12,
      Species == "Grizzly bear" ~ y_base - 0.12,
      TRUE ~ y_base
    )
  ) %>%
  filter(!is.na(y))

plt_combined <- ggplot(plot_df, aes(y = y, color = Species)) +
  geom_segment(aes(x = lower80, xend = upper80, yend = y), linewidth = 3) +
  geom_segment(aes(x = lower95, xend = upper95, yend = y), linewidth = 1.2) +
  geom_point(aes(x = median, shape = Species), size = 6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "blue", linewidth = 1) +
  scale_y_continuous(breaks = seq_along(level_order), labels = level_order) +
  scale_color_manual(values = c("Black bear" = "#1b9e77",
                                "Grizzly bear" = "#7A4988")) +
  scale_shape_manual(values = c("Black bear" = 16,
                                "Grizzly bear" = 17)) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.ticks = element_line(color = "black"),
    axis.line  = element_line(color = "black", linewidth = 0.8)
  ) +
  labs(x = "Posterior estimate", y = NULL,
       color = "Species", shape = "Species")

ggsave(
  file.path(output_dir, "posterior_effects_combined.png"),
  plot = plt_combined,
  width = 10, height = 7, dpi = 300
)

plt_combined

# ===============================================================
#  Diagnostics
#   PP checks, trace plots, Bayesian p-values, Moran's I
# ===============================================================
# ---------------------------------------------------------------
# D0) Libraries
# ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(bayesplot); library(brms); library(posterior)
  library(sf); library(spdep)
})

# ---------------------------------------------------------------
# D1) Paths / folders
# ---------------------------------------------------------------
output_dir   <- "Output"
diag_dir     <- file.path(output_dir, "diagnostics")
pp_dir       <- file.path(diag_dir, "pp_checks")
trace_dir    <- file.path(diag_dir, "trace_plots")
acf_dir      <- file.path(diag_dir, "acf_plots")
acf_site_dir <- file.path(acf_dir,   "by_site")

dir.create(pp_dir,       showWarnings = FALSE, recursive = TRUE)
dir.create(trace_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(acf_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(acf_site_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------
# D2) Inputs (station coords, model file paths, map file paths)
# ---------------------------------------------------------------
coords_df <- read.csv("Data/stations.csv") |>
  rename(Deployment.Location.ID = station_id) |>
  select(Deployment.Location.ID, longitude, latitude)

model_files <- file.path(output_dir, c("black_bear_model.rds", "grizzly_bear_model.rds"))
model_files <- model_files[file.exists(model_files)]

map_paths <- list(
  "black_bear_model.rds"   = file.path(output_dir, "black_bear_model_data_map.rds"),
  "grizzly_bear_model.rds" = file.path(output_dir, "grizzly_bear_model_data_map.rds")
)

# ---------------------------------------------------------------
# D3) Helper functions
# ---------------------------------------------------------------
iso_week_start <- function(week_label) {
  as.Date(paste0(week_label, "-1"), format = "%G-W%V-%u")
}

compute_nb_pearson <- function(model) {
  y_var     <- as.character(model$formula$formula[[2]])
  y_obs     <- model$data[[y_var]]
  mu        <- fitted(model)[, "Estimate"]
  shape_est <- as.numeric(posterior_summary(model, variable = "shape")[,"Estimate"])
  var_mu    <- mu + (mu^2) / shape_est
  (y_obs - mu) / sqrt(var_mu)
}

# Join residuals to map timing to build time series per site & week
build_ts_resids <- function(model, map_path, species_name) {
  if (!file.exists(map_path)) {
    warning("Missing map file for ", species_name, ": ", map_path)
    return(NULL)
  }
  
  map_df <- readRDS(map_path)
  
  resid_vec <- compute_nb_pearson(model)
  stopifnot(length(resid_vec) == nrow(model$data))
  
  md <- model$data %>%
    mutate(
      Residuals = resid_vec,
      Deployment.Location.ID = as.character(Deployment.Location.ID)
    ) %>%
    select(Deployment.Location.ID, Residuals) %>%
    group_by(Deployment.Location.ID) %>%
    mutate(.dup = dplyr::row_number()) %>%
    ungroup()
  
  rd <- map_df %>%
    mutate(
      Deployment.Location.ID = as.character(Deployment.Location.ID),
      week_date = dplyr::case_when(
        "Week_Start" %in% names(map_df) & !all(is.na(week_start)) ~ as.Date(week_start),
        "Week_Label" %in% names(map_df) & !all(is.na(Week_Label)) ~ iso_week_start(Week_Label),
        TRUE ~ as.Date(NA)
      )
    ) %>%
    filter(!is.na(week_date)) %>%
    arrange(Deployment.Location.ID, week_date) %>%
    group_by(Deployment.Location.ID) %>%
    mutate(.dup = dplyr::row_number()) %>%
    ungroup() %>%
    select(Deployment.Location.ID, week_date, .dup)
  
  if (nrow(rd) == 0) {
    warning("No valid week_date in map for ", species_name)
    return(NULL)
  }
  
  counts_md <- md %>% count(Deployment.Location.ID, name = "n_md")
  counts_rd <- rd %>% count(Deployment.Location.ID, name = "n_rd")
  
  counts <- full_join(counts_md, counts_rd, by = "Deployment.Location.ID") %>%
    mutate(n_keep = pmax(0, pmin(coalesce(n_md, 0L), coalesce(n_rd, 0L))))
  
  md_trim <- md %>%
    left_join(counts, by = "Deployment.Location.ID") %>%
    filter(.dup <= n_keep) %>%
    select(-n_md, -n_rd, -n_keep)
  
  rd_trim <- rd %>%
    left_join(counts, by = "Deployment.Location.ID") %>%
    filter(.dup <= n_keep) %>%
    select(-n_md, -n_rd, -n_keep)
  
  inner_join(md_trim, rd_trim, by = c("Deployment.Location.ID", ".dup")) %>%
    select(Deployment.Location.ID, week_date, Residuals) %>%
    filter(is.finite(Residuals), !is.na(week_date))
}

# ---------------------------------------------------------------
# D4) Storage objects
# ---------------------------------------------------------------
pp_plots      <- list()
trace_plots   <- list()
bayes_pvals   <- tibble(Species = character(), Bayes_p_value = numeric())
moran_summary <- tibble(Species = character(), Moran_I = numeric(), P_value = numeric())

# ---------------------------------------------------------------
# D5) Main model loop
# ---------------------------------------------------------------
for (f in model_files) {
  
  model <- readRDS(f)
  fname <- basename(f)
  
  species_name <- if (fname == "black_bear_model.rds") "Black bear" else "Grizzly bear"
  cat("\nRunning diagnostics for:", species_name, "\n")
  
  # -------------------------------------------------------------
  # D5.1) Posterior predictive check (mean)
  # -------------------------------------------------------------
  pp_plot <- pp_check(model, type = "stat", stat = "mean", ndraws = 200) +
    ggtitle(paste("PP check (mean):", species_name))
  
  ggsave(
    file.path(pp_dir, paste0(gsub(" ", "_", tolower(species_name)), "_pp_mean.png")),
    plot = pp_plot, width = 7, height = 5, dpi = 300
  )
  pp_plots[[species_name]] <- pp_plot
  
  # -------------------------------------------------------------
  # D5.2) Trace plots (fixed effects)
  # -------------------------------------------------------------
  posterior_array <- as.array(model)
  b_params <- grep("^b_", dimnames(posterior_array)[[3]], value = TRUE)
  
  if (length(b_params) > 0) {
    trace_plot <- bayesplot::mcmc_trace(posterior_array, pars = b_params) +
      ggtitle(paste("Trace Plots:", species_name)) +
      theme_bw(base_size = 12)
    
    ggsave(
      file.path(trace_dir, paste0(gsub(" ", "_", tolower(species_name)), "_trace.png")),
      plot = trace_plot, width = 10, height = 6, dpi = 300
    )
    trace_plots[[species_name]] <- trace_plot
  }
  
  # -------------------------------------------------------------
  # D5.3) Bayesian p-value (mean)
  # -------------------------------------------------------------
  y_var  <- as.character(model$formula$formula[[2]])
  y_obs  <- model$data[[y_var]]
  y_rep  <- posterior_predict(model, ndraws = 1000)
  
  T_obs  <- mean(y_obs, na.rm = TRUE)
  T_rep  <- apply(y_rep, 1, function(x) mean(x, na.rm = TRUE))
  p_bayes <- mean(T_rep >= T_obs)
  
  bayes_pvals <- bind_rows(
    bayes_pvals,
    tibble(Species = species_name, Bayes_p_value = round(p_bayes, 3))
  )
  
  # -------------------------------------------------------------
  # D5.4) Moran's I (station-mean Pearson residuals)
  # -------------------------------------------------------------
  resid_vec <- compute_nb_pearson(model)
  
  spatial_data <- model$data |>
    mutate(residuals = resid_vec) |>
    group_by(Deployment.Location.ID) |>
    summarise(mean_resid = mean(residuals, na.rm = TRUE), .groups = "drop") |>
    left_join(coords_df, by = "Deployment.Location.ID") |>
    filter(is.finite(mean_resid)) |>
    tidyr::drop_na(longitude, latitude)
  
  if (nrow(spatial_data) >= 5) {
    
    stations_sf <- sf::st_as_sf(
      spatial_data, coords = c("longitude", "latitude"), crs = 4326
    ) |>
      sf::st_transform(3005)
    
    coords_matrix <- sf::st_coordinates(stations_sf)
    k_use <- min(4, nrow(spatial_data) - 1)
    
    nb <- spdep::knn2nb(spdep::knearneigh(coords_matrix, k = k_use))
    lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
    
    vals <- as.numeric(spatial_data$mean_resid)
    
    if (sd(vals, na.rm = TRUE) > 0) {
      moran_result <- spdep::moran.test(vals, lw, zero.policy = TRUE, na.action = na.omit)
      
      moran_summary <- bind_rows(
        moran_summary,
        tibble(
          Species = species_name,
          Moran_I = round(moran_result$estimate[["Moran I statistic"]], 3),
          P_value = round(moran_result$p.value, 4)
        )
      )
      
    } else {
      message("Variance of station-mean residuals ~ 0; skipping Moran's I: ", species_name)
    }
    
  } else {
    cat("Not enough stations for Moran’s I:", species_name, "\n")
  }
}

# ---------------------------------------------------------------
# D6) Save outputs
# ---------------------------------------------------------------
if (length(pp_plots)) {
  ggsave(
    file.path(diag_dir, "Combined_PP_Checks.pdf"),
    wrap_plots(pp_plots, ncol = 2),
    width = 10, height = 8
  )
}

if (length(trace_plots)) {
  ggsave(
    file.path(diag_dir, "Combined_Trace_Plots.pdf"),
    wrap_plots(trace_plots, ncol = 1),
    width = 10, height = 15
  )
}

write.csv(bayes_pvals,   file.path(diag_dir, "bayesian_p_values.csv"), row.names = FALSE)
write.csv(moran_summary, file.path(diag_dir, "moran_results.csv"),     row.names = FALSE)

cat("\n=== Summaries ===\n")
print(bayes_pvals)
print(moran_summary)


# ===============================================================
# ANALYSIS 2: Diel activity patterns 
# ===============================================================

# ---------------------------------------------------------------
# 2.0) Libraries
# ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(lubridate)
  library(stringr)
  library(activity)
  library(overlap)
})

# ---------------------------------------------------------------
# 2.1) Working directory and inputs
# ---------------------------------------------------------------
setwd("/Volumes/ALI'S USB/DIEL_ACTIVITY") #  too large for local computer
base_dir <- getwd()

model_data <- read.csv(file.path(base_dir, "processed_data",
                                 "combined_bear_model_data.csv"))
img        <- read.csv(file.path(base_dir, "raw_data", "BLT_detection_data_urs_hom.csv"))
locs       <- read.csv(file.path(base_dir, "raw_data", "BLT_station_data.csv"))

# Output folder
output_dir <- file.path(base_dir, "diel_output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------------------------------------------------------------
# 2.2) Parse LOCAL timestamps from Image.ID and join station data
# ---------------------------------------------------------------
img <- img %>%
  mutate(
    datetime_string = str_extract(Image.ID, "\\d{4}-\\d{2}-\\d{2}__\\d{2}-\\d{2}-\\d{2}"),
    datetime_string = str_replace_all(
      datetime_string,
      c("__" = " ", "-" = ":", "(\\d{4}-\\d{2}-\\d{2})" = "\\1")
    ),
    Date_Time.Local = ymd_hms(datetime_string, tz = "America/Vancouver"),
    Date_Time.Captured_utc = with_tz(Date_Time.Local, "UTC")
  ) %>%
  select(-datetime_string)

img_locs <- left_join(img, locs, by = "Deployment.Location.ID")

# Station-specific timestamp corrections
img_locs <- img_locs %>%
  mutate(
    Date_Time.Local = case_when(
      Deployment.Location.ID == "BLT38" ~ Date_Time.Local %m+% months(1),
      Deployment.Location.ID == "BLT29" ~ Date_Time.Local + hours(12),
      TRUE ~ Date_Time.Local
    ),
    Date_Time.Captured_utc = with_tz(Date_Time.Local, "UTC")
  )

# ---------------------------------------------------------------
# 2.3) Filter survey windows (LOCAL) and remove missing coordinates
# ---------------------------------------------------------------
img_locs <- img_locs %>%
  filter(
    (Date_Time.Local >= ymd("2023-07-12") & Date_Time.Local <= ymd("2023-11-01")) |
      (Date_Time.Local >= ymd("2024-05-01") & Date_Time.Local <= ymd("2024-10-30")) |
      (Date_Time.Local >= ymd("2025-04-30") & Date_Time.Local <= ymd("2025-06-24"))
  ) %>%
  filter(!is.na(Date_Time.Local),
         !is.na(Latitude),
         !is.na(Longitude))

# ---------------------------------------------------------------
# 2.4) Convert clock time to solar time
# ---------------------------------------------------------------
st <- solartime(
  dat = img_locs$Date_Time.Local,
  lat = img_locs$Latitude,
  lon = img_locs$Longitude,
  tz  = -7
)

img_locs$solar <- st$solar
img_locs$clock <- st$clock

# ---------------------------------------------------------------
# 2.5) Week_Label (non-ISO) and trail_status assignment
# ---------------------------------------------------------------
if (!"Week_Label" %in% names(img_locs)) {
  img_locs <- img_locs %>%
    mutate(
      Week_Label = paste0(
        strftime(Date_Time.Local, "%Y"),
        "-W",
        sprintf("%02d", as.integer(strftime(Date_Time.Local, "%U")))
      )
    )
}

img_locs <- img_locs %>%
  mutate(
    trail_status = case_when(
      Deployment.Location.ID %in% c("BLT01","BLT03","BLT04","BLT07","BLT10","BLT14B") ~ "open",
      
      Deployment.Location.ID == "BLT16" & Week_Label <  "2024-W27" ~ NA_character_,
      Deployment.Location.ID == "BLT16" & Week_Label >= "2024-W27" ~ "closed",
      
      Deployment.Location.ID %in% c("BLT13","BLT15") &
        Week_Label >= "2024-W24" & Week_Label <= "2024-W29" ~ "open",
      Deployment.Location.ID %in% c("BLT13","BLT15") ~ "closed",
      
      Deployment.Location.ID %in% c("BLT17","BLT19","BLT20","BLT21","BLT22","BLT23",
                                    "BLT25","BLT26","BLT28","BLT29","BLT31","BLT33",
                                    "BLT34","BLT36","BLT40","BLT42","BLT44") ~ "closed",
      TRUE ~ NA_character_
    )
  )

# ---------------------------------------------------------------
# 2.6) Solar-time vectors for activity analyses
# ---------------------------------------------------------------
t_human <- img_locs$solar[img_locs$Species == "Homo sapiens"]

t_arctos_open   <- img_locs$solar[img_locs$Species == "Ursus arctos" &
                                    img_locs$trail_status == "open"]
t_arctos_closed <- img_locs$solar[img_locs$Species == "Ursus arctos" &
                                    img_locs$trail_status == "closed"]

t_am_open       <- img_locs$solar[img_locs$Species == "Ursus americanus" &
                                    img_locs$trail_status == "open"]
t_am_closed     <- img_locs$solar[img_locs$Species == "Ursus americanus" &
                                    img_locs$trail_status == "closed"]

# ---------------------------------------------------------------
# 2.7) Kernel density activity models
# ---------------------------------------------------------------
reps_model <- 1000

m_human             <- fitact(t_human,           sample = "model", reps = reps_model)
m_open_arctos       <- fitact(t_arctos_open,     sample = "model", reps = reps_model)
m_closed_arctos     <- fitact(t_arctos_closed,   sample = "model", reps = reps_model)
m_open_americanus   <- fitact(t_am_open,         sample = "model", reps = reps_model)
m_closed_americanus <- fitact(t_am_closed,       sample = "model", reps = reps_model)

# ---------------------------------------------------------------
# 2.8) Wald tests for activity curve differences
# ---------------------------------------------------------------
cat("\nGrizzly (open vs human):\n")
print(compareCkern(m_open_arctos, m_human, reps = 1000))

cat("\nGrizzly (closed vs human):\n")
print(compareCkern(m_closed_arctos, m_human, reps = 1000))

cat("\nGrizzly (open vs closed):\n")
print(compareCkern(m_open_arctos, m_closed_arctos, reps = 1000))

cat("\nBlack bear (open vs human):\n")
print(compareCkern(m_open_americanus, m_human, reps = 1000))

cat("\nBlack bear (closed vs human):\n")
print(compareCkern(m_closed_americanus, m_human, reps = 1000))

cat("\nBlack bear (open vs closed):\n")
print(compareCkern(m_open_americanus, m_closed_americanus, reps = 1000))

cat("\nInter-specific (closed):\n")
print(compareCkern(m_closed_americanus, m_closed_arctos, reps = 1000))

cat("\nInter-specific (open):\n")
print(compareCkern(m_open_americanus, m_open_arctos, reps = 1000))


# ---------------------------------------------------------------
# 2.9) Overlap estimates + 95% CIs
#   - Dhat1 if min(n) < 50, else Dhat4
# ---------------------------------------------------------------
pick_type <- function(x, y) ifelse(min(length(x), length(y)) < 50, "Dhat1", "Dhat4")

pairs <- list(
  arctos_open_vs_human       = list(x = t_arctos_open,   y = t_human),
  arctos_closed_vs_human     = list(x = t_arctos_closed, y = t_human),
  arctos_open_vs_closed      = list(x = t_arctos_open,   y = t_arctos_closed),
  americanus_open_vs_human   = list(x = t_am_open,       y = t_human),
  americanus_closed_vs_human = list(x = t_am_closed,     y = t_human),
  americanus_open_vs_closed  = list(x = t_am_open,       y = t_am_closed),
  intersp_open               = list(x = t_arctos_open,   y = t_am_open),
  intersp_closed             = list(x = t_arctos_closed, y = t_am_closed)
)

set.seed(42)

ov_rows <- lapply(names(pairs), function(nm) {
  x <- pairs[[nm]]$x
  y <- pairs[[nm]]$y
  n1 <- length(x)
  n2 <- length(y)
  
  typ <- pick_type(x, y)
  est <- overlapEst(x, y, type = typ)
  bs  <- bootstrap(x, y, 1000, type = typ) 
  ci  <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)
  
  data.frame(
    comparison = nm,
    type       = typ,
    estimate   = est,
    lower95    = as.numeric(ci[1]),
    upper95    = as.numeric(ci[2]),
    n1         = n1,
    n2         = n2
  )
})

overlap_df <- do.call(rbind, ov_rows)
print(overlap_df)

write.csv(overlap_df,
          "E:/DIEL_ACTIVITY/diel_output/bear_diel_overlap.csv",
          row.names = FALSE)

# ---------------------------------------------------------------
# 2.10) Plots (Open vs Closed vs Humans) 
# ---------------------------------------------------------------
output_dir <- "E:/DIEL_ACTIVITY/diel_output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

col_open   <- "red"
col_closed <- "blue"
col_human  <- "black"
y_max <- 0.20

# ---- Grizzly
png(filename = file.path(output_dir, "Ursus_arctos_open-closed-humans.png"),
    width = 6, height = 4, units = "in", res = 300)

plot(m_open_arctos, yunit = "density", data = "none", las = 1,
     tline = list(col = col_open, lwd = 2), cline = list(lty = 0),
     main = "Ursus arctos: Open vs Closed vs Humans",
     xlab = "Time of day (solar radians)", ylab = "Activity density",
     ylim = c(0, y_max))

plot(m_closed_arctos, yunit = "density", data = "none", add = TRUE,
     tline = list(col = col_closed, lwd = 2), cline = list(lty = 0))

plot(m_human, yunit = "density", data = "none", add = TRUE,
     tline = list(col = col_human, lwd = 2, lty = 2), cline = list(lty = 0))

legend("topright", legend = c("Open", "Closed", "Humans"),
       col = c(col_open, col_closed, col_human),
       lty = c(1, 1, 2), lwd = 2, bty = "n")

dev.off()

# ---- Black bear
png(filename = file.path(output_dir, "Ursus_americanus_open-closed-humans.png"),
    width = 6, height = 4, units = "in", res = 300)

plot(m_open_americanus, yunit = "density", data = "none", las = 1,
     tline = list(col = col_open, lwd = 2), cline = list(lty = 0),
     main = "Ursus americanus: Open vs Closed vs Humans",
     xlab = "Time of day (solar radians)", ylab = "Activity density",
     ylim = c(0, y_max))

plot(m_closed_americanus, yunit = "density", data = "none", add = TRUE,
     tline = list(col = col_closed, lwd = 2), cline = list(lty = 0))

plot(m_human, yunit = "density", data = "none", add = TRUE,
     tline = list(col = col_human, lwd = 2, lty = 2), cline = list(lty = 0))

legend("topright", legend = c("Open", "Closed", "Humans"),
       col = c(col_open, col_closed, col_human),
       lty = c(1, 1, 2), lwd = 2, bty = "n")

dev.off()

# ---------------------------------------------------------------
# 2.11) Species-vs-species overlap within same status (record-keeping)
# ---------------------------------------------------------------
t_griz_closed <- t_arctos_closed
t_bb_closed   <- t_am_closed

type_closed <- ifelse(min(length(t_griz_closed), length(t_bb_closed)) < 50, "Dhat1", "Dhat4")
est_closed  <- overlapEst(t_griz_closed, t_bb_closed, type = type_closed)
bs_closed   <- bootstrap(t_griz_closed, t_bb_closed, 2000, type = type_closed)
ci_closed   <- quantile(bs_closed, c(0.025, 0.975), na.rm = TRUE)

t_griz_open <- t_arctos_open
t_bb_open   <- t_am_open

type_open <- ifelse(min(length(t_griz_open), length(t_bb_open)) < 50, "Dhat1", "Dhat4")
est_open  <- overlapEst(t_griz_open, t_bb_open, type = type_open)
bs_open   <- bootstrap(t_griz_open, t_bb_open, 2000, type = type_open)
ci_open   <- quantile(bs_open, c(0.025, 0.975), na.rm = TRUE)

species_overlap_df <- data.frame(
  comparison = c("Grizzly vs Black bear (closed)", "Grizzly vs Black bear (open)"),
  type       = c(type_closed, type_open),
  estimate   = c(est_closed, est_open),
  lower95    = c(as.numeric(ci_closed[1]), as.numeric(ci_open[1])),
  upper95    = c(as.numeric(ci_closed[2]), as.numeric(ci_open[2])),
  n_griz     = c(length(t_griz_closed), length(t_griz_open)),
  n_bb       = c(length(t_bb_closed), length(t_bb_open))
)

print(species_overlap_df)

write.csv(species_overlap_df,
          "E:/DIEL_ACTIVITY/diel_output/bear_species_overlap_open_closed.csv",
          row.names = FALSE)

# ---------------------------------------------------------------
# 2.12) Plot overlap coefficients + CIs (Open vs Closed; includes species overlap)
#   - Uses:
#       diel_output/bear_diel_overlap.csv  (human–bear overlap)
#       diel_output/bear_species_overlap_open_closed.csv (species overlap)
# ---------------------------------------------------------------

# ---- Load diel human–bear overlap (keeps ONLY open/closed vs human)
hb <- read.csv(
  file.path("/Volumes", "ALI'S USB", "DIEL_ACTIVITY", "diel_output",
            "bear_diel_overlap.csv"),
  stringsAsFactors = FALSE
)

hb <- hb %>%
  mutate(
    comp_l = tolower(comparison),
    x_cat = case_when(
      grepl("americanus", comp_l) ~ "Black bear",
      grepl("arctos", comp_l)     ~ "Grizzly",
      TRUE ~ NA_character_
    ),
    status = case_when(
      grepl("open_vs_human", comp_l)   ~ "Open",
      grepl("closed_vs_human", comp_l) ~ "Closed",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(x_cat), !is.na(status))

# ---- Load species overlap (open vs closed)
sp <- read.csv(
  file.path("/Volumes", "ALI'S USB", "DIEL_ACTIVITY", "diel_output",
            "bear_species_overlap_open_closed.csv"),
  stringsAsFactors = FALSE
)

sp <- sp %>%
  mutate(
    comp_l = tolower(comparison),
    x_cat  = "Species Overlap",
    status = case_when(
      grepl("open", comp_l)   ~ "Open",
      grepl("closed", comp_l) ~ "Closed",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(status))
# ---- Combine data first
plot_df <- bind_rows(
  hb %>% select(x_cat, status, estimate, lower95, upper95),
  sp %>% select(x_cat, status, estimate, lower95, upper95)
)

# ---- Recode + factor AFTER plot_df exists
plot_df <- plot_df %>%
  mutate(
    x_cat = case_when(
      x_cat == "Black bear"       ~ "Black bear–Human",
      x_cat == "Grizzly bear"     ~ "Grizzly bear–Human",
      x_cat == "Species Overlap"  ~ "Black bear–Grizzly bear",
      TRUE ~ x_cat
    ),
    x_cat  = factor(
      x_cat,
      levels = c("Black bear–Human",
                 "Grizzly bear–Human",
                 "Black bear–Grizzly bear")
    ),
    status = factor(status, levels = c("Open", "Closed"))
  )


print(table(plot_df$x_cat, plot_df$status, useNA = "ifany"))

pd <- position_dodge(width = 0.5)

ggplot(plot_df, aes(x = x_cat, y = estimate, colour = status)) +
  geom_point(position = pd, size = 4) +
  geom_errorbar(aes(ymin = lower95, ymax = upper95),
                position = pd, width = 0.25, linewidth = 1.1) +
  scale_colour_manual(values = c("Open" = "red", "Closed" = "blue")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    x = "",
    y = expression("Temporal overlap coefficient (" * Delta * ")"),
    colour = ""
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 18),
    axis.title.y = element_text(size = 20)
  )


# ===============================================================
# ANALYSIS 3 — AAR (ratios + GLMM)
# ===============================================================
# =============================================================
# 3.1. Load libraries
# =============================================================
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(ggdist)
library(tidyr)
library(tidybayes)
library(posterior)
library(brms)
library(ggridges)
library(scales)

# =============================================================
# 3.2 Define trail status
# =============================================================
bear_model_data <- read.csv(
  "Data/combined_bear_model_data.csv"
) %>%
  select(Deployment.Location.ID, Week_Label, trail_status) %>%
  distinct()

# =============================================================
# 3.3 Load detection data (date-range filter replaces Month filter)
# =============================================================
detections <- read.csv(
  "Data/BLT_detection_data_urs_hom.csv",
  stringsAsFactors = FALSE
) %>%
  mutate(
    datetime_str = str_extract(Image.ID, "\\d{4}-\\d{2}-\\d{2}__\\d{2}-\\d{2}-\\d{2}"),
    datetime_str = str_replace(datetime_str, "__", " "),
    datetime_str = str_replace_all(datetime_str, "(\\d{2})-(\\d{2})-(\\d{2})$", "\\1:\\2:\\3"),
    Date_Time.Captured = ymd_hms(datetime_str),
    Year       = year(Date_Time.Captured),
    Week       = sprintf("W%02d", isoweek(Date_Time.Captured)),
    Week_Label = paste0(Year, "-", Week)
  ) %>%
  select(-datetime_str) %>%
  # --- date windows (replaces: filter(Month >= 5 & Month <= 10)) ---
  filter(
    (Date_Time.Captured >= ymd("2023-07-12") & Date_Time.Captured <= ymd("2023-11-01")) |
      (Date_Time.Captured >= ymd("2024-05-01") & Date_Time.Captured <= ymd("2024-10-30")) |
      (Date_Time.Captured >= ymd("2025-04-30") & Date_Time.Captured <= ymd("2025-06-24"))
  ) %>%
  # keep only weeks present in bear_model_data and join status
  filter(Week_Label %in%bear_model_data$Week_Label) %>%
  left_join(bear_model_data, by = c("Deployment.Location.ID", "Week_Label"))

# Check coverage
cat("Detection date range:\n")
print(range(detections$Date_Time.Captured, na.rm = TRUE))
table(detections$trail_status, useNA = "ifany")

# =============================================================
# 3.4.  Initialize results holder
# =============================================================
aar_list <- list()

# =============================================================
# 3.5 Loop over bear species to build AAR1 intervals
# =============================================================
for (bear_species in c("Ursus americanus", "Ursus arctos")) {
  message("Processing: ", bear_species)
  
  results <- data.frame()
  
  bear_detections <- detections %>%
    filter(Species == bear_species) %>%
    arrange(Deployment.Location.ID, Date_Time.Captured) %>%
    group_by(Deployment.Location.ID) %>%
    mutate(
      Time_Diff   = as.numeric(difftime(Date_Time.Captured, lag(Date_Time.Captured), units = "mins")),
      Independent = is.na(Time_Diff) | Time_Diff >= 30
    ) %>%
    filter(Independent) %>%
    ungroup()
  
  human_detections <- detections %>% filter(Species == "Homo sapiens")
  
  for (cam in unique(detections$Deployment.Location.ID)) {
    bears  <- bear_detections  %>% filter(Deployment.Location.ID == cam)
    humans <- human_detections %>% filter(Deployment.Location.ID == cam)
    trail_status <- unique(bears$trail_status)
    
    if (nrow(bears) < 2 | nrow(humans) < 1) next
    
    for (i in 1:(nrow(bears) - 1)) {
      bear1_time    <- bears$Date_Time.Captured[i]
      next_bear_time <- bears$Date_Time.Captured[i + 1]
      
      humans_between <- humans %>%
        filter(Date_Time.Captured > bear1_time, Date_Time.Captured < next_bear_time)
      
      n_humans         <- nrow(humans_between)
      first_human_time <- if (n_humans > 0) humans_between$Date_Time.Captured[1] else NA
      last_human_time  <- if (n_humans > 0) humans_between$Date_Time.Captured[n_humans] else NA
      
      T1 <- if (!is.na(first_human_time)) as.numeric(difftime(first_human_time, bear1_time, units = "mins")) else NA
      T2 <- if (!is.na(last_human_time))  as.numeric(difftime(next_bear_time, last_human_time,  units = "mins")) else NA
      
      
      max_T <- 7 * 30.44 * 24 * 60  # ~7 months in minutes
      # max_T <- 30 * 24 * 60   # 30 days in minutes
      
      AAR1 <- if (!is.na(T1) & !is.na(T2) &
                  T1 > 0 & T1 <= max_T &
                  T2 <= max_T) {
        T2 / T1
      } else {
        NA
      }
      
      
      results <- bind_rows(results, data.frame(
        Bear_Species = bear_species,
        Deployment.Location.ID = cam,
        trail_status = trail_status,
        Bear1 = bear1_time,
        Human_Start = first_human_time,
        Human_End = last_human_time,
        Bear2 = next_bear_time,
        T1 = T1,
        T2 = T2,
        AAR1 = AAR1,
        Num_Humans_Between = n_humans,
        Human_Present = n_humans > 0
      ))
    }
  }
  
  aar_list[[bear_species]] <- results
  message("  --> Intervals generated: ", nrow(results))
}

# =============================================================
# 3.6 Summaries
# =============================================================
black_bear_results   <- aar_list[["Ursus americanus"]]
grizzly_bear_results <- aar_list[["Ursus arctos"]]

indep_summary <- detections %>%
  filter(Species %in% c("Ursus americanus", "Ursus arctos")) %>%
  arrange(Species, Deployment.Location.ID, Date_Time.Captured) %>%
  group_by(Species, Deployment.Location.ID) %>%
  mutate(
    Time_Diff   = as.numeric(difftime(Date_Time.Captured, lag(Date_Time.Captured), units = "mins")),
    Independent = is.na(Time_Diff) | Time_Diff >= 30
  ) %>%
  filter(Independent) %>%
  ungroup() %>%
  group_by(Species, trail_status) %>%
  summarise(Total_Independent_Detections = n(), .groups = "drop")

aar1_summary <- bind_rows(black_bear_results, grizzly_bear_results) %>%
  group_by(Bear_Species, trail_status) %>%
  summarise(Valid_AAR1 = sum(!is.na(AAR1)), .groups = "drop")

final_summary <- indep_summary %>%
  left_join(aar1_summary, by = c("Species" = "Bear_Species", "trail_status"))

print(final_summary)

# =============================================================
# 3.7. AAR1 boxplot: trail status on x, species as colour
# =============================================================
plot_data <- bind_rows(black_bear_results, grizzly_bear_results) %>%
  filter(!is.na(AAR1), AAR1 > 0, !is.na(trail_status)) %>%
  mutate(
    Species_lbl = dplyr::case_when(
      Bear_Species == "Ursus americanus" ~ "Black bear",
      Bear_Species == "Ursus arctos"     ~ "Grizzly bear",
      TRUE ~ as.character(Bear_Species)
    ),
    Species_lbl  = factor(Species_lbl, levels = c("Black bear", "Grizzly bear")),
    trail_status = factor(trail_status,
                          levels = c("closed", "open"),
                          labels = c("Closed", "Open"))
    
  )


aar_plot <- ggplot(plot_data, aes(x = trail_status, y = AAR1, fill = Species_lbl)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85,
               position = position_dodge(width = 0.8)) +
  scale_y_log10(labels = scales::label_number(accuracy = 0.01)) +
  coord_cartesian(ylim = c(0.1, 4)) +
  scale_fill_manual(values = c("Black bear" = "#1b9e77", "Grizzly bear" = "#7A4988")) +
  labs(
    x = "Trail Status",
    y = "AAR (log scale)",
    fill = "Species"
  ) +
  theme_minimal(base_size = 28) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 18),
    
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.6)
  )


print(aar_plot)
ggsave(file.path(output_dir, "aar_plot_combined.png"),
       plot = aar_plot, width = 10, height = 7, dpi = 300)

aar_plot
# =============================================================
# 3.8. Summary table for AAR1
# =============================================================
summary_table <- plot_data %>%
  group_by(Bear_Species, trail_status) %>%
  summarise(
    Valid_AAR1  = sum(!is.na(AAR1)),
    Median_AAR1 = median(AAR1, na.rm = TRUE),
    Mean_AAR1   = mean(AAR1, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_table)

# =============================================================
# 3.9. Data for brms models
# =============================================================
black_bear_model_data <- black_bear_results %>%
  filter(!is.na(AAR1), AAR1 > 0) %>%
  mutate(
    log_AAR1 = log(AAR1),
    trail_status = factor(trail_status),
    Deployment.Location.ID = factor(Deployment.Location.ID)
  )

grizzly_bear_model_data <- grizzly_bear_results %>%
  filter(!is.na(AAR1), AAR1 > 0) %>%
  mutate(
    log_AAR1 = log(AAR1),
    trail_status = factor(trail_status),
    Deployment.Location.ID = factor(Deployment.Location.ID)
  )

# =============================================================
# 3.10. Fit models
# =============================================================
black_bear_aar1_model <- brm(
  log_AAR1 ~ trail_status + (1 | Deployment.Location.ID),
  data = black_bear_model_data,
  family = student(),
  chains = 4, iter = 4000, warmup = 2000, seed = 123,
  control = list(adapt_delta = 0.99)
)

grizzly_bear_aar1_model <- brm(
  log_AAR1 ~ trail_status + (1 | Deployment.Location.ID),
  data = grizzly_bear_model_data,
  family = student(),
  chains = 4, iter = 4000, warmup = 2000, seed = 123,
  control = list(adapt_delta = 0.99)
)

# =============================================================
# 3.11. Posterior draws & plots of trail_status effect
# =============================================================
extract_trail_status_effect <- function(model, species_label) {
  as_draws_df(model) %>%
    select(b_trail_statusopen) %>%
    summarise(
      Species = species_label,
      median  = median(b_trail_statusopen),
      lower95 = quantile(b_trail_statusopen, 0.025),
      upper95 = quantile(b_trail_statusopen, 0.975),
      lower80 = quantile(b_trail_statusopen, 0.10),
      upper80 = quantile(b_trail_statusopen, 0.90)
    )
}

effects_black   <- extract_trail_status_effect(black_bear_aar1_model, "Black bear")
effects_grizzly <- extract_trail_status_effect(grizzly_bear_aar1_model, "Grizzly bear")

plot_df <- bind_rows(effects_black, effects_grizzly) %>%
  mutate(Species = factor(Species, levels = c("Black bear", "Grizzly bear")))

extract_trail_status_effect_draws <- function(model, species_label) {
  as_draws_df(model) %>%
    select(b_trail_statusopen) %>%
    mutate(Species = species_label)
}

draws_black   <- extract_trail_status_effect_draws(black_bear_aar1_model, "Black bear")
draws_grizzly <- extract_trail_status_effect_draws(grizzly_bear_aar1_model, "Grizzly bear")

draws_df <- bind_rows(draws_black, draws_grizzly) %>%
  mutate(Species = factor(Species, levels = c("Black bear", "Grizzly bear")))

ggplot(draws_df, aes(x = b_trail_statusopen, y = Species, fill = Species)) +
  ggridges::geom_density_ridges(alpha = 0.6, scale = 1.1, rel_min_height = 0.01,
                                color = "black", size = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1) +
  theme_minimal(base_size = 16) +
  labs(x = "Posterior of Trail Status Effect (log AAR1)", y = NULL) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.line.x = element_line(linewidth = 0.8)
  )


summaries <- draws_df %>%
  group_by(Species) %>%
  summarise(
    median = median(b_trail_statusopen),
    low80  = quantile(b_trail_statusopen, 0.10),
    high80 = quantile(b_trail_statusopen, 0.90),
    low95  = quantile(b_trail_statusopen, 0.025),
    high95 = quantile(b_trail_statusopen, 0.975),
    .groups = "drop"
  )

# --- Plot: trail-status posterior effects ---
p_trailstatus <- ggplot(draws_df, aes(x = b_trail_statusopen, y = Species, fill = Species)) +
  stat_halfeye(.width = c(0.80, 0.95), alpha = 0.6, adjust = 1.5) +
  
  # 95% interval lines (thin)
  geom_segment(
    data = summaries,
    aes(x = low95, xend = high95, y = Species, yend = Species, color = Species),
    linewidth = 1
  ) +
  # 80% interval lines (thick)
  geom_segment(
    data = summaries,
    aes(x = low80, xend = high80, y = Species, yend = Species, color = Species),
    linewidth = 2.5
  ) +
  # Median points
  geom_point(
    data = summaries,
    aes(x = median, y = Species, color = Species),
    size = 4
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "blue", linewidth = 1) +
  
  scale_fill_manual(values = c("Black bear" = "#1b9e77", "Grizzly bear" = "#7A4988")) +
  scale_color_manual(values = c("Black bear" = "#1b9e77", "Grizzly bear" = "#7A4988")) +
  
  labs(
    x = "Effect of Trail Status: Open (log AAR)",
    y = NULL
  ) +
  theme_minimal(base_size = 22) + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    plot.title = element_text(size = 17, hjust = 0.5),
    axis.line = element_line(linewidth = 1),
    legend.position = "none"
  )

# --- Display plot ---
p_trailstatus

# --- Save high-resolution output ---
ggsave(
  filename = "Output/trail_status_effects.png",
  plot = p_trailstatus,
  width = 10,
  height = 7,
  dpi = 600
)

# =============================================================
# 3.12. Save outputs
# =============================================================

# Combine into one summary file
all_summaries <- list(
  independent_detections = indep_summary,
  aar1_intervals = aar1_summary,
  aar1_stats = summary_table
)

write.csv(indep_summary, "Output/independent_detections.csv", row.names = FALSE)
write.csv(aar1_summary, "Output/aar1_intervals.csv", row.names = FALSE)
write.csv(summary_table, "Output/brms_output/aar1_stats.csv", row.names = FALSE)

saveRDS(all_summaries, "Output/aar_summaries.rds")

summary(black_bear_aar1_model)
summary(grizzly_bear_aar1_model)



# =============================================================
# 3.13. Diagnostics for AAR models: PP checks, trace, Bayes p
# =============================================================
suppressPackageStartupMessages({
  library(bayesplot); library(brms); library(posterior)
  library(ggplot2); library(dplyr); library(readr)
})

# Ensure output folders exist
output_dir   <- "Output"
diag_dir     <- file.path(output_dir, "diagnostics")
pp_dir       <- file.path(diag_dir, "pp_checks")
trace_dir    <- file.path(diag_dir, "trace_plots")
dir.create(pp_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(trace_dir, showWarnings = FALSE, recursive = TRUE)

# Collect models
models <- list(
  "Black_bear"   = black_bear_aar1_model,
  "Grizzly_bear" = grizzly_bear_aar1_model
)

# Set a seed for reproducibility of posterior_predict subsampling
set.seed(1234)

# Container to accumulate Bayesian p-values
bayes_pvals <- tibble(
  Species   = character(),
  p_mean    = numeric(),
  p_sd      = numeric(),
  p_tailpos = numeric()  # Pr( proportion(y>0)_rep >= observed ) on log scale (0 == AAR=1)
)

for (sp in names(models)) {
  m <- models[[sp]]
  
  # -----------------------------------------
  # 1) Trace plots
  # -----------------------------------------
  # Use all unconstrained parameters (draws array)
  draws_arr <- as_draws_array(m)
  p_trace <- bayesplot::mcmc_trace(draws_arr, np = nuts_params(m), facet_args = list(ncol = 2))
  ggsave(filename = file.path(trace_dir, paste0("trace_", sp, ".png")),
         plot = p_trace, width = 12, height = 8, dpi = 300)
  
  # -----------------------------------------
  # 2) Posterior predictive checks
  # -----------------------------------------
  # Density overlay
  p_dens <- pp_check(m, type = "dens_overlay", ndraws = 100)
  ggsave(filename = file.path(pp_dir, paste0("ppc_dens_", sp, ".png")),
         plot = p_dens, width = 9, height = 6, dpi = 300)
  
  # Mean stat
  p_mean <- pp_check(m, type = "stat", stat = "mean", ndraws = 2000)
  ggsave(filename = file.path(pp_dir, paste0("ppc_mean_", sp, ".png")),
         plot = p_mean, width = 9, height = 6, dpi = 300)
  
  # SD stat
  p_sd <- pp_check(m, type = "stat", stat = "sd", ndraws = 2000)
  ggsave(filename = file.path(pp_dir, paste0("ppc_sd_", sp, ".png")),
         plot = p_sd, width = 9, height = 6, dpi = 300)
  
  # Optional: empirical CDF overlay (often informative for Gaussian responses)
  p_ecdf <- pp_check(m, type = "ecdf_overlay", ndraws = 200)
  ggsave(filename = file.path(pp_dir, paste0("ppc_ecdf_", sp, ".png")),
         plot = p_ecdf, width = 9, height = 6, dpi = 300)
  
  # -----------------------------------------
  # 3) Simple Bayesian p-values (discrepancy checks)
  #    Compare replicated vs observed for:
  #    - mean
  #    - sd
  #    - tail-positive proportion: mean(y > 0), noting y is log(AAR1) so 0 ≡ AAR=1
  # -----------------------------------------
  # Observed y and replicated draws
  y_obs  <- model.response(model.frame(m))           # equals log_AAR1 column used in fit
  yrep   <- posterior_predict(m, ndraws = 2000)      # draws x N
  # Stats on observed
  T_mean_obs    <- mean(y_obs, na.rm = TRUE)
  T_sd_obs      <- sd(y_obs, na.rm = TRUE)
  T_tailpos_obs <- mean(y_obs > 0, na.rm = TRUE)
  
  # Stats on replicated
  T_mean_rep    <- apply(yrep, 1, function(x) mean(x, na.rm = TRUE))
  T_sd_rep      <- apply(yrep, 1, function(x) sd(x, na.rm = TRUE))
  T_tailpos_rep <- apply(yrep, 1, function(x) mean(x > 0, na.rm = TRUE))
  
  # Two one-sided tail choices are common; here we report Pr(T_rep >= T_obs)
  pval_mean    <- mean(T_mean_rep    >= T_mean_obs)
  pval_sd      <- mean(T_sd_rep      >= T_sd_obs)
  pval_tailpos <- mean(T_tailpos_rep >= T_tailpos_obs)
  
  bayes_pvals <- bind_rows(
    bayes_pvals,
    tibble(
      Species   = gsub("_", " ", sp),
      p_mean    = pval_mean,
      p_sd      = pval_sd,
      p_tailpos = pval_tailpos
    )
  )
}

# Save p-values table
readr::write_csv(bayes_pvals, file.path(diag_dir, "bayesian_pvalues_aar_models.csv"))

# Print quick summary to console
bayes_pvals


make_trace <- function(model, species_label) {
  
  draws_df <- as_draws_df(model)
  
  pars_to_show <- c("b_Intercept", "b_trail_statusopen")
  
  p <- mcmc_trace(
    draws_df,
    pars = pars_to_show,
    facet_args = list(ncol = 1)  
  ) +
    ggtitle(species_label) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "none"
    )
  
  return(p)
}
p_black   <- make_trace(black_bear_aar1_model,   "Black bear")
p_grizzly <- make_trace(grizzly_bear_aar1_model, "Grizzly bear")
combined_plot <- p_black | p_grizzly
combined_plot





# ===============================================================
# Additional Figures
# ===============================================================

# Monthly Independent Detections (Figure 3A)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(stringr)
library(tibble)

# -------------------------------------------------
# 0) Load data
# -------------------------------------------------
obs <- read.csv("Data/BLT_30min_Independent.csv")
detections <- read.csv("Data/BLT_detection_data_urs_hom.csv", stringsAsFactors = FALSE)

# -------------------------------------------------
# 1) Bears (already 30‑min independent in `obs`)
# -------------------------------------------------
bear_df <- obs %>%
  filter(Species %in% c("Ursus americanus", "Ursus arctos")) %>%
  mutate(
    datetime    = ymd_hms(Date_Time.Captured, quiet = TRUE),
    month       = as.Date(floor_date(datetime, "month")),
    species_std = ifelse(Species == "Ursus americanus", "Black bear", "Grizzly bear")
  ) %>%
  select(month, species_std)

bear_monthly <- bear_df %>%
  count(species_std, month, name = "detections") %>%
  mutate(species_std = as.character(species_std))

# -------------------------------------------------
# 2) Humans → 1‑min independence + open/closed classification
# -------------------------------------------------
hum_1min <- detections %>%
  mutate(
    datetime_str = str_extract(Image.ID, "\\d{4}-\\d{2}-\\d{2}__\\d{2}-\\d{2}-\\d{2}"),
    datetime_str = str_replace(datetime_str, "__", " "),
    datetime_str = str_replace_all(datetime_str, "(\\d{2})-(\\d{2})-(\\d{2})$", "\\1:\\2:\\3"),
    Date_Time.Captured = ymd_hms(datetime_str, quiet = TRUE)
  ) %>%
  arrange(Deployment.Location.ID, Date_Time.Captured) %>%
  filter(Species == "Homo sapiens") %>%
  group_by(Deployment.Location.ID) %>%
  mutate(
    lag_time = lag(Date_Time.Captured),
    gap_sec  = as.numeric(difftime(Date_Time.Captured, lag_time, units = "secs")),
    is_independent = if_else(is.na(gap_sec) | gap_sec >= 60, TRUE, FALSE)
  ) %>%
  ungroup() %>%
  filter(is_independent) %>%
  mutate(
    week_start = floor_date(Date_Time.Captured, unit = "week", week_start = 1),
    Week_Label = paste0(isoyear(week_start), "-W", sprintf("%02d", isoweek(week_start)))
  ) %>%
  mutate(
    trail_status = case_when(
      Deployment.Location.ID %in% c("BLT01","BLT03","BLT04","BLT07","BLT10","BLT14B") ~ "open",
      Deployment.Location.ID == "BLT16" & Week_Label <= "2024-W27" ~ NA_character_,
      Deployment.Location.ID == "BLT16" & Week_Label  > "2024-W27" ~ "closed",
      Deployment.Location.ID %in% c("BLT13","BLT15") &
        Week_Label >= "2024-W24" & Week_Label <= "2024-W29" ~ "open",
      Deployment.Location.ID %in% c("BLT13","BLT15") ~ "closed",
      Deployment.Location.ID %in% c("BLT17","BLT19","BLT20","BLT21","BLT22","BLT23",
                                    "BLT25","BLT26","BLT28","BLT29","BLT31","BLT33",
                                    "BLT34","BLT36","BLT40","BLT42","BLT44") ~ "closed",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(trail_status)) %>%
  mutate(
    month       = as.Date(floor_date(Date_Time.Captured, "month")),
    species_std = if_else(trail_status == "open", "Humans (Open)", "Humans (Closed)")
  ) %>%
  select(month, species_std)


hum_monthly <- hum_1min %>%
  count(species_std, month, name = "detections") %>%
  mutate(species_std = as.character(species_std))

# -------------------------------------------------
# 3) Combine bears + humans, fill missing months
# -------------------------------------------------
x_start <- as.Date("2023-08-01")
x_end   <- as.Date("2025-06-24")

all_monthly <- bind_rows(bear_monthly, hum_monthly) %>%
  group_by(species_std) %>%
  complete(
    month = seq(from = x_start, to = x_end, by = "month"),
    fill = list(detections = 0)
  ) %>%
  ungroup() %>%
  mutate(
    species_std = factor(
      species_std,
      levels = c("Black bear", "Grizzly bear", "Humans (Open)", "Humans (Closed)")
    )
  )


# -------------------------------------------------
# 4) Define months to keep (seasonal windows)
# -------------------------------------------------
keep_months <- bind_rows(
  tibble(start = as.Date("2023-08-01"), end = as.Date("2023-10-31")),
  tibble(start = as.Date("2024-05-01"), end = as.Date("2024-10-31")),
  tibble(start = as.Date("2025-05-01"), end = as.Date("2025-06-24"))
) %>%
  rowwise() %>%
  mutate(month = list(seq(from = start, to = end, by = "month"))) %>%
  ungroup() %>%
  tidyr::unnest(month) %>%
  distinct(month) %>%
  arrange(month) %>%
  mutate(month_idx = row_number())

# -------------------------------------------------
# 5) Build gap‑aware monthly dataset
# -------------------------------------------------
plot_monthly <- all_monthly %>%
  inner_join(keep_months, by = "month") %>%
  group_by(species_std) %>%
  arrange(month) %>%
  mutate(
    gap_days = as.numeric(difftime(month, lag(month), units = "days")),
    new_seg  = if_else(is.na(gap_days) | gap_days <= 45, 0L, 1L),
    seg_id   = cumsum(new_seg) + 1L
  ) %>%
  ungroup()

# -------------------------------------------------
# 6) Build dashed connectors across gaps
# -------------------------------------------------
connector_df <- plot_monthly %>%
  group_by(species_std) %>%
  arrange(month) %>%
  mutate(
    next_month     = lead(month),
    next_month_idx = lead(month_idx),
    next_det       = lead(detections),
    gap_days_next  = as.numeric(difftime(next_month, month, units = "days"))
  ) %>%
  filter(!is.na(gap_days_next) & gap_days_next > 45) %>%
  transmute(
    species_std,
    x    = month_idx,
    xend = next_month_idx,
    y    = detections,
    yend = next_det
  ) %>%
  ungroup()

# -------------------------------------------------
# 7) Plot
# -------------------------------------------------
p2 <- ggplot() +
  geom_line(
    data = plot_monthly,
    aes(x = month_idx, y = detections, colour = species_std,
        group = interaction(species_std, seg_id)),
    linewidth = 1
  ) +
  geom_segment(
    data = connector_df,
    aes(x = x, xend = xend, y = y, yend = yend, colour = species_std),
    linewidth = 1,
    linetype = "dashed"
  ) +
  geom_point(
    data = plot_monthly,
    aes(x = month_idx, y = detections, colour = species_std),
    size = 2
  ) +
  facet_wrap(~ species_std, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(
    "Black bear"      = "#3B9B3B",
    "Grizzly bear"    = "#7A4988",
    "Humans (Open)"   = "black",
    "Humans (Closed)" = "#4D4D4D"
  )) +
  scale_x_continuous(
    breaks = keep_months$month_idx,
    labels = format(keep_months$month, "%b\n%Y"),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(x = "Month", y = "Independent Detections") +
  theme_classic(base_size = 18) +
  theme(
    panel.spacing.y = unit(10, "pt"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.title  = element_text(size = 16),
    strip.text  = element_text(face = "bold", size = 16),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 25, l = 5)
  )

print(p2)

ggsave(
  filename = "Additional_Figures/detection_monthly_lineplot.png",
  plot = p2, width = 9, height = 10, dpi = 300
)



# Figure 3B spatial detections

library(tidyverse)
library(sf)
library(ggspatial)
library(lubridate)
library(patchwork)

# ------------------------------------------------------------------
# Stations (clean IDs + coords)  -> WGS84 then transform later
# ------------------------------------------------------------------
sta <- read.csv("Data/stations.csv")

if (!"station_id" %in% names(sta)) {
  if ("Deployment.Location.ID" %in% names(sta)) {
    sta <- sta %>% rename(station_id = Deployment.Location.ID)
  } else if ("Deployment.Location" %in% names(sta)) {
    sta <- sta %>% rename(station_id = Deployment.Location)
  } else if ("Site" %in% names(sta)) {
    sta <- sta %>% rename(station_id = Site)
  } else {
    stop("No station ID column found on `sta`.")
  }
}

sta <- sta %>%
  mutate(
    station_id = as.character(trimws(station_id)),
    Longitude  = as.numeric(longitude),
    Latitude   = as.numeric(latitude)
  )

sta_sf <- st_as_sf(sta, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

# ------------------------------------------------------------------
# Trail layer (projected)
# ------------------------------------------------------------------
trail <- st_read("Data/GIS_data/Full_BLT_Trail.gpkg", quiet = TRUE) %>%
  st_transform(3857)

# ------------------------------------------------------------------
# On/Off-trail classification
# ------------------------------------------------------------------
on_trail <- c(
  "BLT01","BLT03","BLT04","BLT07","BLT10",
  "BLT13","BLT15","BLT16","BLT17","BLT19",
  "BLT14B","BLT20","BLT21","BLT22","BLT23",
  "BLT25","BLT26","BLT28","BLT29","BLT31",
  "BLT34","BLT33","BLT36","BLT40","BLT42","BLT44"
)

trail_status_df <- tibble(
  station_id   = sta$station_id,
  trail_status = if_else(station_id %in% on_trail, "On-trail", "Off-trail")
)

# ------------------------------------------------------------------
# Weekly detections -> total detections per 100 trap-days
# ------------------------------------------------------------------
bear_model_data <- read.csv("Data/combined_bear_model_data.csv") %>%
  rename(station_id = Deployment.Location.ID) %>%
  mutate(
    station_id = as.character(trimws(station_id)),
    Week_Start = as.Date(Week_Start)
  ) %>%
  filter(Species %in% c("Ursus americanus", "Ursus arctos"))

bear_rate <- bear_model_data %>%
  group_by(station_id, Species) %>%
  summarise(
    total_detections = sum(Bear_Detections, na.rm = TRUE),
    total_effort     = sum(Effort,           na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    detections_per100 = if_else(
      is.finite(total_effort) & total_effort > 0,
      (total_detections / total_effort) * 100,
      NA_real_
    )
  )

# ------------------------------------------------------------------
# Build plot data
#   pts_base: one row per station (for station markers)
#   pts_rate: joins in bear_rate (two rows per station: one per species)
# ------------------------------------------------------------------
pts_base <- sta_sf %>%
  st_transform(3857) %>%
  mutate(
    x = st_coordinates(.)[, 1],
    y = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  transmute(
    station_id = as.character(station_id),
    x, y
  ) %>%
  left_join(trail_status_df, by = "station_id")

pts_rate <- pts_base %>%
  left_join(bear_rate, by = "station_id")

# ------------------------------------------------------------------
# Bounding box (pad 5 km) based on station locations only
# ------------------------------------------------------------------
pts_sf_3857 <- st_as_sf(pts_base, coords = c("x", "y"), crs = 3857, remove = FALSE)
bb <- st_bbox(st_buffer(st_union(pts_sf_3857), 1000))
xlim_use <- c(bb["xmin"], bb["xmax"])
ylim_use <- c(bb["ymin"], bb["ymax"])

# ------------------------------------------------------------------
# Shared size scale (prevents small positive values from disappearing)
# ------------------------------------------------------------------
all_vals <- pts_rate$detections_per100
all_vals <- all_vals[is.finite(all_vals) & all_vals > 0]

global_max <- ifelse(length(all_vals) == 0, 1, max(all_vals))

size_limits <- c(0, global_max)

size_breaks <- c(0.5, 1, 2, 4, 8)
size_breaks <- size_breaks[size_breaks <= global_max]

if (global_max < 1) {
  size_breaks <- global_max
}


# ------------------------------------------------------------------
# Map styling
# ------------------------------------------------------------------
col_black <- "#1b9e77"
col_grizz <- "#7A4988"

station_fill <- c("On-trail" = "black", "Off-trail" = "grey60")

base_map <- function() {
  ggplot() +
    annotation_map_tile(type = "osm", zoomin = -1) +
    geom_sf(data = trail, color = "grey40", size = 0.6) +
    coord_sf(crs = 3857, xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text  = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.box = "vertical"
    )
}


# ------------------------------------------------------------------
# Black bear map
# ------------------------------------------------------------------
p_black <- base_map() +
  
  # 1) detections first 
  geom_point(
    data = dplyr::filter(pts_rate, Species == "Ursus americanus",
                         is.finite(detections_per100), detections_per100 > 0),
    aes(x = x, y = y, size = detections_per100),
    color = col_black, shape = 16, alpha = 0.75
  ) +
  
  # 2) cameras second 
  geom_point(
    data = pts_base,
    aes(x = x, y = y, fill = trail_status),
    shape = 21, color = "black", size = 2.0, stroke = 0.5
  ) +
  
  scale_fill_manual(
    values = c("On-trail" = "black", "Off-trail" = "grey60"),
    name   = "Camera location"
  ) +
  
  scale_size_continuous(
    name   = "Detections per 100 trap-days",
    limits = c(0, global_max),
    breaks = size_breaks[size_breaks > 0],
    labels = function(x) ifelse(x < 1, "<1", x),
    range  = c(4, 14),
    guide  = guide_legend()
  )



ggsave(
  "Additional_Figures/dot_map_black_bears_per100days.png",
  plot = p_black, width = 7, height = 5, dpi = 300
)

# ------------------------------------------------------------------
# Grizzly bear map
# ------------------------------------------------------------------
p_grizz <- base_map() +
  
  # 1) Grizzly detections  
  geom_point(
    data = dplyr::filter(
      pts_rate,
      Species == "Ursus arctos",
      is.finite(detections_per100),
      detections_per100 > 0
    ),
    aes(x = x, y = y, size = detections_per100),
    color = col_grizz,
    shape = 16,
    alpha = 0.75
  ) +
  
  # 2) Camera locations
  geom_point(
    data = pts_base,
    aes(x = x, y = y, fill = trail_status),
    shape  = 21,
    color  = "black",
    size   = 2.0,
    stroke = 0.5
  ) +
  
  scale_fill_manual(
    values = c("On-trail" = "black", "Off-trail" = "grey60"),
    name   = "Camera location"
  ) +
  
  scale_size_continuous(
    name   = "Detections per 100 trap-days",
    limits = c(0, global_max),
    breaks = size_breaks[size_breaks > 0],
    labels = function(x) ifelse(x < 1, "<1", x),
    range  = c(4, 14),
    guide  = guide_legend()
  )


ggsave(
  "Additional_Figures/dot_map_grizzly_bears_per100days.png",
  plot = p_grizz, width = 7, height = 5, dpi = 300
)

# ------------------------------------------------------------------
# Combined figure (shared legend)
# ------------------------------------------------------------------
combined_maps <- p_black + p_grizz +
  plot_layout(ncol = 1, guides = "collect") +
  plot_annotation(theme = theme(legend.position = "right"))

ggsave(
  "Additional_Figures/dot_map_bears_combined.png",
  plot   = combined_maps,
  width  = 8,
  height = 10,
  dpi    = 300
)

combined_maps

