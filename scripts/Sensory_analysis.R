## Package names
required_packages <- c("tidyverse", "janitor", "readxl", "broom", "car")

## Load packages
pacman::p_load(char = required_packages, install = FALSE)

## Load raw data
raw_sensory_data <- read_xlsx("data/sensory_22.xlsx") |> 
  clean_names() 

## Clean data and format
# Make a vector of sensory characteristic names
trait_names <- names(raw_sensory_data) |> 
  keep(~ str_detect(.x, "(ar|fl|tx|at)$"))

# Use these names to select important columns and change variables to factor
sensory_data <- raw_sensory_data |>
  select(subject_code, day, repetition, sample_id, all_of(trait_names)) |>
  pivot_longer(all_of(trait_names),
               names_to = "attribute",
               values_to = "score",
               values_drop_na = TRUE) |>
  filter(!str_detect(attribute, "other")) |> 
  mutate(across(subject_code:attribute, ~ as.factor(.x))) |> 
  # Discard first two replicates
  filter(!repetition %in% c("1","2"))

## Descriptive statistics
# By sample_id and attribute
desc_stats <- sensory_data %>%
  group_by(sample_id, attribute) %>%
  summarise(
    n    = n(),
    min  = min(score, na.rm = TRUE),
    max  = max(score, na.rm = TRUE),
    mean = mean(score, na.rm = TRUE),
    sd   = sd(score, na.rm = TRUE),
    cv   = (sd / mean) * 100,          # coefficient of variation (%)
    sem  = sd / sqrt(n),               # standard error of the mean
    .groups = "drop"
  )

# By attribute overall
desc_stats_overall <- sensory_data %>%
  group_by(attribute) %>%
  summarise(
    n    = n(),
    min  = min(score, na.rm = TRUE),
    max  = max(score, na.rm = TRUE),
    mean = mean(score, na.rm = TRUE),
    sd   = sd(score, na.rm = TRUE),
    cv   = (sd / mean) * 100,
    sem  = sd / sqrt(n),
    .groups = "drop"
  )

## discrimination for 1 panelist
asd <- aov(score ~ attribute, data = sensory_data |> filter(subject_code == "1"))
summary(asd)
zxc <- TukeyHSD(asd)
sum(zxc$attribute[,4] <0.05)

# can't work out repeatability
sdf <- aov(score ~ repetition*sample_id, data = sensory_data |> filter(subject_code == "1"))
summary(sdf)
zxc <- TukeyHSD(sdf)
sum(zxc$sample_id[,4] >0.05)

# no interaction needs to calculate the average score for the attribute across all panelists from the previous section and determine the number of attributes that they scored similarly to the average (p >0.05)


## one row in the table
# This gives the F-values for the effect of sample on aroma intensity scores
asd <- lmerTest::lmer(
  score ~ sample_id +                     # fixed effect
    (1 | subject_code) +                     # random panelists
    (1 | repetition)  +                     # random sessions
    (1 | subject_code:repetition),              # random interaction
  data = sensory_data |> filter(attribute == "aroma_intensity_ar"),
  REML = TRUE
)
anova(asd)

# This gives the F-values for the effect of panellist on aroma intensity scores
asd <- lmerTest::lmer(
  score ~ subject_code +                     # fixed effect
    (1 | sample_id) +                     # random
    (1 | repetition)  +                     # random
    (1 | sample_id:repetition),              # random interaction
  data = sensory_data |> filter(attribute == "aroma_intensity_ar"),
  REML = TRUE
)
anova(asd)

# This gives the F-values for the effect of panellist on aroma intensity scores
asd <- lmerTest::lmer(
  score ~ repetition +                     # fixed effect
    (1 | sample_id) +                     # random 
    (1 | subject_code)  +                     # random 
    (1 | sample_id:subject_code),              # random interaction
  data = sensory_data |> filter(attribute == "aroma_intensity_ar"),
  REML = TRUE
)
anova(asd)



## ============================================================================
## PANELIST PERFORMANCE ANALYSIS
## ============================================================================
## For each panelist and attribute combination, run ANOVA to assess:
## 1. Discrimination power: Can they distinguish between samples? (p < 0.05 for sample)
## 2. Repeatability: Are they consistent across replicates? (p >= 0.05 for replicate)
## 3. No interaction: Do they agree with panel consensus? (p >= 0.05 for panelist:sample)

# Function to run ANOVA for a single panelist × attribute combination
run_panelist_anova <- function(df) {
  # Fit ANOVA model
    model <- aov(score ~ sample_id * repetition, data = df)
    anova_result <- broom::tidy(model)
    
    # Extract p-values for key effects
    result <- tibble(
      term = c("sample_id", "repetition", "sample_id:repetition"),
      p.value = c(
        anova_result$p.value[anova_result$term == "sample_id"][1],
        anova_result$p.value[anova_result$term == "repetition"][1],
        anova_result$p.value[anova_result$term == "sample_id:repetition"][1]
      )
    )
    return(result)
}

# Run ANOVA for each panelist × attribute combination
panelist_performance <- sensory_data %>%
  group_by(subject_code, attribute) %>%
  nest() %>%
  mutate(anova_results = map(data, run_panelist_anova)) %>%
  select(-data) %>%
  unnest(anova_results) %>%
  pivot_wider(names_from = term, values_from = p.value) %>%
  rename(
    p_sample = sample_id,
    p_replicate = repetition,
    p_interaction = `sample_id:repetition`
  )

# Calculate performance metrics for each panelist
panelist_summary <- panelist_performance %>%
  group_by(subject_code) %>%
  summarise(
    # Discrimination: significant sample effect (p < 0.05)
    discrimination = sum(p_sample < 0.05, na.rm = TRUE),
    # Repeatability: non-significant replicate effect (p >= 0.05)
    repeatability = sum(p_replicate >= 0.05, na.rm = TRUE),
    # No interaction: non-significant interaction (p >= 0.05)
    no_interaction = sum(p_interaction >= 0.05, na.rm = TRUE),
    # Total attributes tested
    total = n(),
    .groups = "drop"
  )

## ============================================================================
## TABLE 1: Panelist Performance Summary
## ============================================================================
table1 <- panelist_summary %>%
  mutate(panelist = as.numeric(as.character(subject_code))) %>%
  select(panelist, discrimination, repeatability, no_interaction, total) %>%
  arrange(panelist)

print("Table 1: Panelist Performance")
print(table1)

## ============================================================================
## ATTRIBUTE-LEVEL ANALYSIS
## ============================================================================
## For Table 2, we need to count across panelists for specific attributes

# Count significant effects across panelists for each attribute
attribute_summary <- panelist_performance %>%
  group_by(attribute) %>%
  summarise(
    # Sample: # panelists with significant sample effect
    sample = sum(p_sample < 0.05, na.rm = TRUE),
    # Panellist: total # panelists tested (or significant panelist effect)
    panellist = n(),
    # Replicate: # panelists with significant replicate effect
    replicate = sum(p_replicate < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

## ============================================================================
## TABLE 2: Attribute Summary for Selected Attributes
## ============================================================================
# Select specific attributes for Table 2
selected_attributes <- c("aroma_intensity_ar", "sweet_fruit_ar", "citrus_ar", "floral_ar")

table2 <- attribute_summary %>%
  filter(attribute %in% selected_attributes) %>%
  mutate(sensory_attribute = case_when(
    attribute == "aroma_intensity_ar" ~ "aroma intensity AR",
    attribute == "sweet_fruit_ar" ~ "sweet fruit AR",
    attribute == "citrus_ar" ~ "citrus AR",
    attribute == "floral_ar" ~ "floral AR"
  )) %>%
  select(sensory_attribute, sample, panellist, replicate)

print("Table 2: Attribute Summary for Selected Attributes")
print(table2)

## ============================================================================
## MIXED EFFECTS MODEL (for reference/additional analysis)
## ============================================================================
## This provides an alternative approach using mixed models with:
## - Fixed effect: sample_id (genotype)
## - Random effects: panelist, session, and their interactions

# Load additional package for mixed models
if (!require(lmerTest)) {
  install.packages("lmerTest")
  library(lmerTest)
}

# Function to fit mixed model per attribute
fit_mixed_per_attr <- function(df_attr) {
  tryCatch({
    lmer(
      score ~ sample_id +                      # fixed effect
        (1 | subject_code) +                   # random panelist effect
        (1 | repetition) +                     # random session/replicate effect
        (1 | subject_code:sample_id),          # random panelist × sample interaction
      data = df_attr,
      REML = TRUE
    )
  }, error = function(e) {
    return(NULL)
  })
}

# Fit mixed models for all attributes
mixed_models <- sensory_data %>%
  group_by(attribute) %>%
  nest() %>%
  mutate(model = map(data, fit_mixed_per_attr)) %>%
  filter(!map_lgl(model, is.null))

# Extract fixed effects (sample_id) significance
mixed_model_results <- mixed_models %>%
  mutate(
    anova_table = map(model, ~ {
      if (!is.null(.x)) {
        anova(.x) %>% 
          as.data.frame() %>%
          rownames_to_column("term") %>%
          as_tibble()
      } else {
        NULL
      }
    })
  ) %>%
  select(attribute, anova_table) %>%
  unnest(anova_table) %>%
  filter(term == "sample_id")

print("Mixed Model Results: Sample Effect Significance")
print(mixed_model_results %>% select(attribute, `F value`, `Pr(>F)`))

## ============================================================================
## OVERALL PANEL ANOVA (Fixed Effects)
## ============================================================================
## Traditional ANOVA approach treating all factors as fixed
## This tests overall panel discrimination ability

# Function to run overall ANOVA per attribute
run_overall_anova <- function(df) {
  tryCatch({
    model <- aov(score ~ sample_id + subject_code + repetition + 
                   sample_id:subject_code, data = df)
    broom::tidy(model)
  }, error = function(e) {
    return(NULL)
  })
}

# Run for all attributes
overall_anova_results <- sensory_data %>%
  group_by(attribute) %>%
  nest() %>%
  mutate(anova_result = map(data, run_overall_anova)) %>%
  select(-data) %>%
  unnest(anova_result)

# Summary of sample discrimination across all attributes
sample_discrimination_summary <- overall_anova_results %>%
  filter(term == "sample_id") %>%
  mutate(significant = p.value < 0.05) %>%
  select(attribute, df, statistic, p.value, significant)

print("Overall Panel Discrimination by Attribute")
print(sample_discrimination_summary)

## ============================================================================
## Export results
## ============================================================================
# Save tables as CSV
write_csv(table1, "output/table1_panelist_performance.csv")
write_csv(table2, "output/table2_attribute_summary.csv")
write_csv(desc_stats, "output/descriptive_stats_by_sample.csv")
write_csv(desc_stats_overall, "output/descriptive_stats_overall.csv")

# Save detailed results
write_csv(panelist_performance, "output/panelist_performance_detailed.csv")
write_csv(attribute_summary, "output/attribute_summary_all.csv")
write_csv(sample_discrimination_summary, "output/overall_discrimination.csv")

print("Analysis complete. Results saved to output folder.")




## Package names
required_packages <- c("tidyverse", "janitor", "lme4")

# ## Install packages
# pak::pak(required_packages)

## Load packages
pacman::p_load(char = required_packages, install = FALSE)
## Load raw data
raw_sensory_data <- read_xlsx("data/sensory_22.xlsx") |> 
  clean_names() 

## Clean data and format
# make a vector of sensory characteristic names
trait_names <- names(raw_sensory_data) |> 
  keep(~ str_detect(.x, "(ar|fl|tx|at)$"))
# Use these names to select important columns and change variables to factor
sensory_data <- raw_sensory_data |>
  select(subject_code, day, repetition, sample_id, all_of(trait_names)) |>
  pivot_longer(all_of(trait_names),
               names_to = "attribute",
               values_to = "score",
               values_drop_na = T) |>
  filter(!str_detect(attribute, "other")) |> 
  mutate(across(subject_code:attribute, ~ as.factor(.x))) |> 
  # Discard first two replicates
  filter(!repetition %in% c("1","2"))
  
## descriptive statistics
# by sample_id and attribute
desc_stats <- sensory_data %>%
  group_by(sample_id, attribute) %>%
  summarise(
    n    = n(),
    min  = min(score, na.rm = TRUE),
    max  = max(score, na.rm = TRUE),
    mean = mean(score, na.rm = TRUE),
    sd   = sd(score, na.rm = TRUE),
    cv   = (sd / mean) * 100,          # coefficient of variation (%)
    sem  = sd / sqrt(n),               # standard error of the mean
    .groups = "drop"
  )

# by attribute
desc_stats_overall <- sensory_data %>%
  group_by(attribute) %>%
  summarise(
    n    = n(),
    min  = min(score, na.rm = TRUE),
    max  = max(score, na.rm = TRUE),
    mean = mean(score, na.rm = TRUE),
    sd   = sd(score, na.rm = TRUE),
    cv   = (sd / mean) * 100,
    sem  = sd / sqrt(n),
    .groups = "drop"
  )

## Fixed effects ANOVA
# Function to fit ANOVA per attribute
fit_anova_per_attr <- function(df) {
  aov(
    score ~ sample_id + Error(subject_code/sample_id),
    data = df
  )
}

# Run for all attributes
anova_results <- sensory_data %>%
  group_split(attribute) %>%
  set_names(levels(sensory_data$attribute)) %>%
  map(~ fit_anova_per_attr(.x))

# Extract tidy tables
anova_tables <- anova_results %>%
  map(~ broom::tidy(.x$Anova))

## Mixed effects model
# Mixed model per attribute
fit_mixed_per_attr <- function(df_attr) {
  lmerTest::lmer(
    score ~ sample_id +                     # fixed effect
      (1 | assessor) +                     # random panelists
      (1 | session)  +                     # random sessions
      (1 | assessor:session),              # random interaction
    data = df_attr,
    REML = TRUE
  )
}

mixed_models <- sensory_data %>%
  group_split(attribute) %>%
  set_names(levels(sensory_data$attribute)) %>%
  map(~ fit_mixed_per_attr(.x))

# ANOVA table (Type III via Satterthwaite DF by default in lmerTest)
anova_mixed <- mixed_models %>%
  map(~ anova(.x))        # look at 'sample_id' p-values for discrimination power

# Effect size: partial eta^2 for sample_id
eta2_mixed <- mixed_models %>%
  map(~ effectsize::eta_squared(.x, partial = TRUE))

## Compare sample group
emm_tables <- mixed_models %>%
  map(~ emmeans(.x, ~ sample_id))           # estimated marginal means per sample_id

pairwise_tukey <- emm_tables %>%
  map(~ pairs(.x, adjust = "tukey"))       # Tukey-adjusted pairwise differences

# Compact Letter Display (grouping letters by Tukey)
cld_letters <- emm_tables %>%
  map(~ multcomp::cld(.x, Letters = letters, adjust = "tukey"))

## plot stats
plot_means_sem <- desc_stats %>%
  ggplot(aes(sample_id, mean, ymin = mean - sem, ymax = mean + sem, color = sample_id)) +
  geom_point(size = 2) +
  geom_errorbar(width = 0.15) +
  facet_wrap(~ attribute, scales = "free_y") +
  labs(
    title = "sample_id means ± SEM by attribute",
    x = "sample_id (Sample Name)",
    y = "Score"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

## model diagnostics
library(DHARMa)

diag_plots <- mixed_models %>%
  map(function(m) {
    sim <- DHARMa::simulateResiduals(m)
    plot(sim)  # QQ, residual vs. predicted, etc.
    performance::check_model(m)  # collinearity, normality, homoscedasticity
  })

