# Script to find Odds Ratios

# Load Packages
library(data.table)
library(tidytable)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

INPUT_TALLY <- args[1]
INPUT_PHENOTYPES <- args[2]
OUTPUT_CARRIERS <- args[3]
OUTPUT_ODDS_RATIOS <- args[4]

# Read in data
tally_data <- fread(INPUT_TALLY)
munged_phenotype_data <- fread(INPUT_PHENOTYPES)

# Add age ranges
tally_phenotype_union <- left_join(tally_data, munged_phenotype_data) %>%
    mutate(age_range = cut(censored_age, breaks = c(40, 50, 60, 70, 80, 90),
                           labels = c("40-49", "40-59", "40-69", "40-79", "40-89"),
                           include.lowest = TRUE, right = FALSE))

# Transform and output variant carriers
variant_carrier_output <- tally_phenotype_union %>%
    select(-age_range) %>%
    filter(non_homref_count > 0)

fwrite(variant_carrier_output, OUTPUT_CARRIERS, sep = "\t")

# Calculate odds ratios
## Code could be better
or_table <- tally_phenotype_union %>%

    ## Cumulative totals and odds
    ### Has to be a simpler way of summarising into cumulative totals
    group_by(age_range) %>%
    summarise(
        total_females = sum(sex == "Female"),
        female_variant_carriers = sum(sex == "Female" & non_homref_count > 0),
        total_males = sum(sex == "Male"),
        male_variant_carriers = sum(sex == "Male" & non_homref_count > 0)
    ) %>%
    ungroup() %>%
    arrange(age_range) %>%
    mutate(
        total_females = cumsum(total_females),
        female_variant_carriers = cumsum(female_variant_carriers),
        total_males = cumsum(total_males),
        male_variant_carriers = cumsum(male_variant_carriers)
    ) %>%

    mutate(
        female_pv_odds = female_variant_carriers /
            ( total_females - female_variant_carriers ),
        male_pv_odds = male_variant_carriers /
            ( total_males - male_variant_carriers ),
        or_male_to_female = case_when(
            female_variant_carriers > 0 & male_variant_carriers > 0
            ~ male_pv_odds / female_pv_odds,
            TRUE
            ~ ((male_variant_carriers + 0.5) /
                   ((total_males + 0.5) - (male_variant_carriers + 0.5))) /
                ((female_variant_carriers + 0.5) /
                     ((total_females + 0.5) - (female_variant_carriers + 0.5)))
        )
    ) %>%

    ## CI's
    mutate(
        lower_95_ci = case_when(
            female_variant_carriers > 0 & male_variant_carriers > 0
            ~ exp(log(or_male_to_female) - qnorm(0.975) * sqrt(
                1/male_variant_carriers + 1/(total_males - male_variant_carriers) +
                    1/female_variant_carriers + 1/(total_females - female_variant_carriers)
            )),
            TRUE
            ~ exp(log(or_male_to_female) - qnorm(0.975) * sqrt(
                1/(male_variant_carriers + 0.5) + 1/((total_males - male_variant_carriers) + 0.5) +
                    1/(female_variant_carriers + 0.5) + 1/((total_females - female_variant_carriers) + 0.5)
            ))
        ),
        upper_95_ci = case_when(
            female_variant_carriers > 0 & male_variant_carriers > 0
            ~ exp(log(or_male_to_female) + qnorm(0.975) * sqrt(
                1/male_variant_carriers + 1/(total_males - male_variant_carriers) +
                    1/female_variant_carriers + 1/(total_females - female_variant_carriers)
            )),
            TRUE
            ~ exp(log(or_male_to_female) + qnorm(0.975) * sqrt(
                1/(male_variant_carriers + 0.5) + 1/((total_males - male_variant_carriers) + 0.5) +
                    1/(female_variant_carriers + 0.5) + 1/((total_females - female_variant_carriers) + 0.5)
            ))
        )
    ) %>%
    ## Test for significance
    rowwise() %>%
    mutate(
        fisher_p_value = case_when(
            female_variant_carriers > 0 & male_variant_carriers > 0
            ~ fisher.test(matrix(c(female_variant_carriers, total_females - female_variant_carriers,
                                   male_variant_carriers, total_males - male_variant_carriers), nrow = 2))$p.value,
            TRUE
            ~ fisher.test(matrix(c(female_variant_carriers + 0.5, total_females - female_variant_carriers + 0.5,
                                   male_variant_carriers + 0.5, total_males - male_variant_carriers + 0.5), nrow = 2))$p.value
        )
    ) %>%
    mutate(
        chi_square_p_value = case_when(
            female_variant_carriers > 0 & male_variant_carriers > 0
            ~ chisq.test(matrix(c(female_variant_carriers, total_females - female_variant_carriers,
                                   male_variant_carriers, total_males - male_variant_carriers), nrow = 2))$p.value,
            TRUE
            ~ chisq.test(matrix(c(female_variant_carriers + 0.5, total_females - female_variant_carriers + 0.5,
                                   male_variant_carriers + 0.5, total_males - male_variant_carriers + 0.5), nrow = 2))$p.value
        )
    )

# Output Data
fwrite(or_table, OUTPUT_ODDS_RATIOS, sep = "\t")
