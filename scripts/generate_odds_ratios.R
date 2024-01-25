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
                           labels = c("40-49", "50-59", "60-69", "70-79", "80-89"),
                           include.lowest = TRUE, right = FALSE))

# Transform and output variant carriers
variant_carrier_output <- tally_phenotype_union %>%
    select(-age_range) %>%
    filter(non_homref_count > 0)

fwrite(variant_carrier_output, OUTPUT_CARRIERS, sep = "\t")

# Calculate odds ratios
## TODO Fisher's test
## Could use functions so DNR code
or_table <- tally_phenotype_union %>%
    group_by(age_range) %>%
    summarise(
        total_females = sum(sex == "Female"),
        female_variant_carriers = sum(sex == "Female" & non_homref_count > 0),
        total_males = sum(sex == "Male"),
        male_variant_carriers = sum(sex == "Male" & non_homref_count > 0)
    ) %>%
    mutate(
        female_pv_odds = female_variant_carriers /
            ( total_females - female_variant_carriers ),
        male_pv_odds = male_variant_carriers /
            ( total_males - male_variant_carriers ),
        or_male_to_female = case_when(
            female_variant_carriers > 0 & male_variant_carriers > 0
            ~ male_pv_odds / female_pv_odds,
            ## Haldane correction if zeros
            TRUE
            ~ ((male_variant_carriers + 0.5) /
                   ((total_males + 0.5) - (male_variant_carriers + 0.5))) /
                ((female_variant_carriers + 0.5) /
                     ((total_females + 0.5) - (female_variant_carriers + 0.5)))
        )
    ) %>%
    ## Confidence intervals, including Haldane correction if needed
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

    ## Test significance
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
    ungroup() %>%
    ## Repeat for cumulative totals and odds
    arrange(age_range) %>%
    mutate(
        cumulative_total_females = cumsum(total_females),
        cumulative_female_variant_carriers = cumsum(female_variant_carriers),
        cumulative_total_males = cumsum(total_males),
        cumulative_male_variant_carriers = cumsum(male_variant_carriers)
    ) %>%
    mutate(
        cumulative_female_pv_odds = cumulative_female_variant_carriers /
            ( cumulative_total_females - cumulative_female_variant_carriers ),
        cumulative_male_pv_odds = cumulative_male_variant_carriers /
            ( cumulative_total_males - cumulative_male_variant_carriers ),
        cumulative_or_male_to_female = case_when(
            cumulative_female_variant_carriers > 0 & cumulative_male_variant_carriers > 0
            ~ cumulative_male_pv_odds / cumulative_female_pv_odds,
            TRUE
            ~ ((cumulative_male_variant_carriers + 0.5) /
                   ((cumulative_total_males + 0.5) - (cumulative_male_variant_carriers + 0.5))) /
                ((cumulative_female_variant_carriers + 0.5) /
                     ((cumulative_total_females + 0.5) - (cumulative_female_variant_carriers + 0.5)))
        )
    ) %>%

    ## Repeat CI's
    mutate(
        cumulative_lower_95_ci = case_when(
            cumulative_female_variant_carriers > 0 & cumulative_male_variant_carriers > 0
            ~ exp(log(cumulative_or_male_to_female) - qnorm(0.975) * sqrt(
                1/cumulative_male_variant_carriers + 1/(cumulative_total_males - cumulative_male_variant_carriers) +
                    1/cumulative_female_variant_carriers + 1/(cumulative_total_females - cumulative_female_variant_carriers)
            )),
            TRUE
            ~ exp(log(cumulative_or_male_to_female) - qnorm(0.975) * sqrt(
                1/(cumulative_male_variant_carriers + 0.5) + 1/((cumulative_total_males - cumulative_male_variant_carriers) + 0.5) +
                    1/(cumulative_female_variant_carriers + 0.5) + 1/((cumulative_total_females - cumulative_female_variant_carriers) + 0.5)
            ))
        ),
        cumulative_upper_95_ci = case_when(
            cumulative_female_variant_carriers > 0 & cumulative_male_variant_carriers > 0
            ~ exp(log(cumulative_or_male_to_female) + qnorm(0.975) * sqrt(
                1/cumulative_male_variant_carriers + 1/(cumulative_total_males - cumulative_male_variant_carriers) +
                    1/cumulative_female_variant_carriers + 1/(cumulative_total_females - cumulative_female_variant_carriers)
            )),
            TRUE
            ~ exp(log(cumulative_or_male_to_female) + qnorm(0.975) * sqrt(
                1/(cumulative_male_variant_carriers + 0.5) + 1/((cumulative_total_males - cumulative_male_variant_carriers) + 0.5) +
                    1/(cumulative_female_variant_carriers + 0.5) + 1/((cumulative_total_females - cumulative_female_variant_carriers) + 0.5)
            ))
        )
    ) %>%
    ## Repeat test for significance
    rowwise() %>%
    mutate(
        cumulative_fisher_p_value = case_when(
            cumulative_female_variant_carriers > 0 & cumulative_male_variant_carriers > 0
            ~ fisher.test(matrix(c(cumulative_female_variant_carriers, cumulative_total_females - cumulative_female_variant_carriers,
                                   cumulative_male_variant_carriers, cumulative_total_males - cumulative_male_variant_carriers), nrow = 2))$p.value,
            TRUE
            ~ fisher.test(matrix(c(cumulative_female_variant_carriers + 0.5, cumulative_total_females - cumulative_female_variant_carriers + 0.5,
                                   cumulative_male_variant_carriers + 0.5, cumulative_total_males - cumulative_male_variant_carriers + 0.5), nrow = 2))$p.value
        )
    )


# Output Data
fwrite(or_table, OUTPUT_ODDS_RATIOS, sep = "\t")
