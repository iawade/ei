# SMARCA4 Related

library(data.table)
library(tidytable)
library(tidyverse)

#####

# Map out different variant classes
# LOF/PTV ...42; clinvar 1*; clinvar 1*/2* ...49 (cldn), ...54 review, ...55 clinsig, ...58 var type (may as well have); vars w/ freq 0.01% - take from VCF

# Map out different inds
# Males
# Females
# then each by age
# then females by prev cancer status and age

variant_data <- fread("SMARCA4_ukb_variants_paricipants_filtered_normalised_annotated.tsv.gz") %>%
    select(-FILTER)

print("successfully read in")

# Split the INFO column into a list
variant_data <- variant_data %>%
    mutate(INFO_split = str_split(INFO, "\\|"))

print("successfully performed one column operation")

# Use unnest_wider to create columns dynamically
variant_data <- variant_data %>%
    tidytable::unnest_wider(INFO_split) %>%
    rename(
        loftee = `...42`,
        clndn = `...49`,
        review_status = `...54`,
        clinsig = `...55`,
        clinvar_var_type = `...58`,
        hgvs_c = `...11`,
        mane_transcript = `...24`,
        gene = `...4`,
        vep_var_type = `...2`,
        INFO_1 = `...1`
    ) %>%
    select(-contains("..."))

print("successfully split vars")

# Rename further
variant_data <- variant_data %>%
    mutate(INFO_1_split = str_split(INFO_1, ";")) %>%
    tidytable::unnest_wider(INFO_1_split) %>%
    rename(
        ukb_af_total = `...1`,
        ukb_ac_eur_unrelated = `...3`,
        ukb_an_eur_unrelated = `...4`
    ) %>%
    select(-contains("..."), -INFO, -QUAL, -INFO_1)

# More Columns
## AF in filtered subset
## clinvar review status
## clinvar disease name
variant_data <- variant_data %>%
    mutate(across(starts_with("ukb_"), ~str_extract(.x, "(?<==)\\S+"))) %>%
    mutate(ukb_af_eur_unrelated = as.numeric(ukb_ac_eur_unrelated) / as.numeric(ukb_an_eur_unrelated))

# This needs to be as in input; break in pipeline? Requires some kind of input and only works post annotation?
# Can probably do upstream on unix with awk
relevant_clndn <- c("Small_cell_carcinoma_of_the_ovary", "Hereditary_cancer-predisposing_syndrome",
                    "Rhabdoid_tumor_predisposition_syndrome_2",
                    "bilateral_breast_cancer", "Malignant_tumor_of_breast", "Neoplasm_of_ovary",
                    "_hypercalcemic_type")
# Star status
# no assertion = no stars
# 2 star or greater: "multiple", "expert" "practice"
variant_data <- variant_data %>%
    mutate(clinvar_star = case_when(
        grepl("assertion", review_status, ignore.case = TRUE) ~ 0,
        grepl("multiple", review_status, ignore.case = TRUE) ~ 2,
        grepl("expert", review_status, ignore.case = TRUE) ~ 3,
        grepl("practice", review_status, ignore.case = TRUE) ~ 4,
        is.na(review_status) | review_status == "" ~ NA_real_,
        TRUE ~ 1
    ))

# Give a label for each tranche of variants
# PTV
# PTV + clinvar 1
# PTV + clinvar 2
# Rare (need to think about excluding benign variants)
# Four new columns
variant_data <- variant_data %>%
    mutate(
        ptv_variants = as.integer(loftee == "HC"),
        ptv_clinvar_1 = as.integer(loftee == "HC" | (any(str_detect(clndn, regex(paste(relevant_clndn, collapse = "|"), ignore_case = TRUE))) & (grepl("Pathogenic", clinsig) | grepl("Likely_pathogenic", clinsig)) & clinvar_star > 0)),
        ptv_clinvar_2 = as.integer(loftee == "HC" | (any(str_detect(clndn, regex(paste(relevant_clndn, collapse = "|"), ignore_case = TRUE))) & (grepl("Pathogenic", clinsig) | grepl("Likely_pathogenic", clinsig)) & clinvar_star > 1)),
        rare = as.integer(ukb_af_eur_unrelated < 1e-04 & !grepl("Benign|Likely_benign", clinsig))
    )

print("variant_data successfully generated")

# Separate out variant tranches
# Light QC - remove vars with missingness > 5%
ptv_variants <- variant_data %>%
    filter(ptv_variants == 1) %>%
    select(-ptv_variants,-ptv_clinvar_1,-ptv_clinvar_2,-rare) %>%
    filter(ukb_an_eur_unrelated > 769294)

ptv_clinvar_1 <- variant_data %>%
    filter(ptv_clinvar_1 == 1) %>%
    select(-ptv_variants,-ptv_clinvar_1,-ptv_clinvar_2,-rare) %>%
    filter(ukb_an_eur_unrelated > 769294)

ptv_clinvar_2 <- variant_data %>%
    filter(ptv_clinvar_2 == 1) %>%
    select(-ptv_variants,-ptv_clinvar_1,-ptv_clinvar_2,-rare) %>%
    filter(ukb_an_eur_unrelated > 769294)

rare <- variant_data %>%
    filter(rare == 1) %>%
    select(-ptv_variants,-ptv_clinvar_1,-ptv_clinvar_2,-rare) %>%
    filter(ukb_an_eur_unrelated > 769294)
# Worth also saving these as objects with and without the genotypes for review

print("variants split out")

# QC for frequencies - need to do something with missingness as well as depth at each site
# some kind of bug in my test data - not present in any of vcf's on HPC?
# must not have most up to date one

# @ some point use phenotyping_data id's as columns to go across for freq's

#####

phenotyping_data <- fread("../../reference_data/European_Ancestry_Unrelated_All_UKBB_Participants_404k.tsv")
# Need age and need cancer status (reported, registry or otherwise); need sex
# For cancer data good to know if it's registry or self-reported
# ICD 10 = C56*; 183.0 in ICD 9

# Data dictionary
data_dictionary <- fread("../../reference_data/Data_Dictionary_Showcase.tsv")


# Calculate age
phenotyping_data <- phenotyping_data %>%
    mutate(
        birthdate = make_date(p34, month = match(p52, month.name)),
        current_age = floor(as.numeric(difftime(Sys.Date(), birthdate, units = "days")) / 365.25)
    ) %>%
    select(-birthdate)


phenotyping_data <- phenotyping_data %>%
    unite("p20001_collapsed", starts_with("p20001"), sep = ";", na.rm = TRUE) %>%
    unite("p2453_collapsed", starts_with("p2453"), sep = ";", na.rm = TRUE) %>%
    unite("p84_collapsed", starts_with("p84"), sep = ";", na.rm = TRUE) %>%
    unite("p40021_collapsed", starts_with("p40021"), sep = ";", na.rm = TRUE) %>%
    unite("p40008_collapsed", starts_with("p40008"), sep = ";", na.rm = TRUE) %>%
    unite("p40006_collapsed", starts_with("p40006"), sep = ";", na.rm = TRUE) %>%
    unite("p40013_collapsed", starts_with("p40013"), sep = ";", na.rm = TRUE)

phenotyping_data_updated_cols <- as.data.table(phenotyping_data)

# Update column names using data dictionary
setnames(
    phenotyping_data_updated_cols,
    sapply(names(phenotyping_data_updated_cols), function(col) {
        if (startsWith(col, "p")) {
            field_id <- parse_number(col)
            matching_row <- data_dictionary$FieldID == field_id
            if (any(matching_row)) {
                return(data_dictionary$Field[matching_row])
            }
        }
        return(col)
    })
)

phenotyping_data_updated_cols <- phenotyping_data_updated_cols %>%
    select(-`Year of birth`, -`Month of birth`) %>%
    mutate(eid = as.character(eid))

# Notes to wheel back to - lost to follow up? ~30 ind's have deaths reported - is that relevant - probably not in this case

# MVP frequencies and frequencies up to each age group
# MVP ignore homozygous calls - shouldn't be any - if there are it's likely an artefact for SMARCA4
# QC remove vars where missingness is high
# 404892 - 809784 alleles ; missingness of 5% is < 769295
# Exon 27 and 30 issue

print("phenotyping data generated")

#####
# Frequencies
# Male/Female overalland in ten year age gaps

# New table: transpose/pivot so have person: PTV carrier, rare var carrier, age and sex
numeric_col_names <- names(ptv_variants)[grepl("^\\d+$", names(ptv_variants))]

# Create a new data.table with numeric column names as rows
ind_carrier_status <- data.table(participant = numeric_col_names)

# Check if any row in ptv_variants starts with "0/1" only for specified columns
ind_carrier_status[, ptv_carrier := apply(ptv_variants[, .SD, .SDcols = numeric_col_names], 2, function(x) any(grepl("^0/1", x)))]
ind_carrier_status[, rare_var_carrier := apply(rare[, .SD, .SDcols = numeric_col_names], 2, function(x) any(grepl("^0/1", x)))]

ind_carrier_status <- ind_carrier_status %>%
    left_join.(select(phenotyping_data_updated_cols, eid, current_age, Sex), by = c("participant" = "eid"))

# Convert the Sex column to a factor for proper ordering
ind_carrier_status$Sex <- factor(ind_carrier_status$Sex, levels = c("Male", "Female"))

# Create age bands with more descriptive labels
age_labels <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-100")
ind_carrier_status <- ind_carrier_status %>%
    mutate(age_band = cut(current_age, breaks = seq(0, 100, by = 10), labels = age_labels, include.lowest = TRUE))

# Calculate the total number of rows in ind_carrier_status
total_rows <- nrow(ind_carrier_status)

# Calculate frequencies overall
overall_freq <- ind_carrier_status %>%
    group_by(Sex) %>%
    summarise(
        ptv_carrier_freq = sum(ptv_carrier) / total_rows,
        rare_var_carrier_freq = sum(rare_var_carrier) / total_rows
    )

# Calculate frequencies in ten-year age bands with more descriptive columns
# Should they be cumulative?
age_band_freq <- ind_carrier_status %>%
    group_by(Sex, age_band) %>%
    summarise(
        ptv_carrier_freq_women = sum(ptv_carrier[Sex == "Female"]) / sum(Sex == "Female"),
        ptv_carrier_freq_men = sum(ptv_carrier[Sex == "Male"]) / sum(Sex == "Male"),
        rare_var_carrier_freq_women = sum(rare_var_carrier[Sex == "Female"]) / sum(Sex == "Female"),
        rare_var_carrier_freq_men = sum(rare_var_carrier[Sex == "Male"]) / sum(Sex == "Male")
    ) %>%
    pivot_wider(names_from = Sex, values_from = c(ptv_carrier_freq_women, ptv_carrier_freq_men, rare_var_carrier_freq_women, rare_var_carrier_freq_men)) %>%
    select(-ptv_carrier_freq_women_Male, -ptv_carrier_freq_men_Female, -rare_var_carrier_freq_women_Male, -rare_var_carrier_freq_men_Female) %>%
    rename(
        ptv_carrier_freq_Male = ptv_carrier_freq_men_Male,
        rare_var_carrier_freq__Female = rare_var_carrier_freq_women_Female,
        rare_var_carrier_freq_Male = rare_var_carrier_freq_men_Male,
        ptv_carrier_freq_Female = ptv_carrier_freq_women_Female
    )

# Print the results
print("Overall Frequencies:")
print(overall_freq)

print("\nFrequencies in Ten-Year Age Bands:")
print(age_band_freq)

fwrite(ptv_variants, "ptv_variants.tsv", sep = "\t")
fwrite(rare, "rare_variants.tsv", sep = "\t")
fwrite(phenotyping_data_updated_cols, "ukb_phenotyping_data.tsv", sep = "\t")
fwrite(overall_freq, "overall_frequency.tsv", sep = "\t")
fwrite(age_band_freq, "frequency_over_ages.tsv", sep = "\t")

# Same but by ovarian cancer status
# What about all cancer status?
