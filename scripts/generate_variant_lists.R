# Script to Find Rare and Pathogenic Variants for VEP annotated VCF

library(data.table)
library(tidytable)
library(tidyverse)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

INPUT <- args[1]
GENE_LIST <- args[2]
PTV_CLINVAR_2_OUTPUT <- args[3]

gene_list <- readLines(GENE_LIST)

# # Read in clinvar phenotypes as vector
# relevant_clndn <- readLines(CLNDN)

# Process input data
## Read in Data
variant_data <- fread(INPUT) %>%
    select(-FILTER) %>%

    ## Split INFO column of VCF and rename for final output
    mutate(INFO_split = str_split(INFO, "\\|")) %>%
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
    select(-contains("...")) %>%
    filter(gene %in% gene_list) %>%

    ## Obtain allele frequency in filtered population
    ### Not flexible - Only really for European ancestry and unrelated
    mutate(INFO_1_split = str_split(INFO_1, ";")) %>%
    tidytable::unnest_wider(INFO_1_split) %>%
    rename(
        query_af_total = `...1`,
        query_ac_filtered = `...3`,
        query_an_filtered = `...4`
    ) %>%
    select(-contains("..."), -INFO, -QUAL, -INFO_1) %>%
    mutate(across(starts_with("query_"), ~str_extract(.x, "(?<==)\\S+"))) %>%
    mutate(query_af_filtered =
               as.numeric(query_ac_filtered) / as.numeric(query_an_filtered)) %>%

    ## Identify ClinVar review status
    mutate(clinvar_star = case_when(
        grepl("assertion", review_status, ignore.case = TRUE) ~ 0,
        grepl("multiple", review_status, ignore.case = TRUE) ~ 2,
        grepl("expert", review_status, ignore.case = TRUE) ~ 3,
        grepl("practice", review_status, ignore.case = TRUE) ~ 4,
        is.na(review_status) | review_status == "" ~ NA_real_,
        TRUE ~ 1
    )) %>%

    ## Identify variants for each tranche
    mutate(
        ptv_variants = as.integer(loftee == "HC"),
        ptv_clinvar_1 = as.integer(loftee == "HC" | (grepl("Pathogenic", clinsig) | grepl("Likely_pathogenic", clinsig)) & clinvar_star > 0),
        ptv_clinvar_2 = as.integer(loftee == "HC" | (grepl("Pathogenic", clinsig) | grepl("Likely_pathogenic", clinsig)) & clinvar_star > 1),
        rare = as.integer(query_af_filtered < as.numeric(0.005) & !grepl("Benign|Likely_benign", clinsig) & !grepl("synonymous", vep_var_type))
    )

# Separate out variant tranches
separate_variants <- function(data, variant_tranche) {
    result <- data %>%
        filter(!!sym(variant_tranche) == 1) %>%
        select(-ptv_variants, -ptv_clinvar_1, -ptv_clinvar_2, -rare)
    return(result)
}

# separate_variants_two_cats <- function(data, variant_tranche_1, variant_tranche_2) {
#     result <- data %>%
#         filter(!!sym(variant_tranche_1) == 1  | !!sym(variant_tranche_2) == 1) %>%
#         select(-ptv_variants, -ptv_clinvar_1, -ptv_clinvar_2, -rare)
#     return(result)
# }

# ptv_variants <- separate_variants(variant_data, "ptv_variants")
# ptv_clinvar_1 <- separate_variants(variant_data, "ptv_clinvar_1")
ptv_clinvar_2 <- separate_variants(variant_data, "ptv_clinvar_2")
# rare <- separate_variants(variant_data, "rare")
#
# ptv_clinvar_1_rare <- separate_variants_two_cats(variant_data, "ptv_clinvar_1", "rare")
# ptv_clinvar_2_rare <- separate_variants_two_cats(variant_data, "ptv_clinvar_2", "rare")

# Write output tables
write_data <- function(data, output_name){
    fwrite(data, output_name, sep = "\t")
}

# write_data(ptv_variants, PTV_OUTPUT)
# write_data(ptv_clinvar_1, PTV_CLINVAR_1_OUTPUT)
write_data(ptv_clinvar_2, PTV_CLINVAR_2_OUTPUT)
# write_data(rare, RARE_OUTPUT)
# write_data(ptv_clinvar_1_rare, PTV_CLINVAR_1_RARE_OUTPUT)
# write_data(ptv_clinvar_2_rare, PTV_CLINVAR_2_RARE_OUTPUT)
