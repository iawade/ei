# Script for Generating Relevant Phenotypic Information from UKBB

# Load Packages
library(data.table)
library(tidytable)
library(tidyverse)

# Read in UKBB Phenotying Data
phenotyping_data <- fread("all_participants_updated_fields_participant.tsv")

## Ovarian Cancer
## ICD 10 = C56*; 183.0 in ICD 9

# Read in UKBB Data dictionary
data_dictionary <- fread("Data_Dictionary_Showcase.tsv")

# Collapse Columns
phenotyping_data <- phenotyping_data %>%
    unite("p40008_collapsed", starts_with("p40008"), sep = ";", na.rm = TRUE) %>%
    unite("p40006_collapsed", starts_with("p40006"), sep = ";", na.rm = TRUE) %>%
    unite("p40013_collapsed", starts_with("p40013"), sep = ";", na.rm = TRUE) %>%
    unite("p40011_collapsed", starts_with("p40011"), sep = ";", na.rm = TRUE) %>%
    unite("p40012_collapsed", starts_with("p40012"), sep = ";", na.rm = TRUE) %>%
    unite("p134_collapsed", starts_with("p134"), sep = ";", na.rm = TRUE) %>%
    unite("p20001_collapsed", starts_with("p20001"), sep = ";", na.rm = TRUE) %>%
    unite("p84_collapsed", starts_with("p84"), sep = ";", na.rm = TRUE) %>%
    unite("p20007_collapsed", starts_with("p20007"), sep = ";", na.rm = TRUE) %>%
    unite("p40007_collapsed", starts_with("p40007"), sep = ";", na.rm = TRUE) %>%
    unite("p40001_collapsed", starts_with("p40001"), sep = ";", na.rm = TRUE)

# Convert to datatable for processing/munging
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

# Rename columns further
phenotyping_data_updated_cols <- phenotyping_data_updated_cols %>%
    select(-`Interpolated Age of participant when non-cancer illness first diagnosed`) %>%
    rename(
        eid = 1,
        sex = 2,
        number_cancer_occurrences_in_register = 3,
        age_at_cancer_diagnosis_from_register = 4,
        type_of_cancer_ICD10_from_register = 5,
        type_of_cancer_ICD9_from_register = 6,
        histogology_info_from_register = 7,
        tumour_behaviour_info_from_register =8,
        number_of_self_reported_cancers = 9,
        cancer_code_self_reported = 10,
        cancer_year_age_self_reported = 11,
        age_at_cancer_diagnosis_interpolated_self_reported = 12,
        year_of_birth = 13,
        month_of_birth = 14,
        age_at_death = 15,
        primary_cause_of_death_ICD10 = 16,
        date_lost_to_follow_up = 17,
        reason_lost_to_follow_up = 18
    )

# split age at death and just retain one
phenotyping_data_updated_cols <- phenotyping_data_updated_cols %>%
    mutate(
        age_at_death = ifelse(str_detect(age_at_death, ";"), gsub(".*;\\s*([^;]+)$", "\\1", age_at_death), age_at_death)
    )


# Obtain "censored" age
phenotyping_data_updated_cols <- phenotyping_data_updated_cols %>%
    mutate(
        birthdate = make_date(year_of_birth, month = match(month_of_birth, month.name)),
        current_age = floor(as.numeric(difftime(Sys.Date(), birthdate, units = "days")) / 365.25)
    ) %>%
    mutate(
        age_lost_to_follow_up = floor(as.numeric(difftime(date_lost_to_follow_up, birthdate, units = "days")) / 365.25)
    ) %>%
    select(-birthdate) %>%
    mutate(
        age_at_death = ifelse(str_detect(age_at_death, ";"), gsub(".*;\\s*([^;]+)$", "\\1", age_at_death), age_at_death)
    ) %>%
    mutate(
        censored_age = ifelse(
            !is.na(age_at_death) & age_at_death != 0 & age_at_death != "",
            as.integer(age_at_death),
            ifelse(!is.na(age_lost_to_follow_up), age_lost_to_follow_up, current_age)
        )
    ) %>%
    select(-current_age, -date_lost_to_follow_up) %>%
    mutate(eid = as.character(eid))

# Tidy up cancers and save somewhere on HPC (with date in the name)
## Keeping self-reported data in messy state - reflects nature of acquiring info
phenotyping_data_updated_cols <- phenotyping_data_updated_cols %>%
    mutate(
        number_cancer_occurrences_in_register = coalesce(number_cancer_occurrences_in_register, 0),
        type_of_cancer_ICD10_from_register = pmap_chr(
            list(
                strsplit(type_of_cancer_ICD10_from_register, ";"),
                number_cancer_occurrences_in_register
            ),
            function(x, n) if (n == 0) x[1] else paste(x[1:n], collapse = ";")
        ),
        type_of_cancer_ICD9_from_register = pmap_chr(
            list(
                strsplit(type_of_cancer_ICD9_from_register, ";"),
                number_cancer_occurrences_in_register
            ),
            function(x, n) if (n == 0) x[1] else paste(x[1:n], collapse = ";")
        ),
        histogology_info_from_register = pmap_chr(
            list(
                strsplit(histogology_info_from_register, ";"),
                number_cancer_occurrences_in_register
            ),
            function(x, n) if (n == 0) x[1] else paste(x[1:n], collapse = ";")
        ),
        tumour_behaviour_info_from_register = pmap_chr(
            list(
                strsplit(tumour_behaviour_info_from_register, ";"),
                number_cancer_occurrences_in_register
            ),
            function(x, n) if (n == 0) x[1] else paste(x[1:n], collapse = ";")
        ),
        primary_cause_of_death_ICD10 = str_remove(primary_cause_of_death_ICD10, ";(?=[^;]*$)"),
        cancer_code_self_reported = pmap_chr(
            list(
                number_of_self_reported_cancers,
                cancer_code_self_reported
            ),
            function(n, code) if (n == 0 || all(strsplit(as.character(n), ";")[[1]] %in% "0")) "" else code
        )
    ) %>%
    select(-cancer_year_age_self_reported)

fwrite(phenotyping_data_updated_cols, "ukb_phenotyping_data_munged_all_participants.tsv", sep="\t")
