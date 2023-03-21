library(tidyverse)
library(bigrquery)
library(data.table)
dir.create("./data")
dir.create("./results")

if(!file.exists("./data/dataset_20230309_condition.csv")){

    dataset_20230309_condition_sql <- paste("
        SELECT
            c_occurrence.person_id,
            c_occurrence.condition_concept_id,
            c_standard_concept.concept_name as standard_concept_name,
            c_standard_concept.concept_code as standard_concept_code,
            c_standard_concept.vocabulary_id as standard_vocabulary,
            c_occurrence.condition_start_datetime,
            c_occurrence.condition_end_datetime,
            c_occurrence.condition_type_concept_id,
            c_type.concept_name as condition_type_concept_name,
            c_occurrence.stop_reason,
            c_occurrence.visit_occurrence_id,
            visit.concept_name as visit_occurrence_concept_name,
            c_occurrence.condition_source_value,
            c_occurrence.condition_source_concept_id,
            c_source_concept.concept_name as source_concept_name,
            c_source_concept.concept_code as source_concept_code,
            c_source_concept.vocabulary_id as source_vocabulary,
            c_occurrence.condition_status_source_value,
            c_occurrence.condition_status_concept_id,
            c_status.concept_name as condition_status_concept_name 
        FROM
            `condition_occurrence` c_occurrence 
        LEFT JOIN
            `concept` c_standard_concept 
                ON c_occurrence.condition_concept_id = c_standard_concept.concept_id 
        LEFT JOIN
            `concept` c_type 
                ON c_occurrence.condition_type_concept_id = c_type.concept_id 
        LEFT JOIN
            `visit_occurrence` v 
                ON c_occurrence.visit_occurrence_id = v.visit_occurrence_id 
        LEFT JOIN
            `concept` visit 
                ON v.visit_concept_id = visit.concept_id 
        LEFT JOIN
            `concept` c_source_concept 
                ON c_occurrence.condition_source_concept_id = c_source_concept.concept_id 
        LEFT JOIN
            `concept` c_status 
                ON c_occurrence.condition_status_concept_id = c_status.concept_id 
        WHERE c_source_concept.vocabulary_id IN ('ICD9CM', 'ICD10CM')", sep="")

    # Formulate a Cloud Storage destination path for the data exported from BigQuery.
    # NOTE: By default data exported multiple times on the same day will overwrite older copies.
    #       But data exported on a different days will write to a new location so that historical
    #       copies can be kept as the dataset definition is changed.
    condition_20230309_path <- file.path(
      Sys.getenv("WORKSPACE_BUCKET"),
      "bq_exports",
      Sys.getenv("OWNER_EMAIL"),
      strftime(lubridate::now(), "%Y%m%d"),  # Comment out this line if you want the export to always overwrite.
      "condition_20230309",
      "condition_20230309_*.csv")
    message(str_glue('The data will be written to {condition_20230309_path}. Use this path when reading ',
                     'the data into your notebooks in the future.'))

    # Perform the query and export the dataset to Cloud Storage as CSV files.
    # NOTE: You only need to run `bq_table_save` once. After that, you can
    #       just read data from the CSVs in Cloud Storage.
    bq_table_save(
      bq_dataset_query(Sys.getenv("WORKSPACE_CDR"), dataset_20230309_condition_sql, billing = Sys.getenv("GOOGLE_PROJECT")),
      condition_20230309_path,
      destination_format = "CSV")

    # Copy the data directly from Cloud Storage to the Jupyter disk.
    dir.create("./data/condition",recursive=T)
    system2('gsutil',paste("-m cp ",condition_20230309_path,"./data/condition"))
}



if(!file.exists("./data/dataset_20230309_person.csv")){
    # This query represents dataset for domain "person" and was generated for All of Us Controlled Tier Dataset v6
    dataset_20230309_person_sql <- paste("
        SELECT
            person.person_id,
            person.gender_concept_id,
            p_gender_concept.concept_name as gender,
            person.birth_datetime as date_of_birth,
            person.race_concept_id,
            p_race_concept.concept_name as race,
            person.ethnicity_concept_id,
            p_ethnicity_concept.concept_name as ethnicity,
            person.sex_at_birth_concept_id,
            p_sex_at_birth_concept.concept_name as sex_at_birth 
        FROM
            `person` person 
        LEFT JOIN
            `concept` p_gender_concept 
                ON person.gender_concept_id = p_gender_concept.concept_id 
        LEFT JOIN
            `concept` p_race_concept 
                ON person.race_concept_id = p_race_concept.concept_id 
        LEFT JOIN
            `concept` p_ethnicity_concept 
                ON person.ethnicity_concept_id = p_ethnicity_concept.concept_id 
        LEFT JOIN
            `concept` p_sex_at_birth_concept 
                ON person.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id", sep="")

    # Formulate a Cloud Storage destination path for the data exported from BigQuery.
    # NOTE: By default data exported multiple times on the same day will overwrite older copies.
    #       But data exported on a different days will write to a new location so that historical
    #       copies can be kept as the dataset definition is changed.
    person_20230309_path <- file.path(
      Sys.getenv("WORKSPACE_BUCKET"),
      "bq_exports",
      Sys.getenv("OWNER_EMAIL"),
      strftime(lubridate::now(), "%Y%m%d"),  # Comment out this line if you want the export to always overwrite.
      "person_20230309",
      "person_20230309_*.csv")
    message(str_glue('The data will be written to {person_20230309_path}. Use this path when reading ',
                     'the data into your notebooks in the future.'))

    # Perform the query and export the dataset to Cloud Storage as CSV files.
    # NOTE: You only need to run `bq_table_save` once. After that, you can
    #       just read data from the CSVs in Cloud Storage.
    bq_table_save(
      bq_dataset_query(Sys.getenv("WORKSPACE_CDR"), dataset_20230309_person_sql, billing = Sys.getenv("GOOGLE_PROJECT")),
      person_20230309_path,
      destination_format = "CSV")

    # Copy the data directly from Cloud Storage to the Jupyter disk.
    dir.create("./data/person",recursive=T)
    system2('gsutil',paste("-m cp ",person_20230309_path,"./data/person"))
}



if(!file.exists("./data/dataset_20230309_zip_code_socioeconomic.csv")){

    # This query represents dataset for domain "zip_code_socioeconomic" and was generated for All of Us Controlled Tier Dataset v6
    dataset_20230309_zip_code_socioeconomic_sql <- paste("
        SELECT
            observation.person_id,
            observation.observation_datetime,
            zip_code.zip3_as_string as zip_code,
            zip_code.fraction_assisted_income as assisted_income,
            zip_code.fraction_high_school_edu as high_school_education,
            zip_code.median_income,
            zip_code.fraction_no_health_ins as no_health_insurance,
            zip_code.fraction_poverty as poverty,
            zip_code.fraction_vacant_housing as vacant_housing,
            zip_code.deprivation_index,
            zip_code.acs as american_community_survey_year 
        FROM
            `zip3_ses_map` zip_code 
        JOIN
            `observation` observation 
                ON CAST(SUBSTR(observation.value_as_string,
            0,
            STRPOS(observation.value_as_string,
            '*') - 1) AS INT64) = zip_code.zip3 
            AND observation_source_concept_id = 1585250 
            AND observation.value_as_string NOT LIKE 'Res%'", sep="")

    # Formulate a Cloud Storage destination path for the data exported from BigQuery.
    # NOTE: By default data exported multiple times on the same day will overwrite older copies.
    #       But data exported on a different days will write to a new location so that historical
    #       copies can be kept as the dataset definition is changed.
    zip_code_socioeconomic_20230309_path <- file.path(
      Sys.getenv("WORKSPACE_BUCKET"),
      "bq_exports",
      Sys.getenv("OWNER_EMAIL"),
      strftime(lubridate::now(), "%Y%m%d"),  # Comment out this line if you want the export to always overwrite.
      "zip_code_socioeconomic_20230309",
      "zip_code_socioeconomic_20230309_*.csv")
    message(str_glue('The data will be written to {zip_code_socioeconomic_20230309_path}. Use this path when reading ',
                     'the data into your notebooks in the future.'))

    # Perform the query and export the dataset to Cloud Storage as CSV files.
    # NOTE: You only need to run `bq_table_save` once. After that, you can
    #       just read data from the CSVs in Cloud Storage.
    bq_table_save(
      bq_dataset_query(Sys.getenv("WORKSPACE_CDR"), dataset_20230309_zip_code_socioeconomic_sql, billing = Sys.getenv("GOOGLE_PROJECT")),
      zip_code_socioeconomic_20230309_path,
      destination_format = "CSV")


    # Copy the data directly from Cloud Storage to the Jupyter disk.
    dir.create("./data/zip_code_socioeconomic",recursive=T)
    system2('gsutil',paste("-m cp ",zip_code_socioeconomic_20230309_path,"./data/zip_code_socioeconomic"))
}



if(!file.exists("./data/dataset_20230309_measurement.csv")){
    # This query represents dataset for domain "measurement" and was generated for All of Us Controlled Tier Dataset v6
    dataset_20230309_measurement_sql <- paste("
        SELECT
            measurement.person_id,
            measurement.measurement_concept_id,
            m_standard_concept.concept_name as standard_concept_name,
            m_standard_concept.concept_code as standard_concept_code,
            m_standard_concept.vocabulary_id as standard_vocabulary,
            measurement.measurement_datetime,
            measurement.measurement_type_concept_id,
            m_type.concept_name as measurement_type_concept_name,
            measurement.operator_concept_id,
            m_operator.concept_name as operator_concept_name,
            measurement.value_as_number,
            measurement.value_as_concept_id,
            m_value.concept_name as value_as_concept_name,
            measurement.unit_concept_id,
            m_unit.concept_name as unit_concept_name,
            measurement.range_low,
            measurement.range_high,
            measurement.visit_occurrence_id,
            m_visit.concept_name as visit_occurrence_concept_name,
            measurement.measurement_source_value,
            measurement.measurement_source_concept_id,
            m_source_concept.concept_name as source_concept_name,
            m_source_concept.concept_code as source_concept_code,
            m_source_concept.vocabulary_id as source_vocabulary,
            measurement.unit_source_value,
            measurement.value_source_value 
        FROM
            ( SELECT
                * 
            FROM
                `measurement` measurement 
            WHERE
                (
                    measurement_source_concept_id IN  (
                        SELECT
                            DISTINCT c.concept_id 
                        FROM
                            `cb_criteria` c 
                        JOIN
                            (
                                select
                                    cast(cr.id as string) as id 
                                FROM
                                    `cb_criteria` cr 
                                WHERE
                                    concept_id IN (
                                        903124, 903121, 903133
                                    ) 
                                    AND full_text LIKE '%_rank1]%'
                            ) a 
                                ON (
                                    c.path LIKE CONCAT('%.',
                                a.id,
                                '.%') 
                                OR c.path LIKE CONCAT('%.',
                                a.id) 
                                OR c.path LIKE CONCAT(a.id,
                                '.%') 
                                OR c.path = a.id) 
                            WHERE
                                is_standard = 0 
                                AND is_selectable = 1
                            )
                    )
                ) measurement 
            LEFT JOIN
                `concept` m_standard_concept 
                    ON measurement.measurement_concept_id = m_standard_concept.concept_id 
            LEFT JOIN
                `concept` m_type 
                    ON measurement.measurement_type_concept_id = m_type.concept_id 
            LEFT JOIN
                `concept` m_operator 
                    ON measurement.operator_concept_id = m_operator.concept_id 
            LEFT JOIN
                `concept` m_value 
                    ON measurement.value_as_concept_id = m_value.concept_id 
            LEFT JOIN
                `concept` m_unit 
                    ON measurement.unit_concept_id = m_unit.concept_id 
            LEFT JOIn
                `visit_occurrence` v 
                    ON measurement.visit_occurrence_id = v.visit_occurrence_id 
            LEFT JOIN
                `concept` m_visit 
                    ON v.visit_concept_id = m_visit.concept_id 
            LEFT JOIN
                `concept` m_source_concept 
                    ON measurement.measurement_source_concept_id = m_source_concept.concept_id", sep="")

    # Formulate a Cloud Storage destination path for the data exported from BigQuery.
    # NOTE: By default data exported multiple times on the same day will overwrite older copies.
    #       But data exported on a different days will write to a new location so that historical
    #       copies can be kept as the dataset definition is changed.
    measurement_20230309_path <- file.path(
      Sys.getenv("WORKSPACE_BUCKET"),
      "bq_exports",
      Sys.getenv("OWNER_EMAIL"),
      strftime(lubridate::now(), "%Y%m%d"),  # Comment out this line if you want the export to always overwrite.
      "measurement_20230309",
      "measurement_20230309_*.csv")
    message(str_glue('The data will be written to {measurement_20230309_path}. Use this path when reading ',
                     'the data into your notebooks in the future.'))

    # Perform the query and export the dataset to Cloud Storage as CSV files.
    # NOTE: You only need to run `bq_table_save` once. After that, you can
    #       just read data from the CSVs in Cloud Storage.
    bq_table_save(
      bq_dataset_query(Sys.getenv("WORKSPACE_CDR"), dataset_20230309_measurement_sql, billing = Sys.getenv("GOOGLE_PROJECT")),
      measurement_20230309_path,
      destination_format = "CSV")
    
    # Copy the data directly from Cloud Storage to the Jupyter disk.
    dir.create("./data/measurement",recursive=T)
    system2('gsutil',paste("-m cp ",measurement_20230309_path,"./data/measurement"))

}



domains <- c("condition","measurement","person","zip_code_socioeconomic")

# combine files into large files for easier processing, then delete

for(i in 1:length(domains)){
    ifiles <- list.files(paste0("./data/",domains[i]),".csv$",full.name=T)
    outfile <- paste0("./data/dataset_20230309_",domains[i],".csv")
    cmdLine <- paste0("awk 'FNR==1 && NR!=1{next;}{print}' ./data/",domains[i],"/*.csv >", outfile)
    system(cmdLine)
}



