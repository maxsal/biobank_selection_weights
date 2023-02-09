library(glue)

## extract file paths based on data version
get_files <- function(mgi_version = "20220822", ukb_version = "20221117") {
  
  if (mgi_version == "20210318") {
    mgi_icd9_file  <- "/net/junglebook/magic_data/MGI_Phenome_20210318/MGI_20210318_ICD9_Phecodes_Birthyears.Rsav"
    mgi_icd10_file <- "/net/junglebook/magic_data/MGI_Phenome_20210318/MGI_20210318_ICD10_Phecodes_Birthyears.Rsav"
    mgi_pim0_file  <- "/net/junglebook/magic_data/MGI_Phenome_20210318/MGI_20210318_PEDMASTER_0.txt"
    mgi_demo_file  <- "/net/junglebook/magic_data/ehr_data_20210318/HPI5635_Demographics.txt"
    mgi_smk_file   <- "/net/junglebook/magic_data/ehr_data_20210318/HPI5635_SocialHx_DEIDENTIFIED.txt"
    mgi_bmi_file   <- "/net/junglebook/magic_data/ehr_data_20210318/HPI5635_Anthropometrics_DEIDENTIFIED.txt"
  } else {
    mgi_path              <- glue("/net/junglebook/magic_data/EHRdata/{mgi_version}/")
    mgi_icd9_file         <- glue("{mgi_path}phenomes/UNFILTERED_{mgi_version}/UNFILTERED_{mgi_version}_ICD9_Phecodes_Birthyears.Rsav")
    mgi_icd10_file        <- glue("{mgi_path}phenomes/UNFILTERED_{mgi_version}/UNFILTERED_{mgi_version}_ICD10_Phecodes_Birthyears.Rsav")
    mgi_phecode_dsb_file  <- glue("{mgi_path}phenomes/UNFILTERED_{mgi_version}/UNFILTERED_{mgi_version}_Phecodes_Birthyears.Rsav")
    mgi_pim0_file         <- glue("{mgi_path}phenomes/UNFILTERED_{mgi_version}/UNFILTERED_{mgi_version}_PEDMASTER_0.txt")
    mgi_cov_file          <- glue("{mgi_path}MGI_{mgi_version}.txt")
    mgi_phe_overview_file <- glue("{mgi_path}phenomes/UNFILTERED_{mgi_version}/UNFILTERED_{mgi_version}_Phenotype_Overview_All_Phecodes1plus.Rsav")
    mgi_demo_file         <- glue("/net/junglebook/magic_data/Data_Pulls_from_Data_Office/{mgi_version}/Demographics_{as.Date(mgi_version, '%Y%m%d')}.txt")
    mgi_smk_file          <- glue("/net/junglebook/magic_data/Data_Pulls_from_Data_Office/{mgi_version}/SocialHx_{as.Date(mgi_version, '%Y%m%d')}.txt")
    mgi_bmi_file          <- glue("/net/junglebook/magic_data/Data_Pulls_from_Data_Office/{mgi_version}/Anthropometrics_{as.Date(mgi_version, '%Y%m%d')}.txt")
  }
  
  if (ukb_version == "20221117") {
    ukb_pim0_file        <- "/net/junglebook/home/mmsalva/createUKBphenome/results/UKB_PHENOME_20221117.txt"
    ukb_phecode_dsb_file <- "/net/junglebook/home/mmsalva/createUKBphenome/results/UKB_PHECODE_DSB_MAPPED_20221117.txt"
    ukb_demo_file        <- "/net/junglebook/home/mmsalva/createUKBphenome/results/UKB_SEX_20221117.txt"
  }
  
  mgi_out <- list(
    icd9_file  = mgi_icd9_file,
    icd10_file = mgi_icd10_file,
    demo_file  = mgi_demo_file,
    pim0_file  = mgi_pim0_file,
    smk_file   = mgi_smk_file,
    bmi_file   = mgi_bmi_file
  )
  
  if (mgi_version != "20210318") {
    mgi_out[["phecode_dsb_file"]]  <- mgi_phecode_dsb_file
    mgi_out[["cov_file"]]          <- mgi_cov_file
    mgi_out[["phe_overview_file"]] <- mgi_phe_overview_file
  }
  
  return(list(
    "mgi" = mgi_out,
    "ukb" = list(
      pim0_file        = ukb_pim0_file,
      icd_phecode_file = ukb_phecode_dsb_file,
      demo_file        = ukb_demo_file
    )
  ))
  
}
