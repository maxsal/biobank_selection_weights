## extract file paths based on data version
get_files <- function(mgi_version = "20210318", ukb_version = "20221117") {
  
  if (mgi_version == "20210318") {
    icd9_file  <- "/net/junglebook/magic_data/MGI_Phenome_20210318/MGI_20210318_ICD9_Phecodes_Birthyears.Rsav"
    icd10_file <- "/net/junglebook/magic_data/MGI_Phenome_20210318/MGI_20210318_ICD10_Phecodes_Birthyears.Rsav"
    pim0_file  <- "/net/junglebook/magic_data/MGI_Phenome_20210318/MGI_20210318_PEDMASTER_0.txt"
    demo_file  <- "/net/junglebook/magic_data/ehr_data_20210318/HPI5635_Demographics.txt"
  }
  
  if (ukb_version == "20221117") {
    ukb_pim0_file        <- "/net/junglebook/home/mmsalva/createUKBphenome/results/UKB_PHENOME_20221117.txt"
    ukb_phecode_dsb_file <- "/net/junglebook/home/mmsalva/createUKBphenome/results/UKB_PHECODE_DSB_MAPPED_20221117.txt"
    ukb_demo_file        <- "/net/junglebook/home/mmsalva/createUKBphenome/results/UKB_SEX_20221117.txt"
  }
  
  return(list(
    "mgi" = list(
      icd9_file  = icd9_file,
      icd10_file = icd10_file,
      demo_file  = demo_file,
      pim0_file  = pim0_file
    ),
    "ukb" = list(
      pim0_file        = ukb_pim0_file,
      icd_phecode_file = ukb_phecode_dsb_file,
      demo_file        = ukb_demo_file
    )
  ))
  
}
