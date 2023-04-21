# loads NHANES data directly from CDC specifying wave letter, years,
#   and datasets
# author: max salvatore
# date:   20230110
suppressPackageStartupMessages({
  require(data.table)
  require(haven)
  require(glue)
})

download_nhanes_data <- function(wave_letter = "J",
                                 wave_years = "2017-2018",
                                 datasets = c("DEMO", "BMX", "SMQ", "DIQ", "MCQ", "BPQ")) {
  # get urls
  if (wave_letter == "P") {
    url_paths <- glue(
      "https://wwwn.cdc.gov/Nchs/Nhanes/",
      "{wave_years}/{wave_letter}_{datasets}.xpt"
    )
  } else {
    url_paths <- glue(
      "https://wwwn.cdc.gov/Nchs/Nhanes/",
      "{wave_years}/{datasets}_{wave_letter}.xpt"
    )
  }

  # download
  out <- list()
  pb <- txtProgressBar(max = length(url_paths), width = 50, style = 3)
  for (i in seq_along(url_paths)) {
    out[[i]] <- read_xpt(url_paths[i]) |>
      as.data.table()
    setTxtProgressBar(pb, i)
  }
  close(pb)

  # merge and output
  Reduce(\(x, y) merge.data.table(x, y, all = TRUE, by = "SEQN"), out)
}
