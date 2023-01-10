# loads NHANES data directly from CDC specifying wave letter, years,
#   and datasets
# author: max salvatore
# date:   20230110
require(data.table)
require(haven)
require(glue)
require(cli)

download_nhanes_data <- function(
    wave_letter = "J",
    wave_years = "2017-2018",
    datasets = c("DEMO", "BMX", "SMQ", "DIQ", "MCQ")
) {
  
  # get urls
  if (wave_letter == "P") {
    url_paths <- glue("https://wwwn.cdc.gov/Nchs/Nhanes/",
                      "{wave_years}/{wave_letter}_{datasets}.xpt")
  } else {
    url_paths <- glue("https://wwwn.cdc.gov/Nchs/Nhanes/",
                      "{wave_years}/{datasets}_{wave_letter}.xpt")
  }
  
  # download
  out <- list()
  cli_progress_bar(
    total = length(url_paths),
    format = "{dataset} | {pb_bar} {pb_percent}")
  for (i in seq_along(url_paths)) {
    dataset <- datasets[i]
    cli_progress_update()
    out[[i]] <- read_xpt(url_paths[i]) |>
      as.data.table()
  }
  
  # merge and output
  Reduce(\(x, y) merge.data.table(x, y, all = TRUE, by = "SEQN"), out)
  
}
