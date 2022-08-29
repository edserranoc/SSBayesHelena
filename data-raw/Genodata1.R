library(usethis,readr)

Genodata1 <- read_csv("data-raw/Genodata1.csv",col_names = FALSE)

usethis::use_data(Genodata1, overwrite = TRUE)

