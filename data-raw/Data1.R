library(usethis,readr)
Data1 <- read_csv("data-raw/Data1.csv",col_names = FALSE)

usethis::use_data(Data1, overwrite = TRUE)
