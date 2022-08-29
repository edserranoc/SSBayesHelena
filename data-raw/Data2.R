library(usethis,readr)
Data2 <- read_csv("data-raw/Data2.csv")

usethis::use_data(Data2, overwrite = TRUE)
