library(usethis,readr)
Data2 <- read.csv("data-raw/Data2.csv",sep=" ")

usethis::use_data(Data2, overwrite = TRUE)
