library(tidyverse)
library(readr)
library(dplyr)

sb_range_bouts_cs <- function(d, bouts) {

  data.frame(
    d,
    mean_SB_bout_min = mean(bouts$lengths_min),
    median_SB_bout_min = median(bouts$lengths_min),
    sb_0_14_hr = sum(ifelse(
      bouts$lengths_min < 15, bouts$lengths_min, 0
    )) / 60,
    sb_15_29_hr = sum(ifelse(
      bouts$lengths_min >= 15 & bouts$lengths_min < 30, bouts$lengths_min, 0
    )) / 60,
    sb_30_60_hr = sum(ifelse(
      bouts$lengths_min >= 30 & bouts$lengths_min < 60, bouts$lengths_min, 0
    )) / 60,
    sb_60_Inf_hr = sum(ifelse(
      bouts$lengths_min >=60, bouts$lengths_min, 0
    )) / 60
  ) %T>%
  {stopifnot(isTRUE(all.equal(
    sum(rev(.)[ ,1:4]), sum(bouts$lengths_min) / 60,
    scale = 1, tolerance = 1/60/10
  )))}
}

environment(sb_range_bouts_cs) <- asNamespace('PBpatterns')
assignInNamespace("sb_range_bouts", sb_range_bouts_cs, ns='PBpatterns')

chap_files <- list.files("CHAP_predictions", pattern="*.csv", full.names=TRUE, recursive=FALSE)
epoch_files <- list.files("Epoch_Preprocessed", pattern="*.csv", full.names=TRUE, recursive=FALSE)

