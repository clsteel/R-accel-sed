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

lapply(chap_files, function(x) {
  epoch_file <- read.csv(x, header=TRUE)
  #create a variable containing the participant identifer pulled from the filename
  string <- basename(x)
  string <- substr(string, 1, nchar(string)-4)
  epoch_file$filename <- string
  #pull date from timestamp
  #epoch_file$Date <- as.Date(epoch_file$timestamp,format="%Y/%m/%d")
  #format timestamp as datetime
  epoch_file$timestamp <- strptime(epoch_file$timestamp, format="%Y-%m-%d %H:%M")
  #this is removing any rows with no timestamp--a leftover quirk from some actical merging
  #epoch_file <- epoch_file[!is.na(epoch_file$timestamp),]
  #create logical Actigraph weartime marker based on w/nw classification previously calculated on the data using Choi algorithm
  epoch_file$is_wear_AG <- epoch_file$AG_nonwear == 'w'
  #change this 600 to another threshold if 10 hours of wear is not the desired number. valid_index marks epochs on days
  #where there were at least 10 hours of weartime and therefore we want to use the data
  epoch_file$valid_index_AG <- sapply(
    epoch_file$Date,
    function(i, vec) vec[match(i, names(vec))],
    vec=tapply(epoch_file$is_wear_AG, epoch_file$Date, sum)) >= 6000
  #this separates behavior only into sedentary (sub 100 counts per actigraph) vs non-sedentary, no concern for mvpa
  epoch_file$posture_factor <- factor(epoch_file$prediction)
  #write an output file ready to have sed pattern vars calculated on it.
  write_csv(epoch_file, file=paste0("test_output/epoch/", string ,"_modified.csv"))
  #create and save dataset of all sedentary bouts with filename saved as participant identifier
  sed_bouts <- PBpatterns::analyze_bouts(
    x=epoch_file$prediction,
    target='sitting', 'rle_standard', epoch_length_sec=10, is_wear=epoch_file$is_wear_AG, valid_indices=epoch_file$valid_index_AG
  )
  sed_bouts$ppt_id <- string
  write_csv(sed_bouts, file=paste0("test_output/bouts/", string ,"_boutsCHAP.csv"))

  #create csv of all days with sed bout variables calculated
  sb_bouts_day <- purrr::map_df(
    split(epoch_file, epoch_file$Date),
    ~ PBpatterns::analyze_bouts(
      .x$prediction, "sitting", "SB_summary",
      is_wear = .x$is_wear_AG,
      valid_indices = .x$valid_index_AG,
      epoch_length_sec = 10,
    )
  )
  #now use the non-sb vars to create the vars needed for the period vars
  NONSB_day <- purrr::map_df(
    split(epoch_file, epoch_file$Date),
    ~ PBpatterns::analyze_bouts(
      .x$prediction, "not-sitting", "SB_summary",
      is_wear = .x$is_wear_AG,
      valid_indices = .x$valid_index_AG,
      epoch_length_sec = 10,
    )
  )
  #pull out the mean and median non sed bouts vars and rename to period. also calculate active to sed probability as reciprocal of mean
  NONSB_day$period_mean <- NONSB_day$mean_SB_bout_min
  NONSB_day$period_median <- NONSB_day$median_SB_bout_min
  NONSB_day$act_sed_prob <- 1/NONSB_day$period_mean
  day_subset <- subset(NONSB_day, select = c("period_mean","period_median","act_sed_prob"))
  #concatenate the new columns to the other ones. this is not a merge on date so it could introduce unforseen errors down the line.
  SB_day <- cbind(sb_bouts_day, day_subset)
  #add participant id variable
  SB_day$ppt_id <- string
  #aggregate by day and bind to results since paul's code doesn't retain the date in the daily sed pattern file.
  date_agg <- aggregate(epoch_file, by=list(epoch_file$Date), FUN='first')
  SB_day <- cbind(SB_day, date_agg$Date)
  write_csv(SB_day, file=paste0('test_output/sed_patterns/', string ,'_sed_patternsCHAP.csv'))
})