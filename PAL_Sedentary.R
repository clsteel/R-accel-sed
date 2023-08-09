library(tidyverse)
library(readr)
library(dplyr)

activpal.file.reader <-
  function(file.name.and.path)
  {
    data <- read.csv(file.name.and.path, stringsAsFactors=FALSE)
    data <- data[,(1:17)]
    names(data) <- c("time","time_approx","datacount","event_type","duration","waking_day","cum_steps","methrs", 'abs_sum_diffX','abs_sum_diffY',
                     'abs_sum_diffZ','timeUpright','timeUpsidedown','TimeBackLying','TimeFrontLying','TimeLeftLying','TimeRightLying')

    data$time <- sub("#","",data$time)
    data$time <- sub("#","",data$time)
    data[,3] <- as.numeric(as.character(data[,3]))
    data[,4] <- as.numeric(as.character(data[,4]))
    #event type codes: 0=sed, 1=stand, 2=stepping, 2.1=cycling, 3.1=primary lying, 3.2=secondary lying, 4=nonwear, 5=seated transport
    data[,5] <- as.numeric(as.character(data[,5]))
    data[,6] <- as.numeric(as.character(data[,6]))
    data[,7] <- as.numeric(as.character(data[,7]))
    data[,8] <- as.numeric(as.character(data[,8]))
    data[,9] <- as.numeric(as.character(data[,9]))
    data[,10] <- as.numeric(as.character(data[,10]))
    data[,11] <- as.numeric(as.character(data[,11]))
    data[,12] <- as.numeric(as.character(data[,12]))
    data[,13] <- as.numeric(as.character(data[,13]))
    data[,14] <- as.numeric(as.character(data[,14]))
    data[,15] <- as.numeric(as.character(data[,15]))
    data[,16] <- as.numeric(as.character(data[,16]))
    data[,17] <- as.numeric(as.character(data[,17]))

    t <- dim(data)[1]

    data <- data[!(data[,"time"] == "1899-12-30"),]
    data <- data[!(data[,"time"] == "0"),]
    n <- dim(data)[1]

    if(is.character(data$time)==TRUE&t==n)
    {
      data$time <- as.numeric(data$time)
      data$time <- as.POSIXct(as.Date(data$time,origin="1899-12-30"))
      data$time <- as.POSIXlt(data$time,tz="GMT")
      data$time <- strptime(data$time,format="%Y-%m-%d %H:%M:%S")
    }

    return(data)
  }

second.by.second <-
  function(data)
  {
    sec.by.sec.data <- data.frame(time=NA, date=NA, ap.posture=NA, mets=NA, met.hours=NA, steps=NA)
    sec.by.sec.data <- sec.by.sec.data[-1,]

    data$duration <- as.numeric(data$duration)

    data$methrs <- as.numeric(data$methrs)

    n <- dim(data)[1]
    time.of.each.event <- as.vector(difftime(strptime(data$time[seq_len(n - 1) + 1],format="%Y-%m-%d %H:%M:%S"),
                                             strptime(data$time[seq_len(n - 1)],format="%Y-%m-%d %H:%M:%S"), units="secs"))
    start.time <- strptime(data$time[1],format="%Y-%m-%d %H:%M:%S")

    time.of.each.event <- c(time.of.each.event, round(data[n,"duration"],0))
    te <- length(time.of.each.event)
    time.of.each.event[is.na(time.of.each.event)==T] <- 1
    events <- rep((1:te),time.of.each.event)

    acts <- rep(data$event_type,time.of.each.event)
    n <- length(acts)
    # The met hours per second in the interval.
    met.hours <- data$methrs/data$duration
    met.hours <- rep(met.hours,time.of.each.event)
    # To compute mets per second in the interval, multiply methours by 3600 sec/hour and divide by number of seconds.
    mets <- data$methrs * 3600/data$duration
    mets <- rep(mets,time.of.each.event)
    steps <- rep(data$cum_steps,time.of.each.event)
    # Make 15-sec epoch variable and METs
    times <- start.time+(0:(n-1))
    fifteen.sec.times <- start.time + (15*rep(0:(floor(n/15)),each=15,length=n))
    fifteen.sec.mets <- tapply(mets, fifteen.sec.times, mean)
    fifteen.sec.mets <- rep(fifteen.sec.mets, each=15, length=n)

    # Make 1-min epoch variable and METs
    times <- start.time+(0:(n-1))
    one.min.times <- start.time + (60*rep(0:(floor(n/60)),each=60,length=n))
    one.min.mets <- tapply(mets, one.min.times, mean)
    one.min.mets <- rep(one.min.mets, each=60, length=n)

    date <- substring(format(times),1,10)

    sec.by.sec.data <- merge(sec.by.sec.data, data.frame(time=times, date=date, ap.posture=acts, mets=mets, fifteen.sec.mets=fifteen.sec.mets,
                                                         one.min.mets=one.min.mets, met.hours=met.hours, steps=steps, num.events=events,
                                                         stringsAsFactors=FALSE), all=TRUE)

    sec.by.sec.data$mets <- signif(sec.by.sec.data$mets,3)

    return(sec.by.sec.data)
  }

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

pal_files <- list.files("PAL_events", pattern="*.csv", full.names=TRUE, recursive=FALSE)
lapply(pal_files, function(x) {
  pal_file_full <- activpal.file.reader(x)
  pal_file <- second.by.second(pal_file_full)
  #in pal file nonwear is marked as 4, using that for nonwear
  pal_file$is_wear <- !pal_file$ap.posture == 4
  #making a variable that shows which postures are 'sedentary' postures. (per jordan, sed is 0, 3.2, and 5)
  pal_file <- pal_file %>% mutate(posture_factor = case_when(pal_file$ap.posture==0 ~ 'SED', pal_file$ap.posture==3.2 ~ 'SED',
                                             pal_file$ap.posture == 5 ~ 'SED', TRUE ~ 'OTHER'))
  pal_file$posture_factor <- factor(pal_file$posture_factor)
  pal_file$date <- as.Date(pal_file$date,format="%Y-%m-%d")
  #creating valid day indicator based on 36000 seconds of wear time per day (10 hrs)
  pal_file$valid_index <- sapply(
    pal_file$date,
    function(i, vec) vec[match(i, names(vec))],
    vec=tapply(pal_file$is_wear, pal_file$date, sum)) >= 36000
  #running the code to create the day level sed patterns file (hopefully using my modified code)
  sb_PAL_day <- purrr::map_df(
    split(pal_file, pal_file$date),
    ~PBpatterns::analyze_bouts(.x$posture_factor, 'SED',"SB_summary",
                               is_wear=.x$is_wear,
                               valid_indices=.x$valid_index,
                               epoch_length_sec=1)
  )
  NONSB_day <- purrr::map_df(
    split(pal_file, pal_file$date),
    ~ PBpatterns::analyze_bouts(
      .x$posture_factor, "OTHER", "SB_summary",
      is_wear = .x$is_wear,
      valid_indices = .x$valid_index,
      epoch_length_sec = 1,
    )
  )
  #pull out the mean and median non sed bouts vars and rename to period. also calculate active to sed probability as reciprocal of mean
  NONSB_day$period_mean <- NONSB_day$mean_SB_bout_min
  NONSB_day$period_median <- NONSB_day$median_SB_bout_min
  NONSB_day$act_sed_prob <- 1/NONSB_day$period_mean
  day_subset <- subset(NONSB_day, select = c("period_mean","period_median","act_sed_prob"))
  #concatenate the new columns to the other ones. this is not a merge on date so it could introduce unforseen errors down the line.
  SB_day <- cbind(sb_PAL_day, day_subset)
  #add participant id variable
  SB_day$ppt_id <- basename(x)
  #aggregate by day and bind to results since paul's code doesn't retain the date in the daily sed pattern file.
  date_agg <- aggregate(pal_file, by=list(pal_file$date), FUN='first')
  SB_day <- cbind(SB_day, date_agg$date)

  write_csv(SB_day, file=paste0('test_output/sed_patterns/',basename(x),'_sed_patternsAP.csv'))
  #running the code to create the sed bout overall file
  sed_bouts_PAL <- PBpatterns::analyze_bouts(
    x=pal_file$posture_factor,
    target="SED", 'rle_standard', epoch_length_sec=1
  )
  sed_bouts_PAL$ppt_id <- basename(x)
  write_csv(sed_bouts_PAL, file=paste0("test_output/bouts/",basename(x),"_boutsAP.csv"))
})