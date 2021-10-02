## ----include=FALSE---------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(eval=TRUE, echo=TRUE, message=FALSE, warnings=FALSE)
solution <- TRUE


## ----datamunge, eval=TRUE--------------------------------------------------------------------------------------------------------------
library(activity)
duikers <- "VideoStartTimes_Daytime.txt"
actdata <- read.table(duikers, header=TRUE)
# data wrangling to make a datetime value from separate
# fields for year, month, day, hour and minute
actdata$date <- paste("2016",
                      sprintf("%02i", actdata$month), 
                      sprintf("%02i", actdata$day),
                      sep="/")
actdata$time <- paste(sprintf("%02i", actdata$hour),
                      sprintf("%02i", actdata$minute),
                      sep=":")
actdata$datetime <- paste(actdata$date, actdata$time)


## ----radian, eval=TRUE-----------------------------------------------------------------------------------------------------------------
actdata$rtime <- gettime(actdata$datetime, "%Y/%m/%d %H:%M")


## ----activity, eval=TRUE---------------------------------------------------------------------------------------------------------------
act_result <- fitact(actdata$rtime, sample="data", reps=100)


## ----actplot, eval=TRUE, fig.cap="Fitted smooth to histogram of camera triggering times for Maxwell's duiker data."--------------------
plot(act_result)


## ----thenumber, eval=TRUE--------------------------------------------------------------------------------------------------------------
print(act_result@act)


## ----organise--------------------------------------------------------------------------------------------------------------------------
daytime <- read.csv("DaytimeDistances.txt", header=TRUE, sep="\t")
daytime <- subset(daytime, select=-c(utm.e, utm.n))
daytime <- daytime[, c("Region.Label", "Area", "multiplier",
                       "Sample.Label", "Effort", "distance")]


## ----utility, eval=TRUE----------------------------------------------------------------------------------------------------------------
library(Distance)
viewangle <- 42
field.of.view <- viewangle / 360
camera.operation.per.day <- 11.5
conversion <- convert_units("meter", NULL, "square kilometer")
trunc.list <- list(right=15, left=2)
mybreaks <- c(seq(2,8,1), 10, 12, 15)


## ----dayhr0, eval=TRUE-----------------------------------------------------------------------------------------------------------------
daytime.hr0 <- ds(daytime, transect = "point", key="hr", adjustment = NULL,
                  cutpoints = mybreaks, truncation = trunc.list,
                  convert.units=conversion)
plot(daytime.hr0, pdf=TRUE)


## ----avmultiplier, eval=TRUE-----------------------------------------------------------------------------------------------------------
prop.camera.time <- camera.operation.per.day / 24
avail <- list(creation=data.frame(rate = act_result@act[1]/prop.camera.time,
                                  SE   = act_result@act[2]/prop.camera.time))


## ----dht2, eval=TRUE-------------------------------------------------------------------------------------------------------------------
daytime_dht2 <- dht2(ddf=daytime.hr0, flatfile=daytime, strat_formula=~1,
                    convert_units=conversion, sample_fraction=field.of.view, 
                    multipliers=avail)
print(daytime_dht2, report="density")


## ----mysummary-------------------------------------------------------------------------------------------------------------------------
mysummary <- function(ests, fit){
  return(data.frame(Label = ests$individuals$D$Label,
                    Dhat = ests$individuals$D$Estimate))
}


## ----multifunc, eval=TRUE--------------------------------------------------------------------------------------------------------------
mult <- list(availability= make_activity_fn(actdata$rtime, sample="data",
                                            detector_daily_duration=camera.operation.per.day))


## ---- bootstrap, results='hide', eval=TRUE---------------------------------------------------------------------------------------------
daytime.boot.hr <- bootdht(model=daytime.hr0, flatfile=daytime,
                          resample_transects = TRUE, nboot=500, cores = 11,
                          summary_fun=mysummary, sample_fraction = field.of.view,
                          convert.units = conversion, multipliers=mult)


## ----bootresult, eval=TRUE-------------------------------------------------------------------------------------------------------------
print(summary(daytime.boot.hr))


## ---- fig.width=8, fig.cap="Distribution of density estimates from bootstrap replicates.", eval=TRUE-----------------------------------
hist(daytime.boot.hr$Dhat, breaks = 20, 
     xlab="Estimated density", main="D-hat estimates bootstraps")
abline(v=quantile(daytime.boot.hr$Dhat, probs = c(0.025,0.975), na.rm=TRUE), lwd=2, lty=3)

