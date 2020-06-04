library(downloader)
library(tidyverse)

psmsl <- read.csv("R/analysis/rslr/USStationsLinearSeaLevelTrends_160927.csv") # File that contains gauge id's, start date, end date and percent complete data

output_file_name <- "AllUS_1983to2001_RSLR" # Define the output file name
date_tag <- "200605" # Give a date tag to track your work

start_date <- 1983 # Start Year
end_date <- 2001 # End Year
min_pc <- 66 # What's your tolerance for incomplete gauge data
year_tol <- 1 # What's your tolerance for gauges that start too late or end too early for the analysis?

if (! file.exists(paste(getwd(), "/", output_file_name, "_", date_tag, sep=""))) { dir.create(paste(getwd(), "/", output_file_name, "_", date_tag, sep="")) } # Create Folder for Summary Output

psmsl_download <- function(gauge_id) { # function for checking if gauge files are there and downloading them if they're not
  if (! file.exists(paste(getwd(), "/dict/", sep=""))) { dir.create(paste(getwd(), "/dict/", sep="")) } # Checks to see if there is a folder caled "dict" which stores all of the gauge data
  this_gauge <- toString(gauge_id) # Identifies the code of the gauge we're using.
  gauge_filepath<-(paste(getwd(), "/dict/", this_gauge, ".csv", sep="")) # Creates a variable with the dict string to look for
  if (file.exists(gauge_filepath) == TRUE) { #looks to see if the file is there
    print("This file exists in your folder.")
  } else {
    print("Downloading file:")
    #tgdl is string data that has has the url for the data we're looking for
    tgdl= paste("https://tidesandcurrents.noaa.gov/sltrends/data/", this_gauge, "_meantrend.csv", sep="") 
    print(tgdl)
    save_path = paste(getwd(),"/dict/", this_gauge, ".csv", sep="") #a string that shows where the download goes
    download(tgdl, save_path, mode = "wb")
    print("File successfully downloaded")
  }
}

# function for calculating slr given gauge and time
calc_rslr <- function(gauge_id, start_date, end_date, gauge_name = "") {
  
  print(gauge_name, max.levels=0)
  
  #col_names<-c("gauge_year", "msl", "code1", "code2")
  this_gauge <- toString(gauge_id) #Identifies the code of the gauge we're using.
  gauge_string<-(paste("dict/",  this_gauge, ".csv", sep="")) #creates a variable with the dict string to attach for
  
  df <- read_csv(gauge_string)
  df <-subset(df, ((Year >= start_date) & (Year <= end_date)) )
  #attach(this_record, warn.conflicts=FALSE)
  
  gauge_completeness <-length(df$Monthly_MSL[! is.na(df$Monthly_MSL) == TRUE]) / length(! is.na(df$Monthly_MSL)) * 100
  
  print(toString(gauge_completeness))
  if (gauge_completeness >= min_pc) {
    string_dates <- paste(df$Year, "-", df$Month, "-01", sep="")
    df["timestamp"] <- as.numeric(as.Date(string_dates, format="%Y-%m-%d"))
    
    new_msl <- df$Monthly_MSL * 1000
    new_year <- (df$timestamp / 365.25 ) + 1970
    
    mslr<-lm(new_msl ~ new_year)
    if (! file.exists(paste(getwd(), "/", output_file_name, "_", date_tag, "/pdf/", sep=""))) { dir.create(paste(getwd(), "/", output_file_name, "_", date_tag, "/pdf/", sep="")) }
    file_pdf = paste(getwd(), "/", output_file_name, "_", date_tag, "/pdf/", this_gauge, ".pdf", sep="")
    pdf(file_pdf,width=6,height=6) 
    
    plot(new_year, new_msl, pch=16, xlab = "Time (year)", ylab = "MSL (mm)", main = paste(gauge_name, this_gauge, sep=" "))
    lines(new_year, new_msl)
    abline(mslr, lty=2, col="red")
    
    co_vector<-data.frame(coefficients(mslr))
    attach(co_vector, warn.conflict=FALSE)
    new_mslr<- round(coefficients.mslr.[2], 2)
    new_se <- round(summary(mslr)$coefficients[2,2], 2)
    rslr_pm <- paste(toString(new_mslr), " +/- ", toString(new_se), " mm/y", sep="")
    legend("topleft", c(rslr_pm), lty=2, col="red")
    
    dev.off()
    
    return(data.frame(new_mslr, new_se))
    
  } else {
    
    return(data.frame(new_mslr = NA, new_se = NA))
  }
}

# subset the database
psmsl_df <- data.frame(psmsl)
psmsl_df <- subset(psmsl_df, ((gauge_start <= (start_date + year_tol)) & (gauge_end >= (end_date - year_tol))))

new_rslr <- c()
new_se <- c()
# for each row in the database
for (i in 1:nrow(psmsl_df)) {
  psmsl_download(psmsl_df$gauge_id[i])
  
  temp_rslr <- calc_rslr(gauge_id = psmsl_df$gauge_id[i], start_date = start_date, end_date = end_date, gauge_name = psmsl_df$gauge_name[i]) # calculate 
  new_rslr <- c(new_rslr, temp_rslr$new_mslr)
  new_se <- c(new_se, temp_rslr$new_se)
}

psmsl_df["rslr"] <- new_rslr
psmsl_df["se"] <- new_se
psmsl_df <- subset(psmsl_df, ! is.na(rslr))

write.table(psmsl_df, paste(output_file_name, "_", date_tag, "/", output_file_name, "_", date_tag, ".csv", sep="" ), row.names=F, sep = ",")
