#-------------------------------------------------------------------------------

library(ggbiplot)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(dplyr)
library(hms)
library(corrplot)
library(zoo)
library(tidyr)
library(HMM)
library(depmixS4)
library(flexmix)

#Extracting Data
DataDf <- read.table("Term_Project_Dataset.txt", header=TRUE,
                     sep=",")
DataDf$Date <-as.Date(DataDf$Date,format = "%d/%m/%Y")
DataDf$Time <- parse_hms(DataDf$Time)
DataDf$Weekday <- lubridate::wday(DataDf$Date,label = TRUE)
DataDf$week_num <- strftime(DataDf$Date, format = "%V")


#Scaling the data
colNames <- c(3:9)
Scaled_data <- DataDf
Scaled_data[colNames]  <- scale(Scaled_data[colNames])
Scaled_data<-na.omit(Scaled_data) #omiting the na values

# Define the time window
start_time <- "18:00:00"
end_time <- "21:00:00"
day_of_week <- "Fri"


#pca on the training data: time window of Friday's at 18:00-21:00
Training_data <- Scaled_data
Training_data <- subset(Training_data, Weekday == day_of_week  & 
                          Time >= parse_hms(start_time) &
                          Time <= parse_hms(end_time))

Training_pca <- prcomp(Training_data[c(3:9)], center = TRUE)
print(Training_pca)
summary(Training_pca)


Training_pca.var <- Training_pca$sdev^2
Training_pca.var.per <- round(Training_pca.var/sum(Training_pca.var)*100,1)
barplot(Training_pca.var.per,main = "Screen Plot", xlab="Principal Component", ylab = "Percent Variation")

plot(Training_pca, type = "l",xlim=c(0,7), ylim=c(0,5),main="Fridays: 18:00-21:00 PC Variances")

loading <- Training_pca$rotation

require(ggbiplot)
ggbiplot(Training_pca,ellipse=TRUE,circle=TRUE, obs.scale = 1, var.scale = 1)

#------------------------------------------------------------------------------------------.
#Part 2


Training_data <- Training_data[c(1,2,3,6)]
# first 2 years is training and last year data is testing
train_data<-subset(Training_data, Date < "2009-01-01")
test_data <- subset(Training_data, Date >= "2009-01-01")


#getting nTimes
trainCount<-(train_data %>% group_by(Date) %>% count())
nTimesTrain<- trainCount %>% pull(n)

#getting nTimesTest
testCount <- (test_data %>% group_by(Date) %>% count())
nTimesTest<- testCount %>% pull(n)


states <- vector()
logs <- vector()
BIC <- vector()
fit_models <- vector()

set.seed(9)

#Using training dataset for training a number of multivariate HMMs that each have a
#different number of states across a range from not less than 4 states to not more than 24 states
for (i in c(4:24)){
  if ((i %in% c(4,6,8,12,14,16,20,24))){
   mod <- depmix(response = list(train_data$Global_intensity ~ 1, train_data$Global_active_power ~ 1), data = train_data, nstates = i, ntimes = nTimesTrain,
                 family = list(gaussian(), gaussian()))
    fm <- fit(mod)
    summary(fm)
    fit_models <- append(fit_models,fm)
    logs <- append(logs,logLik(fm))
    BIC <- append(BIC,BIC(fm))
    states <- append(states,i)
  }
}

Plotdata <- data.frame(states,logs,BIC) #frame containing all the plotting variables

ggplot(Plotdata, aes(states)) +
  geom_line(aes(y = logs, color = "red")) +
  geom_line(aes(y = BIC, color = "green")) +
  scale_color_manual(values = c("red", "green"),labels=c("BIC","Log likelihood"),name="Color")+
  ggtitle("BIC and Log Likelihoood") +
  xlab("States") +
  ylab("Values")



test_logs <- c()
# Iterate through each selected model and calculate  log-likelihood of test-data
for (i in seq_along(fit_models)) {
  model <- fit_models[[i]]
  # Create a new model with the same specification
  new_model <- depmix(response = list(test_data$Global_intensity ~ 1, test_data$Global_active_power ~ 1),
                      data = test_data, nstates = model@nstates,ntimes = nTimesTest, family = list(gaussian(), gaussian()))
  params <- getpars(model)
  new_model <- setpars(new_model,params)
  test_logs[i] <-forwardbackward(new_model,test_data,return.all=FALSE, useC = FALSE)$logLik
}

# Calculate the number of observations in the train data
num_train_observations <- nrow(train_data)
# Calculate the number of observations in the test data
num_test_observations <- nrow(test_data)

ratio = num_train_observations / num_test_observations

# Calculate the normalized log-likelihoods of train and test data
norm_train_loglik <- logs / ratio
norm_test_loglik <- test_logs

Plotdata <- data.frame(states,norm_train_loglik,BIC,norm_test_loglik)

# Comparison of normalized log-likelihoods of train and test data
ggplot(Plotdata, aes(states)) +
  geom_line(aes(y = BIC, color = "green")) +
  geom_line(aes(y = norm_train_loglik, color = "blue")) +
  geom_line(aes(y = norm_test_loglik, color = "red")) +
  scale_color_manual(values = c( "green", "red","blue"),labels=c("Norm Log likelihood Test","BIC", "Norm Log likelihood Train"),name="Color")+
  ggtitle("BIC and Normalized Log Likelihoood of test and train") +
  xlab("States") +
  ylab("Values")

selected_model = fit_models[[which.min(BIC)]] #final selected model based on log likelihood and BIC

#------------------------------------------------------------------------------------------.
#Part 3

# Load the anomalous datasets
anomaly_data1 <- read.table("Dataset_with_Anomalies_1.txt", header=TRUE, sep=",")
anomaly_data2 <- read.table("Dataset_with_Anomalies_2.txt", header=TRUE, sep=",")
anomaly_data3 <- read.table("Dataset_with_Anomalies_3.txt", header=TRUE, sep=",")

#scaling anomaly data 1
scaled_data1 <- anomaly_data1
scaled_data1[c(3:9)]  <- scale(scaled_data1[c(3:9)])
scaled_data1$Date <- as.Date(scaled_data1$Date,format = "%d/%m/%Y")
scaled_data1$Weekday <- wday(scaled_data1$Date,label = TRUE)
scaled_data1$Time <- parse_hms(scaled_data1$Time)
scaled_data1<-na.omit(scaled_data1) #omiting the na values
# Extract the time window for normal data
anomaly_data1_tw <- subset(scaled_data1, Weekday == day_of_week &
                           Time >= parse_hms(start_time) &
                           Time <= parse_hms(end_time))
anomaly_data1_tw <- anomaly_data1_tw[c(1,2,3,6,8)]
#computing nTimes
anomaly_data1_count<-(anomaly_data1_tw %>% group_by(Date) %>% count())
ad1count<- anomaly_data1_count %>% pull(n)

#scaling anomaly data 2
scaled_data2 <- anomaly_data2
scaled_data2[c(3:9)]  <- scale(scaled_data2[c(3:9)])
scaled_data2$Date <- as.Date(scaled_data2$Date,format = "%d/%m/%Y")
scaled_data2$Weekday <- wday(scaled_data2$Date, label = TRUE)
scaled_data2$Time <- parse_hms(scaled_data2$Time)
scaled_data2<-na.omit(scaled_data2) #omiting the na values
# Extract the time window for anomaly dataset 2
anomaly_data2_tw <- subset(scaled_data2, Weekday == day_of_week &
                           Time >= parse_hms(start_time) &
                           Time <= parse_hms(end_time))
anomaly_data2_tw <- anomaly_data2_tw[c(1,2,3,6,8)]
anomaly_data2_count<-(anomaly_data2_tw %>% group_by(Date) %>% count())
ad2count<- anomaly_data2_count %>% pull(n)

#scaling anomaly data 3
scaled_data3 <- anomaly_data3
scaled_data3[c(3:9)]  <- scale(scaled_data3[c(3:9)])
scaled_data3$Date <- as.Date(scaled_data3$Date,format = "%d/%m/%Y")
scaled_data3$Weekday <- wday(scaled_data3$Date, label = TRUE)
scaled_data3$Time <- parse_hms(scaled_data3$Time)
scaled_data3<-na.omit(scaled_data3) #omiting the na values
# Extract the time window for anomaly dataset 3
anomaly_data3_tw <- subset(scaled_data3, Weekday == day_of_week &
                             Time >= parse_hms(start_time) &
                             Time <= parse_hms(end_time))
anomaly_data3_tw <- anomaly_data3_tw[c(1,2,3,6,8)]
anomaly_data3_count<-(anomaly_data3_tw %>% group_by(Date) %>% count())
ad3count<- anomaly_data3_count %>% pull(n)

#calculate  log-likelihood of anomalies data with selected trained model
new_model1 <- depmix(response = list(anomaly_data1_tw$Global_intensity ~ 1, anomaly_data1_tw$Global_active_power ~ 1),
                      data = anomaly_data1_tw, nstates = selected_model@nstates,ntimes = ad1count, family = list(gaussian(), gaussian()))
new_model2 <- depmix(response = list(anomaly_data2_tw$Global_intensity ~ 1, anomaly_data2_tw$Global_active_power ~ 1),
                       data = anomaly_data2_tw, nstates = selected_model@nstates,ntimes = ad2count, family = list(gaussian(), gaussian()))
new_model3 <- depmix(response = list(anomaly_data3_tw$Global_intensity ~ 1, anomaly_data3_tw$Global_active_power ~ 1),
                       data = anomaly_data3_tw, nstates = selected_model@nstates,ntimes = ad3count, family = list(gaussian(), gaussian()))


params <- getpars(selected_model)
new_model1 <- setpars(new_model1,params)
new_model2 <- setpars(new_model2,params)
new_model3 <- setpars(new_model3,params)
  
logs_anomoly1 <-forwardbackward(new_model1, anomaly_data1_tw, return.all = FALSE, useC = FALSE)$logLike
logs_anomoly2 <-forwardbackward(new_model2, anomaly_data2_tw, return.all = FALSE, useC = FALSE)$logLike
logs_anomoly3 <-forwardbackward(new_model3, anomaly_data3_tw, return.all = FALSE, useC = FALSE)$logLike

#log likelihoods of all 3 anomalies data
print(paste0("Log Anomalies 1 = ",logs_anomoly1))
print(paste0("Log Anomalies 2 = ",logs_anomoly2))
print(paste0("Log Anomalies 3 = ",logs_anomoly3))



