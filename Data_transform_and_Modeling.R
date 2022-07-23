####
# Prediction of elemental concentration 
# from soil NIR Spectra
# using Random Forest model
####

# Packages
require(prospectr)
require(dplyr)
require(ggridges)
require(ggplot2)
require(tidyr)
require(ggthemes)
require(clhs)
require(doParallel)
require(caret)
require(randomForest)
require(pls)
require(devtools)
print("Packages loaded!")


# Functions ----
# Function to assess models results
goof <- function(observed, predicted){
  # Coefficient of determination
  rLM <- lm(predicted ~ observed)
  R2 <- as.matrix(summary(rLM)$adj.r.squared)
  
  # Standard error of prediction ^2
  SEP2 <- mean((observed - predicted)^2)
  
  # Standard error of prediction
  SEP <- sqrt(SEP2)
  
  #Bias
  bias <- mean(predicted) - mean(observed)
  
  # residual  variance
  SEP2c <- sum(((predicted - bias - observed)^2) / length(observed))
  SEPc <- sqrt(SEP2c)
  
  # ratio of performance to deviation
  RPD <- sd(observed) / SEP
  
  # Ratio of performance to interquartile distance
  IQ <- c(quantile(observed))[4] - c(quantile(observed))[2]
  RPIQ <- IQ / SEP
  
  # Concordance
  mx <- mean(observed)
  my <- mean(predicted)
  s2x <- var(observed)
  s2y <- var(predicted)
  sxy <- mean((observed-mx) * (predicted-my))
  ccc <- 2 * sxy / (s2x + s2y + (mx - my)^2)
  
  gf <- data.frame(R2=R2, concordance=ccc, MSE=SEP2, RMSE=SEP, bias=bias, 
                MSEc=SEP2c,RMSEc=SEPc, RPD=RPD, RPIQ=RPIQ, row.names=NULL)

  return(gf)
}

# Setting fonts for plots
windowsFonts(
  A=windowsFont('Arial'),
  B=windowsFont('Times New Roman')
)

# Setting up the folder and data location ----
dir <- ("C:/NIR_RF_Prediction_R/files/")

# Creating necessary folders
dir.create(paste0(dir, "Tables"))
dir.create(paste0(dir, "Tables/validation"))
dir.create(paste0(dir, "Tables/calibration"))
dir.create(paste0(dir, "figures"))

# Importing dataset and defining variables ----
db = read.csv(paste0(dir, "dataset.csv"))
str(db[,1:30])
colnames(db[,1:30])

print("WARNING! The following variables must be inserted manually")
print("The insertion must be based on the first and last columns with elemental concentrations")
print("The following numbers correspond to the index columns from the example dataset")

i = 4 #first column with elemental concentration
j = 9 #last column with elemental concentration
# i and j is the firt and last columns, respectively, with known elemental concentrations

db[, i:j] = sapply(db[, i:j], as.numeric)
summary(db[, i:j])

# Checking and removing NAs ----
summary(db[, i:j])
row.has.na = apply(db, 1, function(x){any(is.na(x))})
sum(row.has.na)
db.ok = db[complete.cases(db), ]
summary(db.ok[, i:j])

### ------------------ ###

# Filtering samples (optional)----

db.all = db.ok
db.ok = filter(db.ok, uso == 'CANA') # selecting samples from a specific land use

# Resampling spectra from 0.5 to 1 nm interval to reduce number of variables and noise ----

str(db.ok[, 1:30])
colnames(db.ok)[(j + 1):length(db.ok)] = gsub("X", "", names(db.ok)[(j+1):length(db.ok)])

db.res = db.ok[, (j + 1):length(db.ok)]
wav.band = as.numeric(colnames(db.res)) # retrieving wavelengths

class(db.res)
db.res = as.matrix(db.res)
rownames(db.res) = db.ok$id

db.res = resample(db.res, wav.band, seq(1000, 2500, 1), interpol = "spline") # window size of 1 nm
class(db.res)
db.res = as.data.frame(db.res)

db.res$id = rownames(db.res)
db.res = db.res %>%
  select(c("id"), everything())
str(db.res[,1:5])

str(db.ok[,1:30])
db.res$id = as.integer(c(1:(nrow(db.ok))))

db.ok.res = cbind(db.ok[, i:j], db.res[,2:1502])
str(db.ok.res[,1:30])
colnames(db.ok.res[,1:30])
head(db.ok.res[, c(1:30, length(db.ok.res))])

write.csv(db.ok.res, paste0(dir, 'db_resampled.csv'), row.names = F)

# Splitting dataset intotest and validation using conditioned Latin Hypercube Sampling (cLHS) ----
colnames(db.ok.res)

# the splitting will be made at a 80 to 20 calibration/validation ratio
ree.clhs = clhs(db.ok.res[, (j - 1):length(db.ok.res)],
                  size = as.integer(0.8*nrow(db.ok.res)),
                  iter = 1000, progress = FALSE,
                  simple = FALSE)
print("cLHS splitting completed!")

plot(ree.clhs, mode = 'obj')

db.train = db.ok.res[ree.clhs$index_samples, ] # 80%
head(db.train[, c(1:5, length(db.train))])
db.train$id_spec = rownames(db.train)
db.train = db.train %>%
  select('id_spec', everything())
write.csv(db.train, paste0(dir, " db_train.csv"), row.names = F)

db.val = db.ok.res[-ree.clhs$index_samples, ] # 20%
head(db.val[, c(1:5, length(db.val))])
db.val$id_spec = rownames(db.val)
db.val = db.val %>%
  select('id_spec', everything())
write.csv(db.val, paste0(dir, " db_val.csv"), row.names = F)

### ------------------ ###

# Preprocessing spectra ####

## Multiplicative Scatter Correction
### Training data
colnames(db.train[, 1:30])
db.train.msc = prospectr::msc(db.train[, j:length(names(db.train))])
class(db.train.msc)

db.train.msc = as.data.frame(db.train.msc)
db.train.msc$id_spec = rownames(db.train.msc)
head(db.train.msc[, c(1:5, length(db.train.msc))])

db.train.msc = db.train.msc %>%
  select(c("id_spec"), everything())

#Merging DFs
colnames(db.train[, 1:30])
colnames(db.train.msc[, 1:30])
 
db.train.msc = cbind(db.train[,1:(j-1)], db.train.msc[,2:length(db.train.msc)])
head(db.train.msc[, c(1:6, length(db.train.msc))])

### Validation data
colnames(db.val[,1:30])
db.val.msc = prospectr::msc(db.val[, j:length(names(db.val))])
class(db.val.msc)
db.val.msc = as.data.frame(db.val.msc)
db.val.msc$id_spec = rownames(db.val.msc)
head(db.val.msc[, c(1:6, length(db.val.msc))])

db.val.msc = db.val.msc %>%
  select(c("id_spec"), everything())

#Merging DFs
colnames(db.val[, 1:30])
colnames(db.val.msc[, 1:30])

db.val.msc = cbind(db.val[,1:(j-1)], db.val.msc[,2:length(db.val.msc)])
head(db.val.msc[, c(1:6, length(db.val.msc))])

## Savitzky-Golay
### Training data
colnames(db.train[, 1:30])
head(db.val[, c(1:6, length(db.val))])
db.train.sg = savitzkyGolay(db.train[, j:length(names(db.train))], m = 1, p = 1, w=9)
class(db.train.sg)

db.train.sg = as.data.frame(db.train.sg)
db.train.sg$id_spec = rownames(db.train.sg)
head(db.train.sg[, c(1:5, length(db.train.sg))])

db.train.sg = db.train.sg %>%
  select(c("id_spec"), everything())

#Merging DFs
colnames(db.train[, 1:30])
colnames(db.train.sg[, 1:30])

db.train.sg = cbind(db.train[,1:(j-1)], db.train.sg[,2:length(db.train.sg)])
head(db.train.sg[, c(1:6, length(db.train.sg))])

### Validation data
colnames(db.val[,1:30])
db.val.sg = savitzkyGolay(db.val[, j:length(names(db.train))], m = 1, p = 1, w=9)
class(db.val.sg)
db.val.sg = as.data.frame(db.val.sg)
db.val.sg$id_spec = rownames(db.val.sg)
head(db.val.sg[, c(1:6, length(db.val.sg))])

db.val.sg = db.val.sg %>%
  select(c("id_spec"), everything())

#Merging DFs
colnames(db.val[, 1:30])
colnames(db.val.sg[, 1:30])

db.val.sg = cbind(db.val[,1:(j-1)], db.val.sg[,2:length(db.val.sg)])

## Standard Normal Variate
### Training data
colnames(db.train[, 1:30])
db.train.snv = standardNormalVariate(db.train[, j:length(names(db.train))])
class(db.train.snv)

db.train.snv = as.data.frame(db.train.snv)
db.train.snv$id_spec = rownames(db.train.snv)
head(db.train.snv[, c(1:5, length(db.train.snv))])

db.train.snv = db.train.snv %>%
  select(c("id_spec"), everything())

#Merging DFs
colnames(db.train[, 1:30])
colnames(db.train.snv[, 1:30])

db.train.snv = cbind(db.train[,1:(j-1)], db.train.snv[,2:length(db.train.snv)])
head(db.train.snv[, c(1:6, length(db.train.snv))])

### Validation data
colnames(db.val[,1:30])
db.val.snv = standardNormalVariate(db.val[, j:length(names(db.val))])
class(db.val.snv)
db.val.snv = as.data.frame(db.val.snv)
db.val.snv$id_spec = rownames(db.val.snv)
head(db.val.snv[, c(1:6, length(db.val.snv))])

db.val.snv = db.val.snv %>%
  select(c("id_spec"), everything())

#Merging DFs
colnames(db.val[, 1:30])
colnames(db.val.snv[, 1:30])

db.val.snv = cbind(db.val[,1:(j-1)], db.val.snv[,2:length(db.val.snv)])
head(db.val.snv[, c(1:6, length(db.val.snv))])

## Detrend normalization
### Training data
colnames(db.train[, 1:30])
wav.band.det = as.numeric(colnames(db.train[, j:length(names(db.train))]))
db.train.det = detrend(db.train[, j:length(names(db.train))], wav.band.det, p = 2)
class(db.train.det)

db.train.det = as.data.frame(db.train.det)
db.train.det$id_spec = rownames(db.train.det)
head(db.train.det[, c(1:5, length(db.train.det))])

db.train.det = db.train.det %>%
  select(c("id_spec"), everything())

#Merging DFs
colnames(db.train[, 1:30])
colnames(db.train.det[, 1:30])

db.train.det = cbind(db.train[,1:(j-1)], db.train.det[,2:length(db.train.det)])
head(db.train.det[, c(1:6, length(db.train.det))])

### Validation data
colnames(db.val[,1:30])
db.val.det = detrend(db.val[, j:length(names(db.val))], wav.band.det, p = 2)
class(db.val.det)
db.val.det = as.data.frame(db.val.det)
db.val.det$id_spec = rownames(db.val.det)
head(db.val.det[, c(1:6, length(db.val.det))])

db.val.det = db.val.det %>%
  select(c("id_spec"), everything())

#Merging DFs
colnames(db.val[, 1:30])
colnames(db.val.det[, 1:30])

db.val.det = cbind(db.val[,1:(j-1)], db.val.det[,2:length(db.val.det)])
head(db.val.det[, c(1:6, length(db.val.det))])

## Continuum removal
### Training data
colnames(db.train[, 1:30])
wav.band.cr = as.numeric(colnames(db.train[, j:length(names(db.train))]))
db.train.cr = continuumRemoval(db.train[, j:length(names(db.train))], 
                               wav.band.cr, 
                               type = "R", interpol = 'linear', 
                               method = 'division')
class(db.train.cr)

db.train.cr = as.data.frame(db.train.cr)
db.train.cr$id_spec = rownames(db.train.cr)
head(db.train.cr[, c(1:5, length(db.train.cr))])

db.train.cr = db.train.cr %>%
  select(c("id_spec"), everything())

#Merging DFs
colnames(db.train[, 1:30])
colnames(db.train.cr[, 1:30])

db.train.cr = cbind(db.train[,1:(j-1)], db.train.cr[,2:length(db.train.cr)])
head(db.train.cr[, c(1:6, length(db.train.cr))])

### Validation data
colnames(db.val[,1:30])
db.val.cr = continuumRemoval(db.val[, j:length(names(db.val))], 
                             wav.band.cr, 
                             type = "R", interpol = 'linear', 
                             method = 'division')
class(db.val.cr)
db.val.cr = as.data.frame(db.val.cr)
db.val.cr$id_spec = rownames(db.val.cr)
head(db.val.cr[, c(1:6, length(db.val.cr))])

db.val.cr = db.val.cr %>%
  select(c("id_spec"), everything())

#Merging DFs
colnames(db.val[, 1:30])
colnames(db.val.cr[, 1:30])

db.val.cr = cbind(db.val[,1:(j-1)], db.val.cr[,2:length(db.val.cr)])
head(db.val.cr[, c(1:6, length(db.val.cr))])

#_________________________________----
#
# From now on the models are going to be built
#
#_________________________________
# Random Forest modeling and prediction ----
#
ctrl = caret::trainControl(method = 'cv', number = 10, search = 'grid')
cores = parallel::detectCores()-2
#
#_________________________________
#
# ElementCode = Element code equal to the one found in the dataframe
# MVXX = Maximum Value of Graph, based on max value of element column
# DXX = Distance X of legend, based on size of the graph (MVXX)
# DYX1, DYX2, DYX3, DYX4, DYX5 = Distance Y of legend, based on size of the graph (MVXX)
#_________________________________
#
#
for (w in i:j){
  ElementCode = c(colnames(db.ok[w])) 
  
  for (k in 1:length(db.ok)){
    if (names(db.ok[k]) != ElementCode){
      next
    }
    MVXX <- max(db.ok[k])
    ktrain <- k-(((length(db.ok)-(j+1))/2)-(length(db.train))+(j+1))
  }
  
  DXX = MVXX*0.14
  DYX1 = MVXX*0.98
  DYX2 = MVXX*0.95
  DYX3 = MVXX*0.92
  DYX4 = MVXX*0.89
  DYX5 = MVXX*0.86
  
  # Modeling and prediction of element concentration ----
  
  ## Original data (raw spectra)----
  cl = makePSOCKcluster(cores)
  doParallel::registerDoParallel(cl)
  colnames(db.train[, 1:30])
  
  system.time(Orig.Rf.Model <- train(x=db.train[ ,j:length(names(db.train))], 
                                     y=db.train[,ktrain],
                                     method = "rf", 
                                     trControl = ctrl,
                                     importance = T))  
  stopCluster(cl)
  
  dir.create(paste0(dir, "models"))
  saveRDS(Orig.Rf.Model, file = paste0(dir, 'models/', paste(ElementCode, "Orig.Rf.Model.rda", sep=".")))
  
  print(Orig.Rf.Model)
  assign(paste(ElementCode, "Orig.Rf.Model", sep="."), Orig.Rf.Model)
  
  Orig.RF.Model.Variables = varImp(Orig.Rf.Model, scale = T)

  
  orig.dbval = length(db.val)+1
  colnames(db.val[, 1:30])
  db.val[orig.dbval] = stats::predict(Orig.Rf.Model, db.val[, j:length(names(db.val))])
  # Predicting using externaldata
  
  colnames(db.val)[orig.dbval] = c(paste(ElementCode, "Orig", sep="_"))
  
  metric.orig = goof(observed = db.val[,ktrain], 
                     predicted = db.val[,length(db.val)])
  
  jpeg(paste0(dir, 'figures/', paste(ElementCode, '_Orig_Obs_Pred.jpeg', sep="")), 
       width = 10, height = 10, units = 'in', res = 300)
  par(mar = c(6, 6, 6, 6))
  summary(db.val[, c(ktrain, length(db.val))])
  plot(db.val[,orig.dbval], db.val[,ktrain],
       xlab = expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
       ylab = expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
       xlim = c(0, MVXX),
       ylim = c(0, MVXX),
       type = 'p',
       pch = 16,
       cex.lab= 1.1,
       cex.axis=1.2,
       col = rgb(red=.5 , green = .5, blue = .5, alpha = .7),
       family = 'A')
  abline(0, 1, col='black', lty=1, lwd=2)
  abline(lm(db.val[,ktrain] ~ db.val[,orig.dbval]), col = 'red', lty=1, lwd=2)
  text(DXX, DYX1, bquote(RMSE == .(round(metric.orig$RMSE, 2))))
  text(DXX, DYX2, bquote(R['adj']^2 == .(round(metric.orig$R2, 2))))
  text(DXX, DYX3, bquote(CCC == .(round(metric.orig$concordance, 2))))
  text(DXX, DYX4, bquote(bias == .(round(metric.orig$bias, 2))))
  text(DXX, DYX5, bquote(RPIQ == .(round(metric.orig$RPIQ, 2))))
  dev.off()
  
  ## CR data----
  cl = makePSOCKcluster(cores)
  doParallel::registerDoParallel(cl)
  colnames(db.train.cr[, 1:30])
  
  system.time(CR.Rf.Model <- train(x=db.train.cr[ ,j:length(names(db.train.cr))], 
                                   y=db.train.cr[,ktrain],
                                   method = "rf", 
                                   trControl = ctrl,
                                   importance = T))  
  stopCluster(cl)
  
  dir.create(paste0(dir, "models"))
  saveRDS(CR.Rf.Model, file = paste0(dir, 'models/', paste(ElementCode, "CR.Rf.Model.rda", sep=".")))
  
  print(CR.Rf.Model)
  assign(paste(ElementCode, "CR.Rf.Model", sep="."), CR.Rf.Model)
  
  CR.RF.Model.Variables = varImp(CR.Rf.Model, scale = T)
  
  CR.dbval = length(db.val.cr)+1
  colnames(db.val.cr[, 1:30])
  db.val.cr[CR.dbval] = stats::predict(CR.Rf.Model, db.val.cr[, j:length(names(db.val.cr))])
  # Predicting using externaldata
  
  colnames(db.val.cr)[CR.dbval] = c(paste(ElementCode, "CR", sep="_"))
  
  metric.CR = goof(observed = db.val.cr[,ktrain], 
                   predicted = db.val.cr[,length(db.val.cr)])
  
  jpeg(paste0(dir, 'figures/', paste(ElementCode, '_CR_Obs_Pred.jpeg', sep="")), 
       width = 10, height = 10, units = 'in', res = 300)
  par(mar = c(6, 6, 6, 6))
  summary(db.val.cr[, c(ktrain, length(db.val.cr))])
  plot(db.val.cr[,CR.dbval], db.val.cr[,ktrain],
       xlab = expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
       ylab = expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
       xlim = c(0, MVXX),
       ylim = c(0, MVXX),
       type = 'p',
       pch = 16,
       cex.lab= 1.1,
       cex.axis=1.2,
       col = rgb(red=.5 , green = .5, blue = .5, alpha = .7),
       family = 'A')
  abline(0, 1, col='black', lty=1, lwd=2)
  abline(lm(db.val.cr[,ktrain] ~ db.val.cr[,CR.dbval]), col = 'red', lty=1, lwd=2)
  text(DXX, DYX1, bquote(RMSE == .(round(metric.CR$RMSE, 2))))
  text(DXX, DYX2, bquote(R['adj']^2 == .(round(metric.CR$R2, 2))))
  text(DXX, DYX3, bquote(CCC == .(round(metric.CR$concordance, 2))))
  text(DXX, DYX4, bquote(bias == .(round(metric.CR$bias, 2))))
  text(DXX, DYX5, bquote(RPIQ == .(round(metric.CR$RPIQ, 2))))
  dev.off()
  
  ## Savitzky-Golay data----
  cl = makePSOCKcluster(cores)
  doParallel::registerDoParallel(cl)
  colnames(db.train.sg[, 1:30])
  
  system.time(SG.Rf.Model <- train(x=db.train.sg[ ,j:length(names(db.train.sg))], 
                                   y=db.train.sg[,ktrain],
                                   method = "rf", 
                                   trControl = ctrl,
                                   importance = T))  
  stopCluster(cl)
  
  dir.create(paste0(dir, "models"))
  saveRDS(SG.Rf.Model, file = paste0(dir, 'models/', paste(ElementCode, "SG.Rf.Model.rda", sep=".")))
  
  print(SG.Rf.Model)
  assign(paste(ElementCode, "SG.Rf.Model", sep="."), SG.Rf.Model)
  
  SG.RF.Model.Variables = varImp(SG.Rf.Model, scale = T)

  SG.dbval = length(db.val.sg)+1
  colnames(db.val.sg[, 1:30])
  db.val.sg[SG.dbval] = stats::predict(SG.Rf.Model, db.val.sg[, j:length(names(db.val.sg))])
  # Predicting using externaldata
  
  colnames(db.val.sg)[SG.dbval] = c(paste(ElementCode, "SG", sep="_"))
  
  metric.SG = goof(observed = db.val.sg[,ktrain], 
                   predicted = db.val.sg[,length(db.val.sg)])
  
  jpeg(paste0(dir, 'figures/', paste(ElementCode, '_SG_Obs_Pred.jpeg', sep="")), 
       width = 10, height = 10, units = 'in', res = 300)
  par(mar = c(6, 6, 6, 6))
  summary(db.val.sg[, c(ktrain, length(db.val.sg))])
  plot(db.val.sg[,SG.dbval], db.val.sg[,ktrain],
       xlab = expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
       ylab = expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
       xlim = c(0, MVXX),
       ylim = c(0, MVXX),
       type = 'p',
       pch = 16,
       cex.lab= 1.1,
       cex.axis=1.2,
       col = rgb(red=.5 , green = .5, blue = .5, alpha = .7),
       family = 'A')
  abline(0, 1, col='black', lty=1, lwd=2)
  abline(lm(db.val.sg[,ktrain] ~ db.val.sg[,SG.dbval]), col = 'red', lty=1, lwd=2)
  text(DXX, DYX1, bquote(RMSE == .(round(metric.SG$RMSE, 2))))
  text(DXX, DYX2, bquote(R['adj']^2 == .(round(metric.SG$R2, 2))))
  text(DXX, DYX3, bquote(CCC == .(round(metric.SG$concordance, 2))))
  text(DXX, DYX4, bquote(bias == .(round(metric.SG$bias, 2))))
  text(DXX, DYX5, bquote(RPIQ == .(round(metric.SG$RPIQ, 2))))
  dev.off()
  
  ## Standard Normal Variate data----
  cl = makePSOCKcluster(cores)
  doParallel::registerDoParallel(cl)
  colnames(db.train.snv[, 1:30])
  
  system.time(SNV.Rf.Model <- train(x=db.train.snv[ ,j:length(names(db.train.snv))], 
                                    y=db.train.snv[,ktrain],
                                    method = "rf", 
                                    trControl = ctrl,
                                    importance = T))  
  stopCluster(cl)
  
  dir.create(paste0(dir, "models"))
  saveRDS(SNV.Rf.Model, file = paste0(dir, 'models/', paste(ElementCode, "SNV.Rf.Model.rda", sep=".")))
  
  print(SNV.Rf.Model)
  assign(paste(ElementCode, "SNV.Rf.Model", sep="."), SNV.Rf.Model)
  
  SNV.RF.Model.Variables = varImp(SNV.Rf.Model, scale = T)
  
  SNV.dbval = length(db.val.snv)+1
  colnames(db.val.snv[, 1:30])
  db.val.snv[SNV.dbval] = stats::predict(SNV.Rf.Model, db.val.snv[, j:length(names(db.val.snv))])
  # Predicting using externaldata
  
  colnames(db.val.snv)[SNV.dbval] = c(paste(ElementCode, "SNV", sep="_"))
  
  metric.SNV = goof(observed = db.val.snv[,ktrain], 
                    predicted = db.val.snv[,length(db.val.snv)])
  
  jpeg(paste0(dir, 'figures/', paste(ElementCode, '_SNV_Obs_Pred.jpeg', sep="")), 
       width = 10, height = 10, units = 'in', res = 300)
  par(mar = c(6, 6, 6, 6))
  summary(db.val.snv[, c(ktrain, length(db.val.snv))])
  plot(db.val.snv[,SNV.dbval], db.val.snv[,ktrain],
       xlab = expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
       ylab = expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
       xlim = c(0, MVXX),
       ylim = c(0, MVXX),
       type = 'p',
       pch = 16,
       cex.lab= 1.1,
       cex.axis=1.2,
       col = rgb(red=.5 , green = .5, blue = .5, alpha = .7),
       family = 'A')
  abline(0, 1, col='black', lty=1, lwd=2)
  abline(lm(db.val.snv[,ktrain] ~ db.val.snv[,SNV.dbval]), col = 'red', lty=1, lwd=2)
  text(DXX, DYX1, bquote(RMSE == .(round(metric.SNV$RMSE, 2))))
  text(DXX, DYX2, bquote(R['adj']^2 == .(round(metric.SNV$R2, 2))))
  text(DXX, DYX3, bquote(CCC == .(round(metric.SNV$concordance, 2))))
  text(DXX, DYX4, bquote(bias == .(round(metric.SNV$bias, 2))))
  text(DXX, DYX5, bquote(RPIQ == .(round(metric.SNV$RPIQ, 2))))
  dev.off()
  
  ## Detrend normalization data----
  cl = makePSOCKcluster(cores)
  doParallel::registerDoParallel(cl)
  colnames(db.train.det[, 1:30])
  
  system.time(DET.Rf.Model <- train(x=db.train.det[ ,j:length(names(db.train.det))], 
                                    y=db.train.det[,ktrain],
                                    method = "rf", 
                                    trControl = ctrl,
                                    importance = T))  
  stopCluster(cl)
  
  dir.create(paste0(dir, "models"))
  saveRDS(DET.Rf.Model, file = paste0(dir, 'models/', paste(ElementCode, "DET.Rf.Model.rda", sep=".")))
  
  print(DET.Rf.Model)
  assign(paste(ElementCode, "DET.Rf.Model", sep="."), DET.Rf.Model)
  
  DET.RF.Model.Variables = varImp(DET.Rf.Model, scale = T)
  
  DET.dbval = length(db.val.det)+1
  colnames(db.val.det[, 1:30])
  db.val.det[DET.dbval] = stats::predict(DET.Rf.Model, db.val.det[, j:length(names(db.val.det))])
  # Predicting using externaldata
  
  colnames(db.val.det)[DET.dbval] = c(paste(ElementCode, "DET", sep="_"))
  
  metric.DET = goof(observed = db.val.det[,ktrain], 
                    predicted = db.val.det[,length(db.val.det)])
  
  jpeg(paste0(dir, 'figures/', paste(ElementCode, '_DET_Obs_Pred.jpeg', sep="")), 
       width = 10, height = 10, units = 'in', res = 300)
  par(mar = c(6, 6, 6, 6))
  summary(db.val.det[, c(ktrain, length(db.val.det))])
  plot(db.val.det[,DET.dbval], db.val.det[,ktrain],
       xlab = expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
       ylab = expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
       xlim = c(0, MVXX),
       ylim = c(0, MVXX),
       type = 'p',
       pch = 16,
       cex.lab= 1.1,
       cex.axis=1.2,
       col = rgb(red=.5 , green = .5, blue = .5, alpha = .7),
       family = 'A')
  abline(0, 1, col='black', lty=1, lwd=2)
  abline(lm(db.val.det[,ktrain] ~ db.val.det[,DET.dbval]), col = 'red', lty=1, lwd=2)
  text(DXX, DYX1, bquote(RMSE == .(round(metric.DET$RMSE, 2))))
  text(DXX, DYX2, bquote(R['adj']^2 == .(round(metric.DET$R2, 2))))
  text(DXX, DYX3, bquote(CCC == .(round(metric.DET$concordance, 2))))
  text(DXX, DYX4, bquote(bias == .(round(metric.DET$bias, 2))))
  text(DXX, DYX5, bquote(RPIQ == .(round(metric.DET$RPIQ, 2))))
  dev.off()
  
  ## Multiplicative Scatter Correction data----
  cl = makePSOCKcluster(cores)
  doParallel::registerDoParallel(cl)
  colnames(db.train.msc[, 1:30])
  
  system.time(MSC.Rf.Model <- train(x=db.train.msc[ ,j:length(names(db.train.msc))], 
                                    y=db.train.msc[,ktrain],
                                    method = "rf", 
                                    trControl = ctrl,
                                    importance = T))  
  stopCluster(cl)
  
  dir.create(paste0(dir, "models"))
  saveRDS(MSC.Rf.Model, file = paste0(dir, 'models/', paste(ElementCode, "MSC.Rf.Model.rda", sep=".")))
  
  print(MSC.Rf.Model)
  assign(paste(ElementCode, "MSC.Rf.Model", sep="."), MSC.Rf.Model)
  
  MSC.RF.Model.Variables = varImp(MSC.Rf.Model, scale = T)
  
  MSC.dbval = length(db.val.msc)+1
  colnames(db.val.msc[, 1:30])
  db.val.msc[MSC.dbval] = stats::predict(MSC.Rf.Model, db.val.msc[, j:length(names(db.val.msc))])
  # Predicting using externaldata
  
  colnames(db.val.msc)[MSC.dbval] = c(paste(ElementCode, "MSC", sep="_"))
  
  metric.MSC = goof(observed = db.val.msc[,ktrain], 
                    predicted = db.val.msc[,length(db.val.msc)])
  
  jpeg(paste0(dir, 'figures/', paste(ElementCode, '_MSC_Obs_Pred.jpeg', sep="")), 
       width = 10, height = 10, units = 'in', res = 300)
  par(mar = c(6, 6, 6, 6))
  summary(db.val.msc[, c(ktrain, length(db.val.msc))])
  plot(db.val.msc[,MSC.dbval], db.val.msc[,ktrain],
       xlab = expression(paste('Predicted ', '(', 'mg', ' kg'^'-1',')')),
       ylab = expression(paste('Observed ', '(', 'mg', ' kg'^'-1',')')),
       xlim = c(0, MVXX),
       ylim = c(0, MVXX),
       type = 'p',
       pch = 16,
       cex.lab= 1.1,
       cex.axis=1.2,
       col = rgb(red=.5 , green = .5, blue = .5, alpha = .7),
       family = 'A')
  abline(0, 1, col='black', lty=1, lwd=2)
  abline(lm(db.val.msc[,ktrain] ~ db.val.msc[,MSC.dbval]), col = 'red', lty=1, lwd=2)
  text(DXX, DYX1, bquote(RMSE == .(round(metric.MSC$RMSE, 2))))
  text(DXX, DYX2, bquote(R['adj']^2 == .(round(metric.MSC$R2, 2))))
  text(DXX, DYX3, bquote(CCC == .(round(metric.MSC$concordance, 2))))
  text(DXX, DYX4, bquote(bias == .(round(metric.MSC$bias, 2))))
  text(DXX, DYX5, bquote(RPIQ == .(round(metric.MSC$RPIQ, 2))))
  dev.off()
  
  # Combining all model internal calibration results into a table ----
  
  orig = Orig.Rf.Model$results[Orig.Rf.Model$results$mtry %in% Orig.Rf.Model$bestTune, ]
  orig$Variable = paste(ElementCode); orig$Preproc = 'Raw Spectra'
  
  snv = SNV.Rf.Model$results[SNV.Rf.Model$results$mtry %in% SNV.Rf.Model$bestTune, ]
  snv$Variable = paste(ElementCode); snv$Preproc = 'Standard Normal Variate'
  
  cr = CR.Rf.Model$results[CR.Rf.Model$results$mtry %in% CR.Rf.Model$bestTune, ]
  cr$Variable = paste(ElementCode); cr$Preproc = 'Continuum Removal'
  
  det = DET.Rf.Model$results[DET.Rf.Model$results$mtry %in% DET.Rf.Model$bestTune, ]
  det$Variable = paste(ElementCode); det$Preproc = 'Detrend normalization'
  
  sg = SG.Rf.Model$results[SG.Rf.Model$results$mtry %in% SG.Rf.Model$bestTune, ]
  sg$Variable = paste(ElementCode); sg$Preproc = 'Savitzky-Golay'
  
  msc = MSC.Rf.Model$results[MSC.Rf.Model$results$mtry %in% MSC.Rf.Model$bestTune, ]
  msc$Variable = paste(ElementCode); msc$Preproc = 'Multiplicative Scatter Correction'
  
  cal.res = rbind(orig, snv, cr, det, sg, msc)
  
  head(cal.res)
  
  cal.res = cal.res %>%
    select(c('Variable', 'Preproc'), everything())
  
  write.csv(cal.res, paste0(dir, 'Tables/calibration/', paste(ElementCode, "calibration.csv", sep="_")), 
            row.names = FALSE)
  
  # Combining all model external validation results into a table ----
  
  metric.orig$Variable = paste(ElementCode)
  metric.orig$Preproc = 'Raw Spectra'
  metric.SG$Variable = paste(ElementCode)
  metric.SG$Preproc = 'Savitzky-Golay'
  metric.CR$Variable = paste(ElementCode)
  metric.CR$Preproc = 'Continuum Removal'
  metric.DET$Variable = paste(ElementCode)
  metric.DET$Preproc = 'Detrend Normalization'
  metric.SNV$Variable = paste(ElementCode)
  metric.SNV$Preproc = 'Standard Normal Variate'
  metric.MSC$Variable = paste(ElementCode)
  metric.MSC$Preproc = 'Multiplicative Scatter Correction'
  
  res = rbind(metric.orig, metric.SG, metric.CR, metric.DET, metric.SNV, 
              metric.MSC)
  
  head(res, 5)
  res= res %>%
    select(c('Variable', 'Preproc'), everything())
  
  write.csv(res, paste0(dir, 'Tables/validation/', paste(ElementCode, "validation.csv", sep="_")), 
            row.names = FALSE)
  
}
# Merging all csv files into calibration and validation csv files (for more than one element) ----
library(dplyr)
library(readr)
calibration <- list.files(path=paste0(dir, 'Tables/calibration'), full.names = TRUE) %>% 
  lapply(read.csv) %>% 
  bind_rows
write.csv(calibration,paste0(dir, 'calibration.csv'), row.names = FALSE)

validation <- list.files(path=paste0(dir, 'Tables/validation'), full.names = TRUE) %>% 
  lapply(read.csv) %>% 
  bind_rows
write.csv(validation, paste0(dir, 'validation.csv'), row.names = FALSE)
### ------------------ ###
### ------------------ ###