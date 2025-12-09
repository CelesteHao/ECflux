# the processing function for NEE partitioning
# try to be same to REddyPro online tool
# for Hanford 300A site
# Copy from https://github.com/EarthyScience/REddyProc/blob/master/vignettes/useCase.md
# Create date 2025/12/08
# change the u* percentage as c(0.05, 0.16, 0.25, 0.5, 0.75, 0.84, 0.95)



#load libraries used in this vignette
#+++ load libraries used in this vignette
library(REddyProc)
library(dplyr)

#+++ Load data with 1 header and 1 unit row from (tab-delimited) text file
fileName <- "C:/split/Hn1_input.txt"
EddyData <- fLoadTXTIntoDataframe(fileName)

# calculate VPD
EddyData$VPD <- fCalcVPDfromRHandTair(EddyData$rh, EddyData$Tair)
#+++ Replace long runs of equal NEE values by NA
EddyData <- filterLongRuns(EddyData, "NEE")

#+++ Add time stamp in POSIX time format
EddyDataWithPosix <- fConvertTimeToPosix(
  EddyData, 'YDH',Year = 'Year',Day = 'DoY', Hour = 'Hour') %>% 
  filterLongRuns("NEE")

#+++ Initalize R5 reference class sEddyProc for post-processing of eddy data
#+++ with the variables needed for post-processing later
EProc <- sEddyProc$new(
  'DE-Tha', EddyDataWithPosix, c('NEE','H','LE','Rg','Tair','VPD', 'Ustar'))

EProc$sPlotFingerprintY('NEE', Year = 2024)

EProc$sEstimateUstarScenarios(
  nSample = 200L, probs = c(0.05, 0.16, 0.25, 0.5, 0.75, 0.84, 0.95))

EProc$sGetEstimatedUstarThresholdDistribution()
EProc$sGetUstarScenarios()
EProc$sMDSGapFillUStarScens('NEE')

grep("NEE_.*_f$",names(EProc$sExportResults()), value = TRUE)
grep("NEE_.*_fsd$",names(EProc$sExportResults()), value = TRUE)

EProc$sSetLocationInfo(
  LatDeg       = 46.41,    
  LongDeg      = -119.28,  
  TimeZoneHour = -8        
)

EProc$sMDSGapFill('Tair', FillAll = FALSE,  minNWarnRunLength = NA)
EProc$sMDSGapFill('VPD',  FillAll = FALSE,  minNWarnRunLength = NA)
EProc$sMDSGapFill('Rg',  FillAll = FALSE,  minNWarnRunLength = NA)
EProc$sMDSGapFill('LE', FillAll = FALSE, minNWarnRunLength = NA)
EProc$sMDSGapFill('H',  FillAll = FALSE, minNWarnRunLength = NA)

EProc$sFillVPDFromDew() # fill longer gaps still present in VPD_f

EProc$sMRFluxPartitionUStarScens() # nighttime method
EProc$sGLFluxPartitionUStarScens() # daytime method
grep("GPP.*_f$|Reco",names(EProc$sExportResults()), value = TRUE)

FilledEddyData <- EProc$sExportResults()
uStarSuffixes <- colnames(EProc$sGetUstarScenarios())[-1]
#suffix <- uStarSuffixes[2]

GPPAggCO2 <- sapply(uStarSuffixes, function(suffix) {
  GPPHalfHour <- FilledEddyData[[paste0("GPP_",suffix,"_f")]]
  mean(GPPHalfHour, na.rm = TRUE)
})

molarMass <- 12.011
GPPAgg <- GPPAggCO2 * 1e-6 * molarMass * 3600*24*365.25
print(GPPAgg)

# EProc$sEstimateUstarScenarios( 
#  nSample = 200, probs = seq(0.025,0.975,length.out = 39) )

FilledEddyData <- EProc$sExportResults()
# CombinedData <- cbind(EddyDataWithPosix, FilledEddyData)
CombinedData <- cbind(
  EddyData,
  DateTime = EddyDataWithPosix$DateTime,  
  FilledEddyData
)

fWriteDataframeToFile(
  CombinedData,
  'output.txt',
  Dir = "C:/split"
)

