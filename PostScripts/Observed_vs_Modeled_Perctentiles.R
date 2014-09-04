#Troy Masters 2014
#Compares observed temperature trends (from Cowtay and Way 2014) to
#CMIP5 RCP4.5 model runs


#Set working directory to the base directory of the cloned CMIPDataProcessing repository
#NOTE: Must set this to be specific on your computer
setwd('C:/Users/Troy/Documents/GitHub/CMIPDataProcessing') 
source('Process_CMIP_Regions.r') #bring in library functions
models<-read.table('models.txt')[[1]] #read our list of models to use

#read all rcp45 files
model_rcp45<-readAllProcessedFiles('processed/rcp45/tas', models)
#read all historical files
model_historical<-readAllProcessedFiles('processed/historical/tas', models)

#match up with RCP4.5 with its appropriate ancestor historical run
matchRcp45ToHist<-function(model, run)
{
	if (model == 'CESM1-CAM5') {
		if (run == 'r1i1p1') return('r3i1p1')
		if (run == 'r2i1p1') return('r1i1p1')
		if (run == 'r3i1p1') return('r2i1p1')
	} else if (model == 'HadCM3') {
		return('r1i1p1') #all these originate at same one
	} else {
		return(run) #all others the runs match up
	}
}

#merge the RCP4.5 runs to their ancestor runs
rcp45_merged<-list()
for (i in 1:length(model_rcp45))
{
	for (j in 1:min(length(model_rcp45[[i]]),5)) #limit to 5 runs per model
	{
		modelName<-names(model_rcp45)[i]
		runName<-names(model_rcp45[[i]])[j]
		histRun<-matchRcp45ToHist(modelName,runName)
		monthlyTS<-combineTS(model_historical[[modelName]][[histRun]],model_rcp45[[i]][[j]])
		#print(paste(modelName, tsp(monthlyTS)))
		annualTS<-aggregate(window(monthlyTS, start=1950, end=c(2020,12)),nfreq=1)/12
		if (length(annualTS) >= 71)
		{
			rcp45_merged[[length(rcp45_merged)+1]]<-annualTS
			ts.plot(annualTS)
		}
	}
}

#grab krig'd hadCRUTv4 from Cowtay and Way (2014)
#NOTE: may need to navigate here first in your browser
CW_annual<-ts(read.table("http://www-users.york.ac.uk/~kdc3/papers/coverage2013/had4_krig_annual_v2_0_0.txt")[[2]],start=1850)

#calculateTrends - determine percentiles for observed temperatures (CW_annual) within the RCP4.5 models runs (rcp45_merged)
# firstStartYr - first year to use for a start year
# lasEndYr - last year to use for an end year
# minTrendLength - the minimum trend length in years
calculateTrends<-function(firstStartYr, lastEndYr, minTrendLength)
{
	firstEndYr<-firstStartYr+minTrendLength-1
	#2D array containing all relevant trends
	summaryYr<-array(NA, dim=c(lastEndYr-firstStartYr-minTrendLength+2,lastEndYr-firstStartYr-minTrendLength+2)) 
	for (startYr in firstStartYr:(lastEndYr-minTrendLength+1))
	{
		print(startYr)
		for (endYr in (startYr+minTrendLength-1):lastEndYr)
		{
			#calculate trends for this time period for all models
			trends<-numeric(length(rcp45_merged))
			for (i in 1:length(rcp45_merged))
			{
				tas.window<-window(rcp45_merged[[i]],start=startYr,end=endYr)
				trends[i]<-coef(lm(tas.window ~ time(tas.window)))[2]
			}
			#calculate trends per Cowtan and Way (2014) kriging
			cw.window<-window(CW_annual,start=startYr,end=endYr)
			cw_trend<-coef(lm(cw.window ~ time(cw.window)))[2]
			#calculate percentile of CW14 within model runs
			pct<-length(trends[trends<cw_trend])/length(trends)*100.0
			summaryYr[startYr-firstStartYr+1,endYr-(firstStartYr+minTrendLength-1)+1]<-pct
		}
	}
	summaryYr
}

#Now calculate the >=15 yr trends
firstStartYr<-1950
lastEndYr<-2013
minTrendLength<-15
yrTrends15<-calculateTrends(firstStartYr,lastEndYr,minTrendLength)

#plot all trends
library(fields)
image.plot(z=yrTrends15, x=firstStartYr:(lastEndYr-minTrendLength+1), y=(firstStartYr+minTrendLength-1):lastEndYr,
	xlab='Start Year', ylab='End Year', main="Observed Temp Trend (15Yr Min) Percentile relative to CMIP5 Runs",
	sub='Observations from CW14, limit max runs to 5 per model from RCP4.5 scenario', zlim=c(0.0,100.0))

#Now calculate the >=30 yr trends
minTrendLength<-30
yrTrends30<-calculateTrends(firstStartYr,lastEndYr,minTrendLength)
#plot all trends
image.plot(z=yrTrends30, x=firstStartYr:(lastEndYr-minTrendLength+1), y=(firstStartYr+minTrendLength-1):lastEndYr,
	xlab='Start Year', ylab='End Year', main="Observed Temp Trend (30Yr Min)Percentile relative to CMIP5 Runs",
	sub='Observations from CW14, limit max runs to 5 per model from RCP4.5 scenario', zlim=c(0,100.0))