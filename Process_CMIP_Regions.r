# Process_CMIP_Regions.R
# Troy Masters
# Library of function used to average CMIP5 data across dimensions
# Must first install ncdf library prior to use

#include ncdf library
library(ncdf)

#The spatial weighting function for a particular area
getWeighting<-function(lat_bnds, lon_bnds, bounds)
{
	weighting<-array(0,dim=c(length(lon_bnds)/2,length(lat_bnds)/2))
	for (lat in 1:length(lat_bnds)/2)
	{
		for (lon in 1:length(lon_bnds)/2)
		{
			weighting[lon,lat]<-(lon_bnds[1,lon]-lon_bnds[2,lon])*pi/180*(sin(lat_bnds[1,lat]*pi/180)-sin(lat_bnds[2,lat]*pi/180))		
		}
	}
	weighting<-weighting/(4*pi)
	top<-bounds[1]
	btm<-bounds[2]
	left<-bounds[3]
	right<-bounds[4]
	topWeight<-(top - lat_bnds[1,])/(lat_bnds[2,] - lat_bnds[1,])
	topWeight[topWeight>1]<-1
	topWeight[topWeight<0]<-0
	btmWeight<-(lat_bnds[2,] - (btm))/(lat_bnds[2,] - lat_bnds[1,])
	btmWeight[btmWeight>1]<-1
	btmWeight[btmWeight<0]<-0
	latWeight<-topWeight*btmWeight
	rightWeight<-(right - lon_bnds[1,])/(lon_bnds[2,] - lon_bnds[1,])
	rightWeight[rightWeight>1]<-1
	rightWeight[rightWeight<0]<-0
	leftWeight<-(lon_bnds[2,] - left)/(lon_bnds[2,] - lon_bnds[1,])
	leftWeight[leftWeight>1]<-1
	leftWeight[leftWeight<0]<-0
	lonWeight<-rightWeight*leftWeight
	totalWeight<-(lonWeight %*% t(latWeight)) * weighting
	totalWeight<-totalWeight / sum(totalWeight)
	totalWeight
}

#combine time series over different periods
combineTS<-function(ts1,ts2)
{
	if (!is.null(ts1)) {
		p<-ts.union(ts1,ts2)
		ret<-pmax(p[,1],p[,2],na.rm=TRUE)
	} else {
		ret<-ts2
	}
	ret
}

# getSingleModelVariableResult - Get the results for a single experiment / model / variable for the specified area bounds
#								 Currently only supports averaging across horizontal grid 
#		experiment - experiment this is for (e.g. "historical" or "historicalGHG") 
#		model - name of the model (e.g. "GFDL-CM3")
#		variable - name of the variable to average 
#		area - square bounds given in degrees of area to average over (top,bottom,left,right)
#		dirPath - top level path we should begin searching within for NC.  Assumes typical CMIP5 format of experiment/table/variable/model/run
# returns - list of time series for relevant experiment / model / variable for each of the runs
getSingleModelVariableResult<-function(experiment, model, variable, area=c(90,-90,0,360), dirPath=".")
{
	#first get a list of end-point directories (those with no subdirectories)
	allDirs<-list.dirs(path = dirPath)
	dirList<-character(0)
	for (i in 1:length(allDirs))
	{
		if (length(list.dirs(path = allDirs[i], recursive = FALSE)) == 0)
			dirList<-c(dirList, allDirs[i])
	}

	results<-list() #get the individual run results
	
	#cycle through each directory to find those relevant for our experiment, model, and variable
	for (dNum in 1:length(dirList))
	{
		filePre<-gsub(pattern="/", replacement="_", x=substring(text=dirList[dNum],first=3))

		#check if this is directory has files relevant to the experiment and model
		if (length(grep(paste(experiment,'_',sep=''), filePre)) > 0 && 
			length(grep(paste(model,'_',sep=''), filePre)) > 0 && 
			length(grep(paste(variable,'_',sep=''), filePre)) > 0)
		{
			fileList<-list.files(path=dirList[dNum], full.names=TRUE)
			varTs<-NULL
			for (fileNo in 1:length(fileList))
			{
				print(fileList[fileNo])
				#get the start year for this file
				reg<-regexpr("_[0-9]+-",fileList[[fileNo]])
				match<-regmatches(fileList[[fileNo]],reg)[1]
				fileStartYr<-as.numeric(substr(match,2,nchar(match)-3))
				fileStartMo<-as.numeric(substr(match,6,nchar(match)-1))
				#print(c(fileStartYr,fileStartMo))
				#open the file
				file.nc<-open.ncdf(fileList[fileNo])
				latNum<-length(file.nc$dim$lat$vals)
				lonNum<-length(file.nc$dim$lon$vals)						
				lat_bnds<-get.var.ncdf(file.nc, 'lat_bnds')
				lon_bnds<-get.var.ncdf(file.nc, 'lon_bnds')
				weighting<-getWeighting(lat_bnds, lon_bnds, area)
				monthsPerFile<-length(file.nc$dim$time$vals)
				monthlyData<-numeric(monthsPerFile)
				for (mon in 1:monthsPerFile)
				{
					data<-get.var.ncdf(file.nc, variable, start=c(1,1,mon), count=c(lonNum,latNum,1))
					monthlyData[mon]<-sum(weighting*data, na.rm=TRUE)
				}
				#print(ts(monthlyData,start=c(fileStartYr,fileStartMo),freq=12))
				varTs<-combineTS(varTs,ts(monthlyData,start=c(fileStartYr,fileStartMo),freq=12))
				#print(varTs)
				close.ncdf(file.nc)
			}
			#figure out the run this is associated with
			if (!is.null(varTs))
			{
				runReg=regexpr("r[0-9]+i[0-9]+p[0-9]+",fileList[[1]])
				match<-regmatches(fileList[[1]],runReg)[1]
				run_data<-match
				results[[run_data]]<-varTs
			}
		}
	}
	results
}

# bulkProcessModels - Process the experiment / variable / area combination for all models specified in the modelFile
#		experiment - experiment this is for (e.g. "historical" or "historicalGHG") 
#		variable - name of the variable to average 
#		area - square bounds given in degrees of area to average over (top,bottom,left,right)
#		modelFile - path to text file where the first row of the table contains the model names to process
#		dirPath - top level path we should begin searching within for NC files.  Assumes typical CMIP5 format of experiment/table/variable/model/run
#		outPath - output directory for the text files 
#		
# returns - nothing

bulkProcessAllModels<-function(experiment, variable, area=c(90,-90,0,360), modelFile='models.txt', dirPath='.', outPath='.')
{
	models<-read.table(modelFile)[[1]]
	for (i in 1:length(models))
	{
		#get the results
		result<-getSingleModelVariableResult(experiment, models[i], variable, area, dirPath=dirPath)
		#output each of the runs
		for (run in 1:length(result))
		{
			df<-ts.union(time(result[[run]]), result[[run]])
			colnames(df)[1]<-'year'
			colnames(df)[2]<-variable
			fileName<-paste(outPath,"/",sep="")
			fileName<-paste(fileName, paste(experiment,variable,models[i],names(result)[run],area[1],area[2],area[3],area[4],sep='_'), sep='')
			fileName<-paste(fileName, ".txt",sep='')
			print(fileName)
			write.table(df, fileName)
		}
	}
}
