#Troy Masters 2014

#set working directory
setwd('C:/Climate/CombiningEvidence')

#resolution to consider for distribution
resolution<-seq(0,10,by=0.05)

#get density for aldrin
aldrin_data<-read.table('aldrin_digitized.txt',sep=',')
aldrin_int<-approx(x=aldrin_data[[1]], y=aldrin_data[[2]], xout=resolution)$y
aldrin_int[is.na(aldrin_int) | aldrin_int < 0]<-0

#get density for lewis and curry
lc_data<-read.table('lewis_and_curry_digitized.txt',sep=',')
lc_int<-approx(x=lc_data[[1]], y=lc_data[[2]], xout=resolution)$y
lc_int[is.na(lc_int) | lc_int < 0]<-0

#get density for masters
masters_data<-read.table('masters_digitized.txt',sep=',')
masters_int<-approx(x=masters_data[[1]], y=masters_data[[2]], xout=resolution)$y
masters_int[is.na(masters_int) | masters_int < 0]<-0

#get average of instrumental distributions
instrumental_avg<-(aldrin_int+lc_int+masters_int)/3

#get density for paleosens
paleo_data<-read.table('paleosens_digitized.txt',sep=',')
#must switch from K / Wm^2 to 2xCO2 doubling
paleo_int<-approx(x=paleo_data[[1]]*3.7, y=paleo_data[[2]], xout=resolution)$y
paleo_int[is.na(paleo_int) | paleo_int < 0]<-0

#draw main plot
plot(x=resolution,y=aldrin_int/sum(aldrin_int*.05), type='l', ylim=c(0,1.0), xlim=c(0,7),lty=2, 
	ylab="Density", xlab="ECS (K)")
points(x=resolution,y=lc_int/sum(lc_int*.05), type='l',lty=3)
points(x=resolution,y=masters_int/sum(masters_int*.05), type='l',lty=4)
points(x=resolution,y=paleo_int/sum(paleo_int*.05), type='l',lwd=2,col='blue')
points(x=resolution,y=(instrumental_avg)/sum(instrumental_avg*.05), type='l',lwd=2,col='black')
points(x=resolution,y=(instrumental_avg*paleo_int)/sum(instrumental_avg*paleo_int*.05), type='l',lwd=5,col='red')
legend("topright", legend=c("Aldrin et al. (2012)", "Lewis and Curry (2014)", "Masters (2014)", "PALAEOSENS (2012)", "Instrumental Mean", "Combined Palaeo+Instrumental"),
	col=c("black","black","black","blue","black","red"), lty=c(2,3,4,1,1,1), lwd=c(1,1,1,2,2,3))
	

#function to find percentiles within the given distribution
findPercentile<-function(distribution, percentile)
{
	for (i in 1:length(distribution))
	{
		if (sum(distribution[1:i]) >= percentile)
		{
			return(resolution[i])
		}
	}
}

#find confidence ranges and median
low<-findPercentile((instrumental_avg*paleo_int)/sum(instrumental_avg*paleo_int), 0.025)
r_lo<-findPercentile((instrumental_avg*paleo_int)/sum(instrumental_avg*paleo_int), 0.16)
med<-findPercentile((instrumental_avg*paleo_int)/sum(instrumental_avg*paleo_int), 0.5)
r_hi<-findPercentile((instrumental_avg*paleo_int)/sum(instrumental_avg*paleo_int), 0.84)
high<-findPercentile((instrumental_avg*paleo_int)/sum(instrumental_avg*paleo_int), 0.975)
arrows(x0=low, y0=0, x1=high, y1=0, lwd=1, length=0.15, code=3, col="red")
arrows(x0=r_lo, y0=0, x1=r_hi, y1=0, lwd=3, length=0.15, code=3, col="red")
points(x=med,y=0, lwd=5, col="red")