#Data from Blankenship 2011, assumed to be in same directory as this file
source("processing_functions.R")

dt=0.1
#Store data separately by phenotype- can be done in one big data frame, but this is more readable

WT_files=c("Blankenship2011_WT_01","Blankenship2011_WT_02","Blankenship2011_WT_03","Blankenship2011_WT_04","Blankenship2011_WT_05")

cx45_files=c("Blankenship2011_cx45_01","Blankenship2011_cx45_02","Blankenship2011_cx45_03","Blankenship2011_cx45_04")

cx3645_files=c("Blankenship2011_cx36_45_01","Blankenship2011_cx36_45_02","Blankenship2011_cx36_45_03","Blankenship2011_cx36_45_04","Blankenship2011_cx36_45_05","Blankenship2011_cx36_45_06")



#initiate data frames
WT=data.frame(ci=integer(0),sttc=integer(0),dist=integer(0))

cx45=data.frame(ci=integer(0),sttc=integer(0),dist=integer(0))

cx3645=data.frame(ci=integer(0),sttc=integer(0),dist=integer(0))

#read in each data file, calculate measures and extract vectors of measures and electrode separation, merge data files


for(i in 1:length(WT_files)){
	WT=mapply(c,WT,run_measures_on_hdf5(WT_files[i],dt),SIMPLIFY=FALSE)
}

for(i in 1:length(cx45_files)){
	cx45=mapply(c,cx45,run_measures_on_hdf5(cx45_files[i],dt),SIMPLIFY=FALSE)
}


for(i in 1:length(cx3645_files)){
	cx3645=mapply(c,cx3645,run_measures_on_hdf5(cx3645_files[i],dt),SIMPLIFY=FALSE)
}

#convert distance to factor for analysis

WT$dist=as.factor(WT$dist)
cx45$dist=as.factor(cx45$dist)
cx3645$dist=as.factor(cx3645$dist)

#get medians and quartiles

WT_quarts=list()
cx45_quarts=list()
cx3645_quarts=list()

WT_quarts$med=quantiles_process(WT,0.5)
WT_quarts$low=quantiles_process(WT,0.25)
WT_quarts$up=quantiles_process(WT,0.75)

cx45_quarts$med=quantiles_process(cx45,0.5)
cx45_quarts$low=quantiles_process(cx45,0.25)
cx45_quarts$up=quantiles_process(cx45,0.75)


cx3645_quarts$med=quantiles_process(cx3645,0.5)
cx3645_quarts$low=quantiles_process(cx3645,0.25)
cx3645_quarts$up=quantiles_process(cx3645,0.75)

#############################
#Scale cx45 and cx3645 medians to have same initial value as WT for ci
#for figure insert

scaled_3645=cx3645_quarts$med[[1]][,2]
scaled_45=cx45_quarts$med[[1]][,2]

#scaling factor
SF3645=WT_quarts$med[[1]][1,2]/cx3645_quarts$med[[1]][1,2]
SF45=WT_quarts$med[[1]][1,2]/cx45_quarts$med[[1]][1,2]

scaled_3645=scaled_3645*SF3645
scaled_45=scaled_45*SF45

############################
#get ranges in order to set graph limits

ci_range=range(c(WT_quarts$low[[1]][,2],WT_quarts$up[[1]][,2],cx3645_quarts$low[[1]][,2],cx3645_quarts$up[[1]][,2],cx45_quarts$low[[1]][,2],cx45_quarts$up[[1]][,2]),na.rm=TRUE)

sttc_range=range(c(WT_quarts$low[[2]][,2],WT_quarts$up[[2]][,2],cx3645_quarts$low[[2]][,2],cx3645_quarts$up[[2]][,2],cx45_quarts$low[[2]][,2],cx45_quarts$up[[2]][,2]),na.rm=TRUE)

max_x=max(as.vector(WT_quarts$med[[1]][,1]))
max_x=as.numeric(max_x)



##########################
#Read in spike trains for raster plots 

#WT
a1=h5.read.spikes("Blankenship2011_WT_05.h5")


#Cx45ko
a2=h5.read.spikes("Blankenship2011_cx45_03.h5")


#Cx 36/45 ko
a3=h5.read.spikes("Blankenship2011_cx36_45_01.h5")


#########################
#plot graphs

attach(mtcars)

postscript("Figure_7.eps",width=inch(17.6),height=inch(17.6))

layout(matrix(c(1,1,2,2,3,3,4,5,4,5),5,2,byrow=TRUE))
par(mar=c(1.1,5.1,1.1,1.1))
 par(oma=c(1.8,1.8,1.80,1.8))

#plot rasters   
plot(a1,beg=0,end=600,which=c(31:40),main="",ylab="Wild Type",xlab="",xaxt="n",cex.lab=1.5)
segments(0,1.03,60,1.03,lwd=3,lend="butt")
text(x=40,y=0.89,"1 min",cex=1.5)

plot(a2,beg=0,end=600,which=c(65:75),main="",ylab="Cx45 KO",xlab="",xaxt="n",cex.lab=1.5, col.lab="red")



plot(a3,beg=2500,end=3100,which=c(51:60),main="",ylab="Cx36/Cx45 dKO",,xlab="",xaxt="n",cex.lab=1.5,col.lab="blue")


#plot measure vs distance plots
par(mar=c(4.1,5.1,1.1,1.5))

#note the range of y values- more visually appealing than strict range
plot(as.vector(WT_quarts$med[[1]][,1]),as.vector(WT_quarts$med[[1]][,2]),type="l",lwd=2,col="black",ylim=c(0,65),main="",xlab=expression(paste("Intercell distance (",mu,"m)",sep="")),ylab="Correlation Index",xlim=c(-5,max_x+10),bty="n",cex.lab=1.5,yaxt="n")
axis(2,las=2,at=c(0,15,30,45,60))
arrows(as.numeric(as.vector(WT_quarts$low[[1]][,1])),as.numeric(as.vector(WT_quarts$low[[1]][,2])),as.numeric(as.vector(WT_quarts$up[[1]][,1])),as.numeric(as.vector(WT_quarts$up[[1]][,2])),col="black",angle=90,code=3,length=0.03,lwd=0.5)
 
lines(as.vector(cx3645_quarts$med[[1]][,1]),as.vector(cx3645_quarts$med[[1]][,2]),col="blue",lwd=2)
arrows(as.numeric(as.vector(cx3645_quarts$low[[1]][,1])),as.numeric(as.vector(cx3645_quarts$low[[1]][,2])),as.numeric(as.vector(cx3645_quarts$up[[1]][,1])),as.numeric(as.vector(cx3645_quarts$up[[1]][,2])),col="blue",angle=90,code=3,length=0.03,lwd=0.5)

lines(as.vector(cx45_quarts$med[[1]][,1]),as.vector(cx45_quarts$med[[1]][,2]),col="red",lwd=2)
arrows(as.numeric(as.vector(cx45_quarts$low[[1]][,1])),as.numeric(as.vector(cx45_quarts$low[[1]][,2])),as.numeric(as.vector(cx45_quarts$up[[1]][,1])),as.numeric(as.vector(cx45_quarts$up[[1]][,2])),col="red",angle=90,code=3,length=0.03,lwd=0.5)


#second measure vs distance plot
plot(as.vector(WT_quarts$med[[2]][,1]),as.vector(WT_quarts$med[[2]][,2]),type="l",lwd=2,col="black",ylim=c(sttc_range[1],0.6),main="",xlab=expression(paste("Intercell distance (",mu,"m)",sep="")),ylab="STTC",xlim=c(-5,max_x+10),bty="n",cex.lab=1.5,yaxt="n")

axis(2,las=2,at=c(0,0.15,0.3,0.45,0.6))
arrows(as.numeric(as.vector(WT_quarts$low[[2]][,1])),as.numeric(as.vector(WT_quarts$low[[2]][,2])),as.numeric(as.vector(WT_quarts$up[[2]][,1])),as.numeric(as.vector(WT_quarts$up[[2]][,2])),col="black",angle=90,code=3,length=0.03,lwd=0.5)


lines(as.vector(cx3645_quarts$med[[2]][,1]),as.vector(cx3645_quarts$med[[2]][,2]),col="blue",lwd=2)
arrows(as.numeric(as.vector(cx3645_quarts$low[[2]][,1])),as.numeric(as.vector(cx3645_quarts$low[[2]][,2])),as.numeric(as.vector(cx3645_quarts$up[[2]][,1])),as.numeric(as.vector(cx3645_quarts$up[[2]][,2])),col="blue",angle=90,code=3,length=0.03,lwd=0.5)

lines(as.vector(cx45_quarts$med[[2]][,1]),as.vector(cx45_quarts$med[[2]][,2]),col="red",lwd=2)
arrows(as.numeric(as.vector(cx45_quarts$low[[2]][,1])),as.numeric(as.vector(cx45_quarts$low[[2]][,2])),as.numeric(as.vector(cx45_quarts$up[[2]][,1])),as.numeric(as.vector(cx45_quarts$up[[2]][,2])),col="red",angle=90,code=3,length=0.03,lwd=0.5)

legend(x=175,y=0.60,col=c("black","red","blue"),lty=1,lwd=2,legend=c("Wild Type      (0.31 Hz, n=5)","Cx45 ko        (1.06 Hz, n=4)","Cx36/45 dko (1.99 Hz, n=6)"))

#label plots
text(grconvertX(c(0.06, 0.06, 0.5), from='ndc'),
   grconvertY(c(0.97, 0.38, 0.38), from='ndc'),  c('A',  'B',  'C'), xpd=NA, cex=2, font=2)

#plot insert

par(fig = c(.30, .45, .23, .38), mar=c(0,0,0,0), new=TRUE)

plot(as.vector(WT_quarts$med[[1]][,1]),as.vector(WT_quarts$med[[1]][,2]),type="l",lwd=2,col="black",ylim=c(0,50),main="",xlab="",ylab="",xlim=c(-5,max_x+10),bty="n",yaxt="n",xaxt="n")

Axis(side=1, labels=FALSE,tck=FALSE)
Axis(side=2, labels=FALSE,tck=FALSE)

lines(as.vector(cx3645_quarts$med[[1]][,1]),scaled_3645,col="blue",lwd=2)
lines(as.vector(cx45_quarts$med[[1]][,1]),scaled_45,col="red",lwd=2)


dev.off()





