# Data from Kirkby et al 2013 data assumed to be in same directory as this file

source("processing_functions.R")

dt=0.1

#Store data separately by phenotype- can be done in one big data frame, but this is more readable

WT_files=c("Kirkby2013_02_WT_P5","Kirkby2013_03_WT_P5","Kirkby2013_04_WT_P5","Kirkby2013_07_WT_P5","Kirkby2013_08_WT_P5","Kirkby2013_10_WT_P5")

B2P4_files=c("Kirkby2013_01_B2KO_P4","Kirkby2013_04_B2KO_P4")

B2P5_files=c("Kirkby2013_02_B2KO_P5")

B2P6_files=c("Kirkby2013_05_B2KO_P6","Kirkby2013_07_B2KO_P6","Kirkby2013_08_B2KO_P6","Kirkby2013_09_B2KO_P6","Kirkby2013_10_B2KO_P6")

B2P7_files=c("Kirkby2013_03_B2KO_P7","Kirkby2013_06_B2KO_P7","Kirkby2013_11_B2KO_P7")

#####################

#initiate data frames
WT=data.frame(ci=integer(0),sttc=integer(0),dist=integer(0))

B2P4=data.frame(ci=integer(0),sttc=integer(0),dist=integer(0))

B2P5=data.frame(ci=integer(0),sttc=integer(0),dist=integer(0))

B2P6=data.frame(ci=integer(0),sttc=integer(0),dist=integer(0))

B2P7=data.frame(ci=integer(0),sttc=integer(0),dist=integer(0))

#read in each data file, calculate measures and extract vectors of measures and electrode separation, merge data files


for(i in 1:length(WT_files)){
	WT=mapply(c,WT,run_measures_on_hdf5(WT_files[i],dt),SIMPLIFY=FALSE)
}

for(i in 1:length(B2P4_files)){
	B2P4=mapply(c,B2P4,run_measures_on_hdf5(B2P4_files[i],dt),SIMPLIFY=FALSE)
}

for(i in 1:length(B2P5_files)){
	B2P5=mapply(c,B2P5,run_measures_on_hdf5(B2P5_files[i],dt),SIMPLIFY=FALSE)
}

for(i in 1:length(B2P6_files)){
	B2P6=mapply(c,B2P6,run_measures_on_hdf5(B2P6_files[i],dt),SIMPLIFY=FALSE)
}

for(i in 1:length(B2P7_files)){
	B2P7=mapply(c,B2P7,run_measures_on_hdf5(B2P7_files[i],dt),SIMPLIFY=FALSE)
}



#convert distance to factor for analysis

WT$dist=as.factor(WT$dist)
B2P4$dist=as.factor(B2P4$dist)
B2P5$dist=as.factor(B2P5$dist)
B2P6$dist=as.factor(B2P6$dist)
B2P7$dist=as.factor(B2P7$dist)

#get medians and quartiles

WT_quarts=list()
B2P4_quarts=list()
B2P5_quarts=list()
B2P6_quarts=list()
B2P7_quarts=list()

WT_quarts$med=quantiles_process(WT,0.5)
WT_quarts$low=quantiles_process(WT,0.25)
WT_quarts$up=quantiles_process(WT,0.75)

B2P4_quarts$med=quantiles_process(B2P4,0.5)
B2P4_quarts$low=quantiles_process(B2P4,0.25)
B2P4_quarts$up=quantiles_process(B2P4,0.75)

B2P5_quarts$med=quantiles_process(B2P5,0.5)
B2P5_quarts$low=quantiles_process(B2P5,0.25)
B2P5_quarts$up=quantiles_process(B2P5,0.75)

B2P6_quarts$med=quantiles_process(B2P6,0.5)
B2P6_quarts$low=quantiles_process(B2P6,0.25)
B2P6_quarts$up=quantiles_process(B2P6,0.75)

B2P7_quarts$med=quantiles_process(B2P7,0.5)
B2P7_quarts$low=quantiles_process(B2P7,0.25)
B2P7_quarts$up=quantiles_process(B2P7,0.75)


########################################
#get ranges in order to set graph limits

ci_range=range(c(WT_quarts$low[[1]][,2],WT_quarts$up[[1]][,2],B2P4_quarts$low[[1]][,2],B2P4_quarts$up[[1]][,2],B2P5_quarts$low[[1]][,2],B2P5_quarts$up[[1]][,2],B2P6_quarts$low[[1]][,2],B2P6_quarts$up[[1]][,2],B2P7_quarts$low[[1]][,2],B2P7_quarts$up[[1]][,2]),na.rm=TRUE)

sttc_range=range(c(WT_quarts$low[[2]][,2],WT_quarts$up[[2]][,2],B2P4_quarts$low[[2]][,2],B2P4_quarts$up[[2]][,2],B2P5_quarts$low[[2]][,2],B2P5_quarts$up[[2]][,2],B2P6_quarts$low[[2]][,2],B2P6_quarts$up[[2]][,2],B2P7_quarts$low[[2]][,2],B2P7_quarts$up[[2]][,2]),na.rm=TRUE)

max_x=max(as.vector(WT_quarts$med[[1]][,1]))
max_x=as.numeric(max_x)


#########################
#Read in spike trains for raster plots 

#WT
a1=h5.read.spikes("Kirkby2013_02_WT_P5.h5")


#B2P4
a2=h5.read.spikes("Kirkby2013_01_B2KO_P4.h5")


postscript("Figure_10.eps",width=inch(17.6),height=inch(17.6))

layout(matrix(c(1,1,2,2,3,4,3,4),4,2,byrow=TRUE))

par(mar=c(1.1,5.1,1.1,1.1))
 par(oma=c(1,1,1,0))



plot(a1,beg=300,end=900,which=c(20:30),main="",ylab="Wild Type P5",xlab="",xaxt="n",cex.lab=1.5)


segments(300,0.2,360,0.2,lwd=3,lend="butt")
text(x=340,y=0.11,"1 min",cex=1.5)

plot(a2,beg=600,end=1200,which=c(10:20),main="",ylab=expression(paste(beta,"2KO P4",sep="")),xlab="",xaxt="n",cex.lab=1.5,col.lab=rgb(8/255,104/255,172/255))



#plot measure vs distance graphs
par(mar=c(4.1,5.1,1.1,1.5))


plot(as.numeric(as.vector(WT_quarts$med[[1]][,1])),as.numeric(as.vector(WT_quarts$med[[1]][,2])),type="p",pch=20,col="black",ylim=c(0,100),main="",xlab=expression(paste("Intercell distance (",mu,"m)",sep="")),ylab="Correlation Index",xlim=c(-5,max_x+10),bty="n",cex=1,cex.lab=1.5,cex.axis=1.5,yaxt="n")
axis(2,las=2,at=c(0,25,50,75,100),cex.axis=1.5)

points(as.numeric(as.vector(WT_quarts$med[[1]][,1])),as.numeric(as.vector(WT_quarts$med[[1]][,2])),type="l",lwd=2)
arrows(as.numeric(as.vector(WT_quarts$low[[1]][1,1])),as.numeric(as.vector(WT_quarts$low[[1]][1,2])), as.numeric(as.vector(WT_quarts$up[[1]][1,1])),as.numeric(as.vector(WT_quarts$up[[1]][1,2])),angle=90,code=3,length=0.05,lwd=0.5)

points(as.numeric(as.vector(B2P4_quarts$med[[1]][,1])),as.numeric(as.vector(B2P4_quarts$med[[1]][,2])),col=rgb(8/255,104/255,172/255),pch=20,cex=1)
points(as.numeric(as.vector(B2P4_quarts$med[[1]][,1])),as.numeric(as.vector(B2P4_quarts$med[[1]][,2])),col=rgb(8/255,104/255,172/255),type="l",lwd=2)
arrows(as.numeric(as.vector(B2P4_quarts$low[[1]][1,1])),as.numeric(as.vector(B2P4_quarts$low[[1]][1,2])), as.numeric(as.vector(B2P4_quarts$up[[1]][1,1])),as.numeric(as.vector(B2P4_quarts$up[[1]][1,2])),angle=90,code=3,length=0.05,lwd=0.5,col=rgb(8/255,104/255,172/255))


points(as.numeric(as.vector(B2P5_quarts$med[[1]][,1])),as.numeric(as.vector(B2P5_quarts$med[[1]][,2])),col="cyan2",pch=20,cex=1)
points(as.numeric(as.vector(B2P5_quarts$med[[1]][,1])),as.numeric(as.vector(B2P5_quarts$med[[1]][,2])),col="cyan2",type="l",lwd=2)
arrows(as.numeric(as.vector(B2P5_quarts$low[[1]][1,1])),as.numeric(as.vector(B2P5_quarts$low[[1]][1,2])), as.numeric(as.vector(B2P5_quarts$up[[1]][1,1])),as.numeric(as.vector(B2P5_quarts$up[[1]][1,2])),angle=90,code=3,length=0.05,lwd=0.5,col="cyan2")


points(as.numeric(as.vector(B2P6_quarts$med[[1]][,1])),as.numeric(as.vector(B2P6_quarts$med[[1]][,2])),col=rgb(123/255,204/255,196/255),pch=20,cex=1)
points(as.numeric(as.vector(B2P6_quarts$med[[1]][,1])),as.numeric(as.vector(B2P6_quarts$med[[1]][,2])),col=rgb(123/255,204/255,196/255),type="l",lwd=2)
arrows(as.numeric(as.vector(B2P6_quarts$low[[1]][1,1])),as.numeric(as.vector(B2P6_quarts$low[[1]][1,2])), as.numeric(as.vector(B2P6_quarts$up[[1]][1,1])),as.numeric(as.vector(B2P6_quarts$up[[1]][1,2])),angle=90,code=3,length=0.05,lwd=0.5,col=rgb(123/255,204/255,196/255))


points(as.numeric(as.vector(B2P7_quarts$med[[1]][,1])),as.numeric(as.vector(B2P7_quarts$med[[1]][,2])),col=rgb(186/255,228/255,188/255),pch=20,cex=1)
points(as.numeric(as.vector(B2P7_quarts$med[[1]][,1])),as.numeric(as.vector(B2P7_quarts$med[[1]][,2])),col=rgb(186/255,228/255,188/255),type="l",lwd=2)
arrows(as.numeric(as.vector(B2P7_quarts$low[[1]][1,1])),as.numeric(as.vector(B2P7_quarts$low[[1]][1,2])), as.numeric(as.vector(B2P7_quarts$up[[1]][1,1])),as.numeric(as.vector(B2P7_quarts$up[[1]][1,2])),angle=90,code=3,length=0.05,lwd=0.5,rgb(186/255,228/255,188/255))


#second measure vs distance plot

plot(as.vector(WT_quarts$med[[2]][,1]),as.vector(WT_quarts$med[[2]][,2]),type="p",pch=20,col="black",main="",xlab=expression(paste("Intercell distance (",mu,"m)",sep="")),ylab="STTC",xlim=c(-5,max_x+10),bty="n",cex=1.25,cex.lab=1.5,cex.axis=1.5,ylim=sttc_range,yaxt="n")
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8),cex.axis=1.5)

points(as.numeric(as.vector(WT_quarts$med[[2]][,1])),as.numeric(as.vector(WT_quarts$med[[2]][,2])),type="l",lwd=2)
arrows(as.numeric(as.vector(WT_quarts$low[[2]][1,1])),as.numeric(as.vector(WT_quarts$low[[2]][1,2])), as.numeric(as.vector(WT_quarts$up[[2]][1,1])),as.numeric(as.vector(WT_quarts$up[[2]][1,2])),angle=90,code=3,length=0.05,lwd=0.5)


points(as.numeric(as.vector(B2P4_quarts$med[[2]][,1])),as.numeric(as.vector(B2P4_quarts$med[[2]][,2])),col=rgb(8/255,104/255,172/255),pch=20,cex=1)
points(as.numeric(as.vector(B2P4_quarts$med[[2]][,1])),as.numeric(as.vector(B2P4_quarts$med[[2]][,2])),col=rgb(8/255,104/255,172/255),type="l",lwd=2)
arrows(as.numeric(as.vector(B2P4_quarts$low[[2]][1,1])),as.numeric(as.vector(B2P4_quarts$low[[2]][1,2])), as.numeric(as.vector(B2P4_quarts$up[[2]][1,1])),as.numeric(as.vector(B2P4_quarts$up[[2]][1,2])),angle=90,code=3,length=0.05,lwd=0.5,col=rgb(8/255,104/255,172/255))

points(as.numeric(as.vector(B2P5_quarts$med[[2]][,1])),as.numeric(as.vector(B2P5_quarts$med[[2]][,2])),col="cyan2",pch=20,cex=1)
points(as.numeric(as.vector(B2P5_quarts$med[[2]][,1])),as.numeric(as.vector(B2P5_quarts$med[[2]][,2])),col="cyan2",type="l",lwd=2)
arrows(as.numeric(as.vector(B2P5_quarts$low[[2]][1,1])),as.numeric(as.vector(B2P5_quarts$low[[2]][1,2])), as.numeric(as.vector(B2P5_quarts$up[[2]][1,1])),as.numeric(as.vector(B2P5_quarts$up[[2]][1,2])),angle=90,code=3,length=0.05,lwd=0.5,col="cyan2")


points(as.numeric(as.vector(B2P6_quarts$med[[2]][,1])),as.numeric(as.vector(B2P6_quarts$med[[2]][,2])),col=rgb(123/255,204/255,196/255),pch=20,cex=1)
points(as.numeric(as.vector(B2P6_quarts$med[[2]][,1])),as.numeric(as.vector(B2P6_quarts$med[[2]][,2])),col=rgb(123/255,204/255,196/255),type="l",lwd=2)
arrows(as.numeric(as.vector(B2P6_quarts$low[[2]][1,1])),as.numeric(as.vector(B2P6_quarts$low[[2]][1,2])), as.numeric(as.vector(B2P6_quarts$up[[2]][1,1])),as.numeric(as.vector(B2P6_quarts$up[[2]][1,2])),angle=90,code=3,length=0.05,lwd=0.5,col=rgb(123/255,204/255,196/255))


points(as.numeric(as.vector(B2P7_quarts$med[[2]][,1])),as.numeric(as.vector(B2P7_quarts$med[[2]][,2])),col=rgb(186/255,228/255,188/255),pch=20,cex=1)
points(as.numeric(as.vector(B2P7_quarts$med[[2]][,1])),as.numeric(as.vector(B2P7_quarts$med[[2]][,2])),col=rgb(186/255,228/255,188/255),type="l",lwd=2)
arrows(as.numeric(as.vector(B2P7_quarts$low[[2]][1,1])),as.numeric(as.vector(B2P7_quarts$low[[2]][1,2])), as.numeric(as.vector(B2P7_quarts$up[[2]][1,1])),as.numeric(as.vector(B2P7_quarts$up[[2]][1,2])),angle=90,code=3,length=0.05,lwd=0.5,col=rgb(186/255,228/255,188/255))


legend(x=110,y=0.8,col=c("black",rgb(8/255,104/255,172/255),"cyan2",rgb(123/255,204/255,196/255),rgb(186/255,228/255,188/255)),pch=c(20),lty=1,lwd=1.5,pt.cex=c(1),cex=c(1.25),legend=c("WT    P5 (0.51 Hz, n=6)",expression(paste(beta,"2KO P4 (0.03 Hz, n=2)",sep="")),expression(paste(beta,"2KO P5 (0.01 Hz, n=1)",sep="")),expression(paste(beta,"2KO P6 (0.37 Hz, n=5)",sep="")),expression(paste(beta,"2KO P7 (0.60 Hz, n=3)",sep=""))))


#label plots
text(grconvertX(c(0.03, 0.03,0.5), from='ndc'),
   grconvertY(c(0.96, 0.450,0.45), from='ndc'),  c('A',  'B',  'C'), xpd=NA, cex=2.2, font=2)



dev.off()


