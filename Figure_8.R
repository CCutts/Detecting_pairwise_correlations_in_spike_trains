#Data from Wong et al 1993 and Demas et al 2003 assumed to be in the same directory
library(fields)

dt=0.05

source("processing_functions.R")


D_P9_files=c("Demas2003P9_CTRL_MY1_1A","Demas2003P9_CTRL_MY1_2A")

D_P15_files=c("Demas2003P15_CTRL_MW2_2A","Demas2003P15_CTRL_MS5_1A","Demas2003P15_CTRL_MI1_2B","Demas2003P15_CTRL_MEE1_1A")

D_P21_files=c("Demas2003p21_ctrl_mx1","Demas2003p21_ctrl_2a_short_new","Demas2003p21_ctrl_1a_short_new","Demas2003p22_ctrl_1a")

D_6w_files=c("Demas20036WK_CTRL_ME3_1B","Demas20036WK_CTRL_ME2","Demas20036WK_CTRL_ME1","Demas20036WK_CTRL_MD3_1A")

##########
#initiate data frames

W_P0=run_measures_on_hdf5_for_binning("Wong1993P0",dt)

W_P15=run_measures_on_hdf5_for_binning("Wong1993P15b",dt)

W_P21=run_measures_on_hdf5_for_binning("Wong1993P21b",dt)

W_P30=run_measures_on_hdf5_for_binning("Wong1993P30a",dt)


D_P9=data.frame(ci=integer(0),sttc=integer(0),dist=integer(0))

D_P15=data.frame(ci=integer(0),sttc=integer(0),dist=integer(0))

D_P21=data.frame(ci=integer(0),sttc=integer(0),dist=integer(0))

D_6w=data.frame(ci=integer(0),sttc=integer(0),dist=integer(0))




for(i in 1:length(D_P9_files)){
	D_P9=mapply(c,D_P9,run_measures_on_hdf5(D_P9_files[i],dt),SIMPLIFY=FALSE)
}


for(i in 1:length(D_P15_files)){
	D_P15=mapply(c,D_P15,run_measures_on_hdf5(D_P15_files[i],dt),SIMPLIFY=FALSE)
}



for(i in 1:length(D_P21_files)){
	D_P21=mapply(c,D_P21,run_measures_on_hdf5(D_P21_files[i],dt),SIMPLIFY=FALSE)
}


for(i in 1:length(D_6w_files)){
	D_6w=mapply(c,D_6w,run_measures_on_hdf5(D_6w_files[i],dt),SIMPLIFY=FALSE)
}




##########################
#Wong data needs to be binned 

bins=seq(from=0,to=600,by=20)



W_P0_b=bin_distances(W_P0,bins)


W_P0_quarts=list()
W_P0_quarts$med=quantiles_process(W_P0_b,0.5)
W_P0_quarts$low=quantiles_process(W_P0_b,0.25)
W_P0_quarts$up=quantiles_process(W_P0_b,0.75)



W_P15_b=bin_distances(W_P15,bins)


W_P15_quarts=list()
W_P15_quarts$med=quantiles_process(W_P15_b,0.5)
W_P15_quarts$low=quantiles_process(W_P15_b,0.25)
W_P15_quarts$up=quantiles_process(W_P15_b,0.75)

W_P21_b=bin_distances(W_P21,bins)


W_P21_quarts=list()
W_P21_quarts$med=quantiles_process(W_P21_b,0.5)
W_P21_quarts$low=quantiles_process(W_P21_b,0.25)
W_P21_quarts$up=quantiles_process(W_P21_b,0.75)

W_P30_b=bin_distances(W_P30,bins)


W_P30_quarts=list()
W_P30_quarts$med=quantiles_process(W_P30_b,0.5)
W_P30_quarts$low=quantiles_process(W_P30_b,0.25)
W_P30_quarts$up=quantiles_process(W_P30_b,0.75)


W_ci_range=range(c(W_P0_quarts$low[[1]][,2],W_P0_quarts$up[[1]][,2],W_P15_quarts$low[[1]][,2],W_P15_quarts$up[[1]][,2],W_P21_quarts$low[[1]][,2],W_P21_quarts$up[[1]][,2],W_P30_quarts$low[[1]][,2],W_P30_quarts$up[[1]][,2]),na.rm=TRUE)

W_sttc_range=range(c(W_P0_quarts$low[[2]][,2],W_P0_quarts$up[[2]][,2],W_P15_quarts$low[[2]][,2],W_P15_quarts$up[[2]][,2],W_P21_quarts$low[[2]][,2],W_P21_quarts$up[[2]][,2],W_P30_quarts$low[[2]][,2],W_P30_quarts$up[[2]][,2]),na.rm=TRUE)

W_max_x=max(c(as.vector(W_P0_quarts$med[[1]][,1]),as.vector(W_P21_quarts$med[[1]][,1]),as.vector(W_P15_quarts$med[[1]][,1]),as.vector(W_P30_quarts$med[[1]][,1])))
W_max_x=as.numeric(W_max_x)

#############################
#Demas data is not binned


#convert distance to factor for analysis

D_P9$dist=as.factor(D_P9$dist)
D_P15$dist=as.factor(D_P15$dist)
D_P21$dist=as.factor(D_P21$dist)
D_6w$dist=as.factor(D_6w$dist)



#get medians and quartiles

D_P9_quarts=list()
D_P9_quarts$med=quantiles_process(D_P9,0.5)
D_P9_quarts$low=quantiles_process(D_P9,0.25)
D_P9_quarts$up=quantiles_process(D_P9,0.75)

D_P15_quarts=list()
D_P15_quarts$med=quantiles_process(D_P15,0.5)
D_P15_quarts$low=quantiles_process(D_P15,0.25)
D_P15_quarts$up=quantiles_process(D_P15,0.75)

D_P21_quarts=list()
D_P21_quarts$med=quantiles_process(D_P21,0.5)
D_P21_quarts$low=quantiles_process(D_P21,0.25)
D_P21_quarts$up=quantiles_process(D_P21,0.75)

D_6w_quarts=list()
D_6w_quarts$med=quantiles_process(D_6w,0.5)
D_6w_quarts$low=quantiles_process(D_6w,0.25)
D_6w_quarts$up=quantiles_process(D_6w,0.75)

#get ranges for graphs
D_ci_range=range(c(D_P9_quarts$low[[1]][,2],D_P9_quarts$up[[1]][,2],D_P15_quarts$low[[1]][,2],D_P15_quarts$up[[1]][,2],D_P21_quarts$low[[1]][,2],D_P21_quarts$up[[1]][,2],D_6w_quarts$low[[1]][,2],D_6w_quarts$up[[1]][,2]),na.rm=TRUE)

D_sttc_range=range(c(D_P9_quarts$low[[2]][,2],D_P9_quarts$up[[2]][,2],D_P15_quarts$low[[2]][,2],D_P15_quarts$up[[2]][,2],D_P21_quarts$low[[2]][,2],D_P21_quarts$up[[2]][,2],D_6w_quarts$low[[2]][,2],D_6w_quarts$up[[2]][,2]),na.rm=TRUE)

D_max_x=max(c(as.vector(D_P9_quarts$med[[1]][,1]),as.vector(D_P15_quarts$med[[1]][,1]),as.vector(D_P21_quarts$med[[1]][,1]),as.vector(D_6w_quarts$med[[1]][,1])))
D_max_x=as.numeric(D_max_x)


postscript("Figure_8.eps",height=inch(11.6),width=inch(11.6))
par(mfrow=c(2,2))
par(mar=c(5,4,1,1))
par(oma=c(0,0,1.5,0))

#top left
plot(as.numeric(as.vector(W_P0_quarts$med[[1]][-1,1])),as.numeric(as.vector(W_P0_quarts$med[[1]][-1,2])),ylim=c(0,80),pch=19,col="black",ylab="Correlation Index",xlab=expression(paste("Intercell distance (",mu,"m)",sep="")),main="",bty="n",xlim=c(0,600),yaxt="n",xaxt="n",cex.lab=0.75,cex=0.75)

axis(2,las=2,cex.axis=0.75)
axis(1,at=c(0,150,300,450,600),labels=c(0,150,300,450,600),cex.axis=0.75)

legend(x=200,y=80,pch=c(19,1,15,17),col=c("black","blue","red","purple"),legend=c("P0   (0.74 Hz, n=1)","P15 (0.29 Hz, n=1)","P21 (0.20 Hz, n=1)","P30 (0.15 Hz, n=1)"),cex=0.6)

arrows(as.numeric(as.vector(W_P0_quarts$low[[1]][2,1])),as.numeric(as.vector(W_P0_quarts$low[[1]][2,2])),as.numeric(as.vector(W_P0_quarts$up[[1]][2,1])),as.numeric(as.vector(W_P0_quarts$up[[1]][2,2])), angle=90, code=3, length=0.03,lwd=0.5)

points(as.numeric(as.vector(W_P15_quarts$med[[1]][-1,1])),as.numeric(as.vector(W_P15_quarts$med[[1]][-1,2])),pch=1,col="blue",cex=0.75)
arrows(as.numeric(as.vector(W_P15_quarts$low[[1]][2,1])),as.numeric(as.vector(W_P15_quarts$low[[1]][2,2])),as.numeric(as.vector(W_P15_quarts$up[[1]][2,1])),as.numeric(as.vector(W_P15_quarts$up[[1]][2,2])), angle=90, code=3, length=0.03,lwd=0.5,col="blue")

points(as.numeric(as.vector(W_P21_quarts$med[[1]][-1,1])),as.numeric(as.vector(W_P21_quarts$med[[1]][-1,2])),pch=15,col="red",cex=0.75)
arrows(as.numeric(as.vector(W_P21_quarts$low[[1]][2,1])),as.numeric(as.vector(W_P21_quarts$low[[1]][2,2])),as.numeric(as.vector(W_P21_quarts$up[[1]][2,1])),as.numeric(as.vector(W_P21_quarts$up[[1]][2,2])), angle=90, code=3, length=0.03,lwd=0.5,col="red")

points(as.numeric(as.vector(W_P30_quarts$med[[1]][-1,1])),as.numeric(as.vector(W_P30_quarts$med[[1]][-1,2])),pch=17,col="purple",cex=0.75)
arrows(as.numeric(as.vector(W_P30_quarts$low[[1]][2,1])),as.numeric(as.vector(W_P30_quarts$low[[1]][2,2])),as.numeric(as.vector(W_P30_quarts$up[[1]][2,1])),as.numeric(as.vector(W_P30_quarts$up[[1]][2,2])), angle=90, code=3, length=0.03,lwd=0.5,col="purple")

#top right
plot(as.numeric(as.vector(W_P0_quarts$med[[2]][-1,1])),as.numeric(as.vector(W_P0_quarts$med[[2]][-1,2])),ylim=c(-0.01,0.8),pch=19,col="black",ylab="STTC",xlab=expression(paste("Intercell distance (",mu,"m)",sep="")),main="",bty="n",xlim=c(0,600),yaxt="n",xaxt="n",cex.lab=0.75,cex=0.75)

axis(2,las=2,cex.axis=0.75)
axis(1,at=c(0,150,300,450,600),labels=c(0,150,300,450,600),cex.axis=0.75)

arrows(as.numeric(as.vector(W_P0_quarts$low[[2]][2,1])),as.numeric(as.vector(W_P0_quarts$low[[2]][2,2])),as.numeric(as.vector(W_P0_quarts$up[[2]][2,1])),as.numeric(as.vector(W_P0_quarts$up[[2]][2,2])), angle=90, code=3, length=0.03,lwd=0.5)

points(as.numeric(as.vector(W_P15_quarts$med[[2]][-1,1])),as.numeric(as.vector(W_P15_quarts$med[[2]][-1,2])),pch=1,col="blue",cex=0.75)
arrows(as.numeric(as.vector(W_P15_quarts$low[[2]][2,1])),as.numeric(as.vector(W_P15_quarts$low[[2]][2,2])),as.numeric(as.vector(W_P15_quarts$up[[2]][2,1])),as.numeric(as.vector(W_P15_quarts$up[[2]][2,2])), angle=90, code=3, length=0.03,lwd=0.5,col="blue")

points(as.numeric(as.vector(W_P21_quarts$med[[2]][-1,1])),as.numeric(as.vector(W_P21_quarts$med[[2]][-1,2])),pch=15,col="red",cex=0.75)
arrows(as.numeric(as.vector(W_P21_quarts$low[[2]][2,1])),as.numeric(as.vector(W_P21_quarts$low[[2]][2,2])),as.numeric(as.vector(W_P21_quarts$up[[2]][2,1])),as.numeric(as.vector(W_P21_quarts$up[[2]][2,2])), angle=90, code=3, length=0.03,lwd=0.5,col="red")

points(as.numeric(as.vector(W_P30_quarts$med[[2]][-1,1])),as.numeric(as.vector(W_P30_quarts$med[[2]][-1,2])),pch=17,col="purple",cex=0.75)
arrows(as.numeric(as.vector(W_P30_quarts$low[[2]][2,1])),as.numeric(as.vector(W_P30_quarts$low[[2]][2,2])),as.numeric(as.vector(W_P30_quarts$up[[2]][2,1])),as.numeric(as.vector(W_P30_quarts$up[[2]][2,2])), angle=90, code=3, length=0.03,lwd=0.5,col="purple")

#bottom left
plot(as.numeric(as.vector(D_P9_quarts$med[[1]][,1])),as.numeric(as.vector(D_P9_quarts$med[[1]][,2])),type="p",pch=19,col="black",ylim=c(0,100),main="",xlab=expression(paste("Intercell distance (",mu,"m)",sep="")),ylab="Correlation Index",xlim=c(-5,850),bty="n",yaxt="n",cex.axis=0.75,cex.lab=0.75,cex=0.75)
axis(2,las=2,at=c(0,25,50,75,100),labels=c(0,25,50,75,100),cex.axis=0.75)
legend(x=240,y=100,pch=c(19,1,15,17),col=c("black","blue","red","purple"),legend=c("P9        (0.30 Hz, n=2)","P15      (0.96 Hz, n=4)","P21      (0.70 Hz, n=3)","6 week (2.40 Hz, n=4)"),cex=0.6)
arrows(as.numeric(as.vector(D_P9_quarts$low[[1]][1,1])),as.numeric(as.vector(D_P9_quarts$low[[1]][1,2])),as.numeric(as.vector(D_P9_quarts$up[[1]][1,1])),as.numeric(as.vector(D_P9_quarts$up[[1]][1,2])), angle=90, code=3, length=0.03,lwd=0.5)
points(as.numeric(as.vector(D_P15_quarts$med[[1]][,1])),as.numeric(as.vector(D_P15_quarts$med[[1]][,2])),pch=1,col="blue",cex=0.75)
arrows(as.numeric(as.vector(D_P15_quarts$low[[1]][1,1])),as.numeric(as.vector(D_P15_quarts$low[[1]][1,2])),as.numeric(as.vector(D_P15_quarts$up[[1]][1,1])),as.numeric(as.vector(D_P15_quarts$up[[1]][1,2])), angle=90, code=3, length=0.03,lwd=0.5,col="blue")
points(as.numeric(as.vector(D_P21_quarts$med[[1]][,1])),as.numeric(as.vector(D_P21_quarts$med[[1]][,2])),pch=15,col="red",cex=0.75)
arrows(as.numeric(as.vector(D_P21_quarts$low[[1]][1,1])),as.numeric(as.vector(D_P21_quarts$low[[1]][1,2])),as.numeric(as.vector(D_P21_quarts$up[[1]][1,1])),as.numeric(as.vector(D_P21_quarts$up[[1]][1,2])), angle=90, code=3, length=0.03,lwd=0.5,col="red")
points(as.numeric(as.vector(D_6w_quarts$med[[1]][,1])),as.numeric(as.vector(D_6w_quarts$med[[1]][,2])),pch=17,col="purple",cex=0.75)
arrows(as.numeric(as.vector(D_6w_quarts$low[[1]][1,1])),as.numeric(as.vector(D_6w_quarts$low[[1]][1,2])),as.numeric(as.vector(D_6w_quarts$up[[1]][1,1])),as.numeric(as.vector(D_6w_quarts$up[[1]][1,2])), angle=90, code=3, length=0.03,lwd=0.5,col="purple")

#bottom right
plot(as.numeric(as.vector(D_P9_quarts$med[[2]][,1])),as.numeric(as.vector(D_P9_quarts$med[[2]][,2])),type="p",pch=19,col="black",ylim=c(-0.01,0.8),main="",xlab=expression(paste("Intercell distance (",mu,"m)",sep="")),ylab="Correlation Index",xlim=c(-5,850),bty="n",yaxt="n",cex.axis=0.75,cex.lab=0.75,cex=0.75)
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8),labels=c(0,25,50,75,100),cex.axis=0.75)
legend(x=240,y=100,pch=c(19,1,15,17),col=c("black","blue","red","purple"),legend=c("P9        (0.30 Hz, n=2)","P15      (0.96 Hz, n=4)","P21      (0.70 Hz, n=3)","6 week (2.40 Hz, n=4)"),cex=0.6)
arrows(as.numeric(as.vector(D_P9_quarts$low[[2]][1,1])),as.numeric(as.vector(D_P9_quarts$low[[2]][1,2])),as.numeric(as.vector(D_P9_quarts$up[[2]][1,1])),as.numeric(as.vector(D_P9_quarts$up[[2]][1,2])), angle=90, code=3, length=0.03,lwd=0.5)
points(as.numeric(as.vector(D_P15_quarts$med[[2]][,1])),as.numeric(as.vector(D_P15_quarts$med[[2]][,2])),pch=1,col="blue",cex=0.75)
arrows(as.numeric(as.vector(D_P15_quarts$low[[2]][1,1])),as.numeric(as.vector(D_P15_quarts$low[[2]][1,2])),as.numeric(as.vector(D_P15_quarts$up[[2]][1,1])),as.numeric(as.vector(D_P15_quarts$up[[2]][1,2])), angle=90, code=3, length=0.03,lwd=0.5,col="blue")
points(as.numeric(as.vector(D_P21_quarts$med[[2]][,1])),as.numeric(as.vector(D_P21_quarts$med[[2]][,2])),pch=15,col="red",cex=0.75)
arrows(as.numeric(as.vector(D_P21_quarts$low[[2]][1,1])),as.numeric(as.vector(D_P21_quarts$low[[2]][1,2])),as.numeric(as.vector(D_P21_quarts$up[[2]][1,1])),as.numeric(as.vector(D_P21_quarts$up[[2]][1,2])), angle=90, code=3, length=0.03,lwd=0.5,col="red")
points(as.numeric(as.vector(D_6w_quarts$med[[2]][,1])),as.numeric(as.vector(D_6w_quarts$med[[2]][,2])),pch=17,col="purple",cex=0.75)
arrows(as.numeric(as.vector(D_6w_quarts$low[[2]][1,1])),as.numeric(as.vector(D_6w_quarts$low[[2]][1,2])),as.numeric(as.vector(D_6w_quarts$up[[2]][1,1])),as.numeric(as.vector(D_6w_quarts$up[[2]][1,2])), angle=90, code=3, length=0.03,lwd=0.5,col="purple")

text(grconvertX(c(0.02, 0.52, 0.02,0.52), from='ndc'),
   grconvertY(c(0.9, 0.9, 0.42,0.42), from='ndc'),  c('A',  'B',  'C','D'), xpd=NA, cex=1.5, font=2)

text(grconvertX(c(0.5, 0.5), from='ndc'),
   grconvertY(c(0.98, 0.49), from='ndc'),  c('Ferret',  'Mouse'), xpd=NA, cex=1.5, font=2)



dev.off()





