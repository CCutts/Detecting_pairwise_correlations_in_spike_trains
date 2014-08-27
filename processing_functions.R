#This file contains the processing functions which are common to all Figures

require(sjemea)
require(rhdf5)
library(gdata)
dyn.load("correlation_index.so")
dyn.load("spike_time_tiling_coefficient.so")


#This reads in the hdf5 file and calculates the ci and sttc for each pair of electrodes 
run_measures_on_hdf5=function(filename,dt){

	a=h5.read.spikes(paste(filename,".h5",sep=""))

#this stores both the ci and the sttc for each pair of cells
	indices=array(0,dim=c(a$NCells[[1]],a$NCells[[1]],2))
	
#this finds the ci and the sttc for each pair of cells 
	for(i in 1:a$NCells[[1]]){
		for(j in 1:a$NCells[[1]]){

			b=.C("run_ci",as.integer(length(a$spikes[[i]])),as.integer(length(a$spikes[[j]])),as.double(dt),as.double(a$rec.time),index=as.double(1),as.double(as.vector(a$spikes[[i]])),as.double(as.vector(a$spikes[[j]])))

			indices[i,j,1]=b[[5]]

			c=.C("run_sttc",as.integer(length(a$spikes[[i]])),as.integer(length(a$spikes[[j]])),as.double(dt),as.double(a$rec.time),index=as.double(1),as.double(as.vector(a$spikes[[i]])),as.double(as.vector(a$spikes[[j]])))

			indices[i,j,2]=c[[5]]

		}
	}

	dist=as.matrix(dist(a$layout$pos,upper=TRUE,diag=TRUE))


	#extract lower triangles -matrix symmetric-repeated vals
	#diagonals not plotted- v large for ci
	ci=lowerTriangle(indices[,,1])
	sttc=lowerTriangle(indices[,,2])
	dist=lowerTriangle(dist)
	
	results=list(ci=ci,sttc=sttc,dist=dist)

	return(results)

}


run_measures_on_hdf5_for_binning=function(filename,dt){

	a=h5.read.spikes(paste(filename,".h5",sep=""))

#this stores both the ci and the sttc for each pair of cells
	indices=array(0,dim=c(a$NCells[[1]],a$NCells[[1]],2))
	
#this finds the ci and the sttc for each pair of cells 
	for(i in 1:a$NCells[[1]]){
		for(j in 1:a$NCells[[1]]){

			b=.C("run_ci",as.integer(length(a$spikes[[i]])),as.integer(length(a$spikes[[j]])),as.double(dt),as.double(a$rec.time),index=as.double(1),as.double(as.vector(a$spikes[[i]])),as.double(as.vector(a$spikes[[j]])))

			indices[i,j,1]=b[[5]]

			c=.C("run_sttc",as.integer(length(a$spikes[[i]])),as.integer(length(a$spikes[[j]])),as.double(dt),as.double(a$rec.time),index=as.double(1),as.double(as.vector(a$spikes[[i]])),as.double(as.vector(a$spikes[[j]])))

			indices[i,j,2]=c[[5]]

		}
	}

	dist=as.matrix(dist(a$layout$pos,upper=TRUE,diag=TRUE))


	#extract lower triangles -matrix symmetric-repeated vals
	#diagonals not plotted- v large for ci
	ci=lowerTriangle(indices[,,1])
	sttc=lowerTriangle(indices[,,2])
	dist=lowerTriangle(dist)
	
	results=data.frame(dist=dist,ci=ci,sttc=sttc)

	return(results)

}



quantiles_process=function(points,quant){

	a=aggregate(ci~dist,points,function (x) quantile(x, quant))
	b=aggregate(sttc~dist,points,function (x) quantile(x, quant))
#note, quantile(1,0.25)=1, quantile(c(1,2),0.25)=1.25, quantile(1,0.75)=1 etc so don't need to worry about single values

#ensure that list is in increasing distance order 
	a=a[order(as.vector(a[,1])) ,]
	b=b[order(as.vector(b[,1])) ,]

	vec=list(a,b)

} 


#converts cm to inches (R plot sizing in inches)
inch=function(x){
    x/2.54
}


bin_distances=function(data,bins){
l=length(bins)
	mids=rep(0,l-1)
	for(i in 1:(l-1)){
	mids[i]=(bins[i]+bins[i+1])/2
	}

for(i in 1:length(data$dist)){
f=abs(data$dist[i]-mids)
q=which(f==min(f))
if(length(q)>1){
q=q[1]
}
data$dist[i]=mids[q]
}
return(data)
}







