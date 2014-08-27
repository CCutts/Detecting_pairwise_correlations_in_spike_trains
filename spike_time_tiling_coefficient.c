#include <R.h>
#include <stdio.h>
#include <Rmath.h>

//calculates the tiling coefficient

double run_P(int N1,int N2,double dt,double *spike_times_1,double *spike_times_2){

	int i;
	int j;
	int Nab;
		
	Nab=0;
	j=0;
	for(i=0;i<=(N1-1);i++){
		while(j<N2){	
//check every spike in train 1 to see if there's a spike in train 2 within dt  (don't count spike pairs)
// don't need to search all j each iteration
			if(fabs(spike_times_1[i]-spike_times_2[j])<=dt){
				Nab=Nab+1;	
				break;				
			}
			else if(spike_times_2[j]>spike_times_1[i]){			
				break;
			}
			else{
				j=j+1;
			}		
		}
	}
	return Nab;
}



double run_T(int N1v,double dtv,double startv, double endv, double *spike_times_1){
	
	double dt= dtv;
	double start=startv;
	double end=endv;
	int N1= N1v;
	double time_A;
	int i=0;
	double diff;
	
//maximum 
	time_A=2*(double)N1*dt;

// if just one spike in train 
	if(N1==1){
		
	  if((spike_times_1[0]-start)<dt){
	    	time_A=time_A-start+spike_times_1[0]-dt;
	  }
	  else if((spike_times_1[0]+dt)>end){
	   	 time_A=time_A-spike_times_1[0]-dt+end;
	      }
	  
	}
	
//if more than one spike in train
	else{
			

			while(i<(N1-1)){
			
				diff=spike_times_1[i+1]-spike_times_1[i];
				
				if(diff<2*dt){
					//subtract overlap 	
					time_A=time_A-2*dt+diff;

				}
				
				i++;
			}
				
			//check if spikes are within dt of the start and/or end, if so just need to subract
			//overlap of first and/or last spike as all within-train overlaps have been accounted for
			
			
			if((spike_times_1[0]-start)<dt){
				
			  time_A=time_A-start+spike_times_1[0]-dt;
			}


			if((end-spike_times_1[N1-1])<dt){
				
			  time_A=time_A-spike_times_1[N1-1]-dt+end;
			}
              	}
	
	return time_A;	
}
	


void run_sttc(int *N1v,int *N2v,double *dtv,double *Time,double *index,double *spike_times_1,double *spike_times_2){

	double TA;
	double TB;
	double PA;
	double PB;
	int N1= *N1v;
	int N2= *N2v;
	double dt= *dtv;
	double T;

	
	
	if(N1==0 || N2==0){
	*index=R_NaN;
	}
	else{
	T=Time[1]-Time[0];
	TA=run_T(N1,dt,Time[0],Time[1], spike_times_1);
	TA=TA/T;
	TB=run_T(N2,dt,Time[0],Time[1], spike_times_2);
	TB=TB/T;
	PA=run_P(N1,N2,dt, spike_times_1, spike_times_2);
	PA=PA/(double)N1;
	PB=run_P(N2,N1,dt, spike_times_2, spike_times_1);
	PB=PB/(double)N2;
	*index=0.5*(PA-TB)/(1-TB*PA)+0.5*(PB-TA)/(1-TA*PB);
	

	}

}



