#include <R.h>
#include <time.h>
#include <stdio.h>
#include <Rmath.h>

//calculates the correlation index

void run_ci(int *N1v,int *N2v,double *dtv,double *Time,double *index,double *spike_times_1,double *spike_times_2){

	
	int i;
	int j;
	int u;
	double dt= *dtv;
	int N1= *N1v;
	int N2= *N2v;
	double T;
	int Nab;
	T=Time[1]-Time[0];	
	
	Nab=0;	
	j=0;
	for(i=0;i<N1;i++){
		
		while(j<N2){
			
			if((spike_times_1[i]-spike_times_2[j])>dt){
				j=j+1;
			
			}
			else if(fabs(spike_times_1[i]-spike_times_2[j])<=dt){
				Nab=Nab+1;
				u=j+1;
				while(fabs(spike_times_1[i]-spike_times_2[u])<=dt){
				
					Nab=Nab+1;
					u=u+1;
				}
				break;
			}
			else{
				
				break;
				
			}
		}
	}
	
	*index=(((double)Nab)*T)/(((double)N1)*((double)N2)*2*dt);
	
		
		
}
