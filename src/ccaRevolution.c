#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// m is (n x 3)-matrix. The first col includes x-coordinates, the second y-coordinates
// and the third one includes place -holder for the cluster-id
// n is the length of the list (the number of rows of the matrix, to be precise)
// l ist the cluster distance

void ccaRev(double *m, int *n, double *l, int *w)
{	
	w[0] = 0;
	int zeros, i, k, j, mm, cc;
	zeros = 0;
	i=0;
	mm=1;
	cc=1;
	double dist;
	
	while(zeros<*n){
		k = w[i];
		if(m[2* *n+ k]==0){
			m[2* *n + k] = cc;
			zeros++;
		}
		if(k>0)
			for(j=k-1; j>=0 && (m[k]-m[j])<=*l; j--){
				if(m[2* *n + j]==0 && abs(m[*n+k]-m[*n+j])<=*l){
					dist = sqrt(((m[k]-m[j])*(m[k]-m[j]))+((m[*n+k]-m[*n+j])*(m[*n+k]-m[*n+j])));
					if(dist<=(double)*l){
						w[mm] = j;
						m[2* *n+j] = cc;
						mm++;
						zeros++;
					}
						
				}	
			}
		if(k<*n-1)
			for(j=k+1; j<*n && (m[j]-m[k])<=*l; j++){
				if(m[2* *n + j]==0 && abs(m[*n+k]-m[*n+j])<=*l){
					dist = sqrt(((m[k]-m[j])*(m[k]-m[j]))+((m[*n+k]-m[*n+j])*(m[*n+k]-m[*n+j])));
					if(dist<=(double)*l){
						w[mm] = j;
						m[2* *n+j] = cc;
						mm++;
						zeros++;
					}	
				}	
			}
		i++;
		if(zeros==*n){ break;}
		if(w[i]==0){
			cc++;
			k = 0;
			while(m[2* *n+ k]!=0) k++;
			w[i] = k;
			mm = i+1;
		}
	}
}



void ccaSum(double *m, int *m3, double *mm, int *n){
	int i;
	for(i=0; i<*n; i++)
		mm[m3[i]-1] += 1;
}

void ccaSumT(double *m, int *m3, double *mm, int *n){
	int i;
	for(i=0; i<*n; i++)
		mm[m3[i]-1] += cos(m[*n+i]);
}

void ccaRevT(double *m, int *n, double *l, int *step_w, int *step_h, double *res_x, double *res_y, int *w)
{	
	w[0] = 0;
	int zeros, i, k, j, mm, cc;
	zeros = 0;
	i=0;
	mm=1;
	cc=1;
	double dist;
	
	while(zeros<*n){
		k = w[i];
		if(m[2* *n+ k]==0){
			m[2* *n + k] = cc;
			zeros++;
		}
		if(k>0)
			for(j=k-1; j>=0 && abs((m[*n+k]-m[*n+j])/ *res_y)<=*step_h; j--){
				if(m[2* *n + j]==0 && abs((m[k]-m[j])/ *res_x)<=*step_w){
					dist = acos(sin(m[*n+k])*sin(m[*n+j])+cos(m[*n+k])*cos(m[*n+j])*cos(m[k]-m[j]))*6371000;
					if(dist<=(double)*l){
						w[mm] = j;
						m[2* *n+j] = cc;
						mm++;
						zeros++;
					}
						
				}	
			}
		if(k<*n-1)
			for(j=k+1; j<*n && abs((m[*n+j]-m[*n+k])/ *res_y)<=*step_h; j++){
				if(m[2* *n + j]==0 && abs((m[k]-m[j])/ *res_x)<=*step_w){ 
					dist = acos(sin(m[*n+k])*sin(m[*n+j])+cos(m[*n+k])*cos(m[*n+j])*cos(m[k]-m[j]))*6371000;
					if(dist<=(double)*l){
						w[mm] = j;
						m[2* *n+j] = cc;
						mm++;
						zeros++;
					}	
				}	
			}

		i++;
		if(zeros==*n){ break;}
		if(w[i]==0){
			cc++;
			k = 0;
			while(m[2* *n+ k]!=0) k++;
			w[i] = k;
			mm = i+1;
		}
	}
}

double SCMakse(int *m, int *n, double *theta_0, double *xq, double *l)
{
	int k, j, delt=0;
	double dist, zsum=0;
	
	for(k=0; k<*n; k++){
	
		if(k>0)
			for(j=k-1; j>=0 && (m[k]-m[j])<=*l; j--){
				if(m[2* *n + j]==0 && abs(m[*n+k]-m[*n+j])<=*l){
					dist = sqrt(((m[k]-m[j])*(m[k]-m[j]))+((m[*n+k]-m[*n+j])*(m[*n+k]-m[*n+j])));
					if(dist==(double)*l){
						delt++;
						zsum += (m[3* *n+j]-*xq)*(m[3* *n+k]-*xq);
						
					}
						
				}	
			}
		if(k<*n-1)
			for(j=k+1; j<*n && (m[j]-m[k])<=*l; j++){
				if(m[2* *n + j]==0 && abs(m[*n+k]-m[*n+j])<=*l){
					dist = sqrt(((m[k]-m[j])*(m[k]-m[j]))+((m[*n+k]-m[*n+j])*(m[*n+k]-m[*n+j])));
					if(dist==(double)*l){
						delt++;
						zsum += (m[3* *n+j]-*xq)*(m[3* *n+k]-*xq);
					}	
				}	
			}
	}
	return (zsum/(delt* *theta_0));
}



int min(int a, int b)
{
	if(a<=b) return(a);
	else return(b);
}

int max(int a, int b)
{
	if(a>b) return(a);
	else return(b);
}

void ccaBuffED(int *m, int *nr, int *nc, int *sz)
{
	int i,j,k,l,s=*sz;
		for(i=0; i<*nc; i++){
			for(j=0; j<*nr; j++){
					if(m[i**nr+j]==1)
						for(k=max(0,i-s); k<=min(*nc-1,i+s); k++)
							for(l=max(0,j-s); l<=min(*nr-1,j+s); l++)
								if(sqrt(((k-i)*(k-i))+((l-j)*(l-j)))<=s && m[k**nr+l]==0)
									m[k**nr+l] = -1;
			}
		}
}

