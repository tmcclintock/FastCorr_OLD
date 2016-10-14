#include "gsl/gsl_spline.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979

//The power spectrum
double get_P(double x,double R,double*k,double*P,int Nk,gsl_spline*Pspl,gsl_inter_accel*acc){
  double kmin, kmax = k[0],k[Nk-1];
  double alpha,A;
  if (x/R < kmin){
    alpha = log(P[1]/P[0])/log(k[1]/k[0]);
    A = P[0]/pow(k[0],alpha);
    return A*pow(x/R,alpha);
  }else if(x/R > kmax){
    alpha = log(P[N-1]/P[N-2])/log(k[N-1]/k[N-2]);
    A = P[N-1]/pow(k[N-1],alpha);
    return A*pow(x/R,alpha);
  }//Assume power laws at ends
  return gsl_spline_eval(Pspl,x/R,acc);
}

double calc_corr_at_R(double R,double*k,double*P,int Nk,int N,double h){
  double*zeros = (double*)malloc(N*sizeof(double));
  double*phi = (double*)malloc(N*sizeof(double));
  double*x = (double*)malloc(N*sizeof(double));
  double*dphi = (double*)malloc(N*sizeof(double));
  double*w = (double*)malloc(N*sizeof(double));
  double*term = (double*)malloc(N*sizeof(double));

  gsl_spline*Pspl = gsl_spline_alloc(gsl_interp_cspline,Nk);
  gsl_spline_init(spline,k,P,Nk);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();

  double sum = 0;
  int i;
  for(i=0;i<N;i++){
    zeros[i] = i+1;
    phi[i] = h*zeros[i]*tanh(PI*sinh(h*zeros[i])/2.);
    x[i] = phi[i]*PI/h;
    dphi[i] = (PI*x[i]*cosh(x[i])+sinh(PI*sinh(x[i])))/(1+cosh(PI*sinh(x[i])));
    w[i] = cos(x[i])*sin(x[i])/(x[i]*cos(x[i])-sin(x[i]));
    Pi = get_P(x[i],R,k,P,Nk,Pspl,acc);
    term[i] = w[i]*x[i]*x[i]*Pi*dphi[i];
    sum += PI*term[i];
  }

  free(zeros),free(phi),free(x),free(dphi),free(w),free(term);
  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  return sum/(R*R*R*PI*PI);
}

int calc_corr(double*k,double*P,int Nk,double*R,double*xi,int NR,int N, double h){
  //int N = 200; //Arbitrary
  //double h = 0.005; //Arbitrary
  int i;
  for(i=0;i<NR;i++){
    xi[i] = calc_corr_at_R(R[i],k,P,Nk,N,h);
  }
  return 0;
}
