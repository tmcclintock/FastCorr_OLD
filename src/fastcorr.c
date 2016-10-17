#include "gsl/gsl_spline.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define PI 3.14159265358979

/*
  This is the integrand. It is a cubic spline of the power spectrum
  between kmin and kmax. Beyond the bounds a power law is assumed.

  x: location in kR space to evaluate the integrand
  R: tangential radial distance - units in either Mpc/h or Mpc
  k: wavenumber - units in either h/Mpc or Mpc^-1
  P: power spectrum - units in either (h/Mpc)^3 or Mpc^-3
  Nk: number of k and P points
  Psl: a spline of P(k)
  acc: a spline accelerator
 */
double get_P(double x,double R,double*k,double*P,int Nk,
	     gsl_spline*Pspl,gsl_interp_accel*acc){
  double ki = x/R;
  double kmin = k[0];
  double kmax = k[Nk-1];
  double alpha,A;
  if (ki < kmin){
    alpha = log(P[1]/P[0])/log(k[1]/k[0]);
    A = P[0]/pow(k[0],alpha);
    return A*pow(ki,alpha);
  }else if (ki > kmax){
    alpha = log(P[Nk-1]/P[Nk-2])/log(k[Nk-1]/k[Nk-2]);
    A = P[Nk-1]/pow(k[Nk-1],alpha);
    return A*pow(ki,alpha);
  }// Assume power laws at ends
  return gsl_spline_eval(Pspl,ki,acc);
}

/*
  This is the routine where xi(R) is actually calculated.

  R: tangential radial distance - units in either Mpc/h or Mpc
  k: wavenumber - units in either h/Mpc or Mpc^-1
  P: power spectrum - units in either (h/Mpc)^3 or Mpc^-3
  Nk: number of k and P points
  N: number of roots of j_0 to evaluate
  h: step size for the quadrature routine
  x: locations in kR space to evaluate the integrand
  sinx: precalculated sin(x)
  dpsi: precalculated function dpsi(x) (see Ogata et al. 2005)
 */
double calc_corr_at_R(double R,double*k,double*P,
		      int Nk,int N,double h,
		      double*x,double*sinx,double*dpsi,
		      gsl_spline*Pspl,gsl_interp_accel*acc){
  double f;

  double sum = 0;
  int i;
  for(i=0;i<N;i++){
    f = x[i]*get_P(x[i],R,k,P,Nk,Pspl,acc);
    sum += f*sinx[i]*dpsi[i];
  }

  return sum/(R*R*R*PI);
}

/*
  Calculates the domain locations in kR space to evaluate the
  integrand.

  N: number of roots of j_0 to evaluate
  h: step size for the quadrature routine
  x: locations in kR space to evaluate the integrand (output)
  sinx: precalculated sin(x) (output)
  dpsi: precalculated function dpsi(x) (see Ogata et al. 2005) (output)
 */
int calc_domain(int N,double h,double*x,double*sinx,double*dpsi){
  double zero,t,psi,PIsinht;
  const double PI_h = PI/h;
  const double PI_2 = PI/2.;
  int i;
  for(i=0;i<N;i++){
    zero = i+1;
    t = h*zero;
    psi = h*zero*tanh(sinh(h*zero)*PI_2);
    x[i] = psi*PI_h;
    sinx[i]=sin(x[i]);
    PIsinht = PI*sinh(t);
    dpsi[i] = (PI*t*cosh(t)+sinh(PIsinht))/(1+cosh(PIsinht));
    if (dpsi[i]!=dpsi[i]) dpsi[i]=1.0;
  }
  return 0;
}

/*
  A function call to interface with FastCorr.py.

  k: wavenumber - units in either h/Mpc or Mpc^-1
  P: power spectrum - units in either (h/Mpc)^3 or Mpc^-3
  Nk: number of k and P points
  R: tangential radial distances - units in either Mpc/h or Mpc
  NR: number of radial points
  xi: correlation function - unitless (output)
  N: number of roots of j_0 to evaluate
  h: step size for the quadrature routine
 */
int calc_corr(double*k,double*P,int Nk,double*R,double*xi,
	      int NR,int N, double h){
  int i,status;
  gsl_spline*Pspl = gsl_spline_alloc(gsl_interp_cspline,Nk);
  gsl_spline_init(Pspl,k,P,Nk);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();

  double*x=(double*)malloc(N*sizeof(double));
  double*sinx=(double*)malloc(N*sizeof(double));
  double*dpsi=(double*)malloc(N*sizeof(double));
  status = calc_domain(N,h,x,sinx,dpsi);
  for(i=0;i<NR;i++)
    xi[i] = calc_corr_at_R(R[i],k,P,Nk,N,h,x,sinx,dpsi,Pspl,acc);
  
  gsl_spline_free(Pspl),gsl_interp_accel_free(acc);
  return 0;
}
