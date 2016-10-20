#include "gsl/gsl_spline.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.1415926535897

/*
  This is the integrand. It is a cubic spline of the power spectrum
  between kmin and kmax. Beyond the bounds a power law is assumed.

  x: location in kR space to evaluate the integrand
  R: tangential radial distance - units in either Mpc/h or Mpc
  k: wavenumber - units in either h/Mpc or Mpc^-1
  P: power spectrum - units in either (h/Mpc)^3 or Mpc^-3
  Nk: number of k and P points
  alpha: power law indices - first element is low-k, second element is high-k
  A: power law amplitude - first element is low-k, second element is high-k
  Psl: a spline of P(k)
  acc: a spline accelerator
 */
double get_P(double x,double R,double*k,double*P,int Nk,
	     double*alpha,double*A,
	     gsl_spline*Pspl,gsl_interp_accel*acc){
  double ki = x/R;
  double kmin = k[0];
  double kmax = k[Nk-1];
  if (ki < kmin){
    return A[0]*pow(ki,alpha[0]);
  }else if (ki > kmax){
    return A[1]*pow(ki,alpha[1]);
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
  alpha: power law indices - first element is low-k, second element is high-k
  A: power law amplitude - first element is low-k, second element is high-k
  sinx: precalculated sin(x)
  dpsi: precalculated function dpsi(x) (see Ogata et al. 2005)
 */
double calc_corr_at_R(double R,double*k,double*P,
		      int Nk,int N,double h,
		      double*x,double*sinx,double*dpsi,
		      double*alpha,double*A,
		      gsl_spline*Pspl,gsl_interp_accel*acc){
  double f;

  double sum = 0;
  int i;
  for(i=0;i<N;i++){
    f = x[i]*get_P(x[i],R,k,P,Nk,alpha,A,Pspl,acc);
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
  Compute the power law behavior of the power spectrum
  at its ends.

  k: wavenumber - units in either h/Mpc or Mpc^-1
  P: power spectrum - units in either (h/Mpc)^3 or Mpc^-3
  Nk: number of k and P points
  alpha: power law indices - first element is low-k, second element is high-k (output)
  A: power law amplitude - first element is low-k, second element is high-k (output)
 */
int calc_power_law(double*k,double*P,int Nk,double*alpha,double*A){
  alpha[0] = log(P[1]/P[0])/log(k[1]/k[0]);
  A[0] = P[0]/pow(k[0],alpha[0]);
  alpha[1] = log(P[Nk-1]/P[Nk-2])/log(k[Nk-1]/k[Nk-2]);
  A[1] = P[Nk-1]/pow(k[Nk-1],alpha[1]);
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
  double*alpha=(double*)malloc(2*sizeof(double));
  double*A=(double*)malloc(2*sizeof(double));
  status = calc_domain(N,h,x,sinx,dpsi);
  status|= calc_power_law(k,P,Nk,alpha,A);
  for(i=0;i<NR;i++)
    xi[i] = calc_corr_at_R(R[i],k,P,Nk,N,h,x,sinx,dpsi,alpha,A,Pspl,acc);
  
  free(x),free(sinx),free(dpsi),free(alpha),free(A);
  gsl_spline_free(Pspl),gsl_interp_accel_free(acc);
  return 0;
}
