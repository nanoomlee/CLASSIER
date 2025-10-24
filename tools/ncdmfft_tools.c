#include "trigonometric_integrals.h"
#include "ncdmfft_tools.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#define J _Complex_I

#define sqr(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))

int G0_coeff (Gcoefficients * coeff, double xi_i, double xi_m, double xi_f, double G_i, double Gp_i, double G_m, double Gp_m, double G_f, double Gp_f) {

  double DeltaG = G_m - G_i;
  double DeltaGp = Gp_m - Gp_i;
  double dGdx = (G_m - G_i)/(xi_m - xi_i);
  double dxi = xi_m - xi_i;
  double beta = -2.0*dGdx/dxi/dxi + (Gp_m + Gp_i)/dxi/dxi;
  double xbeta = ( -DeltaGp /beta/dxi + xi_m + xi_i)/2.0;
  double a1 = G_i - dGdx*xi_i - beta * xbeta *xi_i*xi_m;
  double b1 = dGdx  + beta *(xbeta*xi_i + xi_i*xi_m + xi_m*xbeta);
  double c1 = -beta * (xi_m + xi_i + xbeta);
  double d1 = beta;

  DeltaG = G_f - G_m;
  DeltaGp = Gp_f - Gp_m;
  dGdx = (G_f - G_m)/(xi_f - xi_m);
  dxi = xi_f - xi_m;
  beta = -2.0*dGdx/dxi/dxi + (Gp_f + Gp_m)/dxi/dxi;
  xbeta = ( -DeltaGp /beta/dxi + xi_f + xi_m)/2.0;
  double a2 = G_m - dGdx*xi_m - beta * xbeta *xi_m*xi_f;
  double b2 = dGdx  + beta *(xbeta*xi_m + xi_m*xi_f + xi_f*xbeta);
  double c2 = -beta * (xi_f + xi_m + xbeta);
  double d2 = beta;
  coeff->a1 = a1; coeff->b1 = b1; coeff->c1 = c1; coeff->d1 = d1;
  coeff->a2 = a2; coeff->b2 = b2; coeff->c2 = c2; coeff->d2 = d2;

  return _SUCCESS_;
}

// Function to compute the six I0 values
int compute_I0_vectorized_sum(double *xi_list,          // xi's
                              int Nxi,                  // # of xi's
                              double xi1,               // xi1
                              double xi0,               // xi0
                              Gcoefficients * Acoeff,   // the coefficients for the A term
                              Gcoefficients * Bcoeff,   // the coefficients for the B term
                              double *I0){
   
  double ****calI;
  int n,l,index_type;
  double xi,c0,s0,ci,si,Si0,Sii,c1,s1,xid,Si1;
  double a, b, c, d, a2, b2, c2, d2;
  ErrorMsg error_message;

  // Using Table II. in 2506.01956

  calI = (double ****)malloc(4 * sizeof(double ***));
  for (index_type = 0; index_type < 2; index_type++) {
    calI[index_type] = (double ***)malloc(3 * sizeof(double **));
    for (l = 0; l < 3; l++) {
      calI[index_type][l] = (double **)malloc(4 * sizeof(double *));
      for (n = 0; n < 4; n++) {
        calI[index_type][l][n] = (double *)malloc(3 * sizeof(double));
      }
    }
  }

  // Now we'll evaluate the integrals at xi=0
  //  l=0
  for (index_type=0;index_type<2;index_type++){
    if (index_type==0){
      for (l=0;l<3;l++){
        if (l==0){
          calI[index_type][l][0][0] = 0.0;  calI[index_type][l][1][0] = -1.0;
          calI[index_type][l][2][0] = 0.0;  calI[index_type][l][3][0] = (2.0);
        }
        else if (l==1){
          //  l=1
          calI[index_type][l][0][0] = -1.0;  calI[index_type][l][1][0] = 0.0;
          calI[index_type][l][2][0] = -2.0;  calI[index_type][l][3][0] = 0.0;
        }
        else if (l==2){
          //  l=2
          calI[index_type][l][0][0] = 0.0;  calI[index_type][l][1][0] = -2.0;
          calI[index_type][l][2][0] = 0.0;  calI[index_type][l][3][0] = (-8.0);
        }
      }
    }
    else if (index_type==1){
      for (l=0;l<3;l++){
        if (l==0){
          //  l=0    this is actually j_0''
          calI[index_type][l][0][0] = 0.0;  calI[index_type][l][1][0] = -1.0;
          //calI[2][0] = 0.0;  calI[3][0] = -10.0;
          calI[index_type][l][2][0] = 0.0;  calI[index_type][l][3][0] = -6.0;
        }
        else if (l==1){
          //  l=1    this is actually j_1''
          calI[index_type][l][0][0] = 1.0/3.0;  calI[index_type][l][1][0] = 0.0;
          calI[index_type][l][2][0] = -2.0; calI[index_type][l][3][0] = 0.0;
        }
        else if (l==2){
          //  l=2    this is actually j_2''
          calI[index_type][l][0][0] = 0.0;
          calI[index_type][l][1][0] = 0.0;
          calI[index_type][l][2][0] =0.0;
          calI[index_type][l][3][0] = -12.0;
        }
      }
    }
  }

  for (int i=0;i<Nxi;i++){
    xi = xi_list[i];
    // First, some quantities
    c0=cos(xi);   s0=sin(xi);
    ci=1.0;    si=0.0;    // cos(0), sin(0)
    sine_integral(xi, &Si0, error_message);  Sii=0.0;  // sin integral
      
    // First we'll evaluate the integrals at xi
    for (index_type=0;index_type<2;index_type++){
      if (index_type==0){
        for (l=0;l<3;l++){
          if (l==0){
            //  l=0
            calI[index_type][l][0][2] = Si0;  calI[index_type][l][1][2] = -c0;  
            calI[index_type][l][2][2] = -xi*c0+s0;  calI[index_type][l][3][2] = (2.0-xi*xi)*c0+2.0*xi*s0;
          }
          else if (l==1){
            //  l=1
            calI[index_type][l][0][2] = -s0/xi;  calI[index_type][l][1][2] = Si0-s0;
            calI[index_type][l][2][2] = -2.0*c0-xi*s0;  calI[index_type][l][3][2] = -3.0*xi*c0+(3.0-xi*xi)*s0;
          }
          else if (l==2){
            //  l=2
            calI[index_type][l][0][2] = (1.5/xi)*(c0-s0/xi)+Si0/2.0;  calI[index_type][l][1][2] = c0-3.0*s0/xi;
            calI[index_type][l][2][2] = xi*c0-4.0*s0+3.0*Si0;  calI[index_type][l][3][2] = (xi*xi-8.0)*c0-5.0*xi*s0;
          }
        }
      }
      else if (index_type==1){
        for (l=0;l<3;l++){
          if (l==0){
            //  l=0    this is actually j_0''
            calI[index_type][l][0][2] = c0/xi-s0/sqr(xi);  calI[index_type][l][1][2] = c0-2.0*s0/xi;
            calI[index_type][l][2][2] = xi*c0-3.0*s0+2.0*Si0;  calI[index_type][l][3][2] = (xi*xi-6.0)*c0-4.0*xi*s0;
          }
          else if (l==1){
            //  l=1    this is actually j_1''
            calI[index_type][l][0][2] = 2.0*c0/sqr(xi)+s0/xi*(1.0-2.0/sqr(xi));  calI[index_type][l][1][2] = 3.0*c0/xi+(xi*xi-3.0)*s0/sqr(xi);
            calI[index_type][l][2][2] =4.0*c0+(xi-6.0/xi)*s0;  calI[index_type][l][3][2] = 5.0*xi*c0-(11.0-xi*xi)*s0+6.0*Si0;
          }
          else if (l==2){
            //  l=2    this is actually j_2''
            calI[index_type][l][0][2] = (-(xi*(-9.0 + xi*xi)*c0) + (-9.0 + 4.0*xi*xi)*s0)/pow(xi,4.0);  
            calI[index_type][l][1][2] = ((12.0-xi*xi)*xi*c0+(5.0*xi*xi-12.0)*s0)/pow(xi,3.0);
            calI[index_type][l][2][2] = (18.0-xi*xi)*c0/xi-(6.0*(3.0-xi*xi)*s0)/sqr(xi)+Si0;
            calI[index_type][l][3][2] = (24.0-xi*xi)*c0+(7.0*xi*xi-36.0)*s0/xi;
          }
        }
      }
    }

    // Now we'll evaluate the integrals at xi-xi1 if needed
    if (xi > xi1) {
      xid = xi - xi1;            
      c1 = cos(xid);
      s1 = sin(xid);
      sine_integral(xid, &Si1, error_message); 

      for (index_type=0;index_type<2;index_type++){
        if (index_type==0){
          for (l=0;l<3;l++){
            if (l==0){
              //  l=0
              calI[index_type][l][0][1] = Si1;
              calI[index_type][l][1][1] = -c1;
              calI[index_type][l][2][1] = -xid * c1 + s1;
              calI[index_type][l][3][1] = (2.0 - xid * xid) * c1 + 2.0 * xid * s1;
            }
            else if (l==1){
              //  l=1
              calI[index_type][l][0][1] = -s1 / xid;
              calI[index_type][l][1][1] = Si1 - s1;
              calI[index_type][l][2][1] = -2.0 * c1 - xid * s1;
              calI[index_type][l][3][1] = -3.0 * xid * c1 + (3.0 - xid * xid) * s1;
            }
            else if (l==2){
              //  l=2
              calI[index_type][l][0][1] = (1.5 / xid) * (c1 - s1 / xid) + Si1 / 2.0;
              calI[index_type][l][1][1] = c1 - 3.0 * s1 / xid;
              calI[index_type][l][2][1] = xid * c1 - 4.0 * s1 + 3.0 * Si1;
              calI[index_type][l][3][1] = (xid * xid - 8.0) * c1 - 5.0 * xid * s1;
            }
          }
        }
        else if (index_type==1){
          for (l=0;l<3;l++){
            if (l==0){
              //  l=0    this is actually j_0''
              calI[index_type][l][0][1] = c1 / xid - s1 / sqr(xid);
              calI[index_type][l][1][1] = c1 - 2.0 * s1 / xid;
              calI[index_type][l][2][1] = xid * c1 - 3.0 * s1 + 2.0 * Si1;
              calI[index_type][l][3][1] = (xid * xid - 6.0) * c1 - 4.0 * xid * s1;
            }
            else if (l==1){
              //  l=1    this is actually j_1''
              calI[index_type][l][0][1] = 2.0 * c1 / sqr(xid) + s1 / xid * (1.0 - 2.0 / sqr(xid));
              calI[index_type][l][1][1] = 3.0 * c1 / xid + (xid * xid - 3.0) * s1 / sqr(xid);
              calI[index_type][l][2][1] = 4.0 * c1 + (xid - 6.0 / xid) * s1;
              calI[index_type][l][3][1] = 5.0 * xid * c1 - (11.0 - xid * xid) * s1 + 6.0 * Si1;
            }
            else if (l==2){
              //  l=2    this is actually j_2''
              calI[index_type][l][0][1] = (-(xid * (-9.0 + xid * xid) * c1) + (-9.0 + 4.0 * xid * xid) * s1) / pow(xid, 4.0);
              calI[index_type][l][1][1] = ((12.0 - xid * xid) * xid * c1 + (5.0 * xid * xid - 12.0) * s1) / pow(xid, 3.0);
              calI[index_type][l][2][1] = (18.0 - xid * xid) * c1/xid - (6.0 * (3.0 - xid * xid) * s1) / sqr(xid) + Si1;
              calI[index_type][l][3][1] = (24.0 - xid * xid) * c1 + (7.0 * xid * xid - 36.0) * s1 / xid;
            }
          }
        }
      }
    }

    // Now we'll evaluate the convolution I_0(xi) for the six values of l
    // we'll first do the A terms
    // we need to do xi<x1 and xi>x1 separately
    for (index_type=0;index_type<2;index_type++){
      if (index_type==0) {
        a=Acoeff->a1; b=Acoeff->b1; c=Acoeff->c1; d=Acoeff->d1; a2=Acoeff->a2; b2=Acoeff->b2; c2=Acoeff->c2; d2=Acoeff->d2;
      }
      else if (index_type==1) {
        a=Bcoeff->a1; b=Bcoeff->b1; c=Bcoeff->c1; d=Bcoeff->d1; a2=Bcoeff->a2; b2=Bcoeff->b2; c2=Bcoeff->c2; d2=Bcoeff->d2;
      }
      for (l=0;l<3;l++){
        if(xi<=xi1) {
          I0[i+l*Nxi] += (a+b*xi+c*xi*xi+d*xi*xi*xi)*(calI[index_type][l][0][2]-calI[index_type][l][0][0]) 
              - (b+2.0*c*xi+3.0*d*xi*xi)*(calI[index_type][l][1][2]-calI[index_type][l][1][0]) + (c+3.0*d*xi)*(calI[index_type][l][2][2]-calI[index_type][l][2][0])
              - d*(calI[index_type][l][3][2]-calI[index_type][l][3][0]);
          if (xi==0) I0[i+l*Nxi]=0;
        }
        else {
          I0[i+l*Nxi] += (a+b*xi+c*xi*xi+d*xi*xi*xi)*(calI[index_type][l][0][2]-calI[index_type][l][0][1]) 
              - (b+2.0*c*xi+3.0*d*xi*xi)*(calI[index_type][l][1][2]-calI[index_type][l][1][1]) + (c+3.0*d*xi)*(calI[index_type][l][2][2]-calI[index_type][l][2][1])
              - d*(calI[index_type][l][3][2]-calI[index_type][l][3][1]);
          I0[i+l*Nxi] += (a2+b2*xi+c2*xi*xi+d2*xi*xi*xi)*(calI[index_type][l][0][1]-calI[index_type][l][0][0]) 
              - (b2+2.0*c2*xi+3.0*d2*xi*xi)*(calI[index_type][l][1][1]-calI[index_type][l][1][0])
              + (c2+3.0*d2*xi)*(calI[index_type][l][2][1]-calI[index_type][l][2][0]) - d2*(calI[index_type][l][3][1]-calI[index_type][l][3][0]);
        }
      }
    }
  }
  // Free allocated memory
  for (index_type=0;index_type<2;index_type++){
    for (l = 0; l < 3; l++) {
      for (n = 0; n < 4; n++) free(calI[index_type][l][n]);
      free(calI[index_type][l]);
    }
    free(calI[index_type]);
  }
  free(calI);
  return _SUCCESS_;
}

// compute all the K arrays at once
int compute_all_jffts(double x0,
                      int N, 
                      double *k_out, 
                      double scale, 
                      double complex *** kernels_fft){ 

  // Helper functions for C_Ei, C_Log, and C_Exp
  //double gsl_sf_Si(double u);     // Sine integral
  //double gsl_sf_Ci(double u);     // Cosine integral
  double complex C_Exp(double x0, double k);   // Exponential difference term
  double complex Eithing,C_Exp_val,y,y1;
  double j0,j1,j2,j3;
  double complex *j0_fft=kernels_fft[0][0], *j1_fft=kernels_fft[0][1],*j2_fft=kernels_fft[0][2];
  double complex *ddj0_fft=kernels_fft[1][0], *ddj1_fft=kernels_fft[1][1],*ddj2_fft=kernels_fft[1][2];
    
  //double dk = kmax / (N - 1); // Step size in k
  //double dk = 2*kmax / (N - 1); // Step size in k

  y1 = cexp(J*x0);
  j0 = sin(x0)/x0;
  j1 = -cos(x0)/x0 + sin(x0)/sqr(x0);
  j2 = (3.0/sqr(x0) - 1.0)*sin(x0)/x0 - 3.0*cos(x0)/sqr(x0);
  j3 = (15.0/cube(x0) - 6.0/x0)*sin(x0)/x0 - (15.0/sqr(x0)-1.0)*cos(x0)/x0;
   
  int N_mid;
  double si_p, ci_p, si_m, ci_m;
  double k;
  ErrorMsg error_message;
  if (N%2==0) N_mid = N/2;
  else N_mid = (N-1)/2;

  for (int i = N_mid; i < N; i++) {
    //k = (i-N_mid) * dk; // Current k value
    k = k_out[i]*_PI_/scale; // Current k value
    cosine_integral ((k+1.0)*x0, &ci_p, error_message);
    cosine_integral (fabs(k-1.0)*x0, &ci_m,error_message);
    sine_integral ((k+1.0)*x0, &si_p,error_message);
    sine_integral ((k-1.0)*x0, &si_m,error_message);

    // Compute C_Ei, C_Log, and C_Exp
    //Eithing = J*((gsl_sf_Si((k+1.0)*x0)) - (gsl_sf_Si((k-1.0)*x0)))
      //  + ( (gsl_sf_Ci((k+1.0)*x0)  - gsl_sf_Ci(fabs((k-1.0))*x0)))  - log((k+1.0)/fabs(k-1.0));
    Eithing = J*(si_p - si_m) + ( (ci_p  - ci_m))  - log((k+1.0)/fabs(k-1.0));
    y = cexp(J*k * x0);
    C_Exp_val = y*(y1 - 1.0/y1);

    // Compute all the fft arrays
    j0_fft[i] = (-J / 2.0 * (Eithing));
    j1_fft[i] = (k / 2.0 * (Eithing) + J / (2.0 * x0) * C_Exp_val + 1.0);
    j2_fft[i] = ((J / 4.0) * (3.0 * k * k - 1.0) * (Eithing)
                + (-3.0 * k * x0 + 3.0 * x0 + 3.0 * J) / (4.0 * x0 * x0) * C_Exp_val
                + (3.0 / (2.0 * x0)) * cexp(J * (k - 1.0) * x0)
                + (3.0 * J * k) / 2.0);
    ddj0_fft[i] = (-k*k * j0_fft[i] - J * k * (y*j0 -1.0)  - j1*y);
    ddj1_fft[i] = (-k*k * j1_fft[i] - J * k * y*j1 + (y * (j0-2.0*j2)/3.0 -1./3.0));
    ddj2_fft[i] = (-k*k * j2_fft[i] - J * k * y*j2 + (y * (2.0*j1-3.0*j3)/5.0 ));
    if (i!=N_mid) {
      j0_fft[2*N_mid-i] = conj(j0_fft[i]);
      j1_fft[2*N_mid-i] = conj(j1_fft[i]);
      j2_fft[2*N_mid-i] = conj(j2_fft[i]);
      ddj0_fft[2*N_mid-i] = conj(ddj0_fft[i]);
      ddj1_fft[2*N_mid-i] = conj(ddj1_fft[i]);
      ddj2_fft[2*N_mid-i] = conj(ddj2_fft[i]);
    }
  }

  if (N%2==0) {
    //k = (N-N_mid)*dk;
    k = k_out[0]*_PI_/scale;
    cosine_integral ((k+1.0)*x0, &ci_p, error_message);
    cosine_integral (fabs(k-1.0)*x0, &ci_m,error_message);
    sine_integral ((k+1.0)*x0, &si_p,error_message);
    sine_integral ((k-1.0)*x0, &si_m,error_message);

    Eithing = I*(si_p - si_m)
              + ( (ci_p  - ci_m))  - log((k+1.0)/fabs(k-1.0));
    y = cexp(I*k * x0);
    C_Exp_val = y*(y1 - 1.0/y1);

    // Compute all the fft arrays
    j0_fft[0] = (-I / 2.0 * (Eithing));
    j1_fft[0] = (k / 2.0 * (Eithing) + I / (2.0 * x0) * C_Exp_val + 1.0);
    j2_fft[0] = ((I / 4.0) * (3.0 * k * k - 1.0) * (Eithing)
                + (-3.0 * k * x0 + 3.0 * x0 + 3.0 * I) / (4.0 * x0 * x0) * C_Exp_val
                + (3.0 / (2.0 * x0)) * cexp(I * (k - 1.0) * x0)
                + (3.0 * I * k) / 2.0);
    ddj0_fft[0] = (-k*k * j0_fft[0] - I * k * (y*j0 -1.0)  - j1*y);
    ddj1_fft[0] = (-k*k * j1_fft[0] - I * k * y*j1 + (y * (j0-2.0*j2)/3.0 -1./3.0));
    ddj2_fft[0] = (-k*k * j2_fft[0] - I * k * y*j2 + (y * (2.0*j1-3.0*j3)/5.0 ));
    j0_fft[0] = conj(j0_fft[0]);
    j1_fft[0] = conj(j1_fft[0]);
    j2_fft[0] = conj(j2_fft[0]);
    ddj0_fft[0] = conj(ddj0_fft[0]);
    ddj1_fft[0] = conj(ddj1_fft[0]);
    ddj2_fft[0] = conj(ddj2_fft[0]);

  }
  /* 
  if (N%2==0) N_mid = N/2;
  else N_mid = (N+1)/2;
  for (int i = 0; i < N_mid; i++) {
    double k = i * dk; // Current k value
    double si_p, ci_p, si_m, ci_m;
    cosine_integral ((k+1.0)*x0, &ci_p, error_message);
    cosine_integral (fabs(k-1.0)*x0, &ci_m,error_message);
    sine_integral ((k+1.0)*x0, &si_p,error_message);
    sine_integral ((k-1.0)*x0, &si_m,error_message);

    // Compute C_Ei, C_Log, and C_Exp
    //Eithing = I*((gsl_sf_Si((k+1.0)*x0)) - (gsl_sf_Si((k-1.0)*x0)))
      //  + ( (gsl_sf_Ci((k+1.0)*x0)  - gsl_sf_Ci(fabs((k-1.0))*x0)))  - log((k+1.0)/fabs(k-1.0));
    Eithing = I*(si_p - si_m)
        + ( (ci_p  - ci_m))  - log((k+1.0)/fabs(k-1.0));
    y = cexp(I*k * x0);
    C_Exp_val = y*(y1 - 1.0/y1);

    // Compute all the fft arrays
    j0_fft[i] = (-I / 2.0 * (Eithing));
    j1_fft[i] = (k / 2.0 * (Eithing) + I / (2.0 * x0) * C_Exp_val + 1.0);
    j2_fft[i] = ((I / 4.0) * (3.0 * k * k - 1.0) * (Eithing)
                + (-3.0 * k * x0 + 3.0 * x0 + 3.0 * I) / (4.0 * x0 * x0) * C_Exp_val
                + (3.0 / (2.0 * x0)) * cexp(I * (k - 1.0) * x0)
                + (3.0 * I * k) / 2.0);
    ddj0_fft[i] = (-k*k * j0_fft[i] - I * k * (y*j0 -1.0)  - j1*y);
    ddj1_fft[i] = (-k*k * j1_fft[i] - I * k * y*j1 + (y * (j0-2.0*j2)/3.0 -1./3.0));
    ddj2_fft[i] = (-k*k * j2_fft[i] - I * k * y*j2 + (y * (2.0*j1-3.0*j3)/5.0 ));
    if (i!=0) {
      j0_fft[N-i] = conj(j0_fft[i]);
      j1_fft[N-i] = conj(j1_fft[i]);
      j2_fft[N-i] = conj(j2_fft[i]);
      ddj0_fft[N-i] = conj(ddj0_fft[i]);
      ddj1_fft[N-i] = conj(ddj1_fft[i]);
      ddj2_fft[N-i] = conj(ddj2_fft[i]);
    }
  }*/
  return _SUCCESS_;
}

/**
 * @brief calculates the spherical besselfunction of degree l at point x.
 * @param l 
 * @param x 
 * @param y 
 * @return int 
 */
int spherical_bessel_function(int l,
                              double x,
                              double *y){ 
  if (x>0.01) {
    switch (l)
    {
    case 0:
      *y = sin(x)/x;
      break;
    case 1:
      *y = sin(x)/x/x - cos(x)/x;
      break;
    case 2:
      *y = sin(x)/x*(3.0/x/x-1) - 3.0 * cos(x)/x/x;
      break;
    case 3:
      *y = sin(x) / pow(x,4) * 15.0 - sin(x) / pow(x,2) * 6.0 - cos(x) / pow(x,3) * 15.0 + cos(x) / x;
      break;
    case 4:
      *y = cos(x) / pow(x,4) * (10.0 *pow(x,2) - 105.0) + sin(x) / pow(x,5) * (pow(x,4) - 45.0 * pow(x,2) + 105.0);
      break;
    default:
      return _FAILURE_;
    }
  }
  else {
    switch (l)
    {
    case 0:
      *y = 1-pow(x,2)/6.0+pow(x,4)/120.0;
      break;
    case 1:
      *y = x/3.0-pow(x,3)/30.0;
      break;
    case 2:
      *y = pow(x,2)/15.0-pow(x,4)/210.0;
      break;
    case 3:
      *y = pow(x,3)/105.0;
      break;
    case 4:
      *y = pow(x,4)/945.0;
      break;
    default:
      return _FAILURE_;
    }
  }
  return _SUCCESS_;
}


/**
 * @brief calculates the n_th derivative of l_th spherical bessel function at point x. 
 * We use the relation that j_l' = j_(l-1)-((l+1)/x) * j_l
 * 
 * @param n 
 * @param l 
 * @param x 
 * @param y 
 * @return int 
 */
int derivative_of_spherical_bessel(int n,
                                  int l,
                                  double x,
                                  double *y){  
  // Some values are not defined at x=0 even though they converge to some fixed value. 
  // These quantities are recovered by a cut-off!
  if (x<1e-6) {
    x=1e-6;
  }
  int i = n;
  if (n==0){
    spherical_bessel_function(l, x, y);
  }
  else if(n==1){
    double a, b;
    spherical_bessel_function(l,x,&a);
    spherical_bessel_function(l+1,x,&b);
    * y = l * a / x - b;
  }
  else if(n==2){
    double a, b;
    spherical_bessel_function(l,x,&a);
    spherical_bessel_function(l+1,x,&b);
    * y = ( (l-1) * l / pow(x,2) - 1 ) * a + b * 2.0 /x;
  }
  else if(n==3){
    double a, b;
    spherical_bessel_function(l,x,&a);
    spherical_bessel_function(l+1,x,&b);
    * y =  ( (l-2.0) * ( (l-1.0) * l - x * x ) * a + b * x * (x*x -6.0 - l*(l+1.0) ) ) / x / x / x;
  }
  else{
      return _FAILURE_;
  }
  return _SUCCESS_;
}


/**
 * @brief calculate the auxilliary function W_ll_prime which is required in the ncdmfft acceleration.
 * See 2201.11129 equation (4-5). It will calculate:
 * W_ll_prime(x) = i^l_prime * P_l_prime(i * d/dx) j_l(x)
 */
int auxiliary_function_ncdmfft(int l,
                               int l_prime,
                               double x,
                               double *y){   
  double a,b;
  switch (l_prime)
  {
  case 0:
      derivative_of_spherical_bessel(0,l,x,&a);
      *y = a;
      break;
  case 1:
      derivative_of_spherical_bessel(1,l,x,&a);
      *y = -a;
      break;
  case 2:
      derivative_of_spherical_bessel(2,l,x,&a);
      derivative_of_spherical_bessel(0,l,x,&b);
      *y = 3.0/2.0*a+1.0/2.0*b;
      break;
  case 3:
      derivative_of_spherical_bessel(3,l,x,&a);
      derivative_of_spherical_bessel(1,l,x,&b);
      *y = -5.0/2.0*a-3.0/2.0*b;
      break;
  default:
      return _FAILURE_;
  }
  return _SUCCESS_;
}

