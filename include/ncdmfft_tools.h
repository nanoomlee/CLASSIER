#include <complex.h>

typedef struct {   // a1,b1,c1,d1 are coefficients for G_1 and a2,b2,c2,d2 for G_2
  double a;
  double b;
  double c;
  double d;
} Gcoefficients2;

typedef struct {   // a1,b1,c1,d1 are coefficients for G_1 and a2,b2,c2,d2 for G_2
  double a1;
  double b1;
  double c1;
  double d1;
  double a2;
  double b2;
  double c2;
  double d2;
} Gcoefficients;

int G0_coeff2 (Gcoefficients2 * coeff, double xi_i, double xi_m, double G_i, double Gp_i, double G_m, double Gp_m);
int G0_coeff (Gcoefficients * coeff, double xi_i, double xi_m, double xi_f, double G_i, double Gp_i, double G_m, double Gp_m, double G_f, double Gp_f);

int compute_I0_vectorized_psi0only(double *xi_list,         //   xi
               int Nxi,
			   double xi1,         //   xi1
			   double xi0,             //   xi0
			   char type,
			   Gcoefficients * coeff,   // the coefficients for the A term
			   double *I0);            // array of six values

int compute_I0_vectorized_sum(double *xi_list,         //   xi
               int Nxi,
			   double xi1,         //   xi1
			   double xi0,             //   xi0
			   Gcoefficients * Acoeff,   // the coefficients for the A term
			   Gcoefficients * Bcoeff,   // the coefficients for the B term
			   double *I0);            // array of six values
int compute_I0_vectorized(double *xi_list,         //   xi
               int Nxi,
			   double xi1,         //   xi1
			   double xi0,             //   xi0
			   char type,
			   Gcoefficients * coeff,   // the coefficients for the A term
			   double *I0);            // array of six values

int compute_I0(double *xi_list,         //   xi
               int Nxi,
			   double xi1,         //   xi1
			   double xi0,             //   xi0
			   int l,
			   char type,
			   Gcoefficients * coeff,   // the coefficients for the A term
			   double *I0);            // array of six values

int compute_j0ffts(double x0,
                      int N, 
                      double *k_out,
                      double scale,
                      double complex ***kenels_fft);
                      
int compute_all_jffts(double x0,
                      int N, 
                      double *k_out,
                      double scale,
                      double complex ***kenels_fft);
                      //double complex *j0_fft,
                      //double complex *j1_fft, 
                      //double complex *j2_fft,
                      //double complex *ddj0_fft,
                      //double complex *ddj1_fft,
                      //double complex *ddj2_fft);

int compute_all_jffts_pos(double x0,
                      int N, 
                      double *k_out,
                      double scale,
                      double complex ***kenels_fft);

int auxiliary_function_ncdmfft(int l,
                              int l_prime,
                              double x,
                              double *y);

int derivative_of_spherical_bessel(int n,
                              int l,
                              double x,
                              double *y);

int spherical_bessel_function(int l,
                              double x,
                              double *y);

