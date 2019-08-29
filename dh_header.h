// LIBRARIES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <fftw3.h>

// MACROS
#define IJ2K(a,b) (a+b*(b+1)/2) // Converts matrix indices
#define ARRAY_REALLOC_FACTOR 2.0 // Factor for realloc
#define MIN(a,b) ( (a < b) ? (a) : (b))
#define MAX(a,b) ( (a > b) ? (a) : (b))
#define ABS(a) ((a > 0) ? (a): (-a) )
#define ALLOC(p,n)  (p)=malloc( (n) * sizeof(*(p))); if( (p) == NULL){printf("Allocation of '%s' failed. Terminate. \n", #p); exit(2); } 
#define FFT_ALLOC(p,n)  (p)=fftw_malloc( (n) * sizeof(*(p))); if( (p) == NULL){printf("Allocation of '%s' failed. Terminate. \n", #p); exit(2); } 
#define REALLOC(p,n)  (p)=realloc( (p) , (n) * sizeof(*(p))); if( (p) == NULL){printf("Allocation of '%s' failed. Terminate. \n", #p); exit(2); } 

// STRUCT
// Complex numbers
typedef struct complex_z
{
	 double r;
	 double i;
} complex_z;

// FUNCTIONS
void initialise( fftw_complex** ,  fftw_complex** , fftw_complex** ,  fftw_complex** ,  double** , complex_z** , long N, gsl_rng**, const gsl_rng_type**, int);
void initialise_trajectory ( double**, long);
void initialise_inverse_correlation_matrix(double***, long );
void write_correlation_exponents(double*, long, double, double);
void write_correlation(fftw_complex*, double*, long);
void write_inverse_correlation_matrix(double **, long, double);
void generate_random_vector(complex_z*, fftw_complex*, fftw_complex*,long);
void set_to_zero(double*, long);
void integrate_noise(double*, fftw_complex*, double, double, long, double);
void find_fpt(double*, double*, double, long);
void fpt_to_zvar(double, double, double);

extern gsl_rng *r;
extern int seed;
