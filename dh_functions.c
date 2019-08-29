#include "dh_header.h"


void initialise( fftw_complex** correlation,  fftw_complex** circulant_eigenvalues,  fftw_complex** rndW,  fftw_complex** fracGN,  double** correlation_exponents, complex_z** randomComplexGaussian, long N, gsl_rng** r,const gsl_rng_type** T, int seed)
{
	// These are the objects that are N long (the increments)
	FFT_ALLOC(*correlation, 2*N);
	FFT_ALLOC(*circulant_eigenvalues, 2*N);
	FFT_ALLOC(*rndW, 2*N);
	FFT_ALLOC(*fracGN, 2*N);

	ALLOC(*correlation_exponents, 2*N);
	set_to_zero(*correlation_exponents, 2*N);
	ALLOC(*randomComplexGaussian, N);
	
	/* Initialises GSL Random Generator */
        gsl_rng_env_setup();
        *T = gsl_rng_default;
        *r = gsl_rng_alloc (*T);
        if(seed==-1) seed = ((int) (((int) clock() ) % 100000));
        printf("# RNG seed: %i\n", seed);
	gsl_rng_set(*r, seed);
}

void initialise_trajectory(  double ** fracbm, long N)
{
	// N increments give N+1 points
	ALLOC( *fracbm, N+1);
}


void write_correlation_exponents(double* correlation_exponents, long N, double invN, double hurst)
{
	int i;
	correlation_exponents[0] = 0.0; // 0^0 = 0.
	for(i=1; i<=N; i++)
	{
		correlation_exponents[i] =pow( (( double) i) * invN, ( (double) 2*hurst ));
	}
}

void write_correlation(fftw_complex* correlation, double* correlation_exponents, long N)
{
	correlation[0][0] = 2*(correlation_exponents[1] );
	correlation[0][1] = 0.0;
	correlation[N][0] = 0.0;
	correlation[N][1] = 0.0;
	
	int i;
	for(i=1; i<=N-1; i++)
	{
		correlation[2*N - i][0] = correlation[i][0] = (correlation_exponents[i+1] + correlation_exponents[i-1] - 2*correlation_exponents[i]);
		correlation[2*N - i][1] = correlation[i][1] = 0;
	}
}
void generate_random_vector(complex_z* randomvector, fftw_complex* rndW, fftw_complex *circulant_eigenvalues, long N)
{
	double sigma = 1.0;
	// This routine creates the random vector used in Davies Harte generation of subgrid
	int i;
	double invN = 1/((double) N);

	for(i = 0; i < N; i++)
		{
			randomvector[i].r= gsl_ran_gaussian_ziggurat(r,sigma);
			randomvector[i].i = gsl_ran_gaussian_ziggurat(r,sigma);
		}

	rndW[0][0] = sqrt(0.5* circulant_eigenvalues[0][0]*invN)*gsl_ran_gaussian_ziggurat(r,sigma);
	rndW[0][1] = 0.0;

	rndW[N][0] = sqrt(0.5*circulant_eigenvalues[N][0]*invN)*gsl_ran_gaussian_ziggurat(r,sigma);
	rndW[N][1] = 0.0;
	
	for(i=1; i < N; i++)
	{
		rndW[i][0] = (sqrt(0.25*circulant_eigenvalues[i][0]*invN) * (randomvector[i].r));
		rndW[i][1] = (sqrt(0.25*circulant_eigenvalues[i][0]*invN) * (randomvector[i].i));
		rndW[2*N-i][0] =rndW[i][0]; 
		rndW[2*N-i][1] = -rndW[i][1]; 
	}
}

void set_to_zero(double* pointer, long length)
{
	// generic function to "re-calloc" pointer
	long i;
	for(i = 0; i < length; i++)
	{
		pointer[i] = 0.0;
	}
}

void integrate_noise(double* fracbm, fftw_complex* fracGN,double lin_drift, double frac_drift, long N, double hurst)
{	
	// Find the first point to jump over the barrier (if exists). Then throw away all points behind. Take the appropiate inverse matrix and pass it on.
	fracbm[0] = 0.0;
	double delta_t = (1/((double) N));

	int i;
	for(i = 1; i <= N; i++)
	{
		fracbm[i] = (fracbm[i-1] + fracGN[i-1][0] + (lin_drift + frac_drift*(pow((((double)i)-0.5)*delta_t, (2*hurst - 1.0))))*delta_t); // Ito integration, fractional drift is evaluated at midpoint --> Check for consequences
	}
}

void find_fpt(double* fracbm, double* first_passage_times, double passage_height, long N)
{
	long i;
	*first_passage_times = 1.0;
	double t_i; // Time prior
	double t_j; // Time after crossing
	double delta_t = (1./((double) N));
	for(i = 0; i <= N; i++)
	{
		if(fracbm[i] > passage_height)
		{
			t_i = (i-1)*delta_t;
			t_j = i*delta_t;
			*first_passage_times = (t_i + (t_j - t_i)/(fracbm[i] - fracbm[i-1])*(passage_height - fracbm[i-1]));
			i = N+1;
			
		}
		/*
		 *
		 * ADAPT CROSSING TIME
		 * FINISH DAVIES HARTE
		 * WRITE TIMING SHELL THAT COMPARES D&H with ABC 
		 * MAKE PLOT
		 * WRITE ABOUT IT
		 * DON'T START CRYING*/}
	}


double time_time_correlation(double ti, double tj, double hurst)
{
	if(hurst == 0.5){return (2*MIN(ti,tj));}
	else{return ( pow(fabs(ti),2*hurst) + pow(fabs(tj),2*hurst) - pow(fabs(ti-tj),2*hurst));}
}


