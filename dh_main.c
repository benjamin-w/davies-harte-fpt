/* This code simulates a fractional Brownian Motion on the interval [0,1] using the Davies Harte method and finds the first passage time.
 * 
 * Authors: Benjamin Walter (Imperial College) , Kay Wiese (ENS Paris)
 * 
 * For any questions/feedback/comments --> b.walter16@imperial.ac.uk
 * 
 */

#include "dh_header.h"

// GSL RNG
const gsl_rng_type *T;
gsl_rng *r;
int seed;



int main(int argc, char *argv[])
{
	setlinebuf(stdout);
	// VARIABLES
	
	// physics
	// Default values, overwritten by getopt
	double frac_drift = 0.0;
	double lin_drift = 0.0;
	double hurst = 0.5; 	// Hurst parameter 0 < h < 1
	int g =8;		// 2^g is number of points
	
	// experiment
	int iteration = 10000;	// Size of ensemble
	int iter;
	

	
	// observables
	double passage_heights = 0.1;
	double first_passage_times = 0.0;
	int seed = -1;

	// ARGUMENTS
	opterr = 0;
	int c = 0;
        while( (c = getopt (argc, argv, "h:g:G:S:I:m:n:E:") ) != -1)
	{                switch(c)
                        {
				case 'm':
					lin_drift = ( - ( double) atof(optarg)); // Sign chosen in accordance with paper. 
					break;
                        	case 'n':
					frac_drift = (- ( double) atof(optarg));
					break;
                                case 'h':
                                        hurst = atof(optarg);
                                        break;
				case 'g':
					g = atoi(optarg);
					break;
				case 'S':
					seed = atoi(optarg);
					break;
				case 'I':
					iteration = atoi(optarg);
					break;
                       		default:
                                exit(EXIT_FAILURE);
                        }
	}// getopt ends


	//int i; // Index
	long N = ((long) pow(2,g));
	double invN = 1/(( double) N);

	fftw_complex *correlation, *circulant_eigenvalues, *rndW, *fracGN;
       	double *correlation_exponents, *fracbm;
	complex_z *randomComplexGaussian;
	fftw_plan p1, p2 ;

	
	// Header
	printf("# First passage times of fractional Brownian Motion with drift using Davies-Harte method\n# H: %g\n# System size: %i\n# Linear drift: %g\n# Fractional drift: %g\n# Iterations: %i\n", hurst, ((int) pow(2,g)), lin_drift, frac_drift, iteration);
	
	// Initialise observables
	initialise(&correlation, &circulant_eigenvalues, &rndW, &fracGN, &correlation_exponents, &randomComplexGaussian, N, &r, &T, seed);
	initialise_trajectory(&fracbm, N);

	// Initialise FFT plans
	p1 = fftw_plan_dft_1d(2*N , correlation, circulant_eigenvalues, FFTW_FORWARD, FFTW_ESTIMATE);
	p2 = fftw_plan_dft_1d(2*N, rndW, fracGN, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	
	// Write correlation of noise
	write_correlation_exponents(correlation_exponents, N, invN, hurst);
	write_correlation(correlation, correlation_exponents, N);
	
	// FFT into circulant eigenvalues
	fftw_execute(p1); 
	
	for(iter = 0; iter < iteration; iter++)
	{
		
		generate_random_vector(randomComplexGaussian, rndW, circulant_eigenvalues, N);
		
		fftw_execute(p2);	
		
		// Reset first passage times
		first_passage_times = 0.0;
		// Integrate fractional Gaussain noise to fBM
		integrate_noise(fracbm, fracGN, lin_drift, frac_drift, N, hurst);
		
		// Find maximum to recursive depth RECURSION_DEPTH
		find_fpt(fracbm, &first_passage_times, passage_heights, N);
		

		// Output
		printf("%g\n",first_passage_times);	
	}// End iteration

	return 0;
}


