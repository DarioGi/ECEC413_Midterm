// Modified by: Darius Remeika (dr542) and Hardip Sahota (hs452)
// Modified on: 11/04/2015

/*  Purpose: Calculate definite integral using trapezoidal rule.
 *
 * Input:   a, b, n
 * Output:  Estimate of integral from a to b of f(x)
 *          using n trapezoids.
 *
 * Compile: gcc -o trap trap.c -lpthread -lm
 * Usage:   ./trap
 *
 * Note:    The function f(x) is hardwired.
 *
 */

#ifdef _WIN32
#define NOMINMAX 
#endif

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>
#include <pthread.h>

#define LEFT_ENDPOINT 5
#define RIGHT_ENDPOINT 1000
#define NUM_TRAPEZOIDS 100000000
#define NUM_THREADS 4

typedef struct trap_thread
{
	double integral;
	int kStart;
	int size;
	float a;
	float h;
	
}TRAP_THREAD;

double compute_using_pthreads(float, float, int, float);
double compute_gold(float, float, int, float);
void* thread_trap_pthread(void *arg);

int main(void) 
{
	int n = NUM_TRAPEZOIDS;
	float a = LEFT_ENDPOINT;
	float b = RIGHT_ENDPOINT;
	float h = (b-a)/(float)n; // Height of each trapezoid
	struct timeval start, stop;
	
	printf("The height of the trapezoid is %f \n", h);

	gettimeofday(&start, NULL);
	double reference = compute_gold(a, b, n, h);
	gettimeofday(&stop, NULL);
	printf("CPU run time = %0.6f s.\n", (stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec) / (float)1000000));
	printf("Reference solution computed on the CPU = %f \n", reference);

	/* Write this function to complete the trapezoidal on the GPU. */
	gettimeofday(&start, NULL);
	double pthread_result = compute_using_pthreads(a, b, n, h);
	gettimeofday(&stop, NULL);
	printf("CPU run time = %0.6f s.\n", (stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec) / (float)1000000));
	printf("Solution computed using pthreads = %f \n", pthread_result);
} 


/*------------------------------------------------------------------
 * Function:    f
 * Purpose:     Compute value of function to be integrated
 * Input args:  x
 * Output: (x+1)/sqrt(x*x + x + 1)

 */
float f(float x) {
		  return (x + 1)/sqrt(x*x + x + 1);
}  /* f */

/*------------------------------------------------------------------
 * Function:    Trap
 * Purpose:     Estimate integral from a to b of f using trap rule and
 *              n trapezoids
 * Input args:  a, b, n, h
 * Return val:  Estimate of the integral 
 */
double compute_gold(float a, float b, int n, float h) {
   double integral;
   int k;

   integral = (f(a) + f(b))/2.0;
   for (k = 1; k <= n-1; k++) {
     integral += f(a+k*h);
   }
   integral = integral*h;

   return integral;
} 


void* thread_trap_pthread(void *arg)
{
	TRAP_THREAD *data = (TRAP_THREAD*) arg;
	double temp = 0.0;	
	int k;
	for ( k = data->kStart; k < data->size; k++ ) 
		temp += f(data->a+k*data->h);
	data->integral = temp;
} 


double compute_using_pthreads(float a, float b, int n, float h)
{
	pthread_t worker_thread[NUM_THREADS];
	TRAP_THREAD data[NUM_THREADS];
	int size = (n)/NUM_THREADS;
	double integral = (f(a) + f(b))/2.0;
	int i;
	int tempStart = 1;
	for ( i = 0; i < NUM_THREADS; i++ )
	{
		data[i].integral = 0;
		data[i].kStart = data[i].size = tempStart;
		if ( i != NUM_THREADS - 1 )
			data[i].size += size;
		else
			data[i].size += n - tempStart;
		data[i].a = a;
		data[i].h = h;
		if ( ( pthread_create(&worker_thread[i], NULL, thread_trap_pthread, (void *) &data[i]) ) != 0 ) 
		{
			printf("Cannot create thread.\n");
			exit(0);
		}
		tempStart += size;
	}
	
	/* Wait for all the worker threads to finish. */
	for ( i = 0; i < NUM_THREADS; i++ )
		pthread_join(worker_thread[i], NULL);
	
	for ( i = 0; i < NUM_THREADS; i++ )
		integral += (data[i].integral);
	
	return integral * h;
}




