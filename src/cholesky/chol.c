/* Cholesky decomposition.
 *
 * Modified by: Darius Remeika (dr542) and Hardip Sahota (hs452)
 * Modified on: 11/04/2015
 *
 * Compile as follows: 
 * gcc -o chol chol.c chol_gold.c -lpthread -fopenmp -lm -std=c99 -O3
 */

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "chol.h"
#include <pthread.h>
#include <semaphore.h>
#include <sys/time.h>


////////////////////////////////////////////////////////////////////////////////
// declarations, forward

Matrix allocate_matrix(int num_rows, int num_columns, int init);
int perform_simple_check(const Matrix M);
void print_matrix(const Matrix M);
extern Matrix create_positive_definite_matrix(unsigned int, unsigned int);
extern int chol_gold(const Matrix, Matrix);
extern int check_chol(const Matrix, const Matrix);
void chol_using_pthreads(const Matrix, Matrix);
void chol_using_openmp(const Matrix, Matrix);

//#define NUM_THREADS 16
#define TRUE 1
#define FALSE !TRUE
int testing = FALSE;
int num_of_threads = 16;

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) 
{
    struct timeval start, stop;
    float timer;
    
    // Check command line arguments
    if ( testing )
    {
        if ( argc == 2 )
        {
            num_of_threads = atoi(argv[1]);
        }
        else
        {
            printf("Usage: Chol <num of threads>\n");
            exit(-1);
        }
    }
    else
    {
        if( argc > 1){
            printf("Error. This program accepts no arguments. \n");
            exit(0);
        }
    }
	 
	// Matrices for the program
	Matrix A; // The N x N input matrix
	Matrix reference; // The upper triangular matrix computed by the CPU
	Matrix U_pthreads; // The upper triangular matrix computed by the pthread implementation
	Matrix U_openmp; // The upper triangular matrix computed by the openmp implementation 
	
	// Initialize the random number generator with a seed value 
	srand(time(NULL));

	// Create the positive definite matrix. May require a few tries if we are unlucky
	int success = 0;
	while(!success){
		A = create_positive_definite_matrix(MATRIX_SIZE, MATRIX_SIZE);
		if(A.elements != NULL)
				  success = 1;
	}
	// print_matrix(A);
	// getchar();


	reference  = allocate_matrix(MATRIX_SIZE, MATRIX_SIZE, 0); // Create a matrix to store the CPU result
	U_pthreads =  allocate_matrix(MATRIX_SIZE, MATRIX_SIZE, 0); // Create a matrix to store the pthread result
	U_openmp =  allocate_matrix(MATRIX_SIZE, MATRIX_SIZE, 0); // Create a matrix to store the openmp result

    
	// compute the Cholesky decomposition on the CPU; single threaded version	
	printf("Performing Cholesky decomposition on the CPU using the single-threaded version. \n");
    gettimeofday(&start, NULL);
	int status = chol_gold(A, reference);
    gettimeofday(&stop, NULL);
    timer = (float)(stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec) / (float)1000000);
    printf("CPU run time for Serial = %0.6f s.\n", timer);
	if(status == 0){
			  printf("Cholesky decomposition failed. The input matrix is not positive definite. \n");
			  exit(0);
	}
	
	
	printf("Double checking for correctness by recovering the original matrix. \n");
	if(check_chol(A, reference) == 0){
		printf("Error performing Cholesky decomposition on the CPU. Try again. Exiting. \n");
		exit(0);
	}
	printf("Cholesky decomposition on the CPU was successful. \n");

	/* MODIFY THIS CODE: Perform the Cholesky decomposition using pthreads. The resulting upper triangular matrix should be returned in 
	 U_pthreads */
    gettimeofday(&start, NULL);
	chol_using_pthreads(A, U_pthreads);
    gettimeofday(&stop, NULL);
    timer = (float)(stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec) / (float)1000000);
    printf("CPU run time for PThreads = %0.6f s.\n", timer);

	/* MODIFY THIS CODE: Perform the Cholesky decomposition using openmp. The resulting upper traingular matrix should be returned in U_openmp */
    gettimeofday(&start, NULL);
	chol_using_openmp(A, U_openmp);
    gettimeofday(&stop, NULL);
    timer = (float)(stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec) / (float)1000000);
    printf("CPU run time for OpenMP = %0.6f s.\n", timer);



	// Check if the pthread and openmp results are equivalent to the expected solution
	if(check_chol(A, U_pthreads) == 0) 
			  printf("Error performing Cholesky decomposition using pthreads. \n");
	else
			  printf("Cholesky decomposition using pthreads was successful. \n");

	if(check_chol(A, U_openmp) == 0) 
			  printf("Error performing Cholesky decomposition using openmp. \n");
	else	
			  printf("Cholesky decomposition using openmp was successful. \n");



	// Free host matrices
	free(A.elements); 	
	free(U_pthreads.elements);	
	free(U_openmp.elements);
	free(reference.elements); 
	return 1;
}

/* Write code to perform Cholesky decopmposition using pthreads. */
void chol_using_pthreads(const Matrix A, Matrix U)
{
}

/* Write code to perform Cholesky decopmposition using openmp. */
void chol_using_openmp(const Matrix A, Matrix U)
{
    unsigned int i, j, k;
    unsigned int num_rows = U.num_rows;
    unsigned int size = A.num_rows * A.num_columns;
    float *u = U.elements, *a = A.elements;
    
#pragma omp parallel private(i, j, k) shared(U, num_rows) num_threads(num_of_threads)
{
    // Copy the contents of the A matrix into the working matrix U
#pragma omp for
    for (i = 0; i < size; i ++)
        u[i] = a[i];

    // Perform the Cholesky decomposition in place on the U matrix
    for(k = 0; k < num_rows; k++){
        
        // Take the square root of the diagonal element
#pragma omp master 
{
        u[k * num_rows + k] = sqrt(u[k * num_rows + k]);
    
        if(u[k * num_rows + k] <= 0){
            printf("Cholesky decomposition failed. \n");
            exit(0);
        }
}
        // Division step
#pragma omp for
        for(j = (k + 1); j < num_rows; j++)
            u[k * num_rows + j] /= u[k * num_rows + k]; // Division step
        
        // Elimination step
#pragma omp for
        for(i = (k + 1); i < num_rows; i++)
            for(j = i; j < num_rows; j++)
                u[i * num_rows + j] -= u[k * num_rows + i] * u[k * num_rows + j];
        
    }
    
    // As the final step, zero out the lower triangular portion of U
#pragma omp for
    for(i = 0; i < num_rows; i++)
        for(j = 0; j < i; j++)
            u[i * num_rows + j] = 0.0;
}
}


// Allocate a matrix of dimensions height*width
//	If init == 0, initialize to all zeroes.  
//	If init == 1, perform random initialization.
Matrix allocate_matrix(int num_rows, int num_columns, int init)
{
    	Matrix M;
    	M.num_columns = M.pitch = num_columns;
    	M.num_rows = num_rows;
    	int size = M.num_rows * M.num_columns;
		
	M.elements = (float *) malloc(size * sizeof(float));
	for(unsigned int i = 0; i < size; i++){
		if(init == 0) M.elements[i] = 0; 
		else
			M.elements[i] = (float)rand()/(float)RAND_MAX;
	}
    return M;
}	




