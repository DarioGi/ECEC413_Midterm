/* Cholesky decomposition.
 * Compile as follows:
 * gcc -fopenmp -o chol chol.c chol_gold.c -lpthread -lm -std=c99
 */

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include "chol.h"
#include <pthread.h>
#include <semaphore.h>

typedef struct chol_thread
{
	Matrix *A;
	Matrix *U;
	int threadID;
	sem_t *counter_sem; /* Protects access to the counter. */
	sem_t *barrier_sem; /* Signals that barrier is safe to cross. */
	int *count;
	
}CHOL_THREAD;

void* thread_chol_pthread(void *arg);
void pthreads_barrier_sync(CHOL_THREAD *data, int index);



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
#define TRUE 1
#define FALSE !TRUE
int testing = TRUE;
int num_of_threads = 16;
////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) 
{	
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
		printf("Error. This program accepts no arguments. \n");
		exit(0);
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
	struct timeval start, stop;

	// compute the Cholesky decomposition on the CPU; single threaded version	
	printf("Performing Cholesky decomposition on the CPU using the single-threaded version. \n");
	gettimeofday(&start, NULL);
	int status = chol_gold(A, reference);
	gettimeofday(&stop, NULL);
	printf("CPU run time = %0.6f s.\n", (stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec) / (float)1000000));
	if(status == 0){
			  printf("Cholesky decomposition failed. The input matrix is not positive definite. \n");
			  exit(0);
	}
	//print_matrix(reference);
	
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
	printf("CPU run time = %0.6f s.\n", (stop.tv_sec - start.tv_sec + (stop.tv_usec - start.tv_usec) / (float)1000000));
	//print_matrix(U_pthreads);
	/* MODIFY THIS CODE: Perform the Cholesky decomposition using openmp. The resulting upper traingular matrix should be returned in U_openmp */
	//chol_using_openmp(A, U_openmp);


	// Check if the pthread and openmp results are equivalent to the expected solution
	if(check_chol(A, U_pthreads) == 0) 
			  printf("Error performing Cholesky decomposition using pthreads. \n");
	else
			  printf("Cholesky decomposition using pthreads was successful. \n");

	//if(check_chol(A, U_openmp) == 0) 
	//		  printf("Error performing Cholesky decomposition using openmp. \n");
	//else	
	//		  printf("Cholesky decomposition using openmp was successful. \n");



	// Free host matrices
	free(A.elements); 	
	free(U_pthreads.elements);	
	free(U_openmp.elements);
	free(reference.elements); 
	return 1;
}

/* Write code to perform Cholesky decopmposition using pthreads. */
void chol_using_pthreads(Matrix A, Matrix U)
{
	pthread_t *worker_thread;
	worker_thread = malloc(sizeof(pthread_t) * num_of_threads);
	CHOL_THREAD *data;
	data = malloc(sizeof(CHOL_THREAD) * num_of_threads);
	int size = NUM_ROWS/num_of_threads;
	int i;
	int count[3] = {0};
	sem_t counter_sem[3]; /* Protects access to the counter. */
	sem_t barrier_sem[3]; /* Signals that barrier is safe to cross. */
	for ( i = 0; i < 3; i++ )
	{
		sem_init(&counter_sem[i], 0, 1); /* Initialize the semaphore protecting the counter to unlocked. */
		sem_init(&barrier_sem[i], 0, 0); /* Initialize the semaphore protecting the barrier to locked. */
	}
	
	for ( i = 0; i < num_of_threads; i++ )
	{
		data[i].A = &A;
		data[i].U = &U;
		data[i].threadID = i;
		data[i].count = &(count[0]);
		data[i].counter_sem = counter_sem;
		data[i].barrier_sem = barrier_sem;
		if ( ( pthread_create(&worker_thread[i], NULL, thread_chol_pthread, (void *) &data[i]) ) != 0 ) 
		{
			printf("Cannot create thread.\n");
			exit(0);
		}
	}
	/* Wait for all the worker threads to finish. */
	for ( i = 0; i < num_of_threads; i++ )
		pthread_join(worker_thread[i], NULL);
	
	//print_matrix(U);
	free(worker_thread);
	free(data);
}

void* thread_chol_pthread(void *arg)
{
	CHOL_THREAD *data = (CHOL_THREAD*) arg;
	unsigned int i, j, k; 
	unsigned int size = data->U->num_rows * data->U->num_columns;
	Matrix *U = data->U;
	Matrix *A = data->A;
	int id = data->threadID;
	// Copy the contents of the A matrix into the working matrix U
	for ( i = 0+id; i < size; i+=num_of_threads )
		U->elements[i] = A->elements[i];
	pthreads_barrier_sync(data, 2);
	// Perform the Cholesky decomposition in place on the U matrix
	for ( k = 0; k < data->U->num_rows; k++ )
	{
		// Take the square root of the diagonal element
		if ( id == 0 )
		{
			U->elements[k * U->num_rows + k] = sqrt(U->elements[k * U->num_rows + k]);
			if ( U->elements[k * U->num_rows + k] <= 0 )
			{
				printf("Cholesky decomposition failed. \n");
			}
		}
		pthreads_barrier_sync(data, 0);
		// Division step
		for ( j = (k + 1 + id); j < U->num_rows; j+=num_of_threads )
			U->elements[k * U->num_rows + j] /= U->elements[k * U->num_rows + k]; // Division step
		
		pthreads_barrier_sync(data, 1);
		// Elimination step
		for ( i = (k + 1 + id); i < U->num_rows; i+=num_of_threads )
		{
			for ( j = i; j < U->num_rows; j++ )
				 U->elements[i * U->num_rows + j] -= U->elements[k * U->num_rows + i] * 
												 U->elements[k * U->num_rows + j];
		}
		pthreads_barrier_sync(data, 2);
	}
	// As the final step, zero out the lower triangular portion of U
	for ( i = 0+id; i < U->num_rows; i+=num_of_threads )
		for ( j = 0; j < i; j++ )
			U->elements[i * U->num_rows + j] = 0.0;
}

void pthreads_barrier_sync(CHOL_THREAD *data, int index)
{
	sem_wait(&(data->counter_sem[index]));
	
	if ( data->count[index] == (num_of_threads - 1) ) 
	{
		data->count[index] = 0;
		for ( int i = 0; i < (num_of_threads - 1); i++ )
			sem_post(&(data->barrier_sem[index]));
		sem_post(&(data->counter_sem[index]));
	}
	else 
	{
		data->count[index] += 1;
		sem_post(&(data->counter_sem[index]));
		sem_wait(&(data->barrier_sem[index])); 
	}
}

/* Write code to perform Cholesky decopmposition using openmp. */
void chol_using_openmp(const Matrix A, Matrix U)
{
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




