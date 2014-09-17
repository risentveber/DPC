#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#define X_NUM 11

int main(char argc, char* argv[])
{
    MPI_Status Status;
	int size, myrank;
    MPI_Init (&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	int t_num;
	double* array[X_NUM];
	for(int i = 0; i < X_NUM; i++)
		array[i] = (double*)malloc(t_num*sizeof(double));
	
	
	for(int i = 0; i < X_NUM; i++)
		free(array[i]);
	
    MPI_Finalize();
	return 0;
}
