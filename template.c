#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>

int main(char argc, char* argv[])
{
    MPI_Status Status;
	int size, myrank;
    MPI_Init (&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
    MPI_Finalize();
	return 0;
}
