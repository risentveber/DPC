#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char*argv[]){
    int i, j, remain, block_size;
	int myrank, size;
    int* array;
	int n = 1000;
	int main_sum = 0;
	int main0_sum = 0;
	int sub_sum = 0;

    MPI_Status Status;
    MPI_Init (&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    remain = n%size;
    block_size = n/size;

    if ( myrank == 0){
        array = (int*)malloc(n * sizeof(int));

        for (i = 0; i < n; i++){
            array[i] = i;
			main0_sum += i;
		}
        for (i = 0; i < remain; i++)
            MPI_Send(array+i*(block_size+1), block_size + 1, MPI_INT, i, i, MPI_COMM_WORLD);
 
        for (i = remain; i < size; i++)
            MPI_Send(array + remain + i * block_size, block_size, MPI_INT, i, i, MPI_COMM_WORLD);
        free(array);
    }

    if (myrank < remain)
        block_size++;

    array = (int*)malloc(block_size * sizeof(int));
    MPI_Recv(array, block_size , MPI_INT, 0, myrank, MPI_COMM_WORLD, &Status);
    for (i = 0; i < block_size; i++)
        sub_sum += array[i];
    printf("Process #%d summed %d elements, result are: %d\n", myrank, block_size, sub_sum);

    MPI_Send(&sub_sum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

    if (myrank == 0){
        for (i = 0; i < size; i++){
            MPI_Recv(&j, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &Status);
            main_sum += j;
        }
        printf("Sequential sum: %d\n", main_sum);
        printf("Parallel sum: %d\n", main0_sum);
    }

    free(array);
    MPI_Finalize();
    return 0;
}
