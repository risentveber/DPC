//#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

const double ROD_LENGTH = 1.0;
const int ROD_PARTS = 10;

const double TIME_TOTAL = 0.1;
const int TIME_PARTS = 200;

const double TEMPERATURE_COEFF = 1.0;

const int SUM_LENGTH = 5;

void send_data(double* data, int size, int myrank, int all_processes)
{
	if (myrank > 0){
		MPI_Send(data + 1, 1, MPI_DOUBLE, myrank - 1, myrank - 1, MPI_COMM_WORLD);
	}
	if (myrank < all_processes - 1){
		MPI_Send(data + size, 1, MPI_DOUBLE, myrank + 1, myrank + 1, MPI_COMM_WORLD);
	}
}

void recv_data(double* data, int size, int myrank, int all_processes, MPI_Status* status)
{
	if (myrank < all_processes - 1){
		MPI_Recv(data + size + 1, 1, MPI_DOUBLE, myrank + 1, myrank, MPI_COMM_WORLD, status);
	}
	if (myrank > 0){
		MPI_Recv(data, 1, MPI_DOUBLE, myrank - 1, myrank, MPI_COMM_WORLD, status);
	}
}

void swap_array(double **a, double **b)
{
	double* temp = *a;
	*a = *b;
	*b = temp;
}

int main(int argc, char* argv[])
{	
	double dt = TIME_TOTAL / TIME_PARTS;
	double dl = ROD_LENGTH / ROD_PARTS;
	
	int myrank, all_processes;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &all_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	double mytime;
	if (myrank == 0){
		mytime = MPI_Wtime();
	}
	
	//Инициализация
	int block_size = (ROD_PARTS	 - 1) / all_processes;
	if (myrank < (ROD_PARTS - 1) % all_processes){
		++block_size;
	}
	double* temperature = (double*)malloc(sizeof(double) * (block_size + 2));
	double* old_temperature = (double*)malloc(sizeof(double) * (block_size + 2));
	
	for (int i = 0; i < block_size + 2; ++i){
		old_temperature[i] = 1.0;
	}
	if (myrank == 0){
		old_temperature[0] = 0.0;
	}
	if (myrank == all_processes - 1){
		old_temperature[block_size + 1] = 0.0;
	}
	//конец инициализации
	
	for (int t = 0; t < TIME_PARTS; t++){
		if (t != 0){
			if (myrank % 2 == 0){
				send_data(old_temperature, block_size, myrank, all_processes);
				recv_data(old_temperature, block_size, myrank, all_processes, MPI_STATUS_IGNORE);
			}
			else{
				recv_data(old_temperature, block_size, myrank, all_processes, MPI_STATUS_IGNORE);
				send_data(old_temperature, block_size, myrank, all_processes);
			}
		}
		for (int i = 1; i <= block_size; i++){
			temperature[i] = 
				old_temperature[i] + 
				dt * TEMPERATURE_COEFF * (old_temperature[i + 1] - 
				2 * old_temperature[i] + 
				old_temperature[i - 1]) / (dl * dl);
		}
		swap_array(&old_temperature, &temperature);
	}
	
	if (myrank != 0){//Отсылаем главному процессу результат
		MPI_Send(old_temperature + 1, block_size, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD);
	}
	else{
		double *results = (double*)malloc(sizeof(double) * (ROD_PARTS + 1));
		//Крайние точки сохраняют температуру
		results[0] = 0;
		results[ROD_PARTS] = 0;

		for (int i = 1; i <= block_size; i++)
			results[i] = old_temperature[i];

		int start = block_size + 1;
		for (int i = 1; i < all_processes; i++)		{
			int recv_size = (ROD_PARTS - 1) / all_processes;
			if (i < (ROD_PARTS - 1) % all_processes)
				recv_size++;
			MPI_Recv(results + start, recv_size, MPI_DOUBLE, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			start += recv_size;
		}
		
		/*//Вывод результата на экран
		for (int i = 0; i <= ROD_PARTS; i++)
			printf("%.6lf\n", results[i]);
		
		printf("-------------------------------------------\n");

		for (int i = 1; i < ROD_PARTS; i++){
			double x = i * dl;
			double sum = 0;
			for (int m = 0; m < SUM_LENGTH; m++){
				double a = exp(-TEMPERATURE_COEFF * M_PI * M_PI * (2 * m + 1) * (2 * m + 1) * TIME_TOTAL)
				 * sin(M_PI * (2 * m + 1) * x / ROD_LENGTH) / (2 * m + 1);
				sum = sum + 4 * 1 * a / M_PI;
			}
			double value = sum;
			printf("%.6f\n delta = %.6f\n", value, fabs(value - results[i]));
		}*/

		free(results);
	}

	if (myrank == 0){
		mytime = MPI_Wtime() - mytime;
		printf("%f\n", mytime);
	}
	
	free(temperature);
	free(old_temperature);
	
	MPI_Finalize(); 
	return 0;
}
