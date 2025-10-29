//
//  main.c
//  MPI_Map
//
//  Created by Serguei Diaz on 04.11.2024.
//

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

// Define the function pointer type for the mapping function
typedef int (*MapFunction)(int);

// Sample map function: squares the input
int square(int x) {
    return x * x;
}

// Another sample map function: doubles the input
int double_value(int x) {
    return x * 2;
}

// Distributed map function that applies `func` to each element in `data`
void distributed_map(int* data, int data_size, MapFunction func) {
    int world_rank, world_size;

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int chunk_size = data_size / world_size;
    int* sub_data = (int*) malloc(chunk_size * sizeof(int));

    // Scatter data from root to all processes
    MPI_Scatter(data, chunk_size, MPI_INT, sub_data, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Apply the map function to each element of the subarray
    for (int i = 0; i < chunk_size; i++) {
        sub_data[i] = func(sub_data[i]);
    }

    // Gather results back to root process
    MPI_Gather(sub_data, chunk_size, MPI_INT, data, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Free allocated memory for subarray
    free(sub_data);
}

int main(int argc, char** argv) {
    int world_rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int data_size = 100; // Total size of data to process
    int* data = NULL;

    if (world_rank == 0) {
        // Initialize data array only on root process
        data = (int*) malloc(data_size * sizeof(int));
        for (int i = 0; i < data_size; i++) {
            data[i] = i + 1; // Fill with some data, e.g., 1, 2, 3, ..., 100
        }
        
        printf("Original data:\n");
        for (int i = 0; i < data_size; i++) {
            printf("%d ", data[i]);
        }
        printf("\n");
    }

    // Call the distributed map function with the chosen map operation
    distributed_map(data, data_size, square); // Can replace `square` with `double_value`

    // Root process prints the result
    if (world_rank == 0) {
        printf("Mapped data:\n");
        for (int i = 0; i < data_size; i++) {
            printf("%d ", data[i]);
        }
        printf("\n");

        // Free allocated memory on root
        free(data);
    }

    MPI_Finalize();
    return 0;
}

