//
//  main.c
//  MPI_Transpose_Matrix
//
//  Created by Serguei Diaz on 05.11.2024.
//

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define N 10  // Define the matrix size (N x N)

// Function to initialize the matrix with sample values
void initialize_matrix(int matrix[N][N]) {
    int value = 1;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = value++;
        }
    }
}

// Function to print the matrix
void print_matrix(int matrix[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%3d ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char** argv) {
    int world_rank, world_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int matrix[N][N];
    int transposed[N][N];

    if (world_rank == 0) {
        // Initialize the original matrix with values
        initialize_matrix(matrix);
        printf("Original matrix:\n");
        print_matrix(matrix);
    }

    // Define a custom MPI datatype for a row of the matrix
    MPI_Datatype row_type;
    MPI_Type_contiguous(N, MPI_INT, &row_type);
    MPI_Type_commit(&row_type);

    // Determine the number of columns per process
    int columns_per_proc = N / world_size;
    int remainder = N % world_size;

    int sendcounts[world_size];
    int displacements[world_size];
    int offset = 0;

    for (int i = 0; i < world_size; i++) {
        sendcounts[i] = columns_per_proc + (i < remainder ? 1 : 0);
        displacements[i] = offset;
        offset += sendcounts[i];
    }

    // Allocate space for each process to receive columns
    int local_cols = sendcounts[world_rank];
    int* local_data = (int*)malloc(N * local_cols * sizeof(int));

    // Scatter columns to each process based on custom sendcounts and displacements
    if (world_rank == 0) {
        // Pack the columns to send to each process
        int column_data[N * N];
        int pos = 0;
        for (int col = 0; col < N; col++) {
            for (int row = 0; row < N; row++) {
                column_data[pos++] = matrix[row][col];
            }
        }
        MPI_Scatterv(column_data, sendcounts, displacements, row_type,
                     local_data, N * local_cols, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
        MPI_Scatterv(NULL, sendcounts, displacements, row_type,
                     local_data, N * local_cols, MPI_INT, 0, MPI_COMM_WORLD);
    }

    // Transpose: each process rearranges its columns into rows
    int* transposed_data = (int*)malloc(N * local_cols * sizeof(int));
    for (int i = 0; i < local_cols; i++) {
        for (int j = 0; j < N; j++) {
            transposed_data[i * N + j] = local_data[j * local_cols + i];
        }
    }

    // Gather the transposed rows back into the transposed matrix on the root process
    MPI_Gatherv(transposed_data, N * local_cols, MPI_INT, &transposed[0][0],
                sendcounts, displacements, row_type, 0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        printf("Transposed matrix:\n");
        print_matrix(transposed);
    }

    // Clean up
    MPI_Type_free(&row_type);
    free(local_data);
    free(transposed_data);

    MPI_Finalize();
    return 0;
}
