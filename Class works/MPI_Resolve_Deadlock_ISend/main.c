//
//  main.c
//  MPI_Resolve_Deadlock_ISend
//
//  Created by Serguei Diaz on 05.11.2024.
//

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 1024

int main(int argc, char **argv)
{
    int myrank, nprocs, len;
    char name[MPI_MAX_PROCESSOR_NAME];
    int *buf;
    MPI_Status st;
    MPI_Request send_request, recv_request;

    buf = (int *)malloc(sizeof(int) * (SIZE * 1024 + 100));

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Get_processor_name(name, &len);

    printf("Hello from processor %s[%d] %d of %d\n", name, len, myrank, nprocs);

    if (myrank % 2 == 0)
    {
        if (myrank < nprocs - 1)
        {
            int i, cl, sz = SIZE;
            double time;

            for (i = 0; i < SIZE * 1024; i++)
                buf[i] = i + 10;

            for (cl = 0; cl < 11; cl++)
            {
                time = MPI_Wtime();
                for (i = 0; i < 100; i++)
                {
                    // Use nonblocking send and receive
                    MPI_Isend(buf, sz, MPI_INT, myrank + 1, 10, MPI_COMM_WORLD, &send_request);
                    MPI_Irecv(buf, sz + 100, MPI_INT, myrank + 1, 20, MPI_COMM_WORLD, &recv_request);

                    // Wait for both send and receive to complete
                    MPI_Wait(&send_request, MPI_STATUS_IGNORE);
                    MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
                }
                time = MPI_Wtime() - time;
                printf("[%d] Time = %lf  Data=%9.0f KByte\n",
                       myrank,
                       time,
                       sz * sizeof(int) * 200.0 / 1024);
                printf("[%d]  Bandwidth[%d] = %lf MByte/sec\n",
                       myrank,
                       cl,
                       sz * sizeof(int) * 200 / (time * 1024 * 1024));
                sz *= 2;
            }
        }
        else
            printf("[%d] Idle\n", myrank);
    }
    else
    {
        int i, cl, sz = SIZE;

        for (i = 0; i < SIZE * 1024; i++)
            buf[i] = i + 100;

        for (cl = 0; cl < 11; cl++)
        {
            for (i = 0; i < 100; i++)
            {
                // Use nonblocking send and receive
                MPI_Isend(buf, sz, MPI_INT, myrank - 1, 20, MPI_COMM_WORLD, &send_request);
                MPI_Irecv(buf, sz + 100, MPI_INT, myrank - 1, 10, MPI_COMM_WORLD, &recv_request);

                // Wait for both send and receive to complete
                MPI_Wait(&send_request, MPI_STATUS_IGNORE);
                MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
            }
            sz *= 2;
        }
    }

    MPI_Finalize();
    printf("--------------\n");

    free(buf);
    return 0;
}

