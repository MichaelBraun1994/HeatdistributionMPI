#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{

    int rank=0, size=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        char helloStr[] = "Hello World";
        MPI_Send(helloStr, strlen(helloStr) + 1, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
    }
    else if (rank == 1)
    {
        char helloStr[12];
        MPI_Recv(helloStr, 12, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Rank 1 received string %s from Rank 0\n", helloStr);
    }

    MPI_Finalize();
    return 0;
}

