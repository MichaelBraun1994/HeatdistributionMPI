#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_handling.h"

int main(int argc, char **argv)
{

    int rank = 0, workerCount = 0;
    int N = 10 + 2;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &workerCount);

    int root = workerCount - 1;

    /* dont user conductor as worker */
    --workerCount;

    int partitionRowCount = ((N - 2) / (workerCount - 1));

    if (rank == root)
    {
        /* prepare data */
        float *dataField = createField(N);
        initDataField(dataField, N, 1.0f, 1.0f, 0.0f, 0.0f);

        printf("Initial datafield.\n Rank %d, NodeCount: %d PartitionRowCount %d\n", rank, workerCount, partitionRowCount);
        fillInnerDataField(dataField, N, 4.2f);
        printDataField(dataField, N);

        /* scatter data */
        MPI_Scatter(&dataField[N], partitionRowCount * N, MPI_FLOAT, &dataField[N], (N - 2) * N, MPI_FLOAT, root, MPI_COMM_WORLD);

        /* process */
        for(int i = 0; i < (N - 2) * N; ++i)
        {
            dataField[N + i] = (float) 6.66f;
        }
        printf("Data wiped\n");
        printDataField(dataField, N);

        /* gather data */
        MPI_Gather(MPI_IN_PLACE, 1, MPI_FLOAT, &dataField[N], partitionRowCount * N, MPI_FLOAT, root, MPI_COMM_WORLD);

        /* save results */
        printf("Data gathered:\n");
        printDataField(dataField, N);

        /* clean */
        free(dataField);
    }
    else
    {
        /* + 2 to store upper and lower border */
        float *dataPartition = (float *) malloc(sizeof(float) * ((partitionRowCount + 2) * N));

        MPI_Scatter(NULL, 0, MPI_FLOAT, &dataPartition[N], (partitionRowCount * N), MPI_FLOAT, root, MPI_COMM_WORLD);
        MPI_Gather(&dataPartition[N], partitionRowCount * N, MPI_FLOAT, NULL, 0, MPI_FLOAT, root, MPI_COMM_WORLD);

        free(dataPartition);
    }
    MPI_Finalize();
    return 0;
}

