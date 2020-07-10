#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_handling.h"

void sendBorders(float *dataPartition, int partitionSize, int N, int rank, int numberOfWorkers)
{
    MPI_Request req;

    if(rank == 1)
    {
        MPI_Isend(&dataPartition[partitionSize * N], N, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &req);
    }
    else if(rank == (numberOfWorkers - 1))
    {
        MPI_Isend(&dataPartition[N], N, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &req);
    }
    else
    {
        MPI_Isend(&dataPartition[partitionSize * N], N, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &req);
        MPI_Isend(&dataPartition[N], N, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &req);
    }
}

void recvBorders(float *dataPartition, int partitionSize, int N, int rank, int numberOfWorkers, int initial)
{
    if(rank == 1)
    {
        if(initial)
        {
            MPI_Recv(dataPartition, N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        MPI_Recv(&dataPartition[(partitionSize + 1) * N], N, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else if(rank == (numberOfWorkers - 1))
    {
        MPI_Recv(dataPartition, N, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(initial)
        {
            MPI_Recv(&dataPartition[(partitionSize + 1) * N], N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    else
    {
        MPI_Recv(dataPartition, N, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&dataPartition[(partitionSize + 1) * N], N, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

int main(int argc, char **argv)
{
    int rank = 0, workerCount = 0;
    int N = 10 + 2;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &workerCount);

    int root = 0;
    int partitionRowCount = ((N - 2) / (workerCount - 1));

    if (rank == root)
    {
        /* prepare data */
        float *dataField = createField(N);
        initDataField(dataField, N, 1.0f, 0.0f, 0.0f, 1.0f);

        printf("Initial datafield.\n Rank %d, NodeCount: %d PartitionRowCount %d\n", rank, workerCount, partitionRowCount);
        printDataField(dataField, N);

        /* scatter data */
        MPI_Scatter(&dataField[N], partitionRowCount * N, MPI_FLOAT, &dataField[N], (N - 2) * N, MPI_FLOAT, root, MPI_COMM_WORLD);

        /* send upper and lower matrix border */
        MPI_Send(dataField, N, MPI_FLOAT, root + 1, 0, MPI_COMM_WORLD);
        MPI_Send(&dataField[(N-1) * N], N, MPI_FLOAT, workerCount - 1, 0, MPI_COMM_WORLD);

        /* process */


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
        float *workPartition = (float *) malloc(sizeof(float) * ((partitionRowCount + 2) * N));
        MPI_Scatter(NULL, 0, MPI_FLOAT, &dataPartition[N], (partitionRowCount * N), MPI_FLOAT, root, MPI_COMM_WORLD);

        sendBorders(dataPartition, partitionRowCount, N, rank, workerCount);
        recvBorders(dataPartition, partitionRowCount, N, rank, workerCount, 1);

        //while(0)
        {

            /* calculate new entries */
            float currentResidium = iterate(&dataPartition, &workPartition, N, partitionRowCount + 2, 1);

            sendBorders(dataPartition, partitionRowCount, N, rank, workerCount);
            recvBorders(dataPartition, partitionRowCount, N, rank, workerCount, 0);
        }

        printf("Rank %d\n", rank);
        printPartition(dataPartition, partitionRowCount + 2, N);
        MPI_Gather(&dataPartition[N], partitionRowCount * N, MPI_FLOAT, NULL, 0, MPI_FLOAT, root, MPI_COMM_WORLD);

        free(dataPartition);
    }
    MPI_Finalize();
    return 0;
}
