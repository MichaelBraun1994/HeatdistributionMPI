#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrix_handling.h"

#define ERROR_FUNCTION_CASE_1 1
#define ERROR_FUNCTION_CASE_2 2

void sendBorders(float* upperBorder, float* lowerBorder, int N, int rank, int numberOfWorkers)
{
    MPI_Request req; 
    if(rank == 1)
    {
        MPI_Isend(lowerBorder, N, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &req);
    }
    else if(rank == numberOfWorkers)
    {
        MPI_Isend(upperBorder, N, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &req);
    }
    else
    {
        MPI_Isend(lowerBorder, N, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &req);
        MPI_Isend(upperBorder, N, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &req);
    }
}

void recvBorders(float* upperBorder, float* lowerBorder, int N, int rank, int numberOfWorkers)
{
    if(rank == 1)
    {
        MPI_Recv(lowerBorder, N, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else if(rank == numberOfWorkers)
    {
        MPI_Recv(upperBorder, N, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
        MPI_Recv(lowerBorder, N, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(upperBorder, N, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

void distributeData(float *dataField, int N, int numberOfWorkers, int rowsPerWorker)
{
    int sendIndex = 1; 

    for(int i = 0; i < numberOfWorkers; ++i)
    {
        MPI_Send(&dataField[i * rowsPerWorker * N], (rowsPerWorker + 2) * N, MPI_FLOAT, sendIndex, 0, MPI_COMM_WORLD);
        sendIndex++;
    }
}

void collectData(float *dataField, int N, int numberOfWorkers, int rowsPerWorker)
{
    int recvIndex = 1;

    for(int i = 0; i < numberOfWorkers; ++i)
    {
        MPI_Recv(&dataField[(i + 1) * N], (rowsPerWorker * N), MPI_FLOAT, recvIndex, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        recvIndex++;
    }
}

int main(int argc, char **argv)
{
    int rank = 0, numberOfWorkers = 0;
    int N = 10 + 2;
    float residiumLimit = 0.0001;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfWorkers);

    int root = 0;
    numberOfWorkers -= 1;
    int rowsPerWorker = ((N - 2) / numberOfWorkers);

    if (rank == root)
    {
        /* prepare data */
        float *dataField = createField(N);
        initDataField(dataField, N, 1.0f, 0.0f, 0.0f, 1.0f);

        printf("Initial datafield.\n Rank %d, NodeCount: %d PartitionRowCount %d\n", rank, numberOfWorkers, rowsPerWorker);
        printDataField(dataField, N);

        /* scatter data */
        distributeData(dataField, N, numberOfWorkers, rowsPerWorker);

        /* process */
        float globalMaxResidium = 0.0f;
        unsigned char continueExecution = 0x1;
        {
            //while(globalMaxResidium > residiumLimit)
            for(int i = 0; i < 4; ++i)
            {
                MPI_Bcast(&continueExecution, 1, MPI_CHAR, root, MPI_COMM_WORLD);
                MPI_Reduce(MPI_IN_PLACE, &globalMaxResidium, 1, MPI_FLOAT, MPI_MAX, root, MPI_COMM_WORLD);
                printf("Maxresidium %f\n", globalMaxResidium);
            }
        }

        /* stop iterating */
        continueExecution = 0x0;
        MPI_Bcast(&continueExecution, 1, MPI_CHAR, root, MPI_COMM_WORLD);

        /* collectData */
        collectData(dataField, N, numberOfWorkers, rowsPerWorker);

        /* save results */
        printf("Data gathered:\n");
        printDataField(dataField, N);

        /* clean */
        free(dataField);
    }
    else
    {
        /* + 2 to store upper and lower border */
        int dataPartitionSize = (rowsPerWorker + 2) * N;
        float *dataPartition = (float *) malloc(sizeof(float) * dataPartitionSize);
        float *workPartition = (float *) malloc(sizeof(float) * dataPartitionSize);

        /* receive initial data */
        MPI_Recv(dataPartition, dataPartitionSize, MPI_FLOAT, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        /* iterate over dataPartition until conductor signals end of execution */
        unsigned char continueExecution = 1;
        while(continueExecution)
        {
            float localMaxResidium = 0.0f;
            //float localMaxResidium = iterate(&dataPartition, &workPartition, rowsPerWorker, N, ERROR_FUNCTION_CASE_1);

            /* send residium to conductor */
            MPI_Reduce(&localMaxResidium, NULL, 1, MPI_FLOAT, MPI_MAX, root, MPI_COMM_WORLD);

            /* check if another iteration is needed */
            MPI_Bcast(&continueExecution, 1, MPI_CHAR, root, MPI_COMM_WORLD);

            /* update neighbour nodes */
            if(continueExecution)
            {
                sendBorders(dataPartition, &dataPartition[dataPartitionSize - N], N, rank, numberOfWorkers);
                recvBorders(dataPartition, &dataPartition[dataPartitionSize -N], N, rank, numberOfWorkers);
            }
        }

        /* send result to conductor */
        MPI_Send(&dataPartition[N], rowsPerWorker * N, MPI_FLOAT, root, 0, MPI_COMM_WORLD);

        free(dataPartition);
        free(workPartition);
    }
    MPI_Finalize();
    return 0;
}
