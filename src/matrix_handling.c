
/*************************************************************/
/* INCLUDES
/*************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#define _USE_MATH_DEFINES
#include <math.h>

/*************************************************************/
/* DEFINES 
/*************************************************************/
#define __USE_STATIC_CFG 0 

#define UPPER_LEFT_CORNER (0)
#define UPPER_RIGHT_CORNER (N - 1)
#define LOWER_LEFT_CORNER ((N - 1) * N)
#define LOWER_RIGHT_CORNER ((N * N) - 1)

/*************************************************************/
/* TYPES 
/*************************************************************/
typedef float (*errorFunctionType)(int x, int y);

/*************************************************************/
/* FUNCTIONS 
/*************************************************************/
float getArray2DAccess(float* arr, int N, int x, int y)
{
    return arr[y * N + x];
}

void setArray2DAccess(float* arr, int N, int x, int y, float value)
{
    arr[y * N + x] = value;
}

float *createField(int N)
{
    return malloc(sizeof(float) * N * N);
 }

float getDelta(int N, float begin, float end)
 {
     return (end - begin) / (N - 1);
 }

 void setDataFieldCorners(float* dataField, int N, float leftUp, float rightUp, float leftDown, float rightDown)
 {
     dataField[UPPER_LEFT_CORNER]  = leftUp;
     dataField[UPPER_RIGHT_CORNER] = rightUp;
     dataField[LOWER_LEFT_CORNER]  = leftDown;
     dataField[LOWER_RIGHT_CORNER] = rightDown;
 }

 void setLeftBorder(float* dataField, int N)
 {
    float begin = dataField[UPPER_LEFT_CORNER];
    float end = dataField[LOWER_LEFT_CORNER];

    float delta = getDelta(N, begin, end);

    for(int i = 1; i < (N - 1); ++i)
    {
        dataField[i * N] = begin + i * delta;
    }
 }
 void setUpperBorder(float* dataField, int N)
 {
    float begin = dataField[UPPER_LEFT_CORNER];
    float end = dataField[UPPER_RIGHT_CORNER];

    float delta = getDelta(N, begin, end);

    for(int i = 1; i < (N - 1); ++i)
    {
        dataField[i] = begin + i * delta;
    }
 }
 void setRightBorder(float* dataField, int N)
 {
    float begin = dataField[UPPER_RIGHT_CORNER];
    float end = dataField[LOWER_RIGHT_CORNER];

    float delta = getDelta(N, begin, end);

    for(int i = 1; i < (N - 1); ++i)
    {
        dataField[(N - 1) + (i * N)] = begin + i * delta;
    }
 }
 void setLowerBorder(float* dataField, int N)
 {
    float begin = dataField[LOWER_LEFT_CORNER];
    float end = dataField[LOWER_RIGHT_CORNER];

    float delta = getDelta(N, begin, end);

    for(int i = 1; i < (N - 1); ++i)
    {
        dataField[(N - 1) * N + i] = begin + i * delta;
    }
 }

 void setDataFieldBorders(float* dataField, int N)
 {
     setUpperBorder(dataField, N);
     setLeftBorder(dataField,  N);
     setRightBorder(dataField, N);
     setLowerBorder(dataField, N);
 }

 void fillInnerDataField(float* dataField, int N, float value)
 {
    // ((N-2) * (N-2)) equals number of inner points
    for(int i = 0; i < ((N-2) * (N-2)); ++i)
    {
        int x = i / (N - 2);
        int y = i % (N - 2);

        dataField[(y + 1) * N + (x + 1)] = value;
    }
 }

 void initDataField(float* dataField, int N, float leftUp, float rightUp, float leftDown, float rightDown)
 {
     setDataFieldCorners(dataField, N, leftUp, rightUp, leftDown, rightDown);
     setDataFieldBorders(dataField, N);
     fillInnerDataField(dataField, N, 0);
 }

 float errorFunction_CaseA(int x, int y)
 {
     return 0;
 }
 float errorFunction_CaseB(int x, int y)
 {
     return (2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y)) * 1.0/(2*2);
 }

 float calculateResidium(float* dataField, int N, int x, int y, errorFunctionType errorFunction)
 {
     return errorFunction(x, y) -4 * getArray2DAccess(dataField, N, x, y) 
        + getArray2DAccess(dataField, N, x + 1, y)
        + getArray2DAccess(dataField, N, x - 1, y)
        + getArray2DAccess(dataField, N, x, y + 1)
        + getArray2DAccess(dataField, N, x, y - 1);
 }

 float updateValue(float* dataField, int N, int x, int y, float residium)
 {
     return getArray2DAccess(dataField, N, x, y) + (0.25f * residium);
 }

 void swapFields(float **first, float **second)
 {
     float *tmp = *first;
     *first = *second;
     *second = tmp;
 }

 float iterate(float** dataField, float** workField, int N, int M, int errorFunctionCase)
 {
     float maxResidium = 0.0f;

    for(int y = 1; y < (M - 1); ++y)
    {
        for(int x = 1; x < (N - 1); ++x)
        {
            float residium;
            if(errorFunctionCase == 1)
            {
                residium = calculateResidium(*dataField, N, x, y, errorFunction_CaseA);
            }
            else if(errorFunctionCase == 2)
            {
                residium = calculateResidium(*dataField, N, x, y, errorFunction_CaseB);
            }
            float updateVal = updateValue(*dataField, N, x, y, residium);
            setArray2DAccess(*workField, N, x, y, updateVal);

            if(residium > maxResidium)
            {
                maxResidium = residium;
            }
        }
    }
    swapFields(dataField, workField);
    return maxResidium;
 }

void printPartition(float *dataPartition, int rows, int N)
{
     for(int i = 0; i < (rows * N); ++i)
    {
        printf("%4.3f, ", dataPartition[i]);
        
        if((i % N) == (N - 1))
        {
            printf("\n");
        }
    }
}

 void printDataField(float* dataField, int N)
 {
     printPartition(dataField, N, N);
 }

 void saveCSV(char fileName[], float* data, int N, float residiumNorm, int iteration)
 {
     FILE *fileHandle = fopen(fileName, "w+");

     fprintf(fileHandle, "# Dimension %d, ResidiumNorm %f, Iterations %d\n", N, residiumNorm, iteration);

     for(int i = 0; i < (N * N); ++i)
     {
        fprintf(fileHandle, "%.5f", data[i]);
        if((i % N) == (N - 1))
        {
            fprintf(fileHandle, "\n");
        }
        else
        {
            fprintf(fileHandle, ",");
        }
     }

     fclose(fileHandle);
 }
