#ifndef __MATRIX_HANDLING
#define __MATRIX_HANDLING

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
float getArray2DAccess(float* arr, int N, int x, int y);
void setArray2DAccess(float* arr, int N, int x, int y, float value);
float *createField(int N);
float getDelta(int N, float begin, float end);
void setDataFieldCorners(float* dataField, int N, float leftUp, float rightUp, float leftDown, float rightDown);
void setLeftBorder(float* dataField, int N);
void setUpperBorder(float* dataField, int N);
void setRightBorder(float* dataField, int N);
void setLowerBorder(float* dataField, int N);
void setDataFieldBorders(float* dataField, int N);
void fillInnerDataField(float* dataField, int N, float value);
void initDataField(float* dataField, int N, float leftUp, float rightUp, float leftDown, float rightDown);
float errorFunction_CaseA(int x, int y);
float errorFunction_CaseB(int x, int y);
float calculateResidium(float* dataField, int N, int x, int y, errorFunctionType errorFunction);
float updateValue(float* dataField, int N, int x, int y, float residium);
void swapFields(float **first, float **second);
float iterate(float** dataField, float** workField, int N, int errorFunctionCase);
void printPartition(float *dataPartition, int rows, int N);
void printDataField(float* dataField, int N);
void saveCSV(char fileName[], float* data, int N, float residiumNorm, int iteration);

#endif
