//////////////////////////////////////////////////////////////////////////////////////
//                                                                                  //
//  Copyright (c) 2016-2017 Leonardo Consoni <consoni_2519@hotmail.com>             //
//                                                                                  //
//  This file is part of Simple Kalman Filter.                                      //
//                                                                                  //
//  Simple Kalman Filter is free software: you can redistribute it and/or modify    //
//  it under the terms of the GNU Lesser General Public License as published        //
//  by the Free Software Foundation, either version 3 of the License, or            //
//  (at your option) any later version.                                             //
//                                                                                  //
//  Simple Kalman Filter is distributed in the hope that it will be useful,         //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of                  //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                    //
//  GNU Lesser General Public License for more details.                             //
//                                                                                  //
//  You should have received a copy of the GNU Lesser General Public License        //
//  along with Simple Kalman Filter. If not, see <http://www.gnu.org/licenses/>.    //
//                                                                                  //
//////////////////////////////////////////////////////////////////////////////////////


#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "matrix/matrix.h"

#include "kalman_filters.h"

struct _KFilterData
{
  Matrix input;                                       // y
  Matrix state;                                       // x
  Matrix error;                                       // e
  Matrix inputModel;                                  // H
  Matrix gain;                                        // K
  Matrix prediction;                                  // F
  Matrix predictionCovariance;                        // P
  Matrix predictionCovarianceNoise;                   // Q
  Matrix errorCovariance;                             // S
  Matrix errorCovarianceNoise;                        // R
};


KFilter Kalman_CreateFilter( size_t dimensionsNumber )
{
  KFilter newFilter = (KFilter) malloc( sizeof(KFilterData) );
  memset( newFilter, 0, sizeof(KFilterData) );
  
  newFilter->state = Mat_Create( NULL, dimensionsNumber, 1 );
  newFilter->error = Mat_Create( NULL, dimensionsNumber, 1 );
  
  newFilter->gain = Mat_CreateSquare( dimensionsNumber, MATRIX_ZERO );
  
  newFilter->prediction = Mat_CreateSquare( dimensionsNumber, MATRIX_IDENTITY );
  newFilter->predictionCovariance = Mat_CreateSquare( dimensionsNumber, MATRIX_ZERO );
  //Mat_Scale( newFilter->predictionCovariance, 1000.0, newFilter->predictionCovariance );
  newFilter->predictionCovarianceNoise = Mat_CreateSquare( dimensionsNumber, MATRIX_IDENTITY );
  
  newFilter->errorCovariance = Mat_CreateSquare( dimensionsNumber, MATRIX_ZERO );
  newFilter->errorCovarianceNoise = Mat_CreateSquare( dimensionsNumber, MATRIX_IDENTITY );

  Kalman_Reset( newFilter );
  
  return newFilter;
}

void Kalman_DiscardFilter( KFilter filter )
{
  if( filter == NULL ) return;
    
  Mat_Discard( filter->input );
  Mat_Discard( filter->state );
  Mat_Discard( filter->error );
  
  Mat_Discard( filter->inputModel );
  Mat_Discard( filter->gain );
  Mat_Discard( filter->prediction );
  Mat_Discard( filter->predictionCovariance );
  Mat_Discard( filter->predictionCovarianceNoise );
  filter->errorCovariance = Mat_Resize( filter->errorCovariance, 0, 0 );
  // Mat_Discard( filter->errorCovariance ); // This causes a weird crash ("free(): invalid next size (fast)")
  Mat_Discard( filter->errorCovarianceNoise );
  
  free( filter );
}

void Kalman_AddInput( KFilter filter, size_t dimensionIndex )
{
  if( filter == NULL ) return;
  
  size_t dimensionsNumber = Mat_GetHeight( filter->state );
  
  if( dimensionIndex >= dimensionsNumber ) return;
  
  size_t newInputIndex = Mat_GetHeight( filter->input );
  size_t newInputsNumber = newInputIndex + 1;
  
  filter->input = Mat_Resize( filter->input, newInputsNumber, 1 );
  
  filter->inputModel = Mat_Resize( filter->inputModel, newInputsNumber, dimensionsNumber );
  for( size_t stateIndex = 0; stateIndex < dimensionsNumber; stateIndex++ )
    Mat_SetElement( filter->inputModel, newInputIndex, stateIndex, 0.0 );
  Mat_SetElement( filter->inputModel, newInputIndex, dimensionIndex, 1.0 );
  
  if( newInputsNumber > dimensionsNumber )
  {
    filter->errorCovariance = Mat_Resize( filter->errorCovarianceNoise, newInputsNumber, newInputsNumber );
    
    filter->errorCovarianceNoise = Mat_Resize( filter->errorCovarianceNoise, newInputsNumber, newInputsNumber );
    Mat_SetElement( filter->errorCovarianceNoise, newInputIndex, newInputIndex, 1.0 );

    filter->gain = Mat_Resize( filter->gain, newInputsNumber, newInputsNumber );
    filter->error = Mat_Resize( filter->error, newInputsNumber, 1 );
  }
}

void Kalman_SetInput( KFilter filter, size_t inputIndex, double value )
{
  if( filter == NULL ) return;
  
  Mat_SetElement( filter->input, inputIndex, 0, value );
}

void Kalman_SetPredictionFactor( KFilter filter, size_t outputIndex, size_t inputIndex, double ratio )
{
  if( filter == NULL ) return;
  
  Mat_SetElement( filter->prediction, outputIndex, inputIndex, ratio );
}

void Kalman_SetInputMaxError( KFilter filter, size_t inputIndex, double maxError )
{
  if( filter == NULL ) return;
  
  Mat_SetElement( filter->errorCovarianceNoise, inputIndex, inputIndex, maxError * maxError );
}

double* Kalman_Predict( KFilter filter, double* result )
{
  if( filter == NULL ) return NULL;
  
  // x = F*x
  Mat_Dot( filter->prediction, MATRIX_KEEP, filter->state, MATRIX_KEEP, filter->state );                                       // F[nxn] * x[nx1] -> x[nx1]
  
  // P = F*P*F' + Q
  Mat_Dot( filter->prediction, MATRIX_KEEP, filter->predictionCovariance, MATRIX_KEEP, filter->predictionCovariance );         // F[nxn] * P[nxn] -> P[nxn]
  Mat_Dot( filter->predictionCovariance, MATRIX_KEEP, filter->prediction, MATRIX_TRANSPOSE, filter->predictionCovariance );    // P[nxn] * F'[nxn] -> P[nxn]
  Mat_Sum( filter->predictionCovariance, 1.0, filter->predictionCovarianceNoise, 1.0, filter->predictionCovariance );          // P[nxn] + Q[nxn] -> P[nxn]
  
  if( result == NULL ) return NULL;
  
  return Mat_GetData( filter->state, result );
}

double* Kalman_Update( KFilter filter, double* inputsList, double* result )
{
  if( filter == NULL ) return NULL;
  
  if( inputsList != NULL ) Mat_SetData( filter->input, inputsList );
  
  // e = y - H*x
  Mat_Dot( filter->inputModel, MATRIX_KEEP, filter->state, MATRIX_KEEP, filter->error );                             // H[mxn] * x[nx1] -> e[mx1]
  Mat_Sum( filter->input, 1.0, filter->error, -1.0, filter->error );                                                 // y[mx1] - e[mx1] -> e[mx1]
  
  // S = H*P*H' + R
  Mat_Dot( filter->inputModel, MATRIX_KEEP, filter->predictionCovariance, MATRIX_KEEP, filter->errorCovariance );    // H[mxn] * P[nxn] -> S[mxn]
  Mat_Dot( filter->errorCovariance, MATRIX_KEEP, filter->inputModel, MATRIX_TRANSPOSE, filter->errorCovariance );    // S[mxn] * H'[nxm] -> S[mxm]
  Mat_Sum( filter->errorCovariance, 1.0, filter->errorCovarianceNoise, 1.0, filter->errorCovariance );               // S[mxm] + R[mxm] -> S[mxm]
  
  // K = P*H' * S^(-1)
  Mat_Dot( filter->predictionCovariance, MATRIX_KEEP, filter->inputModel, MATRIX_TRANSPOSE, filter->gain );          // P[nxn] * H'[nxm] -> K[nxm]
  if( Mat_Inverse( filter->errorCovariance, filter->errorCovariance ) != NULL )                                      // S^(-1)[mxm] -> S[mxm]
  {
    Mat_Dot( filter->gain, MATRIX_KEEP, filter->errorCovariance, MATRIX_KEEP, filter->gain );                          // K[nxm] * S[mxm] -> K[nxm]
    
    // x = x + K*e
    Mat_Dot( filter->gain, MATRIX_KEEP, filter->error, MATRIX_KEEP, filter->error );                                   // K[nxm] * e[mx1] -> e[nx1]
    Mat_Sum( filter->state, 1.0, filter->error, 1.0, filter->state );                                                  // x[nx1] + e[nx1] -> x[nx1]
    
    // P' = P - K*H*P
    Mat_Dot( filter->gain, MATRIX_KEEP, filter->inputModel, MATRIX_KEEP, filter->gain );                               // K[nxm] * H[mxn] -> K[nxn]
    Mat_Dot( filter->gain, MATRIX_KEEP, filter->predictionCovariance, MATRIX_KEEP, filter->gain );                     // K[nxn] * P[nxn] -> K[nxn]
    Mat_Sum( filter->predictionCovariance, 1.0, filter->gain, -1.0, filter->predictionCovariance );                    // P[nxn] - K[nxn] -> P[nxn]
  }

  if( result == NULL ) return NULL;
  
  return Mat_GetData( filter->state, result );
}

void Kalman_Reset( KFilter filter )
{
  if( filter == NULL ) return;
  
  Mat_Clear( filter->input );
  Mat_Clear( filter->state );
  Mat_Clear( filter->error );
  
  Mat_Clear( filter->gain );
  Mat_Clear( filter->predictionCovariance );
  Mat_Clear( filter->errorCovariance );
}
