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


/// @file kalman_filters.h
/// @brief Kalman filter implementation
///
/// Signal filtering and fusion utilities using Kalman algorithm 
/// Based on http://www.bzarg.com/p/how-a-kalman-filter-works-in-pictures/#mjx-eqn-kalupdatefull

#ifndef KALMAN_FILTERS_H
#define KALMAN_FILTERS_H


/// Single Kalman filter instance internal data structure
typedef struct _KFilterData KFilterData;
/// Opaque reference to Kalman filter data structure
typedef KFilterData* KFilter;

                                                                            
/// @brief Creates and initializes internal matrices of Kalman filter data structure                                               
/// @param[in] dimensionsNumber size (in elements) of the internal state vector                                    
/// @return reference/pointer to allocated and initialized Kalman filter data structure
KFilter Kalman_CreateFilter( size_t dimensionsNumber );

/// @brief Deallocates internal data of given filter                              
/// @param[in] filter reference to filter
void Kalman_DiscardFilter( KFilter filter );

/// @brief Appends field to internal input vector
/// @param[in] filter reference to filter
/// @param[in] dimensionIndex index of the state vector variable related to the added input field 
void Kalman_AddInput( KFilter filter, size_t dimensionIndex );
                                                                  
/// @brief Updates specified value of the given filter input vector
/// @param[in] filter reference to filter
/// @param[in] inputIndex index of the internal input vector where value will be updated
/// @param[in] value new input value
void Kalman_SetInput( KFilter filter, size_t inputIndex, double value );

/// @brief Defines correlation between two state variables for prediction phase                             
/// @param[in] filter reference to filter
/// @param[in] outputIndex index (in state vector) of variable updated during predicition                                         
/// @param[in] inputIndex index (in state vector) of variable used to calculate predicition
/// @param[in] ratio output/input ratio desired on prediction
void Kalman_SetPredictionFactor( KFilter filter, size_t outputIndex, size_t inputIndex, double ratio );

/// @brief Sets maximum measure deviation for given input variable          
/// @param[in] filter reference to filter
/// @param[in] inputIndex index (in input vector)
/// @param[in] maxError maximum deviation for error modeling
void Kalman_SetInputMaxError( KFilter filter, size_t inputIndex, double maxError );

/// @brief Runs prediction phase on given Kalman filter                      
/// @param[in] filter reference to filter
/// @param[in] result pointer to array where predicted internal state will be copied (NULL if not required)
/// @return pointer to @a result array containing predicted filter state (NULL on errors)
double* Kalman_Predict( KFilter filter, double* result );

/// @brief Runs update phase on given Kalman filter                              
/// @param[in] filter reference to filter
/// @param[in] inputsList array of input values (size in elements equal to the number of added inputs)                            
/// @param[in] result pointer to array where updated internal state will be copied (NULL if not required)
/// @return pointer to @a result array containing updated filter state (NULL on errors)
double* Kalman_Update( KFilter filter, double* inputsList, double* result );
                                                                        
/// @brief Resets internal filter matrices data
/// @param[in] filter reference to filter
void Kalman_Reset( KFilter filter );


#endif  // KALMAN_FILTERS_H
