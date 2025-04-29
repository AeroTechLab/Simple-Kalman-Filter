//////////////////////////////////////////////////////////////////////////////////////
//                                                                                  //
//  Copyright (c) 2016-2025 Leonardo Consoni <leonardojc@protonmail.com>            //
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
/// @param[in] statesNumber size (in elements) of the internal estimated state vector    
/// @param[in] measurementsNumber size (in elements) of the measurements vector 
/// @param[in] inputsNumber size (in elements) of the inputs vector 
/// @return reference/pointer to allocated and initialized Kalman filter data structure
KFilter Kalman_CreateFilter( size_t statesNumber, size_t measurementsNumber, size_t inputsNumber );

/// @brief Deallocates internal data of given filter                              
/// @param[in] filter reference to filter
void Kalman_DiscardFilter( KFilter filter );
                                                                  
/// @brief Defines correlation between input and state variables for prediction phase 
/// @param[in] filter reference to filter
/// @param[in] stateIndex index of the correspondent state variable in the internal state vector
/// @param[in] inputIndex index of the input variable in the internal input vector
/// @param[in] ratio output/input ratio desired on prediction
void Kalman_SetInputFactor( KFilter filter, size_t stateIndex, size_t inputIndex, double ratio );

/// @brief Defines correlation between two state variables for state transition on prediction phase                             
/// @param[in] filter reference to filter
/// @param[in] newStateIndex index (in state vector) of variable updated during predicition 
/// @param[in] oldStateIndex index (in state vector) of variable used to calculate predicition                                        
/// @param[in] ratio output/input ratio desired on prediction
void Kalman_SetTransitionFactor( KFilter filter, size_t newStateIndex, size_t oldStateIndex, double ratio );

/// @brief Defines impact of a measurement variable for state estimation          
/// @param[in] filter reference to filter
/// @param[in] measureIndex index of the measure variable in the internal measure vector
/// @param[in] stateIndex index of the correspondent state variable in the internal state vector
/// @param[in] maxError maximum deviation for error modeling
void Kalman_SetMeasureWeight( KFilter filter, size_t measureIndex, size_t stateIndex, double maxError );

/// @brief Sets new value for a single measument
/// @param[in] filter reference to filter
/// @param[in] measureIndex index of the measure variable in the internal measure vector
/// @param[in] value new value of the measurement variable
void Kalman_SetMeasure( KFilter filter, size_t measureIndex, double value );

/// @brief Sets new value for a single input
/// @param[in] filter reference to filter
/// @param[in] inputIndex index of the input variable in the internal input vector
/// @param[in] value new value of the input variable
void Kalman_SetInput( KFilter filter, size_t inputIndex, double value );

/// @brief Runs prediction phase on given Kalman filter                      
/// @param[in] filter reference to filter
/// @param[in] inputsList array of input values (NULL if not updated)
/// @param[out] result pointer to array where predicted internal state will be copied (NULL if not required)
/// @return pointer to @a result array containing predicted filter state (NULL on errors)
double* Kalman_Predict( KFilter filter, double* inputsList, double* result );

/// @brief Runs update phase on given Kalman filter                              
/// @param[in] filter reference to filter
/// @param[in] measuresList array of measurement values (NULL if not updated)             
/// @param[out] result pointer to array where updated internal state will be copied (NULL if not required)
/// @return pointer to @a result array containing updated filter state (NULL on errors)
double* Kalman_Update( KFilter filter, double* measuresList, double* result );
                                                                        
/// @brief Resets internal filter matrices data
/// @param[in] filter reference to filter
void Kalman_Reset( KFilter filter );


#endif  // KALMAN_FILTERS_H
