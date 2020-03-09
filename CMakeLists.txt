cmake_minimum_required( VERSION 2.8 )
project( Simple-Kalman-Filter )

set( CMAKE_C_STANDARD 99 )
set( CMAKE_C_STANDARD_REQUIRED ON )

set( LIBRARY_DIR ${CMAKE_CURRENT_LIST_DIR} CACHE PATH "Relative or absolute path to directory where built shared libraries will be placed" )

if( NOT TARGET Matrix )
  include( ${CMAKE_CURRENT_LIST_DIR}/matrix/CMakeLists.txt )
endif()

add_library( KalmanFilter SHARED ${CMAKE_CURRENT_LIST_DIR}/kalman_filters.c )
set_target_properties( KalmanFilter PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_DIR} )
target_include_directories( KalmanFilter PUBLIC ${CMAKE_CURRENT_LIST_DIR} )
target_compile_definitions( KalmanFilter PUBLIC -DDEBUG )
target_link_libraries( KalmanFilter Matrix )
