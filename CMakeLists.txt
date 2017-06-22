cmake_minimum_required( VERSION 2.8 )
project( Simple-Kalman-Filters )

set( CMAKE_C_STANDARD 99 )
set( CMAKE_C_STANDARD_REQUIRED ON )

set( LIBRARY_DIR ${CMAKE_CURRENT_LIST_DIR} CACHE PATH "Relative or absolute path to directory where built shared libraries will be placed" )

include( ${CMAKE_CURRENT_LIST_DIR}/matrix/CMakeLists.txt )

add_library( KalmanFilters SHARED ${CMAKE_CURRENT_LIST_DIR}/kalman_filters.c )
set_target_properties( KalmanFilters PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_DIR} )
target_include_directories( KalmanFilters PUBLIC ${CMAKE_CURRENT_LIST_DIR} )
target_compile_definitions( KalmanFilters PUBLIC -DDEBUG )
target_link_libraries( KalmanFilters Matrix )
