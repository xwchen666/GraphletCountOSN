cmake_minimum_required(VERSION 2.8)

# enable warnings (always good)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -std=c++11 -Wno-format\
        -fopenmp -lstdc++ -ffast-math -funroll-loops -fno-strict-aliasing \
	 -fomit-frame-pointer -fexpensive-optimizations -funroll-loops \
	 -fmove-loop-invariants -fprefetch-loop-arrays -ftree-loop-optimize \
	 -ftree-vect-loop-version -ftree-vectorize")

# add src include path
include_directories("include")

file( GLOB RW_LIB_SOURCES src/*.cpp )
file( GLOB RW_LIB_HEADERS include/*.h )

add_library( RW_Lib ${RW_LIB_SOURCES} ${RW_LIB_HEADERS} )

# number of triangles
add_executable(triangle_compute app/triangle_compute.cpp)
target_link_libraries(triangle_compute RW_Lib)

# graphlet four nodes graphlet compuation
add_executable(graphlet_test app/graphlet_test.cpp)
target_link_libraries(graphlet_test RW_Lib)

# orca five node graphlet computation
add_executable(orca_test app/orca_test.cpp)
target_link_libraries(orca_test RW_Lib)

################counts#################################
add_executable(subgraphCounts app/subgraphCounts.cpp)
target_link_libraries(subgraphCounts RW_Lib)

################comparison#################################
add_executable(subgraphConcentration app/subgraphConcentration.cpp)
target_link_libraries(subgraphConcentration RW_Lib)


