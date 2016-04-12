cmake_minimum_required(VERSION 2.6.0)

include_directories(../..)

set(CMAKE_CXX_FLAGS "-std=c++11 -fopenmp")

add_library(pir_optim STATIC OptimService.cpp OptimVars.cpp)
