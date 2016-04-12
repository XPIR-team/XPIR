cmake_minimum_required(VERSION 2.6.0)

include_directories(..)
include_directories(../..)

add_library(pir_client_optim STATIC PIROptimizer.cpp)
