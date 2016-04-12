cmake_minimum_required(VERSION 2.6.0)

include_directories(.)
include_directories(..)
include_directories(../..)

add_library(pir_server_pir STATIC PIROptimizer.cpp GenericPIRReplyGenerator.cpp PIRReplyGeneratorFactory.cpp PIRReplyGeneratorGMP.cpp PIRReplyGeneratorNFL_internal.cpp PIRReplyGeneratorTrivial.cpp)
