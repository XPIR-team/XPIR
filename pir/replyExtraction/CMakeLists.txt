cmake_minimum_required(VERSION 2.6.0)

include_directories(..)
include_directories(../..)

add_library(pir_reply STATIC  PIRReplyExtraction_internal.cpp PIRReplyWriter.cpp)
