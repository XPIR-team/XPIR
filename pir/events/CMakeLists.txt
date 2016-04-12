cmake_minimum_required(VERSION 2.6.0)

include_directories(..)
include_directories(../..)

add_library(pir_client_events STATIC CatalogEvent.cpp MessageEvent.cpp WriteEvent.cpp)
