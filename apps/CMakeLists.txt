cmake_minimum_required(VERSION 2.6.0)

add_subdirectory("client")
add_subdirectory("optim")
add_subdirectory("server")
add_subdirectory("simplepir")

set(CMAKE_CXX_FLAGS "-std=c++11")
set (CMAKE_EXE_LINKER_FLAGS "-fopenmp")

include_directories(..)
include_directories(../..)

add_custom_target(build-time-make-directory ALL COMMAND ${CMAKE_COMMAND} -E make_directory "server/db")
add_custom_target(build-time-make-directory1 ALL COMMAND ${CMAKE_COMMAND} -E make_directory "server/exp")
add_custom_target(build-time-make-directory2 ALL COMMAND ${CMAKE_COMMAND} -E make_directory "client/reception")
add_custom_target(build-time-make-directory3 ALL COMMAND ${CMAKE_COMMAND} -E make_directory "client/exp")

#add_executable(simple_pir simplePIR.cpp ../pir/libpir.hpp ./server/DBGenerator.cpp ./server/DBDirectoryProcessor.cpp)

#target_link_libraries(simple_pir pir_server_pir pir_query_gen pir_reply pir pir_client_events pir_client_optim pthread ${MPFR_LIBRARIES} ${Boost_LIBRARIES} ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})
