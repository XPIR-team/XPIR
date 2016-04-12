cmake_minimum_required(VERSION 2.6.0)

set (CMAKE_EXE_LINKER_FLAGS "-fopenmp") 
if (APPLE)
	set(CMAKE_EXE_LINKER_FLAGS "-fopenmp -L$ENV{HOME}/Code/xpir/freshxpir/xpir-ng/local/lib/ -lboost_program_options")

endif()
include_directories(..)
include_directories(../..)

add_executable(pir_client DESC.cpp PIRClient.cpp PIRClientLog.cpp PIRController.cpp PIRView.cpp PIRViewCLI.cpp main.cpp)

target_link_libraries(pir_client pir_client_optim pir pthread ${MPFR_LIBRARIES} ${Boost_LIBRARIES} ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})

