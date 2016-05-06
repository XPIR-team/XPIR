cmake_minimum_required(VERSION 2.6.0)
project(XPIR)
include(CMakeCompilers.txt)

enable_language(CXX ASM)

if(APPLE)
    cmake_policy(SET CMP0042 NEW)
endif()

#set(CMAKE_C_COMPILER "gcc")
#set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-std=c++11 -fopenmp -ggdb -Wno-deprecated-declarations")
set(CMAKE_EXE_LINKER_FLAGS "-fopenmp")
set(CMAKE_LIBRARY_PATH "./local")

# Boost 
#  you can give hints to find_package in order to find your favorite boost build
#  if you use the boost coming with XPIR helper_script.sh, set BOOST_ROOT to ./local as in the following line 
set(BOOST_ROOT "./local")
#  otherwise set it to your own boost build, by uncomenting the following one for example
#  but the order in which cmake.find_package searches potential directories is obscure and 
#  we recommend you double-check the output of cmake to see if it picked the proper build
# set(BOOST_ROOT "/usr/local/Cellar/boost/1.60.0_2/")

SET( Boost_USE_STATIC_LIBS FALSE )

find_package(Boost 1.55.0  REQUIRED
  COMPONENTS atomic chrono date_time exception program_options regex system thread program_options)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})
 message(STATUS "Boost_lib=" ${Boost_LIBRARY_DIR})
 message(STATUS "Boost_inc=" ${Boost_INCLUDE_DIR})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})

find_package(GMP 6 REQUIRED)
find_package(MPFR 3.1.2 REQUIRED)

include_directories("/opt/local/include")

option(SEND_CATALOG "Send the catalog to the client (default is send catalog if |catalog|<1000)" ON)
if(SEND_CATALOG)
	message(STATUS "Send the catalog to the client")
	add_definitions(-DSEND_CATALOG)
else()
	message(STATUS "Do not send the catalog to the client")	
endif()
	
option(MULTI_THREAD "Use multi-threading" ON)
if(MULTI_THREAD)
	message(STATUS "Use multi-threading")
	add_definitions(-DMULTI_THREAD)
else()
	message(STATUS "Do not use multi-threading")		
endif()
	
option(PERF_TIMERS "Show performance measurements during execution" ON)
if(PERF_TIMERS)
	message(STATUS "Show performance measurements during execution")
	add_definitions(-DPERF_TIMERS)
else()
	message(STATUS "Do not show performance measurements during execution")
endif()

add_subdirectory("apps")
add_subdirectory("crypto")
add_subdirectory("pir")
	
	

