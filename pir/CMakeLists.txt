cmake_minimum_required(VERSION 2.8)
#project(nfl CXX ASM-ATT)
enable_testing()
if(APPLE)
    cmake_policy(SET CMP0042 NEW)
endif()

find_package(GMP REQUIRED)

include_directories("/opt/local/include" ${CMAKE_SOURCE_DIR}/crypto ${CMAKE_SOURCE_DIR}/pir ${CMAKE_SOURCE_DIR})

add_subdirectory("events")
add_subdirectory("optim")
add_subdirectory("dbhandlers")
add_subdirectory("queryGen")
add_subdirectory("replyExtraction")
add_subdirectory("replyGenerator")

file(GLOB PIRLIB_ASMS ${CMAKE_SOURCE_DIR}/crypto/prng/*.s )
file(GLOB PIRLIB_SRCS 
	${CMAKE_SOURCE_DIR}/crypto/*.cpp 
	${CMAKE_SOURCE_DIR}/crypto/prng/*.cpp 
	${CMAKE_SOURCE_DIR}/pir/*.cpp 
	${CMAKE_SOURCE_DIR}/apps/client/DESC.cpp 
	${CMAKE_SOURCE_DIR}/apps/server/DB*.cpp 
	${CMAKE_SOURCE_DIR}/pir/events/*.cpp 
	${CMAKE_SOURCE_DIR}/pir/optim/*.cpp 
	${CMAKE_SOURCE_DIR}/pir/dbhandlers/*.cpp 
	${CMAKE_SOURCE_DIR}/pir/queryGen/*.cpp 
	${CMAKE_SOURCE_DIR}/pir/replyExtraction/*.cpp 
	${CMAKE_SOURCE_DIR}/pir/replyGenerator/*.cpp )

add_library(pir SHARED ${PIRLIB_ASMS} ${PIRLIB_SRCS})
target_link_libraries(pir  ${MPFR_LIBRARY} ${Boost_LIBRARIES} ${GMP_LIBRARIES} ${GOBJECT_LIBRARIES})
install(TARGETS pir LIBRARY DESTINATION lib)

add_library(pir_static STATIC ${PIRLIB_ASMS} ${PIRLIB_SRCS})
target_link_libraries(pir_static ${MPFR_LIBRARY} ${Boost_LIBRARIES} ${GMP_LIBRARY} ${GOBJECT_LIBRARIES} gmp libgmp )
install(TARGETS pir_static ARCHIVE DESTINATION lib)
