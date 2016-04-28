cmake_minimum_required(VERSION 2.6.0)


include_directories(.)
include_directories(..)
include_directories(../..)

add_custom_target(build-time-make-directory-simple-pir ALL COMMAND ${CMAKE_COMMAND} -E make_directory "db")
file(COPY "../tools/makedb.sh" DESTINATION ".")


add_executable(simple_pir simplePIR.cpp ../../libpir.hpp)
target_link_libraries(simple_pir pir pthread ${MPFR_LIBRARIES} ${Boost_LIBRARIES} ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})


add_custom_command(
    OUTPUT .simplepirdb
    WORKING_DIR ${PROJECT_BINARY_DIR}/apps/simplepir
    COMMAND ./makedb.sh 1024 16 2>/dev/null
    COMMAND touch .simplepirdb
)

add_custom_target(simplepirdb ALL DEPENDS .simplepirdb)
	
