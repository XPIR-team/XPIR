cmake_minimum_required(VERSION 2.6.0)


include_directories(.)
include_directories(..)
include_directories(../..)

add_executable(pir_server PIRServer.cpp PIRSession.cpp ServerService.cpp main.cpp)

target_link_libraries(pir_server pir ${Boost_LIBRARIES} ${GMP_LIBRARIES} pthread )

file(COPY "../tools/check-correctness.sh" DESTINATION ".")
file(COPY "../tools/makedb.sh" DESTINATION ".")
file(COPY "../tools/mkdb-correctness.sh" DESTINATION ".")

add_custom_command(
    OUTPUT .db
    WORKING_DIR ${PROJECT_BINARY_DIR}/apps/server
    COMMAND ./makedb.sh 1024 16 2>/dev/null
    COMMAND touch .db
)

add_custom_target(db ALL DEPENDS .db)
	
add_custom_command(
    OUTPUT .db-correct
    WORKING_DIR ${PROJECT_BINARY_DIR}/apps/server
    COMMAND ./mkdb-correctness.sh
    COMMAND touch .db-correct
)

add_custom_target(dbcorrect DEPENDS .db-correct)

add_custom_command(
    OUTPUT .check
    WORKING_DIR ${PROJECT_BINARY_DIR}/apps/server
    COMMAND ./check-correctness.sh
)

add_custom_target(check DEPENDS .check dbcorrect pir_client pir_server)
