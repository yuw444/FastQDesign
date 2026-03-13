cmake_minimum_required(VERSION 3.16.0)

project(fastF VERSION 1.0)

set(CMAKE_C_STANDARD 99)
set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_BUILD_TYPE Release)

find_package(OpenMP REQUIRED)

# External library paths
set(HTSLIB_PATH /hpc/apps/htslib/1.22.1)

# Include directories
include_directories(${HTSLIB_PATH}/include)
include_directories(${CMAKE_SOURCE_DIR}/src)

# Link directories
link_directories(${HTSLIB_PATH}/lib)

# Source files (all in src/)
set(SOURCES
    ${CMAKE_SOURCE_DIR}/src/main.c
    ${CMAKE_SOURCE_DIR}/src/argparse.c
    ${CMAKE_SOURCE_DIR}/src/mt19937ar.c
    ${CMAKE_SOURCE_DIR}/src/hashtable.c
    ${CMAKE_SOURCE_DIR}/src/filter.c
    ${CMAKE_SOURCE_DIR}/src/bam2db.c
    ${CMAKE_SOURCE_DIR}/src/utils.c
)

add_executable(fastF ${SOURCES})

target_link_libraries(fastF z pthread hts sqlite3)

# Tests
enable_testing()

# Unity framework
set(UNITY_SOURCES
    ${CMAKE_SOURCE_DIR}/tests/unity/unity.c
)

# Test: utils
set(TEST_UTILS_SOURCES
    ${CMAKE_SOURCE_DIR}/tests/c/test_utils.c
    ${CMAKE_SOURCE_DIR}/src/filter.c
    ${CMAKE_SOURCE_DIR}/src/utils.c
    ${CMAKE_SOURCE_DIR}/src/bam2db.c
    ${CMAKE_SOURCE_DIR}/src/hashtable.c
    ${CMAKE_SOURCE_DIR}/src/mt19937ar.c
    ${UNITY_SOURCES}
)
add_executable(test_utils ${TEST_UTILS_SOURCES})
target_include_directories(test_utils PRIVATE 
    ${CMAKE_SOURCE_DIR}/tests/unity
    ${CMAKE_SOURCE_DIR}/tests/c
    ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(test_utils z m hts sqlite3)
add_test(NAME test_utils COMMAND test_utils
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/c)

# Test: filter
set(TEST_FILTER_SOURCES
    ${CMAKE_SOURCE_DIR}/tests/c/test_filter.c
    ${CMAKE_SOURCE_DIR}/src/filter.c
    ${UNITY_SOURCES}
)
add_executable(test_filter ${TEST_FILTER_SOURCES})
target_include_directories(test_filter PRIVATE 
    ${CMAKE_SOURCE_DIR}/tests/unity
    ${CMAKE_SOURCE_DIR}/tests/c
    ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(test_filter z m)
add_test(NAME test_filter COMMAND test_filter
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/c)

# Test: bam2db
set(TEST_BAM2DB_SOURCES
    ${CMAKE_SOURCE_DIR}/tests/c/test_bam2db.c
    ${CMAKE_SOURCE_DIR}/src/bam2db.c
    ${CMAKE_SOURCE_DIR}/src/utils.c
    ${CMAKE_SOURCE_DIR}/src/filter.c
    ${CMAKE_SOURCE_DIR}/src/hashtable.c
    ${CMAKE_SOURCE_DIR}/src/mt19937ar.c
    ${UNITY_SOURCES}
)
add_executable(test_bam2db ${TEST_BAM2DB_SOURCES})
target_include_directories(test_bam2db PRIVATE 
    ${CMAKE_SOURCE_DIR}/tests/unity
    ${CMAKE_SOURCE_DIR}/tests/c
    ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(test_bam2db z m hts sqlite3)
add_test(NAME test_bam2db COMMAND test_bam2db
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/c)

# Test: integration
set(TEST_INTEGRATION_SOURCES
    ${CMAKE_SOURCE_DIR}/tests/c/test_integration.c
    ${UNITY_SOURCES}
)
add_executable(test_integration ${TEST_INTEGRATION_SOURCES})
target_include_directories(test_integration PRIVATE 
    ${CMAKE_SOURCE_DIR}/tests/unity
    ${CMAKE_SOURCE_DIR}/tests/c
)
target_link_libraries(test_integration z m)
add_test(NAME test_integration COMMAND test_integration
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
