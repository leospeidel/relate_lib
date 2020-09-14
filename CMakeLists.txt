project                     (Relate)
cmake_minimum_required      (VERSION 3.1)

set(CMAKE_C_FLAGS "-std=c99")
set(CMAKE_CXX_FLAGS "-O3")        ## Optimize
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} ) 
endif( ZLIB_FOUND )

# Set dir
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/include/src)
set(TEST_DIR ${CMAKE_SOURCE_DIR}/include/test)
set(GZSTREAM_DIR ${CMAKE_SOURCE_DIR}/include/src/gzstream)
set(TSKIT_DIR ${CMAKE_SOURCE_DIR}/include/src/tskit)
set(EXAMPLE_DIR ${CMAKE_SOURCE_DIR}/include/example)
set(CONVERT_DIR ${CMAKE_SOURCE_DIR}/include/convert)
set(INTROGRESSION_DIR ${CMAKE_SOURCE_DIR}/include/introgression)

#add_definitions(-static)
#set(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -static")

# Set includes
include_directories(${PROJECT_SOURCE_DIR} ${TEST_DIR})
include_directories(${PROJECT_SOURCE_DIR} ${GZSTREAM_DIR})
include_directories(${PROJECT_SOURCE_DIR} ${TSKIT_DIR})

# Add subdirectories
add_subdirectory(${PROJECT_SOURCE_DIR})
add_subdirectory(${GZSTREAM_DIR})
add_subdirectory(${TSKIT_DIR})
add_subdirectory(${TEST_DIR})
add_subdirectory(${EXAMPLE_DIR})
add_subdirectory(${CONVERT_DIR})
add_subdirectory(${INTROGRESSION_DIR})

# tests
add_executable(Tests ${TEST_DIR}/Tests.cpp)
target_link_libraries(Tests relateShared UnitTest)

enable_testing()
add_test(UnitTest ${CMAKE_BINARY_DIR}/Tests)
add_custom_command(TARGET Tests COMMENT "Run tests" POST_BUILD WORKING_DIRECTORY ${CMAKE_BINARY_DIR} COMMAND Tests)
