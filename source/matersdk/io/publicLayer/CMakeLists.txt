# 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/matersdk/io/publicLayer)
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/matersdk/io/publicLayer)
file(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/obj/matersdk/io/publicLayer)
file(MAKE_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})


# 5. 
link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

add_library(CpuNeighborList SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/CpuNeighborList.cc)

# 6.
if (BUILD_TEST)
    list(APPEND CMAKE_CXX_FLAGS "-g")

    add_executable(test_CpuNeighborList ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_CpuNeighborList.cc)
    #target_link_directories(test_CpuNeighborList PRIVATE ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    target_link_libraries(test_CpuNeighborList PRIVATE gtest_main gtest CpuNeighborList)

    gtest_discover_tests(test_CpuNeighborList)    
endif()