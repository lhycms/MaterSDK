# 3. 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/matersdk/io/publicLayer)
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/matersdk/io/publicLayer)
file(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/obj/matersdk/io/publicLayer)
file(MAKE_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})


# 5. 
#add_library(CpuNeighborList SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/CpuNeighborList.cc)
#add_library(structure SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/structure.cc)

# 6.
if (BUILD_TEST)
    list(APPEND CMAKE_CXX_FLAGS "-g")

    # add_executable(test_CpuNeighborList ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_CpuNeighborList.cc)
    # target_link_libraries(test_CpuNeighborList PRIVATE gtest_main gtest CpuNeighborList)

    add_executable(test_structure ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_structure.cc)
    target_link_libraries(test_structure gtest_main gtest)

    add_executable(test_binLinkedList ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_binLinkedList.cc)
    target_link_libraries(test_binLinkedList gtest_main gtest)

    add_executable(test_neighborList ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_neighborList.cc)
    target_link_libraries(test_neighborList gtest gtest_main)

    # gtest_discover_tests(test_CpuNeighborList)    
    gtest_discover_tests(test_structure)
    gtest_discover_tests(test_binLinkedList)
    gtest_discover_tests(test_neighborList)
endif()