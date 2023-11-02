# 2. Include package
execute_process(COMMAND python3 -c "import pybind11; print(pybind11.__path__[0])" OUTPUT_VARIABLE PYBIND11_PATH)
string(REPLACE "\n" "/share/cmake" PYBIND11_CMAKE_PREFIX_PATH ${PYBIND11_PATH})
list(APPEND CMAKE_PREFIX_PATH ${PYBIND11_CMAKE_PREFIX_PATH})
find_package(pybind11 CONFIG REQUIRED)

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
###pybind11_add_module(structure ${CMAKE_CURRENT_SOURCE_DIR}/bind/structure_bind.cc)


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

    add_executable(test_neighborList_lmp ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_neighborList_lmp.cc)
    target_link_libraries(test_neighborList_lmp gtest gtest_main)

    # gtest_discover_tests(test_CpuNeighborList)    
    gtest_discover_tests(test_structure)
    gtest_discover_tests(test_binLinkedList)
    gtest_discover_tests(test_neighborList)
    gtest_discover_tests(test_neighborList_lmp)
endif()