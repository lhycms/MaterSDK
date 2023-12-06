# 2. include 
#add_subdirectory("/data/home/liuhanyu/hyliu/code/matersdk/source/ext/googletest-main")
#include(GoogleTest)
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)


# 3. dirs
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/pymatgen4cpp)
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/pymatgen4cpp)
file(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/obj/pymatgen4cpp)
file(MAKE_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})


# 4. add subdirectory

# 5. compile library
add_library(load_file SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/load_file.cc)
target_include_directories(load_file PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(load_file PRIVATE ${PYTHON_INCLUDE_DIRS})
target_link_libraries(load_file PRIVATE ${PYTHON_LIBRARIES} Python3::NumPy)

# 6. compile test
if (BUILD_TEST)
    add_executable(test_structure_pymatgen ${CMAKE_CURRENT_SOURCE_DIR}/test/test_structure.cc)
    target_include_directories(test_structure_pymatgen PRIVATE ${PYTHON_INCLUDE_DIRS})
    target_link_libraries(test_structure_pymatgen PRIVATE gtest gtest_main load_file ${PYTHON_LIBRARIES} Python3::NumPy)

    gtest_discover_tests(test_structure_pymatgen)
endif()