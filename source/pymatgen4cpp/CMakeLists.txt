# 2. include 
# find_package(pybind11 CONFIG REQUIRED), then you can use ${PYTHON_INCLUDE_DIRS}, ${PYTHON_LIBRARIES}


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
target_include_directories(load_file PRIVATE "/data/home/liuhanyu/anaconda3/envs/mlff_1/lib/python3.9/site-packages/numpy/core/include")
target_link_libraries(load_file PRIVATE ${PYTHON_LIBRARIES})

# 6. compile test
if (BUILD_TEST)
    add_executable(test_structure_pymatgen ${CMAKE_CURRENT_SOURCE_DIR}/test/test_structure.cc)
    target_link_libraries(test_structure_pymatgen PRIVATE gtest gtest_main load_file)
    target_include_directories(test_structure_pymatgen PRIVATE ${PYTHON_INCLUDE_DIRS})
    target_include_directories(test_structure_pymatgen PRIVATE "/data/home/liuhanyu/anaconda3/envs/mlff_1/lib/python3.9/site-packages/numpy/core/include")
    #target_link_libraries(test_structure_pymatgen PRIVATE ${PYTHON_LIBRARIES})
    
    gtest_discover_tests(test_structure_pymatgen)

endif()