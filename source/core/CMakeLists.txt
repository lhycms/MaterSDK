#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/core)
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/core)
file(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/obj/core)
file(MAKE_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

# 5. 
if(BUILD_TEST)
    add_executable(test_vecx ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_vecx.cc)
    target_link_libraries(test_vecx PRIVATE gtest_main gtest)

    add_executable(test_AlignedArray ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_AlignedArray.cc)
    target_link_libraries(test_AlignedArray PRIVATE gtest_main gtest)
    
    gtest_discover_tests(test_vecx)
endif()