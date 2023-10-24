# 3. 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/matersdk/feature/deepmd)
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/matersdk/feature/deepmd)
file(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/obj/matersdk/feature/deepmd)
file(MAKE_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})


# 5. add_library



# 6. test
if (BUILD_TEST)
    add_executable(test_se ${CMAKE_CURRENT_SOURCE_DIR}/test/test_se.cc)
    target_link_libraries(test_se gtest gtest_main)

    add_executable(test_se_lmp ${CMAKE_CURRENT_SOURCE_DIR}/test/test_se_lmp.cc)
    target_link_libraries(test_se_lmp gtest gtest_main)

    add_executable(test_se4pw ${CMAKE_CURRENT_SOURCE_DIR}/test/test_se4pw.cc)
    target_link_libraries(test_se4pw gtest gtest_main)

    add_executable(demo ${CMAKE_CURRENT_SOURCE_DIR}/test/demo.cc)

    gtest_discover_tests(test_se)
    gtest_discover_tests(test_se_lmp)
    gtest_discover_tests(test_se4pw)
endif()