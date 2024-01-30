# 3. 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/descriptor/mtp)
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/descriptor/mtp)
file(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/obj/descriptor/mtp)
file(MAKE_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})


# 5. compile
if (TORCH_OP)    
    execute_process(COMMAND python3 -c "import torch;print(torch.__path__[0])" OUTPUT_VARIABLE TORCH_PATH)
    string(REPLACE "\n" "/share/cmake" TORCH_CMAKE_PREFIX_PATH ${TORCH_PATH})
    list(APPEND CMAKE_PREFIX_PATH ${TORCH_CMAKE_PREFIX_PATH})
    find_package(Torch CONFIG REQUIRED)
    list(APPEND CMAKE_CXX_FLAGS ${TORCH_CXX_FLAGS})

    add_library(mtpMAngularOp SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/mtpMAngularOp.cc)
    target_include_directories(mtpMAngularOp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(mtpMAngularOp PRIVATE ${TORCH_LIBRARIES})
endif()


# 6. compile test
if (BUILD_TEST)
    add_executable(test_radial_basis ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_radial_basis.cc)
    target_link_libraries(test_radial_basis gtest gtest_main)

    add_executable(test_level ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_level.cc)
    target_link_libraries(test_level gtest gtest_main)
    
    if (TORCH_OP)
        add_executable(test_mtpMAngularOp ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_mtpMAngularOp.cc)
        target_link_libraries(test_mtpMAngularOp gtest gtest_main mtpMAngularOp ${TORCH_LIBRARIES})
    endif()

    gtest_discover_tests(test_radial_basis)
    gtest_discover_tests(test_level)
    if (TORCH_OP)
        gtest_discover_tests(test_mtpMAngularOp)
    endif()
endif()