# 3. 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/descriptor/mtp)
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/descriptor/mtp)
file(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/obj/descriptor/mtp)
file(MAKE_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})


# 5. compile
add_library(mtpLevel SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/mtpLevel.cc)

if (TORCH_OP)    
    execute_process(COMMAND python3 -c "import torch;print(torch.__path__[0])" OUTPUT_VARIABLE TORCH_PATH)
    string(REPLACE "\n" "/share/cmake" TORCH_CMAKE_PREFIX_PATH ${TORCH_PATH})
    list(APPEND CMAKE_PREFIX_PATH ${TORCH_CMAKE_PREFIX_PATH})
    find_package(Torch CONFIG REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

    add_library(mtpMAngularOp SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/mtpMAngularOp.cc)
    target_include_directories(mtpMAngularOp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(mtpMAngularOp PRIVATE ${TORCH_LIBRARIES})

    add_library(mtpMRadialOp SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/mtpMRadialOp.cc)
    target_include_directories(mtpMRadialOp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(mtpMRadialOp PRIVATE ${TORCH_LIBRARIES})

    add_library(mtpMModule SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/mtpMModule.cc)
    target_include_directories(mtpMModule PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(mtpMModule PRIVATE ${TORCH_LIBRARIES} mtpMRadialOp mtpMAngularOp)
endif()


# 6. compile test
if (BUILD_TEST)
    add_executable(test_mtpMRadial_de ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_mtpMRadial_de.cc)
    target_link_libraries(test_mtpMRadial_de PRIVATE gtest gtest_main)

    add_executable(test_mtpMRadial ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_mtpMRadial.cc)
    target_link_libraries(test_mtpMRadial PRIVATE gtest gtest_main)

    add_executable(test_level ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_level.cc)
    target_link_libraries(test_level PRIVATE gtest gtest_main)

    add_executable(test_mtpLevel ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_mtpLevel.cc)
    target_link_libraries(test_mtpLevel PRIVATE gtest gtest_main mtpLevel)
    
    if (TORCH_OP)
        add_executable(test_mtpMAngularOp ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_mtpMAngularOp.cc)
        target_link_libraries(test_mtpMAngularOp PRIVATE gtest gtest_main ${TORCH_LIBRARIES} mtpMAngularOp)
        if(OpenMP_CXX_FOUND)
            target_link_libraries(test_mtpMAngularOp PRIVATE OpenMP::OpenMP_CXX)
        endif()

        add_executable(test_mtpMRadialOp ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_mtpMRadialOp.cc)
        target_link_libraries(test_mtpMRadialOp PRIVATE gtest gtest_main ${TORCH_LIBRARIES} mtpMRadialOp)

        add_executable(test_mtpMModule ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_mtpMModule.cc)
        target_link_libraries(test_mtpMModule PRIVATE gtest gtest_main ${TORCH_LIBRARIES} mtpMModule)
    endif()

    gtest_discover_tests(test_mtpMRadial_de)
    gtest_discover_tests(test_mtpMRadial)
    gtest_discover_tests(test_level)
    if (TORCH_OP)
        gtest_discover_tests(test_mtpMAngularOp)
        gtest_discover_tests(test_mtpMRadialOp)
    endif()
endif()
