### Usage: Please `include(matersdk.cmake)` in you CMakeList.txt

# 1. Basic
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
project(matersdk)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)

# 2. Include

# 3. Set variable
get_filename_component(MATERSDK_DIR ${CMAKE_CURRENT_LIST_DIR}/.. ABSOLUTE)

set(MATERSDK_INCLUDE_DIRS)
list(APPEND MATERSDK_INCLUDE_DIRS ${MATERSDK_DIR}/source/include;)
list(APPEND MATERSDK_INCLUDE_DIRS ${MATERSDK_DIR}/source/core/include;)
list(APPEND MATERSDK_INCLUDE_DIRS ${MATERSDK_DIR}/source/matersdk/io/publicLayer/include;)
list(APPEND MATERSDK_INCLUDE_DIRS ${MATERSDK_DIR}/source/matersdk/feature/deepmd/include;)