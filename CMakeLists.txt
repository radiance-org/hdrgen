cmake_minimum_required(VERSION 3.10)
project(hdrgen_tools VERSION 1.0 LANGUAGES C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_policy(SET CMP0079 NEW)
cmake_policy(SET CMP0077 NEW)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()


# ExternalProject_Add(panlib SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/panlib)
add_subdirectory(external/panlib)

# Define executables
add_executable(hdrgen src/hdrgenmain.cpp)
add_executable(hdrcvt src/hdrcvtmain.cpp)

# Link libraries
target_link_libraries(hdrgen PRIVATE pan)
target_link_libraries(hdrcvt PRIVATE pan)
# target_include_directories(hdrgen PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/panlib/src)
# target_include_directories(hdrcvt PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/panlib/src)

# Installation
install(TARGETS hdrgen hdrcvt
    RUNTIME DESTINATION bin
)
