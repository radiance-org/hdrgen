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


add_subdirectory(external/panlib)

add_subdirectory(external/radiance/src/common)

add_subdirectory(external/openexr)

set(tiff-tools OFF)
set(tiff-tests OFF)
set(tiff-docs OFF)
set(tiff-contrib OFF)
add_subdirectory(external/libtiff/)

ExternalProject_Add(libjpeg
SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/libjpeg-turbo
CMAKE_ARGS
	-DENABLE_STATIC=TRUE
	-DENABLE_SHARED=FALSE
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
)

# Include directories
# include_directories(
#     ${CMAKE_SOURCE_DIR}/src
# )

# Define executables
add_executable(hdrgen src/hdrgenmain.cpp)
add_executable(hdrcvt src/hdrcvtmain.cpp)
# add_executable(expose2range src/expose2range.c)
# add_executable(warpimage src/warpimage.c)

# Link libraries
target_link_libraries(hdrgen
        PRIVATE
            panlib
            rtrad
            IlmImf
            Imath
            Half
            Iex
            IlmThread
	    jpeg-static
	    tiff
            # ZLIB::ZLIB
)

# Installation
install(TARGETS hdrgen hdrcvt
    RUNTIME DESTINATION bin
)
