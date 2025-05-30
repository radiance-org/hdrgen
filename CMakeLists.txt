cmake_minimum_required(VERSION 3.10)
project(hdrgen_tools VERSION 1.0 LANGUAGES C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_policy(SET CMP0079 NEW)
cmake_policy(SET CMP0077 NEW)

# Build type
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

# ExternalProject_Add(panlib SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/panlib)
add_subdirectory(external/panlib)

set(PANLIB_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/panlib/src)
set(RAD_COMMON_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/panlib/external/radiance/src/common)

add_executable(hdrgen src/hdrgenmain.cpp)
target_include_directories(hdrgen PRIVATE ${PANLIB_SOURCE_DIR} ${RAD_COMMON_DIR})
target_link_libraries(hdrgen PRIVATE pan)
install(TARGETS hdrgen
    RUNTIME DESTINATION bin
)

add_executable(hdrcvt src/hdrcvtmain.cpp)
target_include_directories(hdrcvt PRIVATE ${PANLIB_SOURCE_DIR} ${RAD_COMMON_DIR})
target_link_libraries(hdrcvt PRIVATE pan)
install(TARGETS hdrcvt
    RUNTIME DESTINATION bin
)

add_executable(PQconvert src/PQconvert.cpp src/pq2.h)
target_include_directories(PQconvert PRIVATE ${PANLIB_SOURCE_DIR} ${RAD_COMMON_DIR})
target_link_libraries(PQconvert PRIVATE pan)
install(TARGETS PQconvert
    RUNTIME DESTINATION bin
)

add_executable(expose2range src/expose2range.c)
target_include_directories(expose2range PRIVATE ${PANLIB_SOURCE_DIR} ${RAD_COMMON_DIR})
target_link_libraries(expose2range PRIVATE pan)
install(TARGETS expose2range
    RUNTIME DESTINATION bin
)

add_executable(warpimage src/warpimage.c)
target_include_directories(warpimage PRIVATE ${PANLIB_SOURCE_DIR} ${RAD_COMMON_DIR})
target_link_libraries(warpimage PRIVATE pan)
install(TARGETS warpimage
    RUNTIME DESTINATION bin
)

add_executable(bitmapop src/bitmapop.cpp)
target_include_directories(bitmapop PRIVATE ${PANLIB_SOURCE_DIR} ${RAD_COMMON_DIR})
target_link_libraries(bitmapop PRIVATE pan)
install(TARGETS bitmapop
    RUNTIME DESTINATION bin
)
