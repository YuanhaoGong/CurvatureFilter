cmake_minimum_required(VERSION 2.8)
set(CMAKE_CONFIGURATION_TYPES Release)
project( CF )
find_package( OpenCV REQUIRED )
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "set build type to release")
include_directories( ${OpenCV_INCLUDE_DIRS} )
add_executable( cf main.cpp )
target_link_libraries( cf ${OpenCV_LIBS} )
set(CMAKE_BUILD_TYPE Release)
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")        ## Optimize
endif()
if(CMAKE_COMPILER_IS_GNUC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native")        ## Optimize
endif()
if(MSVC)

endif(MSVC)
