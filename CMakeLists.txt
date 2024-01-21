cmake_minimum_required(VERSION 3.16.3)

project(MockSZ 
    VERSION 0.1.0 
    DESCRIPTION "Mock SZ maps & signals"
    LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(GSL REQUIRED)
add_library(electronstats SHARED src/include/Stats.cpp)
target_include_directories(electronstats PUBLIC src/include)

add_library(signal SHARED src/include/Signal.cpp)
target_include_directories(signal PUBLIC src/include)

add_library(mocksz SHARED src/cpp/InterfaceCPU.cpp)
target_link_libraries(mocksz PRIVATE electronstats GSL::gsl signal)
#target_include_directories(mocksz PRIVATE src/include)

#if(NOT WIN32)
#    target_compile_options(tiempo2 PRIVATE -pthread)
#endif()

#include(CheckLanguage)

#check_language(CUDA)

#if(CMAKE_CUDA_COMPILER)
#    enable_language(CUDA)
#    set(CMAKE_CUDA_STANDARD 11)
#    file(GLOB CUDAfiles src/cuda/*.cu)

#    add_library(cutiempo2 SHARED ${CUDAfiles})
#    target_include_directories(cutiempo2 PRIVATE src/cuda src/include)
#    target_link_libraries(cutiempo2 PRIVATE timer)
#    set_target_properties(cutiempo2 PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#    set_property(TARGET cutiempo2 PROPERTY CUDA_ARCHITECTURES OFF)
#else()
#    message(STATUS "WARNING: no CUDA compiler detected. Not building CUDA libraries")
#endif()
