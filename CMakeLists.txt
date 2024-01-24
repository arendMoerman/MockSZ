cmake_minimum_required(VERSION 3.16.3)

project(MockSZ 
    VERSION 0.1.2 
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

#if(NOT WIN32)
#    target_compile_options(tiempo2 PRIVATE -pthread)
#endif()

