cmake_minimum_required(VERSION 3.16.3)

project(MockSZ
    VERSION 0.2.4
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

# Add GSL to electronstats and signal
if(TARGET GSL::gsl)
    target_link_libraries(electronstats PRIVATE GSL::gsl)
    target_link_libraries(signal PRIVATE GSL::gsl)
else()
    target_include_directories(electronstats PRIVATE ${GSL_INCLUDE_DIRS})
    target_link_libraries(electronstats PRIVATE ${GSL_LIBRARIES})
    target_include_directories(signal PRIVATE ${GSL_INCLUDE_DIRS})
    target_link_libraries(signal PRIVATE ${GSL_LIBRARIES})
endif()

add_library(romb SHARED src/include/Romberg.cpp)
target_include_directories(romb PUBLIC src/include)

add_library(mocksz SHARED src/cpp/InterfaceCPU.cpp)
target_link_libraries(mocksz PRIVATE electronstats GSL::gsl signal romb)

if(NOT WIN32)
    target_compile_options(mocksz PRIVATE)
endif()
