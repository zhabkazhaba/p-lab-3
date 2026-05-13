cmake_minimum_required(VERSION 3.26)
project(rando_lab)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1")

add_executable(rando_lab
        src/main.cpp
        src/CubeCG.hpp
        src/LCGShuffle.hpp
        src/QuickBBS.hpp
        src/util.hpp
        src/tests.hpp
)