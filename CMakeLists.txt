cmake_minimum_required(VERSION 3.17)
project(runge_kutta)

set(CMAKE_CXX_STANDARD 17)

add_compile_options(-O2 -Wall)

add_executable(runge_kutta main.cpp)