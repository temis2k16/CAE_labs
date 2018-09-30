cmake_minimum_required(VERSION 3.8)
project(myMKE)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp methods.cpp)
add_executable(myMKE ${SOURCE_FILES})