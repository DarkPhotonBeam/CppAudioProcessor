cmake_minimum_required(VERSION 3.29)
project(waveAnalysis)

set(CMAKE_CXX_STANDARD 20)

find_package(Eigen3 REQUIRED)

include_directories(waveAnalysis ${EIGEN3_INCLUDE_DIRS})
include_directories(waveAnalysis /usr/include/python3.12)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/sample.wav DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(waveAnalysis main.cpp)

target_link_libraries(waveAnalysis python3.12)