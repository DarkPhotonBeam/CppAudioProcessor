cmake_minimum_required(VERSION 3.29)
project(waveAnalysis)

set(CMAKE_CXX_STANDARD 20)

find_package(Eigen3 REQUIRED)
add_compile_definitions(WITHOUT_NUMPY)

link_directories(/Library/Frameworks/Python.framework/Versions/3.10/lib)

include_directories(waveAnalysis ${EIGEN3_INCLUDE_DIRS})
include_directories(waveAnalysis /Library/Frameworks/Python.framework/Versions/3.10/Headers)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/audio DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(waveAnalysis main.cpp)

target_link_libraries(waveAnalysis python3.10)