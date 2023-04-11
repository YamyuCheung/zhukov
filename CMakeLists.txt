cmake_minimum_required(VERSION 3.20)
project(zhukov)

set(CMAKE_CXX_STANDARD 17)
find_package(MPI REQUIRED)
add_executable(zhukov main.cpp)
target_link_libraries(zhukov PUBLIC MPI::MPI_CXX)