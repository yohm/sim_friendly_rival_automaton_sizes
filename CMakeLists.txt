cmake_minimum_required(VERSION 3.16)
project(mpi_sample)

set(CMAKE_CXX_STANDARD 14)
find_package(MPI)

add_executable(mpi_sample main.cpp)
target_link_libraries(mpi_sample PRIVATE MPI::MPI_CXX)