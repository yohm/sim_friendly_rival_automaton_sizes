cmake_minimum_required(VERSION 3.16)
project(mpi_sample)

set(CMAKE_CXX_STANDARD 14)
find_package(MPI)
find_package(Eigen3 REQUIRED)
find_package(LAPACK)

include_directories(icecream-cpp /usr/local/opt/lapack/include)
link_directories(/usr/local/opt/lapack/lib)

add_executable(main main.cpp StrategyM3.hpp StrategyM3.cpp Action.cpp Action.hpp DirectedGraph.cpp DirectedGraph.hpp UnionFind.hpp)
target_link_libraries(main PRIVATE MPI::MPI_CXX Eigen3::Eigen lapacke)