cmake_minimum_required(VERSION 3.16)
project(mpi_sample)

set(CMAKE_CXX_STANDARD 14)
find_package(MPI)
find_package(Eigen3 REQUIRED)
find_package(LAPACK)

include_directories(icecream-cpp /usr/local/opt/lapack/include)
link_directories(/usr/local/opt/lapack/lib)

set(SOURCE_FILES Action.cpp Action.hpp DirectedGraph.hpp DirectedGraph.cpp
        UnionFind.hpp Partition.hpp StrategyM3.hpp StrategyM3.cpp)

add_executable(main main.cpp ${SOURCE_FILES})
target_link_libraries(main PRIVATE MPI::MPI_CXX Eigen3::Eigen lapacke)

add_executable(test_automaton_size test_automaton_size.cpp ${SOURCE_FILES})
target_link_libraries(test_automaton_size PRIVATE Eigen3::Eigen lapacke)
