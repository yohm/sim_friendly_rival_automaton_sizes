cmake_minimum_required(VERSION 3.16)
project(sim_friendly_rival_automaton_sizes)

set(CMAKE_CXX_STANDARD 14)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(LAPACK REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/icecream-cpp /usr/local/opt/lapack/include ${CMAKE_SOURCE_DIR}/caravan-lib ${CMAKE_SOURCE_DIR}/caravan-lib/json/include)
link_directories(/usr/local/opt/lapack/lib)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE MPI::MPI_CXX OpenMP::OpenMP_CXX)

add_executable(minimum minimum.cpp)

set(SOURCE_FILES Action.cpp Action.hpp DirectedGraph.hpp DirectedGraph.cpp
        UnionFind.hpp Partition.hpp StrategyM3.hpp StrategyM3.cpp)

add_executable(test_automaton_size test_automaton_size.cpp ${SOURCE_FILES})
target_link_libraries(test_automaton_size PRIVATE Eigen3::Eigen lapacke)
