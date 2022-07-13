#include <iostream>
#include <vector>
#include <array>
#include <bitset>
#include <chrono>
#include <cassert>
#include "icecream.hpp"
#include "DFA_translator.hpp"
#include "StrategyM3.hpp"


void PerformanceTest(const char str[64]) {
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  int times = 100'000;
  for (int i = 0; i < times; i++) {
    DFA_translator::MinimizedDFASize(str);
  }
  end = std::chrono::system_clock::now();
  std::cerr << "MinimizedDFASize(): " << std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() / (double)times << " micro sec" << std::endl;

}

void EqualityTest(const char str[64]) {
  auto partition1 = DFA_translator::MinimizedPartition(str);
  IC(partition1);

  StrategyM3 strategy(str);
  auto partition2 = strategy.MinimizeDFA(true).to_map();
  IC(partition2);

  if (partition1 != partition2) {
    std::cerr << "partition1 != partition2" << std::endl;
    throw std::runtime_error("partition1 != partition2");
  }

  auto partition3 = DFA_translator::MinimizedPartitionSimple(str);
  IC(partition3);

  auto partition4 = strategy.MinimizeDFA(false).to_map();
  IC(partition4);

  if (partition3 != partition4) {
    std::cerr << "partition3 != partition4" << std::endl;
    throw std::runtime_error("partition3 != partition4");
  }
}


int main(int argc, char* argv[]) {

  if (argc == 1) {
    const std::string line = "cddcccddcdcdddddcccddddcdddcddcdcccccdddddcddcddddccdddddcdccccd";

    PerformanceTest(line.c_str());
    EqualityTest(line.c_str());
  }
  else if (argc == 2 || argc == 3) {
    bool is_full = true;
    if (argc == 3) {
      is_full = std::stoi(argv[2]) == 0;
    }
    const std::string line = argv[1];
    auto p = is_full ? DFA_translator::MinimizedPartition(line.c_str()) : DFA_translator::MinimizedPartitionSimple(line.c_str());
    for (auto pair: p) {
      for (auto n: pair.second) {
        std::cout << n << ' ';
      }
      std::cout << std::endl;
    }
  }
  else {
    std::cerr << "Usage: " << argv[0] << " [line] [full:0, simplified:1]" << std::endl;
    return 1;
  }

  return 0;
}