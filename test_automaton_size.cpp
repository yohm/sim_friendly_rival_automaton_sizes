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
  auto partition = DFA_translator::MinimizedPartition(str);
  IC(partition);

  StrategyM3 strategy(str);
  auto partition2 = strategy.MinimizeDFA(true).to_map();
  IC(partition2);

  if (partition != partition2) {
    std::cerr << "partition != partition2" << std::endl;
    throw std::runtime_error("partition != partition2");
  }
}


int main(int argc, char* argv[]) {

  if (argc == 1) {
    const std::string line = "cddcccddcdcdddddcccddddcdddcddcdcccccdddddcddcddddccdddddcdccccd";

    PerformanceTest(line.c_str());
    EqualityTest(line.c_str());
  }
  else if (argc == 2) {
    const std::string line = argv[1];
    auto p = DFA_translator::MinimizedPartition(line.c_str());
    IC(p);
    for (auto pair: p) {
      std::cout << pair.first << ",";
    }
    std::cout << std::endl;
  }
  else {
    std::cerr << "Usage: " << argv[0] << " [line]" << std::endl;
    return 1;
  }

  return 0;
}