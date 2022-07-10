#include <iostream>
#include <fstream>
#include <array>
#include "icecream.hpp"
#include "DFA_translator.hpp"


int main(int argc, char* argv[]) {

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " [files...]" << std::endl;
    return 1;
  }

  std::map<size_t, std::vector<std::string>> automatons;

  for (int i = 1; i < argc; i++) {
    std::cerr << "reading " << argv[i] << std::endl;
    std::string filename = argv[i];
    std::ifstream fin(filename);
    while (fin) {
      std::string line;
      std::getline(fin, line);
      if (line.size() > 0) {
        size_t size = DFA_translator::MinimizedDFASize(line.c_str());
        automatons[size].emplace_back(line);
      }
    }
  }

  for (auto& pair : automatons) {
    size_t s = pair.first;
    std::sort(pair.second.begin(), pair.second.end());
    for (auto& line : pair.second) {
      std::cout << s << " " << line << " ";
      auto m = DFA_translator::MinimizedPartition(line.c_str());
      for (auto& p : m) {
        for (auto& n : p.second) {
          std::cout << n << ' ';
        }
        std::cout << "/ ";
      }
      std::cout << std::endl;
    }
  }

  return 0;
}