#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <sstream>
#include <tuple>
#include "icecream.hpp"
#include "DFA_translator.hpp"


template <class T>
std::string Join(const std::vector<T>& vec, const std::string& delim) {
  std::stringstream ss;
  for (size_t i = 0; i < vec.size(); i++) {
    ss << vec[i];
    if (i != vec.size()-1) {
      ss << delim;
    }
  }
  return ss.str();
}

int main(int argc, char* argv[]) {

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " [files...]" << std::endl;
    return 1;
  }

  std::map<size_t, std::vector<std::string>> automatons;

  for (int i = 1; i < argc; i++) {
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

  std::map<std::string, std::pair<std::string,size_t>> serialized_automatons; // {serialized, [strategy, frequency]}

  for (auto& pair : automatons) {
    size_t s = pair.first;
    std::sort(pair.second.begin(), pair.second.end());
    for (auto& line : pair.second) {
      std::cout << s << " " << line << " ";
      auto m = DFA_translator::MinimizedPartition(line.c_str());
      std::vector<std::string> v;
      for (auto& p : m) {
        std::string s = Join(p.second, ",");
        v.emplace_back(s);
      }
      std::cout << Join(v, ";") << " ";

      auto m2 = DFA_translator::MinimizedPartitionSimple(line.c_str());
      std::vector<std::string> v2;
      for (auto& p : m2) {
        std::string s = Join(p.second, ",");
        v2.emplace_back(s);
      }
      std::cout << Join(v2, ";") << std::endl;

      if (m2.size() < 6) {
        std::string serialized = DFA_translator::SerializeSimpleAutom(line.c_str());
        if (serialized_automatons.find(serialized) == serialized_automatons.end()) {
          serialized_automatons[serialized] = std::make_pair(line, 0);
        }
        serialized_automatons[serialized].second++;
      }
    }
  }

  // sort by number of states
  using serialized_strategy_size_t = std::tuple<std::string, std::string, size_t>;
  std::vector<serialized_strategy_size_t> sorted_automatons;
  for (auto& pair : serialized_automatons) {
    sorted_automatons.emplace_back(pair.first, pair.second.first, pair.second.second);
  }
  std::sort(sorted_automatons.begin(), sorted_automatons.end(),
            [](const serialized_strategy_size_t& a, const serialized_strategy_size_t & b) {
              return std::get<2>(a) > std::get<2>(b);
            });
  for (auto& tuple : sorted_automatons) {
    std::cerr << std::get<0>(tuple) << " " << std::get<1>(tuple) << " " << std::get<2>(tuple) << std::endl;
  }
  // IC(serialized_automatons, serialized_automatons.size());

  return 0;
}