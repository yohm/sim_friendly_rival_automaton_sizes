#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <map>
#include <chrono>
#include <bitset>
#include <mpi.h>
#include <omp.h>
#include "icecream.hpp"
#include "caravan.hpp"
#include "DFA_translator.hpp"


using json = nlohmann::json;
using automaton_sizes_t = std::array<size_t, 65>;
size_t output_automaton_size_max = 0;

void ExpandWildcard(const std::string& line, automaton_sizes_t& automaton_sizes, std::vector<std::string>& found) {
  for (int i = 0; i < line.size(); i++) {
    if (line[i] == '*') {
      std::string temp = line;
      temp[i] = 'c';
      ExpandWildcard(temp, automaton_sizes, found);
      temp[i] = 'd';
      ExpandWildcard(temp, automaton_sizes, found);
      return;
    }
  }
  int size = DFA_translator::MinimizedDFASize(line.c_str());
  automaton_sizes[size]++;
  if (size <= output_automaton_size_max) {
    #pragma omp critical
    {
      found.emplace_back(line);
    }
  }
}

automaton_sizes_t AutomatonSizes(const std::string& line, std::vector<std::string>& found) {
  automaton_sizes_t sizes;
  sizes.fill(0);

  ExpandWildcard(line, sizes, found);
  return sizes;
}

void AddExpandedLines(const std::string &line, std::vector<std::string>& lines, int level) {
  if (level == 0) {
    lines.push_back(line);
    return;
  }
  for (int i = 0; i < line.size(); i++) {
    if (line[i] == '*') {
      std::string temp = line;
      temp[i] = 'c';
      AddExpandedLines(temp, lines, level-1);
      temp[i] = 'd';
      AddExpandedLines(temp, lines, level-1);
      return;
    }
  }
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  int my_rank, mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == 0) {
    int num_threads = omp_get_max_threads();
    std::cerr << "mpi_size = " << mpi_size << std::endl;
    std::cerr << "num_threads = " << num_threads << std::endl;
  }

  if (argc != 3) {
    if (my_rank == 0) {
      std::cerr << "Usage: " << argv[0] << " <input_file> <output_automaton_size_max>" << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }

  std::ifstream fin;
  std::ofstream fout;
  if (my_rank == 0) {
    fin.open(argv[1]);
    if (!fin) {
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    fout.open("automatons");
  }
  output_automaton_size_max = std::stoul(argv[2]);

  std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

  size_t total_num_lines_read = 0;
  automaton_sizes_t total_sizes;
  total_sizes.fill(0);

  auto read_lines_and_push_tasks = [&fin,&start,&total_num_lines_read,mpi_size](caravan::Queue& q) {
    while(fin && q.Size() < mpi_size) {
      json lines;
      size_t num_strategies = 0ul;
      const size_t max_line_size = 8192;
      const size_t max_num_strategies = 1ul << 22;

      std::string line;
      while (std::getline(fin, line)) {
        if (line.empty()) break;
        total_num_lines_read += 1;
        if (total_num_lines_read % 1000000 == 0) {
          std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
          std::chrono::duration<double> elapsed_seconds = end - start;
          std::cerr << "total_num_lines_read:" << total_num_lines_read << ' ' << elapsed_seconds.count() << std::endl;
        }

        lines.emplace_back(line);
        size_t n = 1ul;
        for (size_t i = 0; i < line.size(); i++) {
          if (line[i] == '*') { n = n << 1; }
        }
        num_strategies += n;
        if (lines.size() >= max_line_size || num_strategies >= max_num_strategies) {
          break;
        }
      }
      // std::cerr << num_strategies << ' ' << lines.size() << std::endl;
      q.Push(lines);
    }
  };

  // define a pre-process: create json object that contains parameters of tasks
  // This function is called only at the master process.
  std::function<void(caravan::Queue&)> on_init = [&fin,&read_lines_and_push_tasks](caravan::Queue& q) {
    read_lines_and_push_tasks(q);
  };

  // After the task was executed at a worker process, its result is returned to the master process.
  // When the master process receives the result, this callback function is called at the master process.
  std::function<void(int64_t, const json&, const json&, caravan::Queue&)> on_result_receive = [&fin,&read_lines_and_push_tasks,&total_sizes,&fout](int64_t task_id, const json& input, const json& output, caravan::Queue& q) {
    auto v = output["sizes"].get<std::vector<size_t>>();
    for (size_t i = 0; i < v.size(); i++) {
      total_sizes[i] += v[i];
    }
    if (!output["found"].empty()) {
      for (auto& line : output["found"]) {
        fout << line.get<std::string>() << std::endl;
      }
    }
    if (fin && q.Size() == 0) {
      read_lines_and_push_tasks(q);
    }
  };

  // Define the function which is executed at a worker process.
  // The input parameter for the task is given as the argument.
  std::function<json(const json& input)> do_task = [](const json& input) {
    std::vector<std::string> lines;
    for (const json& j : input) {
      std::string line = j.get<std::string>();
      int n_ast = 0;
      for (size_t i = 0; i < line.size(); i++) {
        if (line[i] == '*') { n_ast++; }
      }

      int max = 18;
      if (n_ast > max) {
        AddExpandedLines(line, lines, n_ast-max);
        // IC(line, n_ast);
      }
      else {
        lines.push_back(line);
      }
    }

    automaton_sizes_t sizes;
    sizes.fill(0);
    std::vector<std::string> found;
    #pragma omp parallel for shared(lines,sizes) schedule(dynamic,1)
    for (size_t i = 0; i < lines.size(); i++) {
      automaton_sizes_t r = AutomatonSizes(lines[i], found);
      for (size_t i = 0; i < r.size(); ++i) {
        #pragma omp atomic
        sizes[i] += r[i];
      }
    }

    json output;
    output["sizes"] = sizes;
    output["found"] = found;
    return output;
  };


  caravan::Start(on_init, on_result_receive, do_task, MPI_COMM_WORLD);

  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;

  if (my_rank == 0) {
    size_t sum = 0;
    double total = 0.0;
    for (int i = 0; i < 65; i++) {
      std::cout << i << ' ' << total_sizes[i] << std::endl;
      sum += total_sizes[i];
      total += i * total_sizes[i];
    }
    IC(sum, total/sum, elapsed_seconds.count(), elapsed_seconds.count()/sum);
  }

  MPI_Finalize();
  return 0;
}
