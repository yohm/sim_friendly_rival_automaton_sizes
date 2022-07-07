# sim_friendly_rival_automaton_sizes

A sample to run DFT minimization in parallel using MPI and OpenMP. Clone this repository with submodules.

```
git submodule --recursive git@github.com:yohm/sim_friendly_rival_automaton_sizes.git
```

Compilation and execution procedure:

1. mpicxx -fopenmp -O3 main.cpp -o main.out
2. mpiexec -np 4 main.out strategy_list

The program assumes that the input is given in the file "strategy_list".
