# sim_friendly_rival_automaton_sizes

A sample to run DFT minimization on MPI.

Compilation and execution procedure:

1. mpicxx main.cpp -o main.out
2. mpiexec -np 4 main.out

The program assumes that the input is given in the file "input_file_%04d".

Function `GetAutomatonSize` is not implemented. Please implement by yourself.
