#!/bin/bash -eux

FCCpx -Nclang -std=c++14 -Kfast -DNDEBUG minimum.cpp -o minimum.out
mpiFCCpx -Nclang -std=c++14 -Kfast -Kopenmp -DNDEBUG -Iicecream-cpp -Icaravan-lib -Icaravan-lib/json/include -o main.out main.cpp

