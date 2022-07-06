#!/bin/bash -eux

mpiFCCpx -Nclang -std=c++14 -Kfast -DNDEBUG -Iicecream-cpp -Icaravan-lib -Icaravan-lib/json/include -o main.out main.cpp

