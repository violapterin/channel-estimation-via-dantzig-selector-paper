#!/usr/bin/env bash

g++ -o ./bin/main ./src/main.cpp ./src/functions.cpp ./src/constants.cpp -I/usr/include/python2.7 -lpython2.7
#g++ -o ./bin/main.o ./src/main.cpp -std=c++11 -I /usr/include/python2.7 -l python2.7

