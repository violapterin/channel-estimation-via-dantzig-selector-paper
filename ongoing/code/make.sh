#!/usr/bin/env bash

#OBJS="main.o functions.o constants.o"
#CC="g++"
#CFLAGS="-Wall -g"
#eval "${CC}" "${CFLAGS}" -o "${OBJS}"

g++ -Wall -O3 -o ./bin/main.o ./src/main.cpp ./src/functions.cpp ./src/constants.cpp 
#g++ -Wall -g -o ./bin/main.o ./src/main.cpp ./src/functions.cpp ./src/constants.cpp 


