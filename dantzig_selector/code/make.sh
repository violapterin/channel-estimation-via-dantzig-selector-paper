#!/usr/bin/env bash

OBJS="main.o functions.o constants.o channel.cpp oommpp.cpp"
CPPS="src/main.cpp src/functions.cpp src/const.cpp src/channel.cpp src/oommpp.cpp src/ddss.cpp"
CC="g++"
CFLAGS="-I lib/ -Wall -g -O3 -o"
BIN="main"

eval "${CC}" "${CFLAGS}" "${BIN}" "${CPPS}"
