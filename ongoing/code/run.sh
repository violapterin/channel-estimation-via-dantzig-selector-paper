#!/usr/bin/env bash

#lldb ./bin/main.o ./dat/table.txt
rm -f ./dat/table.txt ./plt/fig.png
./bin/main.o ./dat/table.txt
./src/read_and_draw.py ./dat/table.txt 

