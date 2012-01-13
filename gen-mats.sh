#!/bin/bash

./time-and-log.sh ./src/genmat 10000 300 0.01 &
./time-and-log.sh ./src/genmat 10000 30 0.001 &
./time-and-log.sh ./src/genmat 10000 150 0.005 &
./time-and-log.sh ./src/genmat 40000 1200 0.01 &
./time-and-log.sh ./src/genmat 40000 120 0.001 &
./time-and-log.sh ./src/genmat 40000 600 0.005 &
./time-and-log.sh ./src/genmat 80000 2400 0.01 &
./time-and-log.sh ./src/genmat 80000 240 0.001 &
./time-and-log.sh ./src/genmat 80000 1200 0.005 &
