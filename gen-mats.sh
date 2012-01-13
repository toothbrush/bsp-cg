#!/bin/bash

./time-and-log.sh ./src/genmat 1000  150 0.01 &
./time-and-log.sh ./src/genmat 1000  15 0.001 &
./time-and-log.sh ./src/genmat 1000  75  0.005 &
./time-and-log.sh ./src/genmat 2000  150 0.01 &
./time-and-log.sh ./src/genmat 2000  15 0.001 &
./time-and-log.sh ./src/genmat 2000  75  0.005 &
./time-and-log.sh ./src/genmat 5000  150 0.01 &
./time-and-log.sh ./src/genmat 5000  15 0.001 &
./time-and-log.sh ./src/genmat 5000  75  0.005 &
./time-and-log.sh ./src/genmat 10000 300 0.01 &
./time-and-log.sh ./src/genmat 10000 30 0.001 &
./time-and-log.sh ./src/genmat 10000 150 0.005 &
./time-and-log.sh ./src/genmat 20000 1200 0.01 &
./time-and-log.sh ./src/genmat 20000 120 0.001 &
./time-and-log.sh ./src/genmat 20000 600 0.005 &
