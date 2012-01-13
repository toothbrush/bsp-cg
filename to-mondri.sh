#!/bin/bash

for p in 1 2 4 8;
do
    for i in *.emm;
    do
        echo $i
        Mondriaan $i $p 0.3
    done
done
