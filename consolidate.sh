#!/bin/bash
#
# grab results from output files and write to stdout
grep csv_answer_data * | awk -F ':' ' { print  } '|tr -d '[:blank:]'
