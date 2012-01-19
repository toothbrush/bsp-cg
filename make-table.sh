grep csv_answer_data *-P1.* | awk -F ':' ' {print $3} '| awk -F ',' '{ print $2 "  &     &   " $3 "   &   " $3/($2*$2)    "\\\\" }'|sort -u | sort -n -k 4
