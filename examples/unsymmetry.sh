#!/bin/bash

awk -F ' ' ' BEGIN  { print  "Doing stuff." > "/dev/stderr"}
      /^%/ {print "Comment encountered. ", $0 > "/dev/stderr";
            next}
      NF == 4 { rows = $1 ; cols = $2; count = $3 }
      NF == 3 { if ($1 != $2) {
                 ++count; # adding a record.

                 # make a new entry!
                 # i.e. transpose the entry.
                 print $2, $1, $3;

             }
             }
      END { print "Rows = " , rows ", cols = ", cols, "count = ", count > "/dev/stderr"}' $1
