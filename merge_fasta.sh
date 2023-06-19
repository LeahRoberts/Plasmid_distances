#!/bin/bash

awk '
    /^>/ { 
        # print the first header
        if (c++ == 0) {print; print ""} 
        next
    } 
    /^$/ {next} 
    {printf "%s", $0} 
    END {print ""}
' $1 > $2
