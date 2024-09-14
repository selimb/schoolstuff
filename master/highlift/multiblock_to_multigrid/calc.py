#!/usr/bin python2
import sys
sub = int(sys.argv[1])
lev = 1
print sub
while ((sub - 1) % 2) == 0 and sub > 2:
    sub = (sub - 1)/2 + 1
    lev += 1
    print sub

