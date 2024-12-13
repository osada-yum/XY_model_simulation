#!/bin/bash

sed -e "s/\(Z'.*'\)/INT(\1, INT64)/" \
    -e 's/IAND(tmp,ISHFT/IAND(int(tmp, int64),ISHFT/' \
    -e 's/(b,Z8)/(INT(b, INT64),Z8)/' \
    -e 's/(b,ZC)/(INT(b, INT64),ZC)/' \
    -e 's/(b,ZE)/(INT(b, INT64),ZE)/' \
    msmt19937.f90 > msmt19937_modified.f90

sed -n 1,1156p msmt19937_modified.f90 > src/gf2xe.f90
sed -n 1162,1720p msmt19937_modified.f90 > src/msmt19937.f90
