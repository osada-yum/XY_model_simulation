#!/bin/bash

set -v -e

filename="xy2d_wolff_relaxation.f90"
exe="${filename//.f90}.out"

dir="relaxation_wolff_kbts_data"
mkdir -p "${dir}"

size=500
kbts=(0.5 0.8 1.0 1.2 2.0)

for kbt in ${kbts[@]}
do
    outputfile="${dir}/xy2d_wolff_relaxation_size${size}_kbt${kbt}.dat"
    sed -i \
        -e "s/nx = [0-9]*,/nx = ${size},/" \
        -e "s/kbt = [0-9.]*d0/kbt = ${kbt}d0/" \
        "${filename}"
    gfortran -O3 "${filename}" -o "${exe}"
    ./"${exe}" > "${outputfile}"&
done

wait

echo "done!"
