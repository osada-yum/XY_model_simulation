#!/bin/bash

set -v -e

filename="xy2d_metropolis_overrelaxation.f90"
exe="${filename//.f90}.out"

dir="relaxation_metropolis_overrelaxation_kbts_data"
mkdir -p "${dir}"

size=2000
kbts=(0.5 0.8 1.0 1.2 2.0)

# size=100, theta
# real    2m48.067s
# user    11m25.832s
# sys     0m2.358s
# size=100, hypot
# real    1m49.958s
# user    7m40.327s
# sys     0m2.710s

for kbt in ${kbts[@]}
do
    outputfile="${dir}/xy2d_metropolis_overrelaxation_size${size}_kbt${kbt}.dat"
    sed -i \
        -e "s/nx = [0-9]*,/nx = ${size},/" \
        -e "s/kbt = [0-9.]*d0/kbt = ${kbt}d0/" \
        "${filename}"
    gfortran -O3 "${filename}" -o "${exe}"
    ./"${exe}" > "${outputfile}" &
done

wait

echo "done!"
