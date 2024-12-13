#!/bin/bash

set -x -u -e

FCFLAGS="-O2"
# FCFLAGS="${FCFLAGS} -g -Mbounds -Minfo=accel -gpu=debug"
nx=14
ny=${nx}
nz=${nx}
mcs_discard=10000
mcs_sample=100000
kbt_begin=3.0
kbt_end=0.1
kbt_nsplit=30
iseed=42

root_dir="./root_$$"
srcfile="./src/xy3d_periodic_m.f90"
progfile="./app/xy3d_periodic_metropolis_equilibrium.f90"
execfile="${root_dir}/bin/xy3d_periodic_metropolis_equilibrium"
outputdir="data/3D-XY_equilibrium"
outputfile="${outputdir}/3D-XY_Metropolis_x${nx}_y${ny}_z${nz}_discard${mcs_discard}_sample${mcs_sample}_kbt${kbt_begin}-${kbt_end}_n${kbt_nsplit}_iseed${iseed}_$(date +%Y%m%d_%H%M%S).dat"
mkdir -p ${outputdir}
tmpdir="data/tmp"
mkdir -p "${tmpdir}"
tmpfile="$(mktemp --tmpdir=${tmpdir})"

fpm clean --skip

sed -i -e "s/mcs_discard = [0-9]*_int32/mcs_discard = ${mcs_discard}_int32/" \
    -e "s/mcs_sample = [0-9]*_int32/mcs_sample = ${mcs_sample}_int32/" \
    -e "s/iseed = [0-9]*_int64/iseed = ${iseed}_int64/" \
    "${progfile}"

start=$(date +%s)
for i in $(seq ${kbt_nsplit})
do
    kbt=$(python3 -c "print(((${kbt_nsplit} - ${i}) * ${kbt_begin} + (${i} - 1) * ${kbt_end}) / (${kbt_nsplit} - 1))")
    echo ${kbt}
    sed -i -e "s/nx = [0-9]*_int64/nx = ${nx}_int64/" \
        -e "s/kbt = [0-9.]*d0/kbt = ${kbt}d0/" \
        "${srcfile}"
    time fpm install $(basename ${execfile}) --prefix="${root_dir}" --compiler="mpif90" --verbose --flag="${FCFLAGS}"
    ${execfile} >> "${tmpfile}"
done
end=$(date +%s)
mv -v "${tmpfile}" "${outputfile}"
echo "output >>> '${outputfile}'"
chmod 400 "${outputfile}"

hour=$( echo "($end - $start) / 3600" | bc)
minute=$( echo "(($end - $start) % 3600) / 60" | bc)
second=$( echo "($end - $start) % 60" | bc)
minute=$(printf "%02d" "${minute}")
second=$(printf "%02d" "${second}")
elapsed_time="${hour}h ${minute}m ${second}s"
#     model,»  size,»     sample,»     mcs,»   kbt, time
echo "3D-XY_MET,  ${nx}x${ny}x${nz} == $((nx*ny*nz)), discard ${mcs_discard}, sample ${mcs_sample}, T: (${kbt_begin}-${kbt_end}_${kbt_nsplit}), iseed ${iseed}, time ${elapsed_time}, ${outputfile}" | tee -a xy3d_equilibrium.log

rm -rf "${root_dir:?}"
