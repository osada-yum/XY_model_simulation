program xy3d_periodic_wolff_equilibrium_simulation
  use, intrinsic :: iso_fortran_env
  use gf2xe
  use msmt19937
  use xy3d_periodic_m
  use variance_kahan_m
  implicit none
  integer(int32), parameter :: mcs_discard = 10000_int32, mcs_sample = 100000_int32
  integer(int64), parameter :: iseed = 42_int64
  type(variance_kahan) :: magnes, energies
  integer(int32) :: i, j

  call init_simulation_xy3d()
  call init_genrand(iseed)

  associate(out => [output_unit, error_unit])
    do i = 1, size(out)
       write(out(i), '(a,i0)'    ) "# Nsize: ", nall
       write(out(i), '(3(a, i0))') "# nx: ", nx, " ny: ", ny, " nz: ", nz
       write(out(i), '(2(a, i0))') "# MCS_discard: ", mcs_discard, " MCS_sample: ", mcs_sample
       write(out(i), '(a, g0)' ) "# 温度: ", kbt
       write(out(i), '(a)' ) "# method: Wolff, model: 3D-XY"
       write(out(i), '(a)') "# Nsize, Nsample, kbt, <m>, <e>, <m^2>, <e^2>, χ, C, m'"
    end do
  end associate
  call init_lattice_order()
  write(error_unit, '(a, i0)') "Discard: ", mcs_discard
  do j = 1, mcs_discard
     call update_wolff()
  end do
  write(error_unit, '(a, i0)') "Sampling: ", mcs_sample
  magnes = variance_kahan()
  energies = variance_kahan()
  do j = 1, mcs_sample
     call update_wolff()
     associate(m => calc_magne_abs(), e => calc_energy())
       call magnes%add_data(m)
       call energies%add_data(e)
     end associate
  end do
  write(output_unit, '(*(g0, 1x))') nall, mcs_sample, kbt, &
       & magnes%mean(), energies%mean(), &
       & magnes%square_mean(), energies%square_mean(), &
       & magnes%var(), energies%var(), 0
end program xy3d_periodic_wolff_equilibrium_simulation
