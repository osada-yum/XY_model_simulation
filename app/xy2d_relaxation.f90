program xy2d_simulation
  use, intrinsic :: iso_fortran_env
  use xy2d_m
  use variance_covariance_kahan_m
  implicit none
  integer(int32), parameter :: mcs = 1000, nsample = 10
  integer(int64), parameter :: nx = 101, ny = 100
  real(real64), parameter :: kbt = 0.90_real64, n_inv_r64 = 1 / real(nx * ny, real64)
  type(xy2d) :: system
  type(variance_covariance_kahan) :: order_parameter(mcs)
  integer(int32) :: i, j
  call system%init(nx, ny, kbt)
  write(output_unit, '(a,i0)'    ) "# Nsize: ", system%nall()
  write(output_unit, '(2(a, i0))') "# nx: ", system%nx(), " ny: ", system%ny()
  write(output_unit, '(2(a, i0))') "# MCS: ", mcs, " Nsample: ", nsample
  write(output_unit , '(a, g0)' ) "# 温度: ", system%kbt()
  write(output_unit , '(a)' ) "# method: METROPOLIS"
  do j = 1, nsample
     write(error_unit, '(a, i0)') "sample: ", j
     call system%set_xy2d_allup()
     do i = 1, mcs
        call system%update()
        associate(m => system%calc_magne_summ(), &
             & e => system%calc_energy_summ())
          call order_parameter(i)%add_data(m * n_inv_r64, e * n_inv_r64)
        end associate
     end do
  end do
  write(output_unit, '(a)') "# Nsize, Nsample, mcs, <m>, <e>, χ, C, m'"
  do i = 1, mcs
     write(output_unit, '(*(g0, 1x))') system%nall(), order_parameter(i)%num_sample(), i, &
          & order_parameter(i)%mean1(), order_parameter(i)%mean2(), &
          & order_parameter(i)%square_mean1(), order_parameter(i)%square_mean2(), &
          & system%nall() * order_parameter(i)%var1(), system%nall() * order_parameter(i)%var2(), &
          & system%nall() * order_parameter(i)%cov()
  end do
end program xy2d_simulation
