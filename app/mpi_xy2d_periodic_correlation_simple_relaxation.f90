program mpi_xy2d_periodic_correlation_simple_simulation
  use, intrinsic :: iso_fortran_env
  use mpi
  use gf2xe
  use msmt19937
  use xy2d_periodic_correlation_simple_m
  use variance_kahan_m
  use mpi_variance_kahan_m
  implicit none
  integer(int32), parameter :: mcs = 10000_int32
  integer(int32), parameter :: nsample = 2_int32
  integer(int64), parameter :: iseed = 42_int64
  integer(int32), parameter :: n_skip = 0_int32
  type(variance_kahan) :: order_parameter(mcs)
  type(variance_kahan), allocatable :: correlations(:, :)
  integer(int32) :: expo
  integer(int32) :: i, j
  integer(int32) :: myrank, num_proc, ierr
  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, num_proc, ierr)

  call init_genrand(iseed)
  !> Skip random numbers. (num_proc * n_skip + myrank) * 2^e
  !> 2^e must be larger than (nx * ny * mcs * nsample).
  expo = ceiling(log(real(2 * nx * ny * (mcs + 1) * nsample, real64)) / log(2.0d0)) + 1
  if (num_proc * n_skip + myrank /= 0) &
       & call mt_jumpahead(num_proc * n_skip + myrank, expo)

  call init_simulation_xy2d()

  order_parameter(:) = variance_kahan()
  allocate(correlations(mcs, n_unique_correlation), source = variance_kahan())

  if (myrank == 0) then
     write(output_unit, '(a,i0)'    ) "# Nsize: ", nall
     write(output_unit, '(2(a, i0))') "# nx: ", nx, " ny: ", ny
     write(output_unit, '(2(a, i0))') "# MCS: ", mcs, " Nsample: ", nsample
     write(output_unit, '(a, g0)' ) "# 温度: ", kbt
     write(output_unit, '(a)' ) "# method: Metropolis, model: 2D-XY"
     write(output_unit, '(a, i0)' ) "# the number of processors: ", num_proc
     write(output_unit, '(a, i0)' ) "# the number unique dist_sq: ", n_unique_correlation

     write(error_unit, '(a,i0)'    ) "# Nsize: ", nall
     write(error_unit, '(2(a, i0))') "# nx: ", nx, " ny: ", ny
     write(error_unit, '(2(a, i0))') "# MCS: ", mcs, " Nsample: ", nsample
     write(error_unit, '(a, g0)' ) "# 温度: ", kbt
     write(error_unit, '(a)' ) "# method: Metropolis, model: 2D-XY"
     write(error_unit, '(a, i0)' ) "# the number of processors: ", num_proc
     write(error_unit, '(a, i0)' ) "# the number unique dist_sq: ", n_unique_correlation
  end if
  do j = 1, nsample
     if (myrank == 0) &
          write(error_unit, '(a, i0)') "sample: ", j
     call init_lattice_order()
     do i = 1, mcs
        call update_metropolis()
        associate(m => calc_magne())
          call order_parameter(i)%add_data(m)
        end associate
        block
          real(real64) :: corr(n_unique_correlation)
          integer(int32) :: k
          call calc_correlation(corr)
          do k = 1, n_unique_correlation
             call correlations(i, k)%add_data(corr(k))
          end do
        end block
     end do
  end do
  block
    type(variance_kahan) :: all_order_params(mcs)
    type(variance_kahan) :: all_correlations(mcs, n_unique_correlation)
    call vk_mpi_multi_gather(mcs, order_parameter, all_order_params, 0, myrank, num_proc, ierr)
    do j = 1, n_unique_correlation
       call vk_mpi_multi_gather(mcs, correlations(:, j), all_correlations(:, j), 0, myrank, num_proc, ierr)
    end do
    if (myrank == 0) then
       write(output_unit, '(a)', advance = "no") "# Nsize, Nsample, mcs, <m>, <e>, <m^2>, <e^2>, χ, C, m'"
       do j = 1, n_unique_correlation
          write(output_unit, '(2(a, es10.4, a))', advance = "no")&
               &  ", G(", unique_dist(j), ", t)", &
               & ", G(", unique_dist(j), ", t)^2"
       end do
       write(output_unit, *)
       do i = 1, mcs
          write(output_unit, '(*(g0, 1x))', advance = "no") nall, nsample * num_proc, i, &
               & all_order_params(i)%mean(), 0, &
               & all_order_params(i)%square_mean(), 0, &
               & 0, 0, 0
          do j = 1, n_unique_correlation
             write(output_unit, '(2(1x, g0))', advance = "no") correlations(i, j)%mean(), correlations(i, j)%square_mean()
          end do
          write(output_unit, *)
       end do
    end if
  end block
  call MPI_Finalize(ierr)
end program mpi_xy2d_periodic_correlation_simple_simulation
