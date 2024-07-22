program mpi_xy3d_periodic_metropolis_simulation
  use, intrinsic :: iso_fortran_env
  use mpi
  use gf2xe
  use msmt19937
  use xy3d_periodic_m
  use variance_kahan_m
  use mpi_variance_kahan_m
  implicit none
  integer(int32), parameter :: mcs = 1000_int32
  integer(int32), parameter :: nsample = 2_int32
  integer(int64), parameter :: iseed = 42_int64
  integer(int32), parameter :: n_skip = 0_int32
  type(variance_kahan) :: magnes(mcs), energies(mcs)
  integer(int32) :: expo
  integer(int32) :: i, j
  integer(int32) :: myrank, num_proc, ierr
  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, num_proc, ierr)

  call init_genrand(iseed)
  !> Skip random numbers. (num_proc * n_skip + myrank) * 2^e
  !> 2^e must be larger than (nx * ny * nz * mcs * nsample).
  expo = ceiling(log(real(2 * nx * ny * nz * (mcs + 1) * nsample, real64)) / log(2.0d0)) + 1
  if (num_proc * n_skip + myrank /= 0) &
       & call mt_jumpahead(num_proc * n_skip + myrank, expo)

  call init_simulation_xy3d()

  magnes(:) = variance_kahan()
  energies(:) = variance_kahan()

  if (myrank == 0) then
     associate(out => [output_unit, error_unit])
       do i = 1, size(out)
          write(out(i), '(a,i0)'    ) "# Nsize: ", nall
          write(out(i), '(3(a, i0))') "# nx: ", nx, " ny: ", ny, " nz: ", nz
          write(out(i), '(2(a, i0))') "# MCS: ", mcs, " Nsample: ", nsample
          write(out(i), '(a, g0)' ) "# 温度: ", kbt
          write(out(i), '(a)' ) "# method: Metropolis, model: 3D-XY"
          write(out(i), '(a, i0)' ) "# the number of processors: ", num_proc
       end do
     end associate
  end if
  do j = 1, nsample
     if (myrank == 0) &
          write(error_unit, '(a, i0)') "sample: ", j
     call init_lattice_order()
     do i = 1, mcs
        call update_metropolis()
        associate(m => calc_magne(), e => calc_energy())
          call magnes(i)%add_data(m)
          call energies(i)%add_data(e)
        end associate
     end do
  end do
  block
    type(variance_kahan) :: all_magnes(mcs), all_energies(mcs)
    call vk_mpi_multi_gather(mcs, magnes, all_magnes, 0, myrank, num_proc, ierr)
    call vk_mpi_multi_gather(mcs, energies, all_energies, 0, myrank, num_proc, ierr)
    if (myrank == 0) then
       write(output_unit, '(a)', advance = "no") "# Nsize, Nsample, mcs, <m>, <e>, <m^2>, <e^2>, χ, C, m'"
       write(output_unit, *)
       do i = 1, mcs
          write(output_unit, '(*(g0, 1x))', advance = "no") nall, nsample * num_proc, i, &
               & all_magnes(i)%mean(), all_energies(i)%mean(), &
               & all_magnes(i)%square_mean(), all_energies(i)%square_mean(), &
               & 0, 0, 0
          write(output_unit, *)
       end do
    end if
  end block
  call MPI_Finalize(ierr)
end program mpi_xy3d_periodic_metropolis_simulation
