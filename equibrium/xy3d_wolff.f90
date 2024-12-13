module xy3d_periodic_m
  use, intrinsic :: iso_fortran_env
  use msmt19937
  implicit none
  private
  real(real64), parameter :: pi = 4 * atan(1d0)
  integer(int64), parameter, public :: nx = 14_int64, ny = nx, nz = nx, nall = nx * ny * nz
  real(real64), parameter :: nall_inv = 1.0d0 / nall

  real(real64), allocatable, public, protected :: spins(:, :, :, :) !> spins(sx or sy, position of x, position of y position of z)
  real(real64), allocatable :: rnds(:, :, :, :) !> for Metropolis. rnds(candidate or update, position of x, position of y position of z)
  integer(int64), allocatable :: stack(:, :) !> for Wolff. stack(x or y or z, maximum size)
  logical, allocatable :: cluster(:, :, :) !> for Wolff. cluster(position of x, position of y, position of z)

  real(real64), parameter, public :: kbt = 0.1d0, beta = 1 / kbt

  public :: init_simulation_xy3d, init_lattice_order, update_metropolis, update_wolff, calc_energy, calc_magne, calc_magne_abs
contains
  !> init: Initialize variables. Allocation.
  impure subroutine init_simulation_xy3d()
    allocate(spins(2, nx, ny, nz))
    allocate(rnds(2, nx, ny, nz))
    allocate(stack(1:3, nall))
    allocate(cluster(nx, ny, nz))
  end subroutine init_simulation_xy3d
  !> init_lattice_order: Set the spins to order state.
  impure subroutine init_lattice_order()
    spins(1, :, :, :) = 1.0d0
    spins(2, :, :, :) = 0.0d0
  end subroutine init_lattice_order
  !> update_metropolis: Update the lattice with 1 MCS.
  impure subroutine update_metropolis()
    integer(int64) :: z, y, x
    integer(int64) :: i
    do z = 1, nz
       do y = 1, ny
          do x = 1, nx
             rnds(1, x, y, z) = grnd()
             rnds(2, x, y, z) = grnd()
          end do
       end do
    end do
    !> parity: odd
    do z = 1, nz
       do y = 1, ny
          do x = 1 + iand(y + z, b'1'), nx, 2
             call local_flip(x, y, z)
          end do
       end do
    end do
    !> parity: even
    do z = 1, nz
       do y = 1, ny
          do x = 2 - iand(y + z, b'1'), nx, 2
             call local_flip(x, y, z)
          end do
       end do
    end do
  contains
    !> local_flip: Update a spin in lattice.
    !> @param x An index on the x-axis of the spin.
    !> @param y An index on the y-axis of the spin.
    !> @param z An index on the z-axis of the spin.
    impure subroutine local_flip(x, y, z)
      integer(int64), intent(in) :: x, y, z
      real(real64) :: summ(1:2)
      real(real64) :: diff_center(1:2)
      real(real64) :: new_theta, new_spin(1:2)
      real(real64) :: delta_e
      integer(int64) :: rx, lx, uy, dy, tz, bz
      summ(1:2) = 0.0d0
      rx = x + 1; if (rx > nx) rx = 1
      uy = y + 1; if (uy > ny) uy = 1
      tz = z + 1; if (tz > nz) tz = 1
      lx = x - 1; if (lx <  1) lx = nx
      dy = y - 1; if (dy <  1) dy = ny
      bz = z - 1; if (bz <  1) bz = nz

      !> right
      summ(1:2) = summ(1:2) + spins(1:2, rx, y, z)
      !> left
      summ(1:2) = summ(1:2) + spins(1:2, lx, y, z)
      !> up
      summ(1:2) = summ(1:2) + spins(1:2, x, uy, z)
      !> down
      summ(1:2) = summ(1:2) + spins(1:2, x, dy, z)
      !> top
      summ(1:2) = summ(1:2) + spins(1:2, x, y, tz)
      !> bottom
      summ(1:2) = summ(1:2) + spins(1:2, x, y, bz)

      new_theta = 2 * pi * rnds(1, x, y, z)
      new_spin(1:2) = [cos(new_theta), sin(new_theta)]
      diff_center(1:2) = new_spin(1:2) - spins(1:2, x, y, z)
      delta_e = - sum(diff_center(1:2) * summ(1:2))
      if (delta_e > 0) then
         if (.not. (rnds(2, x, y, z) < exp(- beta * delta_e))) return
      end if
      !> delta_e <= 0 or rnds(2, x, y, z) < exp(- beta * delta_e).
      spins(1:2, x, y, z) = new_spin(1:2)
    end subroutine local_flip
  end subroutine update_metropolis

  !> update_wolff: Update the spin cluster by Wolff algorithm
  impure subroutine update_wolff()
    integer(int64) :: x, y, z
    integer(int64) :: near_x, near_y, near_z
    integer(int64) :: cxyz(1:3), s
    real(real64) :: reflection_theta, reflection_spin(1:2)
    real(real64) :: r_sigma_c
    x = nx * grnd() + 1
    y = ny * grnd() + 1
    z = nz * grnd() + 1
    cluster(:, :, :) = .false.
    s = 1
    stack(1:3, s) = [x, y, z]
    cluster(x, y, z) = .true.
    reflection_theta = 2 * pi * grnd()
    reflection_spin(1:2) = [cos(reflection_theta), sin(reflection_theta)]
    spins(1:2, x, y, z) = spins(1:2, x, y, z) - 2 * sum(spins(1:2, x, y, z) * reflection_spin(1:2)) * reflection_spin(1:2)
    do while (s > 0_int64)
       cxyz(1:3) = stack(1:3, s)
       s = s - 1
       r_sigma_c = sum(reflection_spin(1:2) * spins(1:2, cxyz(1), cxyz(2), cxyz(3)))
       !> right
       near_x = cxyz(1) + 1
       if (near_x > nx) near_x = 1
       if (.not. cluster(near_x, cxyz(2), cxyz(3))) then
          call wolff_local_update(stack, spins, cluster, s, reflection_spin, r_sigma_c, near_x, cxyz(2), cxyz(3))
       end if
       !> left
       near_x = cxyz(1) - 1
       if (near_x < 1) near_x = nx
       if (.not. cluster(near_x, cxyz(2), cxyz(3))) then
          call wolff_local_update(stack, spins, cluster, s, reflection_spin, r_sigma_c, near_x, cxyz(2), cxyz(3))
       end if
       !> up
       near_y = cxyz(2) + 1
       if (near_y > ny) near_y = 1
       if (.not. cluster(cxyz(1), near_y, cxyz(3))) then
          call wolff_local_update(stack, spins, cluster, s, reflection_spin, r_sigma_c, cxyz(1), near_y, cxyz(3))
       end if
       !> down
       near_y = cxyz(2) - 1
       if (near_y < 1) near_y = ny
       if (.not. cluster(cxyz(1), near_y, cxyz(3))) then
          call wolff_local_update(stack, spins, cluster, s, reflection_spin, r_sigma_c, cxyz(1), near_y, cxyz(3))
       end if
       !> top
       near_z = cxyz(3) + 1
       if (near_z > nz) near_z = 1
       if (.not. cluster(cxyz(1), cxyz(2), near_z)) then
          call wolff_local_update(stack, spins, cluster, s, reflection_spin, r_sigma_c, cxyz(1), cxyz(2), near_z)
       end if
       !> bottom
       near_z = cxyz(3) - 1
       if (near_z < 1) near_z = nz
       if (.not. cluster(cxyz(1), cxyz(2), near_z)) then
          call wolff_local_update(stack, spins, cluster, s, reflection_spin, r_sigma_c, cxyz(1), cxyz(2), near_z)
       end if
    end do
  contains
    impure subroutine wolff_local_update(stack, spins, cluster, s, reflection_spin, r_sigma_c, x, y, z)
      integer(int64), intent(inout) :: stack(1:3, nall)
      logical, intent(inout) :: cluster(nx, ny, nz)
      real(real64), intent(inout) :: spins(1:2, nx, ny, nz)
      integer(int64), intent(inout) :: s
      real(real64), intent(in) :: reflection_spin(1:2), r_sigma_c
      integer(int64), intent(in) :: x, y, z
      real(real64) :: r_sigma_near, beta_sigma_sigma, prob
      r_sigma_near = sum(reflection_spin(1:2) * spins(1:2, x, y, z))
      beta_sigma_sigma = 2 * beta * r_sigma_c * r_sigma_near
      if (beta_sigma_sigma >= 0) return
      prob = 1 - exp(beta_sigma_sigma)
      if (grnd() >= prob) return
      s = s + 1
      stack(1:3, s) = [x, y, z]
      cluster(x, y, z) = .true.
      spins(1:2, x, y, z) = &
           & spins(1:2, x, y, z) - 2 * sum(spins(1:2, x, y, z) * reflection_spin(1:2)) * reflection_spin(1:2)
    end subroutine wolff_local_update
  end subroutine update_wolff

  !> calc_energy: Calculate the magnetization for the lattice.
  pure real(real64) function calc_energy() result(res)
    integer(int64) :: x, y, z
    integer(int64) :: near_x, near_y, near_z
    real(real64) :: summs(1:2)
    res = 0.0d0
    do z = 1, nz
       near_z = z + 1
       if (near_z > nz) near_z = 1
       do y = 1, ny
          near_y = y + 1
          if (near_y > ny) near_y = 1
          do x = 1, nx
             near_x = x + 1
             if (near_x > nx) near_x = 1
             summs(1:2) = spins(1:2, near_x, y, z) + spins(1:2, x, near_y, z) + spins(1:2, x, y, near_z)
             res = res - sum(spins(1:2, x, y, z) * summs(1:2))
          end do
       end do
    end do
    res = res * nall_inv
  end function calc_energy
  !> calc_magne: Calculate the magnetization for the lattice.
  pure real(real64) function calc_magne() result(res)
    integer(int64) :: x, y, z
    res = 0.0d0
    do z = 1, nz
       do y = 1, ny
          do x = 1, nx
             res = res + spins(1, x, y, z)
          end do
       end do
    end do
    res = res * nall_inv
  end function calc_magne
  !> calc_magne_abs: Calculate the absolution magnetization for the lattice.
  pure real(real64) function calc_magne_abs() result(res)
    integer(int64) :: x, y, z
    real(real64) :: summs(1:2)
    summs(1:2) = 0d0
    do z = 1, nz
       do y = 1, ny
          do x = 1, nx
             summs(1:2) = summs(1:2) + spins(1:2, x, y, z)
          end do
       end do
    end do
    summs(:) = summs(:) * nall_inv
    res = hypot(summs(1), summs(2))
  end function calc_magne_abs
end module xy3d_periodic_m

program xy3d_periodic_wolff_equilibrium_simulation
  use, intrinsic :: iso_fortran_env
  use gf2xe
  use msmt19937
  use xy3d_periodic_m
  implicit none
  integer(int32), parameter :: mcs_discard = 100000_int32, mcs_sample = 100000_int32
  integer(int64), parameter :: iseed = 42_int64
  real(real64) :: magne, magne_square, magne_sus, energy, energy_square, specific_heat
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
       write(out(i), '(a)') "# Nsize, Nsample, kbt, <m>, <e>, <m^2>, <e^2>, χ, C"
    end do
  end associate
  call init_lattice_order()
  write(error_unit, '(a, i0)') "Discard: ", mcs_discard
  do j = 1, mcs_discard
     call update_wolff()
  end do
  write(error_unit, '(a, i0)') "Sampling: ", mcs_sample
  magne = 0d0
  magne_square = 0d0
  energy = 0d0
  energy_square = 0d0
  do j = 1, mcs_sample
     call update_wolff()
     associate(m => calc_magne_abs(), e => calc_energy())
       magne = magne + m
       magne_square = magne_square + m ** 2
       energy = energy + e
       energy_square = energy_square + e ** 2
     end associate
  end do
  magne = magne / mcs_sample
  magne_square = magne_square / mcs_sample
  energy = energy / mcs_sample
  energy_square = energy_square / mcs_sample
  magne_sus = magne_square - magne ** 2
  specific_heat = energy_square - energy ** 2
  write(output_unit, '(*(g0, 1x))') nall, mcs_sample, kbt, &
       & magne, energy, &
       & magne_square, energy_square, &
       & magne_sus, specific_heat

  !> output the state of lattice.
  block
    character(len=*), parameter :: filename = "spin_configuration.txt"
    integer(int32) :: myunit
    real(real64) :: theta(nx, ny, nz)
    integer(int64) :: i, j, k
    open(file = filename, newunit = myunit)
    do i = 1, nz
       do j = 1, ny
          do k = 1, nx
             theta(k, j, i) = atan2(spins(2, k, j, i), spins(1, k, j, i))
          end do
       end do
    end do
    write(myunit, '(*(g0, 1x))') theta(:, :, :)
    close(myunit)
  end block
end program xy3d_periodic_wolff_equilibrium_simulation
