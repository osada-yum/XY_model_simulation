module xy2d_dual_lattice_m
  use, intrinsic :: iso_fortran_env
  use msmt19937
  implicit none
  private
  real(real64), parameter :: pi = 4 * atan(1d0)
  integer(int64), parameter, public :: nx = 101, ny = nx - 1, nall = nx * ny
  real(real64), parameter :: nall_inv = 1.0d0 / nall
  real(real64), parameter, public :: kbt = 0.90d0, beta = 1/ kbt
  real(real64), allocatable, private :: spins_odd(:, :)
  real(real64), allocatable, private :: spins_even(:, :)
  public :: init_simulation_dual_lattice_xy2d, init_lattice_order, update_metropolis, calc_energy, calc_magne
contains
  !> init_simulation_dual_lattice_xy2d: Initialize the lattice and the random number genrator.
  !> Skip `n_skip` * 2^(log_2(`nall` * `mcs`)).
  impure subroutine init_simulation_dual_lattice_xy2d()
    allocate(spins_odd(2, 1 - nx / 2 : nall / 2 + (nx + 1) / 2))
    allocate(spins_even(2, 1 - (nx + 1) / 2 : nall / 2 + nx / 2))
  end subroutine init_simulation_dual_lattice_xy2d
  !> init_lattice_order: Initialize spins on the lattice with the all-aligned state.
  impure subroutine init_lattice_order()
    integer(int64) :: i
    do i = 1 - nx / 2, nall / 2 + (nx + 1) / 2
       spins_odd(:, i) = [1.0d0, 0.0d0]
    end do
    do i = 1 - (nx + 1) / 2, nall / 2 + nx / 2
       spins_even(:, i) = [1.0d0, 0.0d0]
    end do
  end subroutine init_lattice_order
  !> update_norishiro: Update the norishiro in the lattice.
  impure subroutine update_norishiro()
    integer(int64) :: i
    do i = 1, (nx + 1) / 2
       spins_odd(:, i + nall / 2) = spins_odd(:, i)
    end do
    do i = 1, nx / 2
       spins_even(:, i + nall / 2) = spins_even(:, i)
    end do
    do i = nall / 2 - nx / 2 + 1, nall / 2
       spins_odd(:, i - nall / 2) = spins_odd(:, i)
    end do
    do i = nall / 2 - (nx + 1) / 2 + 1, nall / 2
       spins_even(:, i - nall / 2) = spins_even(:, i)
    end do
  end subroutine update_norishiro
  !> update_metropolis: Update the lattice with 1 MCS.
  impure subroutine update_metropolis()
    integer(int64), parameter :: norishiro_odd_down_upper = (nx + 1) / 2, norishiro_odd_up_lower = nall / 2 - nx / 2 + 1
    integer(int64), parameter :: norishiro_even_down_upper = nx / 2, norishiro_even_up_lower = nall / 2 - (nx + 1) / 2 + 1
    integer(int64) :: i
    !> odd.
    do i = 1, norishiro_odd_down_upper
       call local_flip_odd(i)
       spins_odd(:, i + nall / 2) = spins_odd(:, i)
    end do
    do i = norishiro_odd_down_upper + 1, norishiro_odd_up_lower - 1
       call local_flip_odd(i)
    end do
    do i = norishiro_odd_up_lower, nall / 2
       call local_flip_odd(i)
       spins_odd(:, i - nall / 2) = spins_odd(:, i)
    end do
    !> even.
    do i = 1, norishiro_even_down_upper
       call local_flip_even(i)
       spins_even(:, i + nall / 2) = spins_even(:, i)
    end do
    do i = norishiro_even_down_upper + 1, norishiro_even_up_lower - 1
       call local_flip_even(i)
    end do
    do i = norishiro_even_up_lower, nall / 2
       call local_flip_even(i)
       spins_even(:, i - nall / 2) = spins_even(:, i)
    end do
  end subroutine update_metropolis
  !> local_flip_odd: Update the lattice with 1 MCS.
  !> @param idx An index of the lattice.
  impure subroutine local_flip_odd(idx)
    integer(int64), intent(in) :: idx
    integer(int32), parameter :: nd = 4
    integer(int64), parameter :: dx(nd) = [-1_int64, 0_int64, nx / 2, - (nx + 1) / 2]
    integer(int64) :: nidx
    real(real64) :: summ(1:2)
    real(real64) :: diff_center(1:2)
    real(real64) :: new_theta, new_spin(1:2)
    real(real64) :: delta_e
    integer(int32) :: d
    summ(1:2) = 0.0d0
    do d = 1, 4
       nidx = idx + dx(d)
       summ(1:2) = summ(1:2) + spins_even(1:2, nidx)
    end do
    new_theta = 2 * pi * grnd()
    new_spin(1:2) = [cos(new_theta), sin(new_theta)]
    diff_center(1:2) = new_spin(1:2) - spins_odd(1:2, idx)
    delta_e = - sum(diff_center(1:2) * summ(1:2))
    if (delta_e <= 0) then
       spins_odd(1:2, idx) = new_spin(1:2)
       return
    end if
    !> delta_e > 0.
    if (grnd() < exp(- beta * delta_e)) &
         & spins_odd(1:2, idx) = new_spin(1:2)
  end subroutine local_flip_odd
  !> local_flip_even: Update the lattice with 1 MCS.
  !> @param idx An index of the lattice.
  impure subroutine local_flip_even(idx)
    integer(int64), intent(in) :: idx
    integer(int32), parameter :: nd = 4
    integer(int64), parameter :: dx(nd) = [0_int64, 1_int64, (nx + 1) / 2, -nx / 2]
    integer(int64) :: nidx
    real(real64) :: summ(1:2)
    real(real64) :: diff_center(1:2)
    real(real64) :: new_theta, new_spin(1:2)
    real(real64) :: delta_e
    integer(int32) :: d
    summ(1:2) = 0.0d0
    do d = 1, 4
       nidx = idx + dx(d)
       summ(1:2) = summ(1:2) + spins_odd(1:2, nidx)
    end do
    new_theta = 2 * pi * grnd()
    new_spin(1:2) = [cos(new_theta), sin(new_theta)]
    diff_center(1:2) = new_spin(1:2) - spins_even(1:2, idx)
    delta_e = - sum(diff_center(1:2) * summ(1:2))
    if (delta_e <= 0) then
       spins_even(1:2, idx) = new_spin(1:2)
       return
    end if
    !> delta_e > 0.
    if (grnd() < exp(- beta * delta_e)) &
         & spins_even(1:2, idx) = new_spin(1:2)
  end subroutine local_flip_even
  !> calc_energy: Calculate the energy for the lattice.
  pure real(real64) function calc_energy() result(res)
    integer(int64) :: i
    res = 0.0d0
    do i = 1, nall / 2
       res = res &
            & - sum(spins_odd(:, i) * (spins_even(:, i) + spins_even(:, i + nx / 2))) &
            & - sum(spins_even(:, i) * (spins_odd(:, i + 1) + spins_odd(:, i + (nx + 1) / 2)))
    end do
    res = res * nall_inv
  end function calc_energy
  !> calc_magne: Calculate the magnetization for the lattice.
  pure real(real64) function calc_magne() result(res)
    integer(int64) :: i
    res = 0.0d0
    do i = 1, nall / 2
       res = res + spins_odd(1, i) + spins_even(1, i)
    end do
    res = res * nall_inv
  end function calc_magne
end module xy2d_dual_lattice_m
