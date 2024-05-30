module xy2d_simple_m
  use, intrinsic :: iso_fortran_env
  use msmt19937
  implicit none
  private
  real(real64), parameter :: pi = 4 * atan(1d0)
  integer(int64), parameter, public :: nx = 101, ny = nx - 1, nall = nx * ny
  real(real64), parameter :: nall_inv = 1.0d0 / nall
  real(real64), parameter, public :: kbt = 0.90d0, beta = 1/ kbt
  real(real64), allocatable, private :: spins(:, :)
  public :: init_simulation_xy2d, init_lattice_order, update_metropolis, calc_energy, calc_magne
contains
  !> init_simulation_xy2d: Initialize the lattice and the random number genrator.
  !> Skip `n_skip` * 2^(log_2(`nall` * `mcs`)).
  impure subroutine init_simulation_xy2d()
    allocate(spins(2, -nx : nall + nx))
  end subroutine init_simulation_xy2d
  !> init_lattice_order: Initialize spins on the lattice with the all-aligned state.
  impure subroutine init_lattice_order()
    integer(int64) :: i
    do i = - nx, nall + nx
       spins(:, i) = [1.0d0, 0.0d0]
    end do
  end subroutine init_lattice_order
  !> update_norishiro: Update the norishiro in the lattice.
  impure subroutine update_norishiro()
    integer(int64) :: i
    do i = 1, nx
       spins(:, i + nall) = spins(:, i)
    end do
    do i = nall - nx + 1, nall
       spins(:, i - nall) = spins(:, i)
    end do
  end subroutine update_norishiro
  !> update_metropolis: Update the lattice with 1 MCS.
  impure subroutine update_metropolis()
    integer(int64) :: i
    !> 奇数.
    do i = 1, nx, 2
       call local_flip(i)
       spins(:, i + nall) = spins(:, i)
    end do
    do i = nx + 2, nall - nx, 2
       call local_flip(i)
    end do
    do i = nall - nx + 2, nall, 2
       call local_flip(i)
       spins(:, i - nall) = spins(:, i)
    end do
    !> 偶数.
    do i = 2, nx, 2
       call local_flip(i)
       spins(:, i + nall) = spins(:, i)
    end do
    do i = nx + 1, nall - nx, 2
       call local_flip(i)
    end do
    do i = nall - nx + 1, nall, 2
       call local_flip(i)
       spins(:, i - nall) = spins(:, i)
    end do
  end subroutine update_metropolis
  !> local_flip: Update the lattice with 1 MCS.
  !> @param idx An index of the lattice.
  impure subroutine local_flip(idx)
    integer(int64), intent(in) :: idx
    integer(int32), parameter :: nd = 4
    integer(int64), parameter :: dx(nd) = [1_int64, nx, -1_int64, -nx]
    integer(int64) :: nidx
    real(real64) :: summ(1:2)
    real(real64) :: diff_center(1:2)
    real(real64) :: new_theta, new_spin(1:2)
    real(real64) :: delta_e
    integer(int32) :: d
    summ(1:2) = 0.0d0
    do d = 1, 4
       nidx = idx + dx(d)
       summ(1:2) = summ(1:2) + spins(1:2, nidx)
    end do
    new_theta = 2 * pi * grnd()
    new_spin(1:2) = [cos(new_theta), sin(new_theta)]
    diff_center(1:2) = new_spin(1:2) - spins(1:2, idx)
    delta_e = - sum(diff_center(1:2) * summ(1:2))
    if (delta_e <= 0) then
       spins(1:2, idx) = new_spin(1:2)
       return
    end if
    !> delta_e > 0.
    if (grnd() < exp(- beta * delta_e)) &
         & spins(1:2, idx) = new_spin(1:2)
  end subroutine local_flip
  !> calc_energy: Calculate the energy for the lattice.
  pure real(real64) function calc_energy() result(res)
    integer(int64) :: i
    res = 0.0d0
    do i = 1, nall
       res = res - sum(spins(:, i) * (spins(:, i + 1) + spins(:, i + nx)))
    end do
    res = res * nall_inv
  end function calc_energy
  !> calc_magne: Calculate the magnetization for the lattice.
  pure real(real64) function calc_magne() result(res)
    integer(int64) :: i
    res = 0.0d0
    do i = 1, nall
       res = res + spins(1, i)
    end do
    res = res * nall_inv
  end function calc_magne
end module xy2d_simple_m
