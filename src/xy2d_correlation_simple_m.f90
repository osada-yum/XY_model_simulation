module xy2d_correlation_simple_m
  use, intrinsic :: iso_fortran_env
  use msmt19937
  implicit none
  private
  real(real64), parameter :: pi = 4 * atan(1d0)
  integer(int64), parameter, public :: nx = 1001, ny = nx - 1, nall = nx * ny
  real(real64), parameter :: nall_inv = 1.0d0 / nall
  real(real64), parameter, public :: kbt = 0.88d0, beta = 1/ kbt
  real(real64), allocatable, private :: spins(:, :)
  real(real64), allocatable, private :: rnds(:, :)

  integer(int64), public, parameter :: n_unique_correlation = ny / 2

  public :: init_simulation_xy2d, init_lattice_order, update_metropolis, calc_energy, calc_magne
  public :: calc_correlation
contains
  !> init_simulation_xy2d: Initialize the lattice and the random number genrator.
  !> Skip `n_skip` * 2^(log_2(`nall` * `mcs`)).
  impure subroutine init_simulation_xy2d()
    integer(int64) :: x, y
    integer(int64) :: i
    allocate(spins(2, -nx : nall + nx))
    allocate(rnds(2, 1:nall))
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
    do i = 1, nall
       rnds(1, i) = grnd()
       rnds(2, i) = grnd()
    end do
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
    new_theta = 2 * pi * rnds(1, idx)
    new_spin(1:2) = [cos(new_theta), sin(new_theta)]
    diff_center(1:2) = new_spin(1:2) - spins(1:2, idx)
    delta_e = - sum(diff_center(1:2) * summ(1:2))
    if (delta_e <= 0) then
       spins(1:2, idx) = new_spin(1:2)
       return
    end if
    !> delta_e > 0.
    if (rnds(2, idx) < exp(- beta * delta_e)) &
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
  !> calc_correlation: Calculate the correlation function for the lattice.
  pure subroutine calc_correlation(corr)
    real(real64), intent(inout) :: corr(n_unique_correlation)
    integer(int64) :: up_idx
    integer(int64) :: i
    integer(int64) :: r
    corr(:) = 0.0d0
    do i = 1, nall
       up_idx = i
       do r = 1, ny / 2
          up_idx = up_idx + nx
          if (up_idx > nall) up_idx = up_idx - nall
          corr(r) = corr(r) + sum(spins(1:2, i) * spins(1:2, up_idx))
       end do
    end do
    corr(:) = corr(:) / nall
    corr(ny / 2) = corr(ny / 2) / 2
  end subroutine calc_correlation
end module xy2d_correlation_simple_m
