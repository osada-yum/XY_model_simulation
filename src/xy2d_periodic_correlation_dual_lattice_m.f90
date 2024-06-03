module xy2d_periodic_correlation_dual_lattice_m
  use, intrinsic :: iso_fortran_env
  use msmt19937
  implicit none
  private
  real(real64), parameter :: pi = 4 * atan(1d0)
  integer(int64), parameter, public :: nx = 500, ny = nx, nall = nx * ny
  real(real64), parameter :: nall_inv = 1.0d0 / nall
  real(real64), parameter, public :: kbt = 0.88d0, beta = 1/ kbt
  real(real64), allocatable, private :: spins_even(:, :, :) !> Black in checker board, forall (x, y) x + y is even.
  real(real64), allocatable, private :: spins_odd(:, :, :) !> White in checker board.
  real(real64), allocatable, private :: rnds(:, :, :)

  integer(int32), parameter :: nd = 4, up = 1, right = 2 !> 1: up, 2: right, 3: down, 4: left.
  integer(int64), parameter :: dy(nd) = [1_int64, 0_int64, -1_int64, 0_int64]
  !>                (even, y even), (odd, y odd)           (odd, y even), (even, y odd)
  integer(int64), parameter :: dx(nd, 0:1) = &
       & reshape([[0_int64, 1_int64, 0_int64, 0_int64], [0_int64, 0_int64, 0_int64, -1_int64]], shape = [4, 2])

  integer(int64), public, parameter :: n_unique_correlation = 2
  real(real64), public, protected :: unique_dist(n_unique_correlation)

  public :: init_simulation_xy2d, init_lattice_order, update_metropolis, calc_energy, calc_magne
  public :: calc_correlation
contains
  !> init_simulation_xy2d: Initialize the lattice and the random number genrator.
  impure subroutine init_simulation_xy2d()
    integer(int64) :: i
    allocate(spins_even(2, nx / 2, ny))
    allocate(spins_odd(2, nx / 2, ny))
    allocate(rnds(2, nx, ny))
    unique_dist(1) = hypot(real(nx / 2 - 1, real64), real(ny / 2, real64))
    unique_dist(n_unique_correlation) = hypot(real(nx / 2, real64), real(ny / 2, real64))
  end subroutine init_simulation_xy2d
  !> init_lattice_order: Initialize spins on the lattice with the all-aligned state.
  impure subroutine init_lattice_order()
    integer(int64) :: y, x
    do y = 1, ny
       do x = 1, nx / 2
          spins_even(:, x, y) = [1.0d0, 0.0d0]
          spins_odd(:, x, y) = [1.0d0, 0.0d0]
       end do
    end do
  end subroutine init_lattice_order
  !> update_metropolis: Update the lattice with 1 MCS.
  impure subroutine update_metropolis()
    integer(int64) :: y, x
    integer(int64) :: i
    do y = 1, ny
       do x = 1, nx
          rnds(1, x, y) = grnd()
          rnds(2, x, y) = grnd()
       end do
    end do
    !> even, y = 1, 2.
    do y = 1, 2
       do x = 1, nx / 2
          call local_flip_even(x, y)
       end do
    end do
    !> even, y = 3, ny.
    !> odd, y = 2, ny - 1.
    do y = 3, ny
       do x = 1, nx / 2
          call local_flip_even(x, y)
          call local_flip_odd(x, y - 1)
       end do
    end do
    do x = 1, nx / 2
       call local_flip_odd(x, ny)
    end do
    do x = 1, nx / 2
       call local_flip_odd(x, 1_int64)
    end do
  end subroutine update_metropolis
  !> local_flip_even: Update a spin in lattice.
  !> @param x An index of x-axis of the lattice.
  !> @param y An index of y-axis of the lattice.
  impure subroutine local_flip_even(x, y)
    integer(int64), intent(in) :: x, y
    real(real64) :: summ(1:2)
    real(real64) :: diff_center(1:2)
    real(real64) :: new_theta, new_spin(1:2)
    real(real64) :: delta_e
    integer(int32) :: d
    integer(int64) :: near_y, near_x
    summ(1:2) = 0.0d0
    do d = 1, nd
       near_y = y + dy(d)
       if (near_y > ny) near_y = 1_int64
       if (near_y < 1) near_y = ny
       near_x = x + dx(d, iand(y, b'1'))
       if (near_x > nx / 2) near_x = 1_int64
       if (near_x < 1) near_x = nx / 2
       summ(1:2) = summ(1:2) + spins_odd(1:2, near_x, near_y)
    end do
    new_theta = 2 * pi * rnds(1, 2 * x - 1 + iand(y + 1, b'1'), y)
    new_spin(1:2) = [cos(new_theta), sin(new_theta)]
    diff_center(1:2) = new_spin(1:2) - spins_even(1:2, x, y)
    delta_e = - sum(diff_center(1:2) * summ(1:2))
    if (delta_e <= 0) then
       spins_even(1:2, x, y) = new_spin(1:2)
       return
    end if
    !> delta_e > 0.
    if (rnds(2, 2 * x - 1 + iand(y + 1, b'1'), y) < exp(- beta * delta_e)) &
         & spins_even(1:2, x, y) = new_spin(1:2)
  end subroutine local_flip_even
  !> local_flip_odd: Update a spin in lattice.
  !> @param x An index of x-axis of the lattice.
  !> @param y An index of y-axis of the lattice.
  impure subroutine local_flip_odd(x, y)
    integer(int64), intent(in) :: x, y
    real(real64) :: summ(1:2)
    real(real64) :: diff_center(1:2)
    real(real64) :: new_theta, new_spin(1:2)
    real(real64) :: delta_e
    integer(int32) :: d
    integer(int64) :: near_y, near_x
    summ(1:2) = 0.0d0
    do d = 1, nd
       near_y = y + dy(d)
       if (near_y > ny) near_y = 1_int64
       if (near_y < 1) near_y = ny
       near_x = x + dx(d, iand(y + 1, b'1'))
       if (near_x > nx / 2) near_x = 1_int64
       if (near_x < 1) near_x = nx / 2
       summ(1:2) = summ(1:2) + spins_even(1:2, near_x, near_y)
    end do
    new_theta = 2 * pi * rnds(1, 2 * x - 1 + iand(y, b'1'), y)
    new_spin(1:2) = [cos(new_theta), sin(new_theta)]
    diff_center(1:2) = new_spin(1:2) - spins_odd(1:2, x, y)
    delta_e = - sum(diff_center(1:2) * summ(1:2))
    if (delta_e <= 0) then
       spins_odd(1:2, x, y) = new_spin(1:2)
       return
    end if
    !> delta_e > 0.
    if (rnds(2, 2 * x - 1 + iand(y, b'1'), y) < exp(- beta * delta_e)) &
         & spins_odd(1:2, x, y) = new_spin(1:2)
  end subroutine local_flip_odd
  !> calc_energy: Calculate the energy for the lattice.
  pure real(real64) function calc_energy() result(res)
    integer(int64) :: x, y
    integer(int64) :: uy, rx
    res = 0.0d0
    do y = 1, ny
       do x = 1, nx / 2
          uy = y + dy(up)
          if (uy > ny) uy = 1_int64
          !> even.
          rx = x + dx(right, iand(y, b'1'))
          if (rx > nx / 2) rx = 1_int64
          res = res - sum(spins_even(:, x, y) * (spins_odd(:, x, uy) + spins_odd(:, rx, y)))

          !> odd.
          rx = x + dx(right, iand(y + 1, b'1'))
          if (rx > nx / 2) rx = 1_int64
          res = res - sum(spins_odd(:, x, y) * (spins_even(:, x, uy) + spins_even(:, rx, y)))
       end do
    end do
    res = res * nall_inv
  end function calc_energy
  !> calc_magne: Calculate the magnetization for the lattice.
  pure real(real64) function calc_magne() result(res)
    integer(int64) :: x, y
    res = 0.0d0
    do y = 1, ny
       do x = 1, nx / 2
          res = res + spins_even(1, x, y) + spins_odd(1, x, y)
       end do
    end do
    res = res * nall_inv
  end function calc_magne
  !> calc_correlation_dxdy: Calculate the correlation function for `Σ_yΣ_x S(x, y) * S(x + dx, y + dy)`.
  pure real(real64) function calc_correlation_dxdy(dx, dy) result(res)
    integer(int64), intent(in) :: dx, dy
    integer(int64) :: rx, ry
    integer(int64) :: x, y
    integer(int64) :: actual_x
    res = 0.0d0
    if (iand(dx + dy, b'1') == 0) then !> even -> even or odd -> odd.
       do y = 1, ny
          do x = 1, nx / 2
             ry = y + dy
             if (ry > ny) ry = ry - ny
             if (ry < 1) ry = ry + ny
             !> even.
             actual_x = 2 * x - 1 + iand(y + 1, b'1')
             rx = actual_x + dx
             if (rx > nx) rx = rx - nx
             if (rx < 1) rx = rx + nx
             res = res + sum(spins_even(1:2, x, y) * spins_even(1:2, (rx + 1) / 2, ry))
             !> odd.
             actual_x = 2 * x - 1 + iand(y, b'1')
             rx = actual_x + dx
             if (rx > nx) rx = rx - nx
             if (rx < 1) rx = rx + nx
             res = res + sum(spins_odd(1:2, x, y) * spins_odd(1:2, (rx + 1) / 2, ry))
          end do
       end do
    else !> even -> odd or odd -> even.
       do y = 1, ny
          do x = 1, nx / 2
             ry = y + dy
             if (ry > ny) ry = ry - ny
             if (ry < 1) ry = ry + ny
             !> even.
             actual_x = 2 * x - 1 + iand(y + 1, b'1')
             rx = actual_x + dx
             if (rx > nx) rx = rx - nx
             if (rx < 1) rx = rx + nx
             res = res + sum(spins_even(1:2, x, y) * spins_odd(1:2, (rx + 1) / 2, ry))
             !> odd.
             actual_x = 2 * x - 1 + iand(y, b'1')
             rx = actual_x + dx
             if (rx > nx) rx = rx - nx
             if (rx < 1) rx = rx + nx
             res = res + sum(spins_odd(1:2, x, y) * spins_even(1:2, (rx + 1) / 2, ry))
          end do
       end do
    end if
    res = res * nall_inv
  end function calc_correlation_dxdy
  !> calc_correlation: Calculate the correlation function for the lattice.
  pure subroutine calc_correlation(corr)
    real(real64), intent(inout) :: corr(n_unique_correlation)
    integer(int64) :: r_idx
    integer(int64) :: x, y
    integer(int64) :: r
    corr(1) = calc_correlation_dxdy(nx / 2 - 1, ny / 2) + calc_correlation_dxdy(nx / 2, ny / 2 - 1)
    corr(1) = corr(1) / 2
    corr(n_unique_correlation) = calc_correlation_dxdy(nx / 2, ny / 2)
  end subroutine calc_correlation
end module xy2d_periodic_correlation_dual_lattice_m
