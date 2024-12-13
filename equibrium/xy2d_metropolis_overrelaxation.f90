module xy2d_periodic_correlation_simple_m
  use, intrinsic :: iso_fortran_env
  use msmt19937
  implicit none
  private
  real(real64), parameter :: pi = 4 * atan(1d0)
  integer(int64), parameter, public :: nx = 2000, ny = nx, nall = nx * ny
  real(real64), parameter :: nall_inv = 1d0 / nall

  real(real64), parameter, public :: kbt = 2.0d0, beta = 1 / kbt
  real(real64), allocatable, public, protected :: spins(:, :, :)
  real(real64), allocatable, private :: rnds(:, :, :)

  integer(int64), allocatable :: stack(:, :) !> for Wolff. stack(x or y or z, maximum size)
  logical, allocatable :: cluster(:, :) !> for Wolff. cluster(position of x, position of y, position of z)

  integer(int32), parameter :: nd = 4
  integer(int64), parameter :: dy(nd) = [1_int64, 0_int64, -1_int64, 0_int64]
  integer(int64), parameter :: dx(nd) = [0_int64, 1_int64, 0_int64, -1_int64]

  integer(int64), public, parameter :: n_unique_correlation = 2
  real(real64), public, protected :: unique_dist(n_unique_correlation)

  public :: init_simulation_xy2d, init_lattice_order, init_lattice_disorder
  public :: calc_energy, calc_magne, calc_magne_abs
  public :: update_metropolis, update_overrelaxation
  public :: update_wolff
  public :: calc_correlation
contains
  !> init_simulation_xy2d: Initialize the lattice and the random number genrator.
  impure subroutine init_simulation_xy2d()
    integer(int64) :: i
    allocate(spins(2, nx, ny))
    allocate(rnds(2, nx, ny))
    unique_dist(1) = hypot(real(nx / 2 - 1, real64), real(ny / 2, real64))
    unique_dist(n_unique_correlation) = hypot(real(nx / 2, real64), real(ny / 2, real64))
    allocate(stack(1:2, nall))
    allocate(cluster(nx, ny))
  end subroutine init_simulation_xy2d
  !> init_lattice_order: Initialize spins on the lattice with the all-aligned state.
  impure subroutine init_lattice_order()
    integer(int64) :: y, x
    do y = 1, ny
       do x = 1, nx
          spins(:, x, y) = [1.0d0, 0.0d0]
       end do
    end do
  end subroutine init_lattice_order
  !> init_lattice_disorder: Initialize spins on the lattice with the disorder state.
  impure subroutine init_lattice_disorder()
    integer(int64) :: y, x
    real(real64) :: theta
    do y = 1, ny
       do x = 1, nx
          theta = 2 * pi * grnd()
          spins(:, x, y) = [cos(theta), sin(theta)]
       end do
    end do
  end subroutine init_lattice_disorder

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
    !> even.
    do y = 1, ny
       do x = 1 + iand(y - 1, b'1'), nx, 2
          call local_flip(x, y)
       end do
    end do
    !> odd.
    do y = 1, ny
       do x = 1 + iand(y, b'1'), nx, 2
          call local_flip(x, y)
       end do
    end do
  end subroutine update_metropolis
  !> local_flip: Update the lattice with 1 MCS.
  !> @param idx An index of the lattice.
  impure subroutine local_flip(x, y)
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
       near_x = x + dx(d)
       if (near_x > nx) near_x = 1_int64
       if (near_x < 1) near_x = nx
       summ(1:2) = summ(1:2) + spins(1:2, near_x, near_y)
    end do
    new_theta = 2 * pi * rnds(1, x, y)
    new_spin(1:2) = [cos(new_theta), sin(new_theta)]
    diff_center(1:2) = new_spin(1:2) - spins(1:2, x, y)
    delta_e = - sum(diff_center(1:2) * summ(1:2))
    if (delta_e <= 0) then
       spins(1:2, x, y) = new_spin(1:2)
       return
    end if
    !> delta_e > 0.
    if (rnds(2, x, y) < exp(- beta * delta_e)) &
         & spins(1:2, x, y) = new_spin(1:2)
  end subroutine local_flip

  !> update_wolff: Update the spin cluster by Wolff algorithm
  impure subroutine update_wolff()
    integer(int64) :: x, y
    integer(int64) :: near_x, near_y, near_z
    integer(int64) :: cxyz(1:2), s
    real(real64) :: reflection_theta, reflection_spin(1:2)
    real(real64) :: r_sigma_c
    x = int(nx * grnd() + 1, int64)
    y = int(ny * grnd() + 1, int64)
    cluster(:, :) = .false.
    s = 1
    stack(1:2, s) = [x, y]
    cluster(x, y) = .true.
    reflection_theta = 2 * pi * grnd()
    reflection_spin(1:2) = [cos(reflection_theta), sin(reflection_theta)]
    spins(1:2, x, y) = spins(1:2, x, y) - 2 * sum(spins(1:2, x, y) * reflection_spin(1:2)) * reflection_spin(1:2)
    do while (s > 0_int64)
       cxyz(1:2) = stack(1:2, s)
       s = s - 1
       r_sigma_c = sum(reflection_spin(1:2) * spins(1:2, cxyz(1), cxyz(2)))
       !> right
       near_x = cxyz(1) + 1
       if (near_x > nx) near_x = 1
       if (.not. cluster(near_x, cxyz(2))) then
          call wolff_local_update(stack, spins, cluster, s, reflection_spin, r_sigma_c, near_x, cxyz(2))
       end if
       !> left
       near_x = cxyz(1) - 1
       if (near_x < 1) near_x = nx
       if (.not. cluster(near_x, cxyz(2))) then
          call wolff_local_update(stack, spins, cluster, s, reflection_spin, r_sigma_c, near_x, cxyz(2))
       end if
       !> up
       near_y = cxyz(2) + 1
       if (near_y > ny) near_y = 1
       if (.not. cluster(cxyz(1), near_y)) then
          call wolff_local_update(stack, spins, cluster, s, reflection_spin, r_sigma_c, cxyz(1), near_y)
       end if
       !> down
       near_y = cxyz(2) - 1
       if (near_y < 1) near_y = ny
       if (.not. cluster(cxyz(1), near_y)) then
          call wolff_local_update(stack, spins, cluster, s, reflection_spin, r_sigma_c, cxyz(1), near_y)
       end if
    end do
  contains
    impure subroutine wolff_local_update(stack, spins, cluster, s, reflection_spin, r_sigma_c, x, y)
      integer(int64), intent(inout) :: stack(1:2, nall)
      logical, intent(inout) :: cluster(nx, ny)
      real(real64), intent(inout) :: spins(1:2, nx, ny)
      integer(int64), intent(inout) :: s
      real(real64), intent(in) :: reflection_spin(1:2), r_sigma_c
      integer(int64), intent(in) :: x, y
      real(real64) :: r_sigma_near, beta_sigma_sigma, prob
      r_sigma_near = sum(reflection_spin(1:2) * spins(1:2, x, y))
      beta_sigma_sigma = 2 * beta * r_sigma_c * r_sigma_near
      if (beta_sigma_sigma >= 0) return
      prob = 1 - exp(beta_sigma_sigma)
      if (grnd() >= prob) return
      s = s + 1
      stack(1:2, s) = [x, y]
      cluster(x, y) = .true.
      spins(1:2, x, y) = &
           & spins(1:2, x, y) - 2 * sum(spins(1:2, x, y) * reflection_spin(1:2)) * reflection_spin(1:2)
    end subroutine wolff_local_update
  end subroutine update_wolff

  !> update_overrelaxation: Update the lattice with overrelaxation algorithm.
  impure subroutine update_overrelaxation()
    integer(int64) :: y, x
    integer(int64) :: i
    !> even.
    do y = 1, ny
       do x = 1 + iand(y - 1, b'1'), nx, 2
          call local_overrelaxation(spins, x, y)
       end do
    end do
    !> odd.
    do y = 1, ny
       do x = 1 + iand(y, b'1'), nx, 2
          call local_overrelaxation(spins, x, y)
       end do
    end do
  contains
    pure subroutine local_overrelaxation(spins, x, y)
      real(real64), intent(inout) :: spins(1:2, 1:nx, 1:ny)
      integer(int64), intent(in) :: x, y
      real(real64) :: local_field(1:2)
      integer(int64) :: lx, rx, dy, uy
      lx = x - 1
      if (lx < 1) lx = nx
      rx = x + 1
      if (rx > nx) rx = 1
      dy = y - 1
      if (dy < 1) dy = ny
      uy = y + 1
      if (uy > ny) uy = 1

      local_field(:) = 0d0
      local_field(1:2) = local_field(1:2) + spins(1:2, lx, y)
      local_field(1:2) = local_field(1:2) + spins(1:2, rx, y)
      local_field(1:2) = local_field(1:2) + spins(1:2, x, dy)
      local_field(1:2) = local_field(1:2) + spins(1:2, x, uy)

      !> local_field 方向へ reflection

      !> 内積を取る, hypotで規格化利用.
      block
        real(real64) :: abs_local_field_inv
        abs_local_field_inv = 1 / hypot(local_field(1), local_field(2))
        local_field(1:2) = local_field(1:2) * abs_local_field_inv
        spins(1:2, x, y) = &
             & (2 * sum(local_field(1:2) * spins(1:2, x, y))) * local_field(1:2) - spins(1:2, x, y)
      end block

      !> 内積を取る, thetaで規格化利用.
      ! block
      !   real(real64) :: local_field_theta
      !   local_field_theta = atan2(local_field(2), local_field(1))
      !   local_field(1:2) = [cos(local_field_theta), sin(local_field_theta)]
      !   spins(1:2, x, y) = &
      !        & (2 * sum(local_field(1:2) * spins(1:2, x, y))) * local_field(1:2) - spins(1:2, x, y)
      ! end block

      !> 内積を取ると, local_field(:) の大きさが小さい場合に不具合が起こる可能性.
      !> 角度に変換して reflection.
      ! block
      !   real(real64) :: theta_local_field, theta_spins, new_theta_spins
      !   theta_local_field = atan2(local_field(2), local_field(1))
      !   theta_spins = atan2(spins(2, x, y), spins(1, x, y))
      !   new_theta_spins = - (theta_spins - theta_local_field) + theta_local_field
      !   spins(1:2, x, y) = [cos(new_theta_spins), sin(new_theta_spins)]
      ! end block
    end subroutine local_overrelaxation
  end subroutine update_overrelaxation

  !> calc_energy: Calculate the energy for the lattice.
  pure real(real64) function calc_energy() result(res)
    integer(int64) :: x, y
    integer(int64) :: uy, rx
    res = 0.0d0
    do y = 1, ny
       do x = 1, nx
          uy = y + 1
          if (uy > ny) uy = 1_int64
          rx = x + 1
          if (rx > nx) rx = 1_int64
          res = res - sum(spins(:, x, y) * (spins(:, x, uy) + spins(:, rx, y)))
       end do
    end do
    res = res * nall_inv
  end function calc_energy
  !> calc_magne: Calculate the magnetization for the lattice.
  pure real(real64) function calc_magne() result(res)
    integer(int64) :: x, y
    res = 0.0d0
    do y = 1, ny
       do x = 1, nx
          res = res + spins(1, x, y)
       end do
    end do
    res = res * nall_inv
  end function calc_magne
  !> calc_magne_abs: Calculate the absolution of the magnetization for the lattice.
  pure real(real64) function calc_magne_abs() result(res)
    integer(int64) :: x, y
    real(real64) :: summ_x, summ_y
    summ_x = 0.0d0
    summ_y = 0.0d0
    do y = 1, ny
       do x = 1, nx
          summ_x = summ_x + spins(1, x, y)
          summ_y = summ_y + spins(2, x, y)
       end do
    end do
    summ_x = summ_x * nall_inv
    summ_y = summ_y * nall_inv
    res = hypot(summ_x, summ_y)
  end function calc_magne_abs
  !> calc_correlation_dxdy: Calculate the correlation function for `Σ_yΣ_x S(x, y) * S(x + dx, y + dy)`.
  pure real(real64) function calc_correlation_dxdy(dx, dy) result(res)
    integer(int64), intent(in) :: dx, dy
    integer(int64) :: rx, ry
    integer(int64) :: x, y
    res = 0.0d0
    do y = 1, ny
       do x = 1, nx
          rx = x + dx
          if (rx > nx) rx = rx - nx
          if (rx < 1) rx = rx + nx
          ry = y + dy
          if (ry > ny) ry = ry - ny
          if (ry < 1) ry = ry + ny
          res = res + sum(spins(1:2, x, y) * spins(1:2, rx, ry))
       end do
    end do
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
end module xy2d_periodic_correlation_simple_m

program xy2d_periodic_metropolis_equilibrium_simulation
  use, intrinsic :: iso_fortran_env
  use gf2xe
  use msmt19937
  use xy2d_periodic_correlation_simple_m
  implicit none
  integer(int32), parameter :: mcs = 100000_int32
  integer(int32), parameter :: overrelaxation_per_mcs = 1
  integer(int64), parameter :: iseed = 42_int64
  integer(int32) :: i, j

  call init_simulation_xy2d()
  call init_genrand(iseed)

  associate(out => [output_unit, error_unit])
    do i = 1, size(out)
       write(out(i), '(a,i0)'    ) "# Nsize: ", nall
       write(out(i), '(2(a, i0))') "# nx: ", nx, " ny: ", ny
       write(out(i), '(a, i0)') "# MCS: ", mcs
       write(out(i), '(a, i0)') "# overrelaxation / 1MCS: ", overrelaxation_per_mcs
       write(out(i), '(a, *(g0, 1x))' ) "# 温度: ", kbt
       write(out(i), '(a)' ) "# method: Metropolis, model: 2D-XY"
       write(out(i), '(a)') "# Nsize, kbt, time, <m>, <e>"
    end do
  end associate
  ! call init_lattice_order()
  call init_lattice_disorder()

  do i = 1, mcs
     call update_metropolis()
     do j = 1, overrelaxation_per_mcs
        call update_overrelaxation()
     end do
     write(error_unit, '(*(g0, 1x))') "MCS: ", i, " T:", kbt
     associate(m => calc_magne_abs(), e => calc_energy())
       write(output_unit, '(*(g0, 1x))') nall, kbt, i, m, e
     end associate
  end do

  ! !> output the state of lattice.
  ! block
  !   character(len=*), parameter :: filename = "spin_configuration.txt"
  !   integer(int32) :: myunit
  !   real(real64) :: theta(nx, ny)
  !   integer(int64) :: i, j, k
  !   open(file = filename, newunit = myunit)
  !   do j = 1, ny
  !      do k = 1, nx
  !         theta(k, j) = atan2(spins(2, k, j), spins(1, k, j))
  !      end do
  !   end do
  !   write(myunit, '(*(g0, 1x))') theta(:, :)
  !   close(myunit)
  ! end block
end program xy2d_periodic_metropolis_equilibrium_simulation
