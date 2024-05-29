module xy2d_simple_m
  use, intrinsic :: iso_fortran_env
  use msmt19937
  implicit none
  integer(int64), parameter :: iseed = 42_int64
  real(real64), parameter :: pi = 4 * atan(1d0)
  integer(int64), parameter :: nx = 101, ny = nx - 1, nall = nx * ny
  real(real64), parameter :: kbt = 0.90d0, beta = 1/ kbt
  real(real64), allocatable :: spins(:, :)

contains
  impure subroutine init(n_skip)
    integer(int64), intent(in) :: n_skip
    integer(int64) :: i
    allocate(spins(2, nall))
    do i = 1_int64, nall
       spins(:, i) = [1.0d0, 0.0d0]
    end do
  end subroutine init
  impure subroutine update_metropolis()
    integer(int32) :: j
    integer(int64) :: i
    do j = 0, 1
       do i = 1 + j, nall, 2
          call local_flip(i)
          if (i <= nx) spins(:, i + nall) = spins(:, i)
          if (i > nall - nx) spins(:, i - nall) = spins(:, i)
       end do
    end do
  end subroutine update_metropolis
  impure subroutine local_flip(idx)
    integer(int64), intent(in) :: idx
    integer(int32), parameter :: nd = 4
    integer(int32), parameter :: dx(nd) = [1, nx, -1, -nx]
    integer(int64) :: nidx
    real(real64) :: summ_x, summ_y
    real(real64) :: diff_center
    integer(int32) :: d
    summ_x = 0.0d0; summ_y = 0.0d0
    do d = 1, 4
       nidx = idx + dx(d)
       summ_x = summ_x + spins(1, nidx)
       summ_y = summ_y + spins(2, nidx)
    end do
  end subroutine local_flip
end module xy2d_simple_m
