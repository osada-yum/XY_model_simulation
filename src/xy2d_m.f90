module xy2d_m
  use, intrinsic :: iso_fortran_env
  use kahan_summation_m
  implicit none
  private
  real(real64), parameter :: pi = 4 * atan(1.0_real64)
  public :: xy2d
  type :: xy2d
     private
     integer(int64) :: nx_, ny_, nall_
     real(real64) :: beta_
     real(real64), allocatable :: spins_(:, :)
   contains
     !> initializer.
     procedure, pass :: init => init_xy2d
     !> setter.
     procedure, pass :: set_xy_allup => set_xy_allup_xy2d
     procedure, pass :: set_xy_random => set_xy_random_xy2d
     procedure, pass :: set_kbt => set_kbt_xy2d
     procedure, pass :: set_beta => set_beta_xy2d
     !> updater.
     procedure, pass :: update => update_xy2d
     procedure, pass, private :: update_onesite => update_onesite_xy2d
     procedure, pass, private :: update_norishiro => update_norishiro_xy2d
     !> calculator.
     procedure, pass :: calc_energy_summ => calc_total_energy_xy2d
     procedure, pass :: calc_magne_summ => calc_total_magne_xy2d
     !> getter.
     procedure, pass :: nx => nx_xy2d
     procedure, pass :: ny => ny_xy2d
     procedure, pass :: nall => nall_xy2d
     procedure, pass :: kbt => kbt_xy2d
     procedure, pass :: beta => beta_xy2d
     procedure, pass, private :: norishiro_begin => norishiro_begin_xy2d
     procedure, pass, private :: norishiro_end => norishiro_end_xy2d
  end type xy2d
contains
  !> init_xy2d: Initialize xy2d once.
  pure subroutine init_xy2d(this, nx, ny, kbt)
    class(xy2d), intent(inout) :: this
    integer(int64), intent(in) :: nx, ny
    real(real64), intent(in) :: kbt
    this%nx_ = nx
    this%ny_ = ny
    this%nall_ = nx * ny
    allocate(this%spins_(this%norishiro_begin() : this%norishiro_end(), 1:2))
    call this%set_xy_allup()
    call this%set_kbt(kbt)
  end subroutine init_xy2d
  !> set_xy_allup_xy2d: Set spins `(1, 0)`.
  pure subroutine set_xy_allup_xy2d(this)
    class(xy2d), intent(inout) :: this
    this%spins_(:, 1) = 1.0_real64
    this%spins_(:, 2) = 0.0_real64
  end subroutine set_xy_allup_xy2d
  !> set_xy_random_xy2d: Set spins `1` or `-1` randomly.
  impure subroutine set_xy_random_xy2d(this)
    class(xy2d), intent(inout) :: this
    real(real64), allocatable :: r(:)
    allocate(r(1:this%nall_))
    call random_number(r)
    this%spins_(1:this%nall_, 1) = cos(2 * pi * r(:))
    this%spins_(1:this%nall_, 2) = sin(2 * pi * r(:))
    call this%update_norishiro()
  end subroutine set_xy_random_xy2d
  !> set_kbt_xy2d: Set parameter `beta` as `1 / kbt`.
  pure subroutine set_kbt_xy2d(this, kbt)
    class(xy2d), intent(inout) :: this
    real(real64), intent(in) :: kbt
    call this%set_beta(1 / kbt)
  end subroutine set_kbt_xy2d
  !> set_beta_xy2d: Set parameter `beta`.
  pure subroutine set_beta_xy2d(this, beta)
    class(xy2d), intent(inout) :: this
    real(real64), intent(in) :: beta
    this%beta_ = beta
  end subroutine set_beta_xy2d

  !> update_xy2d: Update the system by Metropolis method.
  impure subroutine update_xy2d(this)
    class(xy2d), intent(inout) :: this
    real(real64), allocatable :: r(:), candidates(:)
    integer(int64) :: i, j
    allocate(r(this%nall_), candidates(this%nall_))
    call random_number(r)
    call random_number(candidates)
    do j = 1, 2
       do i = j, this%nall_, 2
          call this%update_onesite(i, r(i), candidates(i))
       end do
       call this%update_norishiro()
    end do
  end subroutine update_xy2d
  !> update_onesite_xy2d: Update a spin of the system.
  pure subroutine update_onesite_xy2d(this, idx, r, candidate)
    class(xy2d), intent(inout) :: this
    integer(int64), intent(in) :: idx
    real(real64), intent(in) :: r, candidate
    real(real64) :: delta_e, spin(1:2)
    spin(1:2) = [cos(2 * pi * candidate), sin(2 * pi * candidate)]
    delta_e = sum(- (spin(:) - this%spins_(idx, :)) * (&
         & this%spins_(idx + 1, :) + this%spins_(idx + this%nx_, :) + &
         & this%spins_(idx - 1, :) + this%spins_(idx - this%nx_, :)) )
    if (r < exp(- this%beta() * delta_e)) &
         this%spins_(idx, :) = spin(:)
  end subroutine update_onesite_xy2d
  !> update_norishiro_xy2d: Update norishiro.
  pure subroutine update_norishiro_xy2d(this)
    class(xy2d), intent(inout) :: this
    integer(int64) :: i, j
    do j = 1, 2
       do i = 1_int64, this%nx_
          this%spins_(this%norishiro_begin() + i - 1, j) = this%spins_(this%nall_ - this%nx_ + i, j)
          this%spins_(this%norishiro_end() - this%nx_ + i, j) = this%spins_(i, j)
       end do
    end do
  end subroutine update_norishiro_xy2d

  !> calc_total_energy_xy2d: Calculate the total energy.
  pure real(real64) function calc_total_energy_xy2d(this) result(res)
    class(xy2d), intent(in) :: this
    type(kahan_summation) :: ksum
    integer(int64) :: i
    ksum = kahan_summation(0.0_real64)
    do i = 1_int64, this%nall_
       ksum = ksum + sum(- this%spins_(i, :) * (this%spins_(i + 1, :) + this%spins_(i + this%nx_, :)))
    end do
    res = ksum%val()
  end function calc_total_energy_xy2d
  !> calc_total_magne_xy2d: Calculate the total magne.
  pure real(real64) function calc_total_magne_xy2d(this) result(res)
    class(xy2d), intent(in) :: this
    type(kahan_summation) :: ksum
    integer(int64) :: i
    ksum = kahan_summation(0.0_real64)
    do i = 1_int64, this%nall_
       ksum = ksum + this%spins_(i, 1)
    end do
    res = ksum%val()
  end function calc_total_magne_xy2d

  !> nx_xy2d: Return size of `x` of the system.
  pure integer(int64) function nx_xy2d(this) result(res)
    class(xy2d), intent(in) :: this
    res = this%nx_
  end function nx_xy2d
  !> ny_xy2d: Return size of `y` of the system.
  pure integer(int64) function ny_xy2d(this) result(res)
    class(xy2d), intent(in) :: this
    res = this%ny_
  end function ny_xy2d
  !> nall_xy2d: Return size of the system.
  pure integer(int64) function nall_xy2d(this) result(res)
    class(xy2d), intent(in) :: this
    res = this%nall_
  end function nall_xy2d
  !> kbt_xy2d: Return temperature of the system.
  pure real(real64) function kbt_xy2d(this) result(res)
    class(xy2d), intent(in) :: this
    res = 1 / this%beta_
  end function kbt_xy2d
  !> beta_xy2d: Return inverse temperature of the system.
  pure real(real64) function beta_xy2d(this) result(res)
    class(xy2d), intent(in) :: this
    res = this%beta_
  end function beta_xy2d
  !> norishiro_begin_xy2d: Return start index of `this%spins_(:)`.
  pure integer(int64) function norishiro_begin_xy2d(this) result(res)
    class(xy2d), intent(in) :: this
    res = 1 - this%nx_
  end function norishiro_begin_xy2d
  !> norishiro_end_xy2d: Return end index of `this%spins_(:)`.
  pure integer(int64) function norishiro_end_xy2d(this) result(res)
    class(xy2d), intent(in) :: this
    res = this%nall_ + this%nx_
  end function norishiro_end_xy2d
end module xy2d_m
