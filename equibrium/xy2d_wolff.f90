module gf2xe
  !===============================================================================
  ! Fortran 90/95 Module for GF(2)[x] computation
  !===============================================================================
  use iso_fortran_env

  implicit none
  private
  public :: gf2x_obj
  public :: gf2x_prime_obj
  public :: new,delete
  public :: print_bit,print_hex
  public :: get_deg
  public :: set_coef, set_prime
  public :: assign
  public ::  add,  add_assign
  public :: mult, mult_assign
  public ::  pow, square
  public :: div, rem, divrem
  public :: mult_by_x, div_by_x, mod_by_x
  public :: shift
  public :: deg_i32, mult_i32, square_i32, shift_i32
  public :: mult_i32_old
  public :: pow_mod
  public :: pow_pow_2

  integer(INT32), parameter :: MAX_KARA  = 64

  type gf2x_obj
     integer(INT32), pointer :: c(:) => NULL()
     integer(INT32) :: deg  = -1
     integer(INT32) :: size = -1
  end type gf2x_obj

  type gf2x_prime_obj
     type(gf2x_obj) :: prime_poly
     type(gf2x_obj) :: barrett_poly
     integer(INT32) :: deg
  end type gf2x_prime_obj

  interface new
     module procedure gf2x_new
     module procedure gf2x_delete_prime
  end interface new

  interface delete
     module procedure gf2x_delete
     module procedure gf2x_delete_prime
  end interface delete

  interface print_bit
     module procedure gf2x_print_bit
  end interface print_bit

  interface print_hex
     module procedure gf2x_print_hex
  end interface print_hex

  interface set_coef
     module procedure gf2x_set_coef
  end interface set_coef

  interface set_prime
     module procedure gf2x_set_prime
  end interface set_prime

  interface assign
     module procedure gf2x_assign
  end interface assign

  interface add
     module procedure gf2x_add
  end interface add

  interface add_assign
     module procedure gf2x_add_assign
  end interface add_assign

  interface mult
     module procedure gf2x_mult_kara
  end interface mult

  interface mult_assign
     module procedure gf2x_mult_assign_kara
  end interface mult_assign

  interface pow
     module procedure gf2x_pow
  end interface pow

  interface pow_mod
     !    module procedure gf2x_pow
     module procedure gf2x_pow_mod
  end interface pow_mod

  interface square
     module procedure gf2x_square
  end interface square

  interface mult_by_x
     module procedure gf2x_mult_by_x
  end interface mult_by_x

  interface mod_by_x
     module procedure gf2x_mod_by_x
  end interface mod_by_x

  interface div_by_x
     module procedure gf2x_div_by_x
  end interface div_by_x

  interface div
     module procedure gf2x_div
  end interface div

  interface rem
     module procedure gf2x_rem
     module procedure gf2x_rem_barrett
  end interface rem

  interface divrem
     module procedure gf2x_divrem
  end interface divrem

  interface shift
     module procedure gf2x_shift
  end interface shift

  interface pow_pow_2
     module procedure gf2x_pow_pow_2
  end interface pow_pow_2

contains

  !!DEC$ ATTRIBUTES FORCEINLINE :: get_size
  function get_size(deg) result(size)
    integer(INT32) :: deg,size
    size = CEILING(real(deg+1,kind=REAL64)/32.0_REAL64)
    return
  end function get_size

  subroutine gf2x_new(this,deg)
    type(gf2x_obj), intent(inout) :: this
    integer(INT32), intent(in) :: deg
    integer(INT32) :: isize
    intrinsic :: SIZE
    if (deg < 0) then
       this%deg  = -1
       this%size = -1
       return
    endif
    isize = get_size(deg)
    this%size = isize
    this%deg  = deg
    if (.not.associated(this%c)) then
       allocate(this%c(0:isize-1))
    else
       if (SIZE(this%c) < this%size) then
          deallocate(this%c)
          NULLIFY(this%c)
          allocate(this%c(0:isize-1))
       endif
    endif
    this%c(:) = 0
    return
  end subroutine gf2x_new

  subroutine gf2x_delete(this)
    type(gf2x_obj), intent(inout) :: this
    integer(INT32) :: ierr
    if (associated(this%c)) then
       deallocate(this%c,STAT=ierr)
    endif
    NULLIFY(this%c)
    this%deg  = -1
    this%size = -1
    return
  end subroutine gf2x_delete

  subroutine gf2x_print_bit(this)
    type(gf2x_obj), intent(in) :: this
    integer(INT32) :: i,ib,iw,deg
    deg = get_deg(this)
    if (deg < 0) then
       write(*,'("0")')
       return
    endif
    do i=deg,0,-1
       ib = mod(i,32)
       iw = i/32
       if (BTEST(this%c(iw),ib)) then
          write(*,'("1",$)')
       else
          write(*,'("0",$)')
       endif
    enddo
    write(*,'("")')
    return
  end subroutine gf2x_print_bit

  subroutine gf2x_print_hex(this)
    type(gf2x_obj), intent(in) :: this
    integer(INT32) :: i,ib,iw,isize
    character(9) :: str
    if (is_zero(this)) then
       write(*,'("0")')
       return
    endif
    isize = get_size(this%deg)
    i = isize-1
    write(str,'(Z8)')this%c(i)
    write(*,'(A,$)')TRIM(ADJUSTL(str))
    do i=isize-2,0,-1
       write(*,'(Z8.8,$)')this%c(i)
    enddo
    write(*,'("")')
    return
  end subroutine gf2x_print_hex

  subroutine gf2x_assign(c,a)
    type(gf2x_obj), intent(inout) :: c  ! c := a
    type(gf2x_obj), intent(in)    :: a
    integer(INT32) :: ia,isa,i

    call delete(c)
    if (is_zero(a)) then
       return
    endif

    ia = get_deg(a)
    isa = get_size(ia)
    call new(c,ia)
    do i=0,isa-1
       c%c(i) = a%c(i)
    enddo

    return
  end subroutine gf2x_assign

  function is_zero(a) result(is)
    type(gf2x_obj), intent(in) :: a
    logical :: is
    integer(INT32) :: deg
    deg = get_deg(a)
    if (deg==-1) then
       is = .true.
    else
       is = .false.
    endif
    return
  end function is_zero

  !!DEC$ ATTRIBUTES FORCEINLINE :: get_deg
  function get_deg(a) result(deg)
    type(gf2x_obj), intent(in) :: a
    integer(INT32) :: deg
    integer(INT32) :: isize,i,top_deg
    intrinsic :: SIZE
    deg=-1
    if (.not.associated(a%c)) return
    isize = SIZE(a%c)
    do i=isize-1,0,-1
       if (a%c(i) /= 0) then
          top_deg = deg_i32(a%c(i))
          deg = 32*i + top_deg
          return
       endif
    enddo
    return
  end function get_deg

  subroutine gf2x_set_coef(a,i)
    type(gf2x_obj), intent(inout) :: a
    integer(INT32), intent(in) :: i
    type(gf2x_obj), pointer :: w
    integer(INT32) :: ib,iw
    NULLIFY(w)
    if (is_zero(a)) then
       call new(a,i)
    endif
    allocate(w)
    call new(w,i)
    iw =     i/32
    ib = mod(i,32)
    w%c(iw) = ibset(w%c(iw),ib)
    call add_assign(w,a)  ! w := w + a
    call assign(a,w)      ! a := w
    call delete(w)
    deallocate(w)
    NULLIFY(w)
    a%deg = get_deg(a)
    return
  end subroutine gf2x_set_coef

  subroutine gf2x_add_assign(c,a)
    type(gf2x_obj), intent(inout) :: c  ! c := c + a
    type(gf2x_obj), intent(in)    :: a
    type(gf2x_obj), pointer :: w
    integer(INT32) :: ia,ic
    integer(INT32) :: isa,isc,i
    if (is_zero(a)) then
       return
    endif
    if (is_zero(c)) then
       call assign(c,a)
       return
    endif
    ia = a%deg
    ic = c%deg
    isa = a%size
    isc = c%size
    if (isc < isa) then
       NULLIFY(w)
       allocate(w)
       call new(w,MAX(ia,ic))
       do i=0,isc-1
          w%c(i) = IEOR(c%c(i),a%c(i))
       enddo
       do i=isc,isa-1
          w%c(i) = a%c(i)
       enddo
       call assign(c,w)
       call delete(w)
       deallocate(w)
       NULLIFY(w)
    else
       do i=0,isa-1
          c%c(i) = IEOR(c%c(i),a%c(i))
       enddo
       c%deg = get_deg(c)
       c%size = get_size(c%deg)
    endif
    return
  end subroutine gf2x_add_assign

  subroutine gf2x_add(c,a,b)
    type(gf2x_obj), intent(inout) :: c   ! c := a + b
    type(gf2x_obj), intent(in)    :: a,b
    integer(INT32) :: ia,ib,ic
    integer(INT32) :: isa,isb,isc,i
    if (is_zero(a) .and. is_zero(b)) then
       return
    endif
    if (is_zero(a)) then
       call assign(c,b)
       return
    endif
    if (is_zero(b)) then
       call assign(c,a)
       return
    endif
    ia = get_deg(a)
    ib = get_deg(b)
    isa = get_size(ia)
    isb = get_size(ib)
    if (c%deg < MAX(ia,ib)) call new(c,MAX(ia,ib))
    if (isa < isb) then
       do i=0,isa-1
          c%c(i) = IEOR(a%c(i),b%c(i))
       enddo
       do i=isa,isb-1
          c%c(i) = b%c(i)
       enddo
    else
       do i=0,isb-1
          c%c(i) = IEOR(a%c(i),b%c(i))
       enddo
       do i=isb,isa-1
          c%c(i) = a%c(i)
       enddo
    endif
    c%deg  = get_deg(c)
    c%size = get_size(c%deg)
    return
  end subroutine gf2x_add

  subroutine gf2x_pow(c,a,e)
    type(gf2x_obj), intent(inout) :: c ! c = a**e
    type(gf2x_obj), intent(in)    :: a
    integer(INT32), intent(in) :: e
    type(gf2x_obj), pointer :: w
    integer(INT32) :: ch,cl
    integer(INT32) :: i,deg
    NULLIFY(w)
    call delete(c)
    if (e==1) then
       call assign(c,a)
       return
    endif
    if (e==0) then
       call set_coef(c,0)
       return
    endif
    if (e<0) then
       write(*,*)"pow: c = a^e : exponent should be e>=0."
       stop
    endif
    if (is_zero(a)) return

    deg = deg_i32(e)

    allocate(w)
    call set_coef(c,0)
    do i=deg,0,-1
       call square(w,c)        ! w := c**2
       if (BTEST(e,i)) then
          call mult(c,w,a)      ! c := w * a
       else
          call assign(c,w)      ! c := w
       endif
    enddo
    call delete(w)
    deallocate(w)
    NULLIFY(w)

    return
  end subroutine gf2x_pow

  subroutine gf2x_square(c,a)
    type(gf2x_obj), intent(inout) :: c ! c := a**2
    type(gf2x_obj), intent(in)    :: a
    integer(INT32) :: ch,cl
    integer(INT32) :: i,deg
    call delete(c)
    if (is_zero(a)) return
    deg = a%deg*2
    call new(c,deg)
    do i=0,a%size-1
       if (a%c(i) == 0) cycle
       call square_i32(a%c(i),ch,cl)
       if (cl /= 0) c%c(2*i)   = IEOR(c%c(2*i),  cl)
       if (ch /= 0) c%c(2*i+1) = IEOR(c%c(2*i+1),ch)
    enddo
    c%deg = get_deg(c)
    c%size = get_size(c%deg)
    return
  end subroutine gf2x_square

  recursive subroutine gf2x_mult_kara(c,a,b)
    !
    ! multiply 2 polyomials using Karatsuba algorithm
    !
    type(gf2x_obj), intent(inout) :: c    ! c := a * b
    type(gf2x_obj), intent(in)    :: a,b
    type(gf2x_obj), pointer :: ah,al,bh,bl,ahbh,albl,ahl,bhl,ahlbhl
    integer(INT32) :: isa,isb,isc
    integer(INT32) :: i,j,deg

    NULLIFY(ah,al,bh,bl,ahbh,albl,ahl,bhl,ahlbhl)
    call delete(c)
    if (is_zero(a)) return
    if (is_zero(b)) return

    isa = a%size
    isb = b%size
    isc = MAX(isa,isb)
    if (isc < MAX_KARA) then
       call gf2x_mult_normal(c,a,b)
       return
    endif

    if (mod(isc,2)/=0) then
       isc = isc + 1
    endif

    allocate(ah,al,bh,bl,ahbh,albl,ahl,bhl,ahlbhl)
    deg = 32*(isc/2)-1
    call new(al,deg)
    call new(bl,deg)
    call new(ah,deg)
    call new(bh,deg)

    do i=0,MIN(isc/2-1,isa-1)
       al%c(i) = a%c(i)
    enddo
    do i=0,MIN(isc/2-1,isb-1)
       bl%c(i) = b%c(i)
    enddo
    do i=isc/2,isa-1
       ah%c(i-isc/2) = a%c(i)
    enddo
    do i=isc/2,isb-1
       bh%c(i-isc/2) = b%c(i)
    enddo
    ah%deg = get_deg(ah)
    al%deg = get_deg(al)
    bh%deg = get_deg(bh)
    bl%deg = get_deg(bl)
    ah%size = get_size(ah%deg)
    al%size = get_size(al%deg)
    bh%size = get_size(bh%deg)
    bl%size = get_size(bl%deg)

    !===================================

    call add(ahl,ah,al)
    call add(bhl,bh,bl)
    call gf2x_mult_kara(ahlbhl,ahl,bhl)
    call delete(ahl)
    call delete(bhl)
    deallocate(ahl,bhl)

    !===================================

    call gf2x_mult_kara(ahbh,ah,bh)
    call delete(ah)
    call delete(bh)
    deallocate(ah,bh)

    call add_assign(ahlbhl,ahbh)

    !===================================

    call gf2x_mult_kara(albl,al,bl)
    call delete(al)
    call delete(bl)
    deallocate(al,bl)

    call add_assign(ahlbhl,albl)

    !===================================
    deg = a%deg + b%deg
    call new(c,deg)

    do i=0,MIN(c%size,albl%size)-1
       c%c(i) = albl%c(i)
    enddo
    call delete(albl)

    if (.not. is_zero(ahlbhl)) then
       do i=isc/2,MIN(c%size,isc/2+ahlbhl%size)-1
          c%c(i) = IEOR(c%c(i),ahlbhl%c(i-isc/2))
       enddo
    endif
    call delete(ahlbhl)

    if (.not. is_zero(ahbh)) then
       do i=isc,MIN(c%size,isc+ahbh%size)-1
          c%c(i) = IEOR(c%c(i),ahbh%c(i-isc))
       enddo
    endif
    call delete(ahbh)
    deallocate(ahbh,albl,ahlbhl)
    NULLIFY(ah,al,bh,bl,ahbh,albl,ahl,bhl,ahlbhl)
    c%deg  = get_deg(c)
    c%size = get_size(c%deg)
    return
  end subroutine gf2x_mult_kara

  subroutine gf2x_mult_assign_kara(a,b)
    type(gf2x_obj), intent(inout) :: a  ! a := a * b
    type(gf2x_obj), intent(in)    :: b
    type(gf2x_obj), pointer :: w
    NULLIFY(w)
    if (is_zero(a)) then
       call delete(a)
       return
    endif
    if (is_zero(b)) then
       call delete(a)
       return
    endif
    allocate(w)
    call gf2x_mult_kara(w,a,b)
    call assign(a,w)
    call delete(w)
    deallocate(w)
    NULLIFY(w)
    return
  end subroutine gf2x_mult_assign_kara

  subroutine gf2x_mult_assign_normal(a,b)
    type(gf2x_obj), intent(inout) :: a  ! a := a * b
    type(gf2x_obj), intent(in)    :: b
    type(gf2x_obj), pointer :: w
    integer(INT32) :: ch,cl
    integer(INT32) :: i,j,deg
    NULLIFY(w)
    allocate(w)
    deg = a%deg + b%deg
    call new(w,deg)
    call gf2x_mult_normal(w,a,b)
    call assign(a,w)  ! a := w
    call delete(w)
    deallocate(w)
    NULLIFY(w)
    return
  end subroutine gf2x_mult_assign_normal

  subroutine gf2x_mult_normal(c,a,b)
    type(gf2x_obj), intent(inout) :: c    ! c := a * b
    type(gf2x_obj), intent(in)    :: a,b
    integer(INT32) :: ch,cl
    integer(INT32) :: i,j,ij,deg,kk,mm
    integer(INT32), allocatable :: hi(:,:),lo(:,:)

    call delete(c)
    if (is_zero(a) .or. is_zero(b) ) then
       return
    endif

    deg = a%deg + b%deg
    call new(c,deg)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do j=0,c%size-2
       kk = MIN(j,  a%size-1)
       mm = MAX(0,j-b%size+1)
       do i=mm,kk
          call mult_i32(a%c(i),b%c(j-i),ch,cl)
          c%c(j)   = IEOR(c%c(j),  cl)
          c%c(j+1) = IEOR(c%c(j+1),ch)
       enddo
    enddo
    j=c%size-1
    kk = a%size-1
    mm = c%size-b%size
    do i=mm,kk
       call mult_i32(a%c(i),b%c(j-i),ch,cl)
       c%c(j)   = IEOR(c%c(j),  cl)
    enddo

    !  do j=0,b%size-1
    !  if (b%c(j) == 0) cycle
    !  do i=0,a%size-1
    !  if (a%c(i) == 0) cycle

    !    ij = i + j
    !    call mult_i32(a%c(i),b%c(j),ch,cl)
    !                     c%c(ij)   = IEOR(c%c(ij),  cl)
    !    if (ij+1<c%size) c%c(ij+1) = IEOR(c%c(ij+1),ch)

    !  enddo
    !  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    c%deg = get_deg(c)
    c%size = get_size(c%deg)

    return
  end subroutine gf2x_mult_normal

  subroutine gf2x_shift(c,a,i)
    type(gf2x_obj), intent(inout) :: c  ! c := shift(a,i)
    type(gf2x_obj), intent(in)    :: a
    integer(INT32), intent(in) :: i
    integer(INT32) :: j,isn,iw,ib,ida,isa,ch,cm,cl
    if (i==0) then
       call assign(c,a)
       return
    endif
    ida = get_deg(a)
    isa = get_size(ida)
    if (ida + i < 0) then
       call delete(c)
       return
    endif
    iw = abs(i)/32
    ib = mod(abs(i),32)
    call delete(c)
    call new(c,ida+i)
    if (i > 0) then
       do j=0,isa-1
          call shift_i32(a%c(j),+ib,ch,cm,cl)
          if (ch /= 0) c%c(j+iw+1) = IEOR(c%c(j+iw+1),ch)
          if (cm /= 0) c%c(j+iw)   = IEOR(c%c(j+iw)  ,cm)
       enddo
    else
       call shift_i32(a%c(iw),-ib,ch,cm,cl)
       if (cm /= 0) c%c(0)   = IEOR(c%c(0),cm)
       do j=iw+1,isa-1
          call shift_i32(a%c(j),-ib,ch,cm,cl)
          if (cm /= 0) c%c(j-iw)   = IEOR(c%c(j-iw)  ,cm)
          if (cl /= 0) c%c(j-iw-1) = IEOR(c%c(j-iw-1),cl)
       enddo
    endif
    c%deg = get_deg(c)
    c%size = get_size(c%deg)
    return
  end subroutine gf2x_shift

  subroutine gf2x_divrem(q,r,a,b)
    ! a =: q * b + r
    type(gf2x_obj), intent(inout) :: q  ! q := a div b
    type(gf2x_obj), intent(inout) :: r  ! r := a mod b
    type(gf2x_obj), intent(in)    :: a,b
    type(gf2x_obj), pointer :: w,t,s
    integer(INT32) :: ida,idb,idw
    call delete(q)
    call delete(r)
    ida = a%deg
    idb = b%deg
    if (ida < idb) then
       call assign(r,a)
       return
    endif
    NULLIFY(w,t,s)
    allocate(w,t,s)
    call assign(w,a)
    idw = w%deg
    do
       call mult_by_x(t,b,idw-idb)  ! t := b * x^(deg(w)-deg(b))
       call set_coef(s,idw-idb)     ! s := s + x^(deg(w)-deg(b))
       call add_assign(w,t)         ! w := w + t
       call delete(t)
       idw = w%deg
       if (idw < idb) exit
    enddo
    call assign(r,w)
    call delete(w)
    call assign(q,s)
    call delete(s)
    deallocate(w,t,s)
    NULLIFY(w,t,s)
    return
  end subroutine gf2x_divrem

  subroutine gf2x_div(q,a,b)
    ! a =: q * b + r
    type(gf2x_obj), intent(inout) :: q  ! q := a div b
    type(gf2x_obj), intent(in)    :: a,b
    type(gf2x_obj), pointer :: w,t,s
    integer(INT32) :: ida,idb,idw
    call delete(q)
    ida = a%deg
    idb = b%deg
    if (ida < idb) then
       return
    endif
    NULLIFY(w,t,s)
    allocate(w,t,s)
    call assign(w,a)
    idw = w%deg
    do
       call mult_by_x(t,b,idw-idb)  ! t := b * x^(deg(w)-deg(b))
       call set_coef(s,idw-idb)     ! s := s + x^(deg(w)-deg(b))
       call add_assign(w,t)         ! w := w + t
       call delete(t)
       idw = w%deg
       if (idw < idb) exit
    enddo
    call delete(w)
    call assign(q,s)
    call delete(s)
    deallocate(w,t,s)
    NULLIFY(w,t,s)
    return
  end subroutine gf2x_div

  subroutine gf2x_rem(r,a,b)
    type(gf2x_obj), intent(inout) :: r   ! r := a mod b
    type(gf2x_obj), intent(in)    :: a,b
    type(gf2x_obj), pointer :: w,t
    integer(INT32) :: ida,idb,idw
    call delete(r)
    ida = a%deg
    idb = b%deg
    if (ida < idb) then
       call assign(r,a)
       return
    endif
    NULLIFY(w,t)
    allocate(w,t)
    call assign(w,a)
    idw = w%deg
    do
       call mult_by_x(t,b,idw-idb)  ! t := b * x^(deg(w)-deg(b))
       call add_assign(w,t)         ! w := w + t
       call delete(t)
       idw = w%deg
       if (idw < idb) exit
    enddo
    call assign(r,w)
    call delete(w)
    deallocate(w,t)
    NULLIFY(w,t)
    return
  end subroutine gf2x_rem

  subroutine gf2x_set_prime(mp,m)
    !
    ! Set a primitive polynomial to the cotainer.
    ! the container contains the prime poly and precomputed polynomial for Barrett reduciont.
    ! This routine does not check the primitivity.
    !  mp : container
    !   m : primitive polynomial
    !
    type(gf2x_prime_obj), intent(inout) :: mp
    type(gf2x_obj),       intent(in)    :: m
    type(gf2x_obj), pointer :: xx
    integer(INT32) :: deg
    call delete(mp)
    call assign(mp%prime_poly,m)
    deg = get_deg(m)
    mp%deg = deg
    NULLIFY(xx)
    allocate(xx)
    call set_coef(xx,2*deg)
    call div(mp%barrett_poly,xx,m)
    call delete(xx)
    deallocate(xx)
    NULLIFY(xx)
    return
  end subroutine gf2x_set_prime

  subroutine gf2x_delete_prime(mp)
    type(gf2x_prime_obj), intent(inout) :: mp
    call delete(mp%prime_poly)
    call delete(mp%barrett_poly)
    return
  end subroutine gf2x_delete_prime

  subroutine gf2x_rem_barrett(r,a,m)
    !
    ! compute  r := a mod m using Barrett algorithm
    !
    type(gf2x_obj), intent(inout) :: r     ! r := a mod m
    type(gf2x_obj), intent(in)    :: a
    type(gf2x_prime_obj), intent(in) :: m  ! precomputed polynomial for Barrett algorithm
    type(gf2x_obj), pointer :: q,p
    integer(INT32) :: deg

    call delete(r)
    deg = m%deg
    if (a%deg < deg) then
       call assign(r,a)
       return
    endif

    NULLIFY(q,p)
    allocate(q,p)

    call div_by_x(q,a,deg)              ! q = a  /  x**deg
    call mod_by_x(r,a,deg)              ! r = a mod x**deg

    call mult_assign(q,m%barrett_poly)  ! q = q  *  mu
    call div_by_x(p,q,deg)              ! p = q  /  x**deg

    call mult_assign(p,m%prime_poly)    ! p = p  *  m
    call mod_by_x(q,p,deg)              ! q = p mod x**deg

    call add_assign(r,q)                ! r = r + q

    call delete(p)
    call delete(q)

    deallocate(p,q)
    NULLIFY(p,q)

    return
  end subroutine gf2x_rem_barrett

  subroutine gf2x_mod_by_x(c,a,i)
    type(gf2x_obj), intent(inout) :: c  ! c := a mod x^i
    type(gf2x_obj), intent(in)    :: a
    integer(INT32), intent(in) :: i
    type(gf2x_obj), pointer :: w
    integer(INT32) :: iw,ib,j
    call delete(c)
    if (a%deg < i) then
       call assign(c,a)
       return
    endif
    if (i == 0) then
       call delete(c)
       return
    endif
    if (i < 0) then
       write(*,'("mod_by_x: error, negative i:",I10)')i
       stop
    endif
    iw = i/32
    ib = mod(i,32)
    NULLIFY(w)
    allocate(w)
    call new(w,i)
    do j=0,w%size-1
       w%c(j) = a%c(j)
    enddo
    w%c(w%size-1) = IAND(w%c(w%size-1),2**ib-1)
    call assign(c,w)
    call delete(w)
    deallocate(w)
    NULLIFY(w)
    return
  end subroutine gf2x_mod_by_x

  subroutine gf2x_mult_by_x(c,a,i)
    type(gf2x_obj), intent(inout) :: c  ! c := a * x^i
    type(gf2x_obj), intent(in)    :: a
    integer(INT32), intent(in) :: i
    if (i < 0) then
       write(*,'("mult_by_x: error, negative i:",I10)')i
       stop
    endif
    if (i == 0) then
       call assign(c,a)
       return
    endif
    call shift(c,a,i)
    return
  end subroutine gf2x_mult_by_x

  subroutine gf2x_div_by_x(c,a,i)
    type(gf2x_obj), intent(inout) :: c  ! c := a div x^i
    type(gf2x_obj), intent(in)    :: a
    integer(INT32), intent(in) :: i
    if (i < 0) then
       write(*,'("div_by_x: error, negative i:",I10)')i
       stop
    endif
    if (i == 0) then
       call assign(c,a)
       return
    endif
    call shift(c,a,-i)
    return
  end subroutine gf2x_div_by_x


  subroutine gf2x_pow_pow_2(c,e,m)
    type(gf2x_obj), intent(inout) :: c  ! c := x**(2**e) mod m
    integer(INT32), intent(in)           :: e
    type(gf2x_prime_obj), intent(in) :: m  ! precomputed polynomial for Barrett algorithm
    integer(INT32) :: i,ee
    type(gf2x_obj), pointer :: w,s

    ee = CEILING(log(REAL(m%deg))/log(2.0))
    call delete(c)
    if (ee > e) then
       call set_coef(c,2**e)
       return
    endif

    NULLIFY(w,s)
    allocate(w,s)
    call set_coef(w,2**ee)
    call rem(s,w,m)      ! s = w mod m
    do i=ee+1,e
       call square(w,s)   ! w = s**2
       call rem(s,w,m)    ! s = w mod m
    enddo
    call assign(c,s)
    call delete(w)
    call delete(s)
    deallocate(w,s)
    NULLIFY(w,s)
    return
  end subroutine gf2x_pow_pow_2

  subroutine gf2x_pow_mod(c,a,e,m)
    type(gf2x_obj), intent(inout) :: c  ! c := a**e mod m
    type(gf2x_obj), intent(in)    :: a
    integer(INT32), intent(in)           :: e
    type(gf2x_prime_obj), intent(in) :: m  ! precomputed polynomial for Barrett algorithm
    type(gf2x_obj), pointer :: w
    integer(INT32) :: i,deg
    NULLIFY(w)
    call delete(c)
    if (e==1) then
       if (a%deg >= m%deg) then
          call rem(c,a,m)
          return
       else
          call assign(c,a)
          return
       endif
    endif
    if (e==0) then
       call set_coef(c,0)
       return
    endif
    if (e<0) then
       write(*,*)"pow: c = a^e mod m : exponent should be e>=0."
       stop
    endif
    if (is_zero(a)) return

    deg = deg_i32(e)

    allocate(w)
    call set_coef(c,0)
    do i=deg,0,-1
       call square(w,c)        ! c := c**2 mod m
       call rem(c,w,m)
       if (BTEST(e,i)) then
          call mult(w,c,a)      ! c := c * a mod m
          call rem(c,w,m)
       endif
    enddo
    call delete(w)
    deallocate(w)
    NULLIFY(w)
    return
  end subroutine gf2x_pow_mod

  !========================================================================
  function deg_i32(a) result(d)
    integer(INT32) :: a,d,i
    d=-1
    do i=31,0,-1
       if (BTEST(a,i)) then
          d=i
          exit
       endif
    enddo
    return
  end function deg_i32

  function deg_i64(a) result(d)
    integer(INT64) :: a
    integer(INT32) :: d,i
    do i=63,0,-1
       if (BTEST(a,i)) then
          d=i
          exit
       endif
    enddo
    return
  end function deg_i64

  subroutine square_i32(a,ch,cl)
    integer(INT32), intent(in) :: a
    integer(INT32), intent(out) :: ch,cl   ! (ch,cl) = a**2
    integer(INT32) :: ia,i
    integer(INT64) :: da,dc
    da = a
    if (da < 0) da = da + 2_8**32 ! convert to unsigned
    dc = int(Z'0', int64)
    ia = deg_i32(a)
    do i = 0,ia
       if (BTEST(a,i)) then
          dc = ibset(dc,i*2)
       endif
    enddo
    ch = ISHFT(dc,-32)
    cl = dc
    return
  end subroutine square_i32

  !DEC$ ATTRIBUTES FORCEINLINE :: mult_i32
  subroutine mult_i32(a,b,ch,cl)

    integer(INT32), intent(in) :: a,b
    integer(INT32), intent(out) :: ch,cl  ! (ch,cl) = a*b
    integer(INT32) :: tmp,u(0:3)
    integer(INT64), parameter :: ZE = int(Z'eeeeeeee', int64)
    integer(INT64), parameter :: ZC = int(Z'cccccccc', int64)
    integer(INT64), parameter :: Z8 = int(Z'88888888', int64)

    if (a==0 .or. b ==0) then
       ch = 0
       cl = 0
       return
    endif

    u(0) = 0
    u(1) = a
    u(2) = ISHFT(u(1),+1)
    u(3) =  IEOR(u(2),a)

    cl =                  IEOR(ISHFT(u(     ISHFT(b,-30)   ),2),u(IAND(ISHFT(b,-28),3)))
    ch =                  ISHFT(cl,-28)
    cl = IEOR(ISHFT(cl,4),IEOR(ISHFT(u(IAND(ISHFT(b,-26),3)),2),u(IAND(ISHFT(b,-24),3))))
    ch =  IOR(ISHFT(ch,4),ISHFT(cl,-28))
    cl = IEOR(ISHFT(cl,4),IEOR(ISHFT(u(IAND(ISHFT(b,-22),3)),2),u(IAND(ISHFT(b,-20),3))))
    ch =  IOR(ISHFT(ch,4),ISHFT(cl,-28))
    cl = IEOR(ISHFT(cl,4),IEOR(ISHFT(u(IAND(ISHFT(b,-18),3)),2),u(IAND(ISHFT(b,-16),3))))
    ch =  IOR(ISHFT(ch,4),ISHFT(cl,-28))
    cl = IEOR(ISHFT(cl,4),IEOR(ISHFT(u(IAND(ISHFT(b,-14),3)),2),u(IAND(ISHFT(b,-12),3))))
    ch =  IOR(ISHFT(ch,4),ISHFT(cl,-28))
    cl = IEOR(ISHFT(cl,4),IEOR(ISHFT(u(IAND(ISHFT(b,-10),3)),2),u(IAND(ISHFT(b, -8),3))))
    ch =  IOR(ISHFT(ch,4),ISHFT(cl,-28))
    cl = IEOR(ISHFT(cl,4),IEOR(ISHFT(u(IAND(ISHFT(b, -6),3)),2),u(IAND(ISHFT(b, -4),3))))
    ch =  IOR(ISHFT(ch,4),ISHFT(cl,-28))
    cl = IEOR(ISHFT(cl,4),IEOR(ISHFT(u(IAND(ISHFT(b, -2),3)),2),u(IAND(      b     ,3))))

    tmp = -IAND(ISHFT(a,-31),1)
    tmp =  IAND(INT(tmp, INT64),ISHFT(IAND(INT(b, INT64),ZE),-1))
    ch  =  IEOR(ch,tmp)

    tmp = -IAND(ISHFT(a,-30),1)
    tmp =  IAND(INT(tmp, INT64),ISHFT(IAND(INT(b, INT64),ZC),-2))
    ch  =  IEOR(ch,tmp)

    tmp = -IAND(ISHFT(a,-29),1)
    tmp =  IAND(INT(tmp, INT64),ISHFT(IAND(INT(b, INT64),Z8),-3))
    ch  =  IEOR(ch,tmp)

    return
  end subroutine mult_i32


  subroutine mult_i32_old(a,b,ch,cl)
    integer(INT32), intent(in) :: a,b
    integer(INT32), intent(out) :: ch,cl  ! (ch,cl) = a*b
    integer(INT32) :: ia,ib,i
    integer(INT64) :: da,db,dc
    da = a
    db = b
    if (da < 0) da = da + 2_8**32 ! convert to unsigned
    if (db < 0) db = db + 2_8**32 ! convert to unsigned
    ia = deg_i32(a)
    ib = deg_i32(b)
    dc = int(Z'0', int64)
    do i = 0,ia
       if (BTEST(a,i)) then
          dc = IEOR(dc,db)
       endif
       dc = ISHFTC(dc,-1)
    enddo
    dc = ISHFTC(dc,ia+1)
    ch = ISHFT(dc,-32)
    cl = dc
    !  write(*,'(B64.64)')dc
    !  write(*,'(B32.32)')ch
    !  write(*,'(B64.64)')cl
    return
  end subroutine mult_i32_old

  subroutine shift_i32(a,i,ch,cm,cl)
    integer(INT32), intent(in) :: a
    integer(INT32), intent(in) :: i
    integer(INT32), intent(out) :: ch,cm,cl  ! (ch,cm,cl) = shift(a,i)
    integer(INT64) :: dc
    if (abs(i) >= 32) then
       write(*,*)"shift_int32: error i=",i
       stop
    endif
    select case (i)
    case (0)
       ch = 0; cm = a; cl = 0
       return
    case (1:31)
       dc = a
       if (dc < 0) dc = dc + 2_8**32 ! convert to unsigned
       dc = ISHFT(dc,i)
       ch = ISHFT(dc,-32)
       cm = dc
       cl = 0
       return
    case (-31:-1)
       dc = a
       if (dc < 0) dc = dc + 2_8**32 ! convert to unsigned
       dc = ISHFT(dc,i+32)
       ch = 0
       cm = ISHFT(dc,-32)
       cl = dc
       return
    end select
    return
  end subroutine shift_i32

end module

module msmt19937

  !-------------------------------------------------------------------------------
  !   This is a Fortran translation from C-program for MT19937-64
  !   (2004/9/29 version)
  !   originally coded by Takuji Nishimura and Makoto Matsumoto
  !   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html
  !   done originally by Rﾃｩmi Piatek. However, this code have the
  !   following modifications:
  !
  ! * Save state to restart PRNG from where it was interrupted
  ! * JUMP AHEAD prng, allowing to split single random number series
  !   from 64bit Mersenne Twister (MT19937_64) into (almost) multiple independent
  !   streams. Therefore, splitted streams can be used on parallel simulations
  !   with MPI. Code for jump ahead is based on mt_stream.f90 module
  !   originally coded by Ken-Ichi Ishikawa
  !   http://theo.phys.sci.hiroshima-u.ac.jp/~ishikawa/PRNG/mt_stream_en.html

  !   Before generate any pseudorandom number, initialize the state by using
  !       call init_genrand(seed)
  !   or
  !       call init_by_array(init_key)
  !   To jump ahead pseudorandom sequence by id*2^jp steps
  !       call mt_jumpahead(id,jp)
  !
  ! Tayroni Alves
  !    Physics Department
  !    Piauﾃｭ State Federal University
  !    email: tay @ ufpi.edu.br
  !-------------------------------------------------------------------------------
  !   This is a Fortran translation of the 64-bit version of
  !   the Mersenne Twister pseudorandom number generator
  !   Before using, initialize the state by using
  !       call init_genrand64(seed)
  !   or
  !       call init_by_array64(init_key)
  !   Translated from C-program for MT19937-64 (2004/9/29 version)
  !   originally coded by Takuji Nishimura and Makoto Matsumoto
  !   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html
  !   Fortran translation by Rﾃｩmi Piatek
  !   The University of Copenhagen
  !   Department of Economics
  !   email: {first}.{last}@econ.ku.dk
  !-------------------------------------------------------------------------------
  !   A C-program for MT19937-64 (2004/9/29 version).
  !   Coded by Takuji Nishimura and Makoto Matsumoto.
  !   This is a 64-bit version of Mersenne Twister pseudorandom number
  !   generator.
  !   Before using, initialize the state by using init_genrand64(seed)
  !   or init_by_array64(init_key, key_length).
  !   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
  !   All rights reserved.
  !   Redistribution and use in source and binary forms, with or without
  !   modification, are permitted provided that the following conditions
  !   are met:
  !     1. Redistributions of source code must retain the above copyright
  !        notice, this list of conditions and the following disclaimer.
  !     2. Redistributions in binary form must reproduce the above copyright
  !        notice, this list of conditions and the following disclaimer in the
  !        documentation and/or other materials provided with the distribution.
  !     3. The names of its contributors may not be used to endorse or promote
  !        products derived from this software without specific prior written
  !        permission.
  !   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  !   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  !   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  !   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
  !   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
  !   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
  !   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
  !   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
  !   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
  !   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  !   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  !   References:
  !   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
  !     ACM Transactions on Modeling and
  !     Computer Simulation 10. (2000) 348--357.
  !   M. Matsumoto and T. Nishimura,
  !     ``Mersenne Twister: a 623-dimensionally equidistributed
  !       uniform pseudorandom number generator''
  !     ACM Transactions on Modeling and
  !     Computer Simulation 8. (Jan. 1998) 3--30.
  !   Any feedback is very welcome.
  !   http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
  !   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
  !-------------------------------------------------------------------------------

  use gf2xe
  use iso_fortran_env

  implicit none

  private

  public :: mt_jumpahead
  public :: savestate_genrand
  public :: readstate_genrand
  public :: erasestate_genrand
  public :: init_genrand
  public :: init_by_array
  public :: grnd
  public :: grnd1
  public :: grnd2
  public :: grndint

  integer, parameter :: n = 312
  integer, parameter :: m = 156

  integer*8, parameter :: seed_def = 5489_int64
  integer*8, parameter :: matrix_a = -5403634167711393303_int64
  integer*8, parameter :: umask    = -2147483648_int64 ! most significant 33 bits
  integer*8, parameter :: lmask    =  2147483647_int64 ! least significant 31 bits
  integer*8, parameter :: maska    = 6148914691236517205_int64
  integer*8, parameter :: maskb    = 8202884508482404352_int64
  integer*8, parameter :: maskc    = -2270628950310912_int64

  integer, parameter :: shift0   = -29_int64
  integer, parameter :: shiftb   = 17_int64
  integer, parameter :: shiftc   = 37_int64
  integer, parameter :: shift1   = -43_int64

  real*8, parameter :: pi253_1  = 1.d0/(2.d0**53 - 1.d0)
  real*8, parameter :: pi253    = 1.d0/(2.d0**53)
  real*8, parameter :: pi252    = 1.d0/(2.d0**52)

  integer(int64), public, protected :: mt(n)       ! array for the state vector
  integer(int32), public, protected :: mti = n+1   ! mti==n+1 means mt is not initialized

  integer, parameter :: unitnumber1 = 100001 !Save mt PRNG state on unit=unitnumber1
  integer, parameter :: unitnumber2 = 100002 !Save mti PRNG state on unit=unitnumber2

contains

  !-----------------------------------------------------------------------------
  ! Jump ahead by id*2^jp steps
  !-----------------------------------------------------------------------------

  subroutine mt_jumpahead(id,jp)

    implicit none

    integer :: jp ! exponent (jump step = id*2^jp)
    integer :: id ! id       (jump step = id*2^jp)
    integer*8 p(n) !jumpahead polynomial coefficients
    integer np
    integer i,iwp,ibp
    integer*8 :: v(n),w(n),s(n)

    ! compute jump ahead polynomial
    ! p(x) coefficients for this MT parameter
    call f_get_coeff(n,m,31,32,matrix_a,id,jp,np,p)

    ! multiply p(B) on a state vector v
    !  p(x) : jump ahead polynomial
    !    B  : transition matrix
    ! with simple Horner's method
    !       w = p(B) v = B^(2^jp) v
    v = mt
    iwp = (np-1)/32
    ibp = mod(np-1,32)

    if (btest(p(iwp+1),ibp)) then
       w = v
    end if

    do i = np-2,0,-1
       iwp = i/32
       ibp = mod(i,32)
       call mt_matvec(w,s)! s = B w
       if (btest(p(iwp+1),ibp)) then
          w = ieor(v,s)   ! w = 1 v + s
       else
          w = s           ! w = 0 v + s
       end if
    end do

    mt = w
    mti = n

    return

  end subroutine mt_jumpahead

  !-----------------------------------------------------------------------------
  ! Multiply transition matrix on a state vector v
  !-----------------------------------------------------------------------------

  subroutine mt_matvec(old,new)

    !  new = B old
    !  this : MT parameters(transition matrix)
    !     old : input vector
    !     new : output vector

    implicit none

    integer*8 old(n)
    integer*8 new(n)
    integer*8 :: mag01(0:1) = (/0_int64, matrix_a/)
    integer*8 x
    integer i

    new(1) = iand(old(2),umask)
    new(2:n-1) = old(3:n)
    x = ior(iand(old(1), umask), iand(old(2), lmask))
    new(n) = ieor(ieor(old(m+1), ishft(x, -1_int64)), mag01(iand(x, 1_int64)))

    return

  end subroutine mt_matvec

  !-----------------------------------------------------------------------------
  ! Compute MT jump ahead polynomial coefficients (depends on the gf2xe module)
  !-----------------------------------------------------------------------------

  subroutine f_get_coeff(nn,mm,rr,ww,avec,id,jp,np,pp)

    implicit none

    integer id,jp
    integer :: nn,mm,rr,ww
    integer*8 :: avec
    integer np
    integer*8 pp(nn)
    type(gf2x_obj) :: af,bf,ff,f1,f2
    type(gf2x_prime_obj) :: fp
    integer :: i,ib,nws

    !MT characteristic polynomial
    !ff : MT char poly.
    call set_coef(af,nn)
    call set_coef(af,mm)   ! af = x^nn + x^mm
    call set_coef(bf,nn-1)
    call set_coef(bf,mm-1) ! bf = x^(nn-1) + x^(mm-1)

    call pow(f1,af,ww-rr)  ! f1 = af^(ww-rr)
    call pow(f2,bf,rr)     ! f2 = bf^(rr)
    call mult(ff,f1,f2)    ! ff = f1*f2

    do i=0,rr-1
       ib = mod(i,ww)
       if (btest(avec,ib)) then
          call pow(f2,bf,rr-1-i)
          call mult_assign(f2,f1)
          call add_assign(ff,f2)
       endif
    enddo

    do i=rr,ww-1
       ib = mod(i,ww)
       if (btest(avec,ib)) then
          call pow(f1,af,ww-1-i)
          call add_assign(ff,f1)
       endif
    enddo

    call delete(af)
    call delete(bf)
    call delete(f1)
    call delete(f2)

    !set ff for Barrett reduction
    call set_prime(fp,ff)  ! fp = ff
    call delete(ff)
    call delete(f1)
    call delete(f2)

    !jump ahead

    !long jump
    call pow_pow_2(f1,jp,fp)  ! f1 = x**(2**jp) mod fp

    !short jump
    call pow_mod(ff,f1,id,fp) ! ff = f1**id mod fp

    pp(:) = 0
    np = get_deg(ff)+1
    nws = ceiling(dfloat(np)/32.d0)
    pp(1:nws) = ff%c(0:nws-1)

    call delete(f1)
    call delete(f2)
    call delete(ff)
    call delete(fp)

    return

  end subroutine f_get_coeff

  !-----------------------------------------------------------------------------
  ! Save a copy of PRNG state
  !-----------------------------------------------------------------------------

  subroutine savestate_genrand(statenumber)

    implicit none

    integer statenumber !PRNG state copy number
    character*11, parameter :: formatstatenumber='(I4.4)'
    character(len=4) stringstatenumber

    write(stringstatenumber,formatstatenumber) statenumber

    open(unit=unitnumber1,file='statemt_'//stringstatenumber//'.dat')
    open(unit=unitnumber2,file='statemti_'//stringstatenumber//'.dat')

    write(unitnumber1,*) mt
    write(unitnumber2,*) mti

    close(unitnumber1)
    close(unitnumber2)

    return

  end subroutine savestate_genrand

  !-----------------------------------------------------------------------------
  ! Read a copy of PRNG state
  !-----------------------------------------------------------------------------

  subroutine readstate_genrand(statenumber)

    implicit none

    integer statenumber !PRNG state copy number
    character*11, parameter :: formatstatenumber='(I4.4)'
    character(len=4) stringstatenumber
    logical existingfile1, existingfile2

    write(stringstatenumber,formatstatenumber) statenumber

    !See if a 'statenumber' copy is saved
    inquire(file='statemt_'//stringstatenumber//'.dat' , exist=existingfile1)
    inquire(file='statemti_'//stringstatenumber//'.dat', exist=existingfile2)

    !If state files exists, then they will be read
    if (existingfile1 .and. existingfile2) then

       open(unit=unitnumber1,file='statemt_'//stringstatenumber//'.dat' )
       open(unit=unitnumber2,file='statemti_'//stringstatenumber//'.dat')

       read(unitnumber1,*) mt
       read(unitnumber2,*) mti

       close(unitnumber1)
       close(unitnumber2)

    end if

    return

  end subroutine readstate_genrand

  !-----------------------------------------------------------------------------
  ! Erase a copy of PRNG state
  !-----------------------------------------------------------------------------

  subroutine erasestate_genrand(statenumber)

    implicit none

    integer statenumber
    character*11, parameter :: formatstatenumber='(I4.4)'
    character(len=4) stringstatenumber
    integer stat
    logical existingfile1, existingfile2

    write(stringstatenumber,formatstatenumber) statenumber

    !See if a 'statenumber' copy is saved
    inquire(file='statemt_'//stringstatenumber//'.dat' , exist=existingfile1)
    inquire(file='statemti_'//stringstatenumber//'.dat', exist=existingfile2)

    !If state files exists, erase then
    if (existingfile1 .and. existingfile2) then

       open(unit=unitnumber1, iostat=stat, file='statemt_'//stringstatenumber//'.dat' , status='old')
       if (stat .eq. 0) close(unitnumber1, status='delete')

       open(unit=unitnumber2, iostat=stat, file='statemti_'//stringstatenumber//'.dat', status='old')
       if (stat .eq. 0) close(unitnumber2, status='delete')

    end if

    return

  end subroutine erasestate_genrand

  !-----------------------------------------------------------------------------
  ! Initializes mt(n) with a seed
  !-----------------------------------------------------------------------------

  subroutine init_genrand(seed)

    implicit none

    integer*8, intent(in) :: seed
    integer               :: i

    mt(1) = seed
    do i = 1, n-1
       mt(i+1) = 6364136223846793005_int64 * ieor(mt(i), ishft(mt(i), -62_int64)) + i
    end do

    mti = n

    return

  end subroutine init_genrand

  !-----------------------------------------------------------------------------
  ! Initializes by an array with array-length - init_key is the array for initializing keys
  !-----------------------------------------------------------------------------

  subroutine init_by_array(init_key)

    implicit none

    integer*8, intent(in) :: init_key(:)
    integer*8, parameter  :: c1 = 3935559000370003845_int64
    integer*8, parameter  :: c2 = 2862933555777941757_int64
    integer*8             :: i, j, k, kk, key_length

    call init_genrand(19650218_int64)

    key_length = size(init_key)
    i = 1; j = 0_int64
    k = max(n, key_length)

    do kk = 1, k

       mt(i+1) = ieor(mt(i+1), c1 * ieor(mt(i), ishft(mt(i), -62_int64))) &
            & + init_key(j+1) + j

       i = i+1; j = j+1_int64

       if (i >= n) then
          mt(1) = mt(n)
          i = 1
       end if

       if(j >= key_length) j = 0_int64

    end do

    do kk = 1, n-1

       mt(i+1) = ieor(mt(i+1), c2 * ieor(mt(i), ishft(mt(i), -62_int64))) - i

       i = i+1

       if (i >= n) then
          mt(1) = mt(n)
          i = 1
       end if

    end do

    mt(1) = ishft(1_int64, 63_int64)  ! MSB is 1; assuring non-zero initial array

    return

  end subroutine init_by_array

  !-----------------------------------------------------------------------------
  ! Generates a random number on [-2^63, 2^63-1]-interval
  !-----------------------------------------------------------------------------

  function grndint() result(x)

    implicit none

    integer*8 :: mag01(0:1) = (/0_int64, matrix_a/)
    integer*8 :: x
    integer   :: i

    if (mti >= n) then ! generate nn words at one time

       ! if init_genrand() has not been called, a default initial seed is used
       if(mti == n+1) call init_genrand(seed_def)

       do i = 1, n-m
          x = ior(iand(mt(i),umask), iand(mt(i+1), lmask))
          mt(i) = ieor(ieor(mt(i+m), ishft(x, -1_int64)), mag01(iand(x, 1_int64)))
       end do

       do i = n-m+1, n-1
          x = ior(iand(mt(i), umask), iand(mt(i+1), lmask))
          mt(i) = ieor(ieor(mt(i+m-n), ishft(x, -1_int64)), mag01(iand(x, 1_int64)))
       end do

       x = ior(iand(mt(n), umask), iand(mt(1), lmask))
       mt(n) = ieor(ieor(mt(m), ishft(x, -1_int64)), mag01(iand(x, 1_int64)))

       mti = 0

    end if

    mti = mti + 1
    x = mt(mti)

    x = ieor(x, iand(ishft(x, shift0), maska))
    x = ieor(x, iand(ishft(x, shiftb), maskb))
    x = ieor(x, iand(ishft(x, shiftc), maskc))
    x = ieor(x, ishft(x, shift1))

    return

  end function grndint

  !-----------------------------------------------------------------------------
  ! Generates a random number on [0,1) real interval
  !-----------------------------------------------------------------------------

  function grnd() result(r)

    implicit none

    real*8 :: r

    r = dfloat(ishft(grndint(), -11_int64)) * pi253

  end function grnd

  !-----------------------------------------------------------------------------
  ! Generates a random number on [0,1] real interval
  !-----------------------------------------------------------------------------

  function grnd1() result(r)

    implicit none

    real*8 :: r

    r = dfloat(ishft(grndint(), -11_int64)) * pi253_1

    return

  end function grnd1

  !-----------------------------------------------------------------------------
  ! Generates a random number on (0,1) real interval
  !-----------------------------------------------------------------------------

  function grnd2() result(r)

    implicit none

    real*8 :: r

    r = dfloat(ishft(grndint(), -12_int64))
    r = (r + 0.5d0) * pi252

    return

  end function grnd2

end module msmt19937

module xy2d_periodic_correlation_simple_m
  use, intrinsic :: iso_fortran_env
  use msmt19937
  implicit none
  private
  real(real64), parameter :: pi = 4 * atan(1d0)
  integer(int64), parameter, public :: nx = 5000, ny = nx, nall = nx * ny
  real(real64), parameter :: nall_inv = 1.0d0 / nall
  real(real64), parameter, public :: kbt = 0.88d0, beta = 1/ kbt
  real(real64), allocatable, public, protected :: spins(:, :, :)
  real(real64), allocatable, private :: rnds(:, :, :)

  integer(int64), allocatable :: stack(:, :) !> for Wolff. stack(x or y or z, maximum size)
  logical, allocatable :: cluster(:, :) !> for Wolff. cluster(position of x, position of y, position of z)

  integer(int32), parameter :: nd = 4
  integer(int64), parameter :: dy(nd) = [1_int64, 0_int64, -1_int64, 0_int64]
  integer(int64), parameter :: dx(nd) = [0_int64, 1_int64, 0_int64, -1_int64]

  integer(int64), public, parameter :: n_unique_correlation = 2
  real(real64), public, protected :: unique_dist(n_unique_correlation)

  public :: init_simulation_xy2d, init_lattice_order, update_metropolis, calc_energy, calc_magne, calc_magne_abs
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
    integer(int64) :: near_x, near_y
    integer(int64) :: cxy(1:2), s
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
       cxy(1:2) = stack(1:2, s)
       s = s - 1
       r_sigma_c = sum(reflection_spin(1:2) * spins(1:2, cxy(1), cxy(2)))
       !> right
       near_x = cxy(1) + 1
       if (near_x > nx) near_x = 1
       if (.not. cluster(near_x, cxy(2))) then
          call wolff_local_update(stack, spins, cluster, s, reflection_spin, r_sigma_c, near_x, cxy(2))
       end if
       !> left
       near_x = cxy(1) - 1
       if (near_x < 1) near_x = nx
       if (.not. cluster(near_x, cxy(2))) then
          call wolff_local_update(stack, spins, cluster, s, reflection_spin, r_sigma_c, near_x, cxy(2))
       end if
       !> up
       near_y = cxy(2) + 1
       if (near_y > ny) near_y = 1
       if (.not. cluster(cxy(1), near_y)) then
          call wolff_local_update(stack, spins, cluster, s, reflection_spin, r_sigma_c, cxy(1), near_y)
       end if
       !> down
       near_y = cxy(2) - 1
       if (near_y < 1) near_y = ny
       if (.not. cluster(cxy(1), near_y)) then
          call wolff_local_update(stack, spins, cluster, s, reflection_spin, r_sigma_c, cxy(1), near_y)
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

program xy2d_periodic_wolff_equilibrium_simulation
  use, intrinsic :: iso_fortran_env
  use gf2xe
  use msmt19937
  use xy2d_periodic_correlation_simple_m
  implicit none
  integer(int32), parameter :: mcs_discard = 100000_int32, mcs_sample = 100000_int32
  integer(int64), parameter :: iseed = 42_int64
  real(real64) :: magne, magne_square, magne_sus, energy, energy_square, specific_heat
  integer(int32) :: i, j

  call init_simulation_xy2d()
  call init_genrand(iseed)

  associate(out => [output_unit, error_unit])
    do i = 1, size(out)
       write(out(i), '(a,i0)'    ) "# Nsize: ", nall
       write(out(i), '(2(a, i0))') "# nx: ", nx, " ny: ", ny
       write(out(i), '(2(a, i0))') "# MCS_discard: ", mcs_discard, " MCS_sample: ", mcs_sample
       write(out(i), '(a, g0)' ) "# 温度: ", kbt
       write(out(i), '(a)' ) "# method: Wolff, model: 2D-XY"
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
    real(real64) :: theta(nx, ny)
    integer(int64) :: i, j, k
    open(file = filename, newunit = myunit)
    do j = 1, ny
       do k = 1, nx
          theta(k, j) = atan2(spins(2, k, j), spins(1, k, j))
       end do
    end do
    write(myunit, '(*(g0, 1x))') theta(:, :)
    close(myunit)
  end block
end program xy2d_periodic_wolff_equilibrium_simulation
