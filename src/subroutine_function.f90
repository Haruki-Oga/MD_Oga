module constant_mod
    double precision, parameter :: pi = 4d0*atan(1d0)
    double precision, parameter :: bk = 1d0
end module constant_mod

subroutine error_msg(crrchar)
    character(32) crrchar
    write(*,*) crrchar
    stop
end subroutine error_msg

function char_cmt(char1, cmtchar, ncmtchar)
    character(128) char_cmt, char1
    integer ncmtchar
    character(ncmtchar) cmtchar
    if(index(char1,cmtchar)==0)then
        char_cmt = trim(char1)
    else
        char_cmt = trim(char1(1:index(char1,cmtchar)-1))
    end if
    return
end function char_cmt

subroutine random_seed_clock
    implicit none
    integer :: seedsize,i
    integer,allocatable :: seed(:)

    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    do i=1,seedsize
        call system_clock(count=seed(i))
    end do
    call random_seed(put=seed(:))
    !write(*,*) 'seedsize=',seedsize
    !write(*,*) 'seed=',seed(:)
end subroutine random_seed_clock


function random(sig,ave) ! Box-Muller's method
    use constant_mod
    implicit none
    double precision random
    double precision,intent(in) :: sig,ave
    integer i
    double precision x,y
    call random_number(x)
    call random_number(y)
    random = dsqrt(-2d0*log(x))*cos(2d0*pi*y)*sig + ave
end function random

subroutine random3(random,sig,ave)! Box-Muller's method
    implicit none
    real(8) random(3),sig,ave
    integer i
    real(8) x,y,a,b,pi
    pi = 4.0d0*atan(1.0d0)
    call random_number(x)
    call random_number(y)
    call random_number(a)
    call random_number(b)
    random(1) = dsqrt(-2.0d0*log(x))*cos(2.0d0*pi*y)*sig + ave
    random(2) = dsqrt(-2.0d0*log(x))*sin(2.0d0*pi*y)*sig + ave
    random(3) = dsqrt(-2.0d0*log(a))*cos(2.0d0*pi*b)*sig + ave
end subroutine random3

!function random3(sig,ave) ! Box-Muller's method
!    use constant_mod
!    implicit none
!    double precision random3(3)
!    double precision,intent(in) :: sig,ave
!    integer i
!    double precision x,y,a,b
!    call random_number(x)
!    call random_number(y)
!    call random_number(a)
!    call random_number(b)
!    random3(1) = dsqrt(-2.0d0*log(x))*cos(2.0d0*pi*y)*sig + ave
!    random3(2) = dsqrt(-2.0d0*log(x))*sin(2.0d0*pi*y)*sig + ave
!    random3(3) = dsqrt(-2.0d0*log(a))*cos(2.0d0*pi*b)*sig + ave
!end function random3

function char_dim(idim)
    implicit none
    character(1) char_dim
    integer idim
    if(idim==1) char_dim = "x"
    if(idim==2) char_dim = "y"
    if(idim==3) char_dim = "z"
end function char_dim

function dim_true(n,idim)
    implicit none
    logical(1) dim_true
    integer n,idim
    if(idim==1)then
        if(n/100==1)then
            dim_true = .true.
        else
            dim_true = .false.
        end if
    else if(idim==2)then
        if(mod(n,100)/10==1)then
            dim_true = .true.
        else
            dim_true = .false.
        end if
    else if(idim==3)then
        if(mod(n,10)==1)then
            dim_true = .true.
        else
            dim_true = .false.
        end if
    end if
    return
end function

subroutine switchint(a, b)
    implicit none
    integer a,b,c
    c = a
    a = b
    b = c
end subroutine switchint
