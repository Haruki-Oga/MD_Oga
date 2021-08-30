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

function dim_char(cdim)
    implicit none
    character(1) cdim
    integer dim_char
    dim_char = 0
    if(cdim=="x") dim_char = 1
    if(cdim=="y") dim_char = 2
    if(cdim=="z") dim_char = 3
    if(dim_char==0) call error_msg("error@dim_char")
end function dim_char

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

subroutine inverse(a,c,n)
    !============================================================
    ! Inverse matrix
    ! Method: Based on Doolittle LU factorization for Ax=b
    ! Alex G. December 2009
    !-----------------------------------------------------------
    ! input ...
    ! a(n,n) - array of coefficients for matrix A
    ! n      - dimension
    ! output ...
    ! c(n,n) - inverse matrix of A
    ! comments ...
    ! the original matrix a(n,n) will be destroyed
    ! during the calculation
    !===========================================================
    implicit none
    integer n
    double precision a(n,n), c(n,n)
    double precision L(n,n), U(n,n), b(n), d(n), x(n)
    double precision coeff
    integer i, j, k

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L=0.0
    U=0.0
    b=0.0

    ! step 1: forward elimination
    do k=1, n-1
        do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
                a(i,j) = a(i,j)-coeff*a(k,j)
            end do
        end do
    end do

    ! Step 2: prepare L and U matrices
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
        L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
        do i=1,j
            U(i,j) = a(i,j)
        end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
        b(k)=1.0
        d(1) = b(1)
        ! Step 3a: Solve Ld=b using the forward substitution
        do i=2,n
            d(i)=b(i)
            do j=1,i-1
                d(i) = d(i) - L(i,j)*d(j)
            end do
        end do
        ! Step 3b: Solve Ux=d using the back substitution
        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
            x(i) = d(i)
            do j=n,i+1,-1
                x(i)=x(i)-U(i,j)*x(j)
            end do
            x(i) = x(i)/u(i,i)
        end do
        ! Step 3c: fill the solutions x(n) into column k of C
        do i=1,n
            c(i,k) = x(i)
        end do
        b(k)=0.0
    end do
end subroutine inverse
