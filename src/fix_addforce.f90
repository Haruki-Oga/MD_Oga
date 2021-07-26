! fix atom (dim, atype)
module fix_addforce_mod
    integer,allocatable :: myfc(:)
    integer :: nlist = 0
    integer,allocatable :: atlist(:)
    integer(1),allocatable :: dimlist(:,:)
    integer,parameter :: ndata = 1
    double precision,allocatable :: fadd(:,:)
    !
contains
    subroutine init_fix_addforce
        use update_mod
        use commn
        use file_mod
        use fix_compute_mod
        implicit none
        integer i,iatype,ndim,int3(1:3)
        integer,allocatable :: a(:)
        double precision,allocatable :: b(:)
        logical(1) dim_true
        double precision fadd_in(3)
        !
        ! === read conf file ===
        read(ifilecalc,*) iatype, ndim, fadd_in(:)
        write(*,'(A,2I,3F)') __FILE__, iatype, ndim, fadd_in(:)
        int3(:) = 0
        do i=1,3
            if(dim_true(ndim,i)) int3(i) = 1
        end do
        !
        atlist = [atlist, iatype]
        a = reshape([dimlist(:,:)], [3*nlist])
        dimlist = reshape( [a,int3], [3,nlist+1] )
        !
        b = reshape([fadd(:,:)], [3*nlist])
        fadd_in(:) = fadd_in(:)*int3(:)
        fadd = reshape( [b,fadd_in], [3,nlist+1] )
        !
        nlist = nlist + 1
        ! --- set fixflag ---
        do i=1,3
            if(dimlist(i,nlist)==1)then
                if(pfixflag(i,atlist(nlist))) call error_msg("pfixflag true @"//__FILE__)
                if(vfixflag(i,atlist(nlist))) call error_msg("vfixflag true @"//__FILE__)
                ffixflag(i,atlist(nlist)) = .true.
            end if
        end do
        ! --- return ---
        myfc = [myfc, nfc]
        nfcdata = [nfcdata, ndata]
        ! ===
    end subroutine init_fix_addforce
    !
    subroutine fix_addforce
        use update_mod
        use commn
        implicit none
        integer i,j,imol
        !
        if( .not. fflag ) return
        !
        do i=1,nlist
            do j=1,natseq(atlist(i))
                imol = imolatseq(i0atseq(atlist(i)) + j)
                ! --- add force ---
                ! dimlist = 1 or 0
                f(:,imol) = f(:,imol) + dimlist(:,i)*fadd(:,i)
                ! ---
            end do
        end do
    end subroutine fix_addforce
    !
    subroutine fix_addforce_compute
        use update_mod
        use commn
        use fix_compute_mod
        implicit none
        double precision work_now(1:nlist)
        integer i,j,imol
        !
        work_now(:) = 0d0
        do i=1,nlist
            do j=1,natseq(atlist(i))
                imol = imolatseq(i0atseq(atlist(i)) + j)
                ! --- calculate work ---
                ! dimlist = 1 or 0
                work_now(i) = work_now(i) + sum( fadd(:,i)*dp(:,imol) )
                ! ---
            end do
        end do
        !
        eadd = eadd + sum(work_now(:))
        ! --- return ---
        do i=1,nlist
            fcdata(i0fcdata(myfc(i))+1) = work_now(i)/dt
        end do
    end subroutine fix_addforce_compute
    !
end module fix_addforce_mod

