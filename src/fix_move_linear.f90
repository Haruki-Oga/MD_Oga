! fix atom (dim, atype)
module fix_move_linear_mod
    integer,allocatable :: myfc(:)
    integer,parameter :: ndata = 1
    integer :: nlist = 0
    integer,allocatable :: atlist(:)
    integer(1),allocatable :: dimlist(:,:)
    double precision,allocatable :: vmove(:,:)
    !
contains
    subroutine init_fix_move_linear
        use update_mod
        use commn
        use file_mod
        use fix_compute_mod
        implicit none
        integer i,iatype,ndim,int3(1:3)
        integer,allocatable :: a(:)
        double precision,allocatable :: b(:)
        logical(1) dim_true
        double precision vmove_in(3)
        !
        ! === read conf file ===
        read(ifilecalc,*) iatype, ndim, vmove_in(:)
        write(*,'(A,2I,3F)') __FILE__,iatype, ndim, vmove_in(:)
        int3(:) = 0
        do i=1,3
            if(dim_true(ndim,i)) int3(i) = 1
        end do
        !
        atlist = [atlist, iatype]
        a = reshape([dimlist(:,:)], [3*nlist])
        dimlist = reshape( [a,int3], [3,nlist+1] )
        !
        b = reshape([vmove(:,:)], [3*nlist])
        vmove = reshape( [b,vmove_in], [3,nlist+1] )
        !
        nlist = nlist + 1
        ! --- set fixflag ---
        do i=1,3
            if(dimlist(i,nlist)==1)then
                if(pfixflag(i,atlist(nlist))) call error_msg("pfixflag true @"//__FILE__)
                if(ffixflag(i,atlist(nlist))) call error_msg("ffixflag true @"//__FILE__)
                if(vfixflag(i,atlist(nlist))) call error_msg("vfixflag true @"//__FILE__)
                pfixflag(i,atlist(nlist)) = .true.
                vfixflag(i,atlist(nlist)) = .true.
            end if
        end do
        ! --- return ---
        myfc = [myfc, nfc]
        nfcdata = [nfcdata, ndata]
        ! ===
    end subroutine init_fix_move_linear
    !
    subroutine fix_move_linear
        use update_mod
        use commn
        implicit none
        integer i,j,imol
        !
        if( fflag ) return
        !
        do i=1,nlist
            do j=1,natseq(atlist(i))
                imol = imolatseq(i0atseq(atlist(i)) + j)
                ! --- freeze atom ---
                ! dimlist = 1 or 0
                if( pflag ) dp(:,imol) = dimlist(:,i)*vmove(:,i)*dt + (1-dimlist(:,i))*dp(:,imol)
                if( vflag ) v(:,imol) = dimlist(:,i)*vmove(:,i) + (1-dimlist(:,i))*v(:,imol)
                ! ---
            end do
        end do
    end subroutine fix_move_linear
    ! --- compute culumtive enerty ---
    subroutine fix_move_linear_compute
        use update_mod
        use commn
        use constant_mod
        use fix_compute_mod
        implicit none
        integer j,i, idim, inum,imol
        double precision work_now(1:nlist)
        !
        !
        work_now(:) = 0d0
        do i=1,nlist
            do j=1,natseq(atlist(i))
                imol = imolatseq(i0atseq(atlist(i)) + j)
                ! --- calculate culumtive energy ---
                work_now(i) = work_now(i) - sum( dimlist(:,i)*0.5*(fbuf(:,imol) + f(:,imol))*dp(:,imol) )
            end do
        end do
        eadd = eadd + sum(work_now(:))
        ! --- return ---
        do i=1,nlist
            fcdata(i0fcdata(myfc(i))+1) = work_now(i)/dt
        end do
    end subroutine fix_move_linear_compute
    !
end module fix_move_linear_mod

