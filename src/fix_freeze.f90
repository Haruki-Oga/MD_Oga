! fix atom (dim, atype)
module fix_freeze_mod
    integer,allocatable :: myfc(:)
    integer :: nlist = 0
    integer,allocatable :: atlist(:)
    integer(1),allocatable :: dimlist(:,:)
    integer,parameter :: ndata = 0
    !
contains
    subroutine init_fix_freeze
        use update_mod
        use commn
        use file_mod
        use fix_compute_mod
        implicit none
        integer ndata,i,iatype,ndim,int3(1:3)
        integer,allocatable :: a(:)
        logical(1) dim_true
        !
        ! === read conf file ===
        read(ifilecalc,*) iatype, ndim
        write(*,'(A,I,I)') __FILE__, iatype, ndim
        int3(:) = 0
        do i=1,3
            if(dim_true(ndim,i)) int3(i) = 1
        end do
        !
        atlist = [atlist, iatype]
        a = reshape([dimlist(:,:)], [3*nlist])
        dimlist = reshape( [a,int3], [3,nlist+1] )
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
    end subroutine init_fix_freeze

    subroutine fix_freeze
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
                if( pflag ) dp(:,imol) = dimlist(:,i)*0d0 + (1-dimlist(:,i))*dp(:,imol)
                if( vflag ) v(:,imol) = dimlist(:,i)*0d0 + (1-dimlist(:,i))*v(:,imol)
                ! ---
            end do
        end do
    end subroutine fix_freeze
end module fix_freeze_mod

