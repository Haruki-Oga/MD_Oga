! langevin thermostat (dim, atype)
module fix_NoseHoover_mod
    integer,allocatable :: myfc(:)
    integer,parameter :: ndata = 1
    integer :: nlist = 0
    integer,allocatable :: atlist(:)
    integer(1),allocatable :: dimlist(:,:)
    double precision,allocatable :: tsetNH(:), QNH(:), zeta(:)
    double precision,allocatable ::  fNH(:,:),fNHbuf(:,:)
    !
contains
    subroutine init_fix_NoseHoover
        use update_mod
        use commn
        use file_mod
        use constant_mod
        use fix_compute_mod
        implicit none
        integer i,iatype,ndim,int3(1:3)
        integer,allocatable :: a(:)
        logical(1) dim_true
        double precision tsetNH_in, Q_in
        !
        ! === read conf file ===
        read(ifilecalc,*) iatype, ndim, tsetNH_in, Q_in
        if(Q_in<=0d0) call error_msg("Q<0@"//__FILE__)
        write(*,'(A,2I,2F)') __FILE__,iatype, ndim, tsetNH_in, Q_in
        int3(:) = 0
        do i=1,3
            if(dim_true(ndim,i)) int3(i) = 1
        end do
        !
        atlist = [atlist, iatype]
        a = reshape([dimlist(:,:)], [3*nlist])
        dimlist = reshape( [a,int3], [3,nlist+1] )
        !
        tsetNH = [tsetNH,tsetNH_in]
        QNH = [QNH, Q_in]
        zeta = [zeta, 0d0]
        if(nlist==0)then
            allocate(fNH(1:3,1:nmol),fNHbuf(1:3,1:nmol))
            fNH(:,:) = 0d0
            fNHbuf(:,:) = 0d0
        end if
        !
        nlist = nlist + 1
        ! --- set fixflag ---
        do i=1,3
            if(dimlist(i,nlist)==1)then
                if(pfixflag(i,atlist(nlist))) call error_msg("pfixflag true @"//__FILE__)
                if(ffixflag(i,atlist(nlist))) call error_msg("ffixflag true @"//__FILE__)
                if(vfixflag(i,atlist(nlist))) call error_msg("vfixflag true @"//__FILE__)
                ffixflag(i,atlist(nlist)) = .true.
            end if
        end do
        ! --- return ---
        myfc = [myfc, nfc]
        nfcdata = [nfcdata, ndata]
        ! ===
    end subroutine init_fix_NoseHoover
    !
    ! --- langevin thermostat ---
    subroutine fix_NoseHoover
        use update_mod
        use commn
        use constant_mod
        implicit none
        integer j,i, idim, inum,imol
        double precision work_now(1:nlist), r83(3), temp_now
        !
        if( .not. fflag ) return
        !
        fNHbuf(:,:) = fNH(:,:)
        do j=1,nlist
            temp_now = 0d0
            do i=1,natseq(atlist(j))
                imol = imolatseq(i0atseq(atlist(j)) + i)
                temp_now = temp_now + sum(v(:,imol)**2d0)
            end do
            temp_now = temp_now*wm(atlist(j))/3d0/bk
            do i=1,natseq(atlist(j))
                imol = imolatseq(i0atseq(atlist(j)) + i)
                ! --- langevin thermostat ---
                zeta(j) = zeta(j) + 3d0*natseq(atlist(j))*bk/2d0/QNH(j)*(temp_now - tsetNH(j))*dt
                fNH(:,imol) =  - dimlist(:,j)*wm(atlist(j))*zeta(j)*v(:,imol)
                f(:,imol) = f(:,imol) + fNH(:,imol)
                ! ---
            end do
        end do
    end subroutine fix_NoseHoover
    !
    ! --- compute culumtive enerty ---
    subroutine fix_NoseHoover_compute
        use update_mod
        use commn
        use constant_mod
        use fix_compute_mod
        implicit none
        integer j,i, idim, inum,imol
        double precision random, work_now(1:nlist)
        !
        !
        work_now(:) = 0d0
        do j=1,nlist
            do i=1,natseq(atlist(j))
                imol = imolatseq(i0atseq(atlist(j)) + i)
                ! --- calculate culumtive energy ---
                !work_now(j) = work_now(j) +&
                !    fNH(idim,imol)*v(idim,imol)*dt
                work_now(j) = work_now(j) + &
                    sum( 0.5*(fNH(:,imol) + fNHbuf(:,imol))*dp(:,imol) )
            end do
        end do
        eadd = eadd + sum(work_now(:))
        ! --- return ---
        do i=1,nlist
            fcdata(i0fcdata(myfc(i))+1) = work_now(i)/dt
        end do
    end subroutine fix_NoseHoover_compute
end module fix_NoseHoover_mod

