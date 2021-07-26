! langevin thermostat (dim, atype)
module fix_langevin_mod
    integer,allocatable :: myfc(:)
    integer,parameter :: ndata = 1
    integer :: nlist = 0
    integer,allocatable :: atlist(:)
    integer(1),allocatable :: dimlist(:,:)
    double precision,allocatable :: tsetL(:), gamma(:), sigR(:)
    double precision,allocatable ::  fL(:,:),fLbuf(:,:)
    !
contains
    subroutine init_fix_langevin
        use update_mod
        use commn
        use file_mod
        use constant_mod
        use fix_compute_mod
        implicit none
        integer i,iatype,ndim,int3(1:3)
        integer,allocatable :: a(:)
        logical(1) dim_true
        double precision tsetL_in, gamma_in, sigR_in
        !
        ! === read conf file ===
        read(ifilecalc,*) iatype, ndim, tsetL_in, gamma_in
        write(*,'(A,2I,2F)') __FILE__,iatype, ndim, tsetL_in, gamma_in
        int3(:) = 0
        do i=1,3
            if(dim_true(ndim,i)) int3(i) = 1
        end do
        !
        atlist = [atlist, iatype]
        a = reshape([dimlist(:,:)], [3*nlist])
        dimlist = reshape( [a,int3], [3,nlist+1] )
        !
        gamma_in = wm(iatype)/gamma_in
        sigR_in = sqrt(2d0*gamma_in*bk*tsetL_in/dt)
        tsetL = [tsetL,tsetL_in]
        gamma = [gamma,gamma_in]
        sigR = [sigR,sigR_in]
        if(nlist==0)then
            allocate(fL(1:3,1:nmol),fLbuf(1:3,1:nmol))
            fL(:,:) = 0d0
            fLbuf(:,:) = 0d0
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
    end subroutine init_fix_langevin
    !
    ! --- langevin thermostat ---
    subroutine fix_langevin
        use update_mod
        use commn
        implicit none
        integer j,i, idim, inum,imol
        double precision work_now(1:nlist), r83(3)
        !
        if( .not. fflag ) return
        !
        fLbuf(:,:) = fL(:,:)
        do j=1,nlist
            do i=1,natseq(atlist(j))
                imol = imolatseq(i0atseq(atlist(j)) + i)
                ! --- langevin thermostat ---
                call random3(r83,sigR(j),0.0d0)
                fL(:,imol) =  dimlist(:,j)*&
                    ( - gamma(j)*v(:,imol) + r83(:) )
                f(:,imol) = f(:,imol) + fL(:,imol)
                ! ---
            end do
        end do
    end subroutine fix_langevin
    !
    ! --- compute culumtive enerty ---
    subroutine fix_langevin_compute
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
                !    fL(idim,imol)*v(idim,imol)*dt
                work_now(j) = work_now(j) + &
                    sum( 0.5*(fL(:,imol) + fLbuf(:,imol))*dp(:,imol) )
            end do
        end do
        eadd = eadd + sum(work_now(:))
        ! --- return ---
        do i=1,nlist
            fcdata(i0fcdata(myfc(i))+1) = work_now(i)/dt
        end do
    end subroutine fix_langevin_compute
end module fix_langevin_mod

