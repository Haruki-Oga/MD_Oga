! langevin thermostat (dim, atype)
module fix_vel_scaling_mod
    integer,allocatable :: myfc(:)
    integer,parameter :: ndata = 1
    integer :: nlist = 0
    integer,allocatable :: atlist(:)
    integer(1),allocatable :: dimlist(:,:)
    double precision,allocatable :: tsetL(:), temp_now(:)
    !
contains
    subroutine init_fix_vel_scaling
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
        temp_now = [temp_now,0d0]
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
    end subroutine init_fix_vel_scaling
    !
    ! --- langevin thermostat ---
    subroutine fix_vel_scaling
        use update_mod
        use commn
        use constant_mod
        implicit none
        integer j,i, idim, inum,imol
        double precision work_now(1:nlist)
        !
        if( .not. vflag ) return
        !
        ! calculate temperature
        do j=1,nlist
            temp_now(j) = 0d0
            do i=1,natseq(atlist(j))
                imol = imolatseq(i0atseq(atlist(j)) + i)
                temp_now(j) = temp_now(j) + sum(v(:,imol)**2.0)
            end do
            temp_now(j) = temp_now(j)*wm(atlist(j))/(3.0*bk)/natseq(atlist(j))
            ! velocity scaling
            do i=1,natseq(atlist(j))
                imol = imolatseq(i0atseq(atlist(j)) + i)
                ! --- velocity scaling ---
                v(:,imol) = dimlist(:,j)*sqrt(tsetL(j)/temp_now(j))*v(:,imol) + (1-dimlist(:,j))*v(:,imol)
                ! ---
            end do
        end do
    end subroutine fix_vel_scaling
    !
    ! --- compute culumtive enerty ---
    subroutine fix_vel_scaling_compute
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
                work_now(j) = work_now(j) + 0.5d0*wm(atlist(j))*sum(dimlist(:,j)*v(:,imol)**2d0*(1d0 - sqrt(temp_now(j)/tsetL(j))))
            end do
        end do
        eadd = eadd + sum(work_now(:))
        ! --- return ---
        do i=1,nlist
            fcdata(i0fcdata(myfc(i))+1) = work_now(i)/dt
        end do
    end subroutine fix_vel_scaling_compute
end module fix_vel_scaling_mod

