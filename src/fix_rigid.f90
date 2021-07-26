! fix atom (dim, atype)
module fix_rigid_mod
    integer,allocatable :: myfc(:)
    integer,parameter :: ndata = 1
    integer :: nlist = 0
    integer,allocatable :: atlist(:)
    integer(1),allocatable :: dimlist(:,:)
    double precision,allocatable :: fsum(:,:), fsumbuf(:,:), dpcom(:,:), vcom(:,:)
    double precision,allocatable :: mass(:)
    !
contains
    subroutine init_fix_rigid
        use update_mod
        use commn
        use file_mod
        use fix_compute_mod
        implicit none
        integer i,iatype,ndim,int3(1:3)
        integer,allocatable :: a(:)
        double precision,allocatable :: b(:)
        double precision vcomin(3)
        integer imol
        logical(1) dim_true
        !
        ! === read conf file ===
        read(ifilecalc,*) iatype, ndim
        write(*,*) __FILE__, iatype, ndim
        int3(:) = 0
        do i=1,3
            if(dim_true(ndim,i)) int3(i) = 1
        end do
        !
        atlist = [atlist, iatype]
        a = reshape([dimlist(:,:)], [3*nlist])
        dimlist = reshape( [a,int3], [3,nlist+1] )
        !
        b = reshape([fsum(:,:)], [3*nlist])
        fsum = reshape( [b,0d0,0d0,0d0], [3,nlist+1] )
        fsumbuf = reshape( [b,0d0,0d0,0d0], [3,nlist+1] )
        dpcom = reshape( [b,0d0,0d0,0d0], [3,nlist+1] )
        vcom = reshape( [b,0d0,0d0,0d0], [3,nlist+1] )
        mass = [mass, 0d0]
        nlist = nlist + 1
        do i=1,natseq(iatype)
            imol = imolatseq(i0atseq(iatype) + i)
            vcom(:,nlist) = vcom(:,nlist) + v(:,imol)
        end do
        vcom(:,nlist) = vcom(:,nlist)/natseq(iatype)
        mass(nlist) = natseq(iatype)*wm(iatype)
        !
        ! --- set fixflag ---
        do i=1,3
            if(dimlist(i,nlist)==1)then
                if(pfixflag(i,atlist(nlist))) call error_msg("pfixflag true @"//__FILE__)
                if(vfixflag(i,atlist(nlist))) call error_msg("vfixflag true @"//__FILE__)
                pfixflag(i,atlist(nlist)) = .true.
                vfixflag(i,atlist(nlist)) = .true.
            end if
        end do
        ! --- return ---
        myfc = [myfc, nfc]
        nfcdata = [nfcdata, ndata]
        ! ===
    end subroutine init_fix_rigid
    !
    subroutine fix_rigid
        use update_mod
        use commn
        implicit none
        integer i,j,imol
        !
        if( pflag )then
            do i=1,nlist
                dpcom(:,i) = dt*(vcom(:,i) + fsum(:,i)*dt/mass(i)*0.5d0)
                !
                do j=1,natseq(atlist(i))
                    imol = imolatseq(i0atseq(atlist(i)) + j)
                    dp(:,imol) = dimlist(:,i)*dpcom(:,i) + (1-dimlist(:,i))*dp(:,imol)
                end do
            end do
        end if
        !
        if( fflag )then
            fsumbuf(:,:) = fsum(:,:)
            fsum(:,:) = 0d0
            do i=1,nlist
                do j=1,natseq(atlist(i))
                    imol = imolatseq(i0atseq(atlist(i)) + j)
                    fsum(:,i) = fsum(:,i) + f(:,imol)
                end do
                !
                do j=1,natseq(atlist(i))
                    imol = imolatseq(i0atseq(atlist(i)) + j)
                    f(:,imol) = dimlist(:,i)*fsum(:,i)/natseq(atlist(i)) + (1-dimlist(:,i))*f(:,imol)
                end do
            end do
        end if
        !
        if( vflag )then
            do i=1,nlist
                vcom(:,i) = vcom(:,i) + (fsumbuf(:,i)+fsum(:,i))*dt/mass(i)*0.5d0
                !
                do j=1,natseq(atlist(i))
                    imol = imolatseq(i0atseq(atlist(i)) + j)
                    v(:,imol) = dimlist(:,i)*vcom(:,i) + (1-dimlist(:,i))*v(:,imol)
                end do
            end do
        end if
    end subroutine fix_rigid
end module fix_rigid_mod

