! compute temperature
module compute_chunk_1d_mod
    integer,allocatable :: myfc(:)
    integer ndata
    integer :: nlist = 0
    integer,allocatable :: nsamp(:), atlist(:)
    double precision,allocatable :: dx(:)
    double precision,allocatable :: ctemp(:), cvel(:,:)
    integer,allocatable :: cdens(:)
    integer,allocatable :: nx(:), i0x(:)
    integer,allocatable :: idir(:)
    !
contains
    subroutine init_compute_chunk_1d
        use update_mod
        use commn
        use file_mod
        use constant_mod
        use fix_compute_mod
        implicit none
        integer i, ios, iatype, n, j
        double precision,allocatable :: a(:)
        double precision d
        !
        nfc = nfc + 1
        ! === read conf file ===
        read(ifilecalc,*) iatype, i, j, d
        write(*,'(A,3I10, F15.5)') __FILE__, iatype, i, j, d
        !
        nsamp = [nsamp, i]
        atlist = [atlist, iatype]
        idir = [idir, j]
        dx = [dx, d]
        !
        n = int(vl(j)/d) + 1
        nx = [nx,n]
        if(nlist==0) i0x = [0]
        i0x = [i0x,i0x(nlist+1)+n]
        !
        do i=1,n
            ctemp = [ctemp,0d0]
            cdens = [cdens,0]
        end do
        cvel = reshape( [cvel, 0d0,0d0,0d0], [3,i0x(nlist+1)+n] ,pad=[0d0])
        !
        nlist = nlist + 1
        ! --- return ---
        fcnum = [fcnum, ifc]
        myfc = [myfc, nfc]
        ndata = n*5
        nfcdata = [nfcdata, ndata]
        ! --- calc fcnumi ---
        if(size(fcnumi) < ifc)then
            do i=1,ifc-size(fcnumi)
                fcnumi = [fcnumi, -1]
            end do
        end if
        fcnumi(ifc) = nfc
        ! ===
    end subroutine init_compute_chunk_1d
    !
    !
    ! --- compute temperature ---
    subroutine compute_chunk_1d
        use update_mod
        use commn
        use constant_mod
        use fix_compute_mod
        use file_mod
        implicit none
        integer j,i,imol, ix
        double precision random, work_now(1:nlist)
        character(128) densfile
        !
        ! --- calculate temperature ---
        do i=1,nlist
            if(mod(irep,nsamp(i)) /= 0) cycle
            do j=1,natseq(atlist(i))
                imol = imolatseq(i0atseq(atlist(i)) + j)
                ! calculate temperature
                ix = int((p(idir(i),imol) - lo(idir(i)))/dx(i)) + 1
                ix = ix + i0x(i)
                !
                ctemp(ix) = ctemp(ix) + sum(v(:,imol)**2d0)
                cdens(ix) = cdens(ix) + 1
                cvel(:,ix) = cvel(:,ix) + v(:,imol)
            end do
        end do
        ! --- return ---
        if(irep/=nrep) return
        do i=1,nlist
            do j=1,nx(i)
                ix = i0x(i) + j
                if( cdens(ix) == 0 ) cycle
                !
                cvel(:,ix) = cvel(:,ix)/cdens(ix)
                ctemp(ix) = ctemp(ix)/cdens(ix)/3d0/bk*wm(atlist(i))
            end do
            !
            write(densfile,'(I0,"_"I0)') atlist(i), idir(i)
            densfile = trim(resultdir)//"chunk"//trim(densfile)//".dat"
            open(10,file=densfile)
            write(10,'(6A15)') "# position", "density", "temperature", "vx", "vy", "vz"
            do j=1,nx(i)
                ix = i0x(i) + j
                write(10,'(6E15.5e3)') (j-0.5d0)*dx(i), &
                    dble(cdens(ix))/nrep*wm(atlist(i))/dx(i)/va(idir(i)), ctemp(ix), cvel(:,ix)
            end do
            close(10)
        end do
        !
    end subroutine compute_chunk_1d
end module compute_chunk_1d_mod

