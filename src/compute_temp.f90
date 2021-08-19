! compute temperature
module compute_temp_mod
    integer,allocatable :: myfc(:)
    integer,parameter :: ndata = 4
    integer :: nlist = 0
    integer,allocatable :: nsamp(:), atlist(:)
    double precision,allocatable :: ctemp(:,:)
    !
contains
    subroutine init_compute_temp
        use update_mod
        use commn
        use file_mod
        use constant_mod
        use fix_compute_mod
        use mpivar
        implicit none
        integer i, ios, iatype
        double precision,allocatable :: a(:)
        !
        if(Iam/=master)then
            read(ifilecalc,*)
            return
        end if
        nfc = nfc + 1
        ! === read conf file ===
        read(ifilecalc,*) iatype, i
        write(*,'(A,2I)') __FILE__, iatype, i
        !
        nsamp = [nsamp, i]
        atlist = [atlist, iatype]
        !
        ctemp = reshape( [0d0], [4,nlist+1], pad=[0d0] )
        nlist = nlist + 1
        ! --- return ---
        fcnum = [fcnum, ifc]
        myfc = [myfc, nfc]
        nfcdata = [nfcdata, ndata]
        ! --- calc fcnumi ---
        if(size(fcnumi) < ifc)then
            do i=1,ifc-size(fcnumi)
                fcnumi = [fcnumi, -1]
            end do
        end if
        fcnumi(ifc) = nfc
        ! ===
    end subroutine init_compute_temp
    !
    !
    ! --- compute temperature ---
    subroutine compute_temp
        use update_mod
        use commn
        use constant_mod
        use fix_compute_mod
        implicit none
        integer j,i,imol
        double precision random, work_now(1:nlist)
        !
        ! --- calculate temperature ---
        ctemp(:,:) = 0d0
        do i=1,nlist
            if(mod(irep,nsamp(i)) /= 0) cycle
            do j=1,natseq(atlist(i))
                imol = imolatseq(i0atseq(atlist(i)) + j)
                ! calculate temperature
                ctemp(1:3,i) = ctemp(1:3,i)&
                    + v(:,imol)**2d0*wm(atlist(i))/bk
            end do
            ctemp(1:3,i) = ctemp(1:3,i)/dble(natseq(atlist(i)))
            ctemp(4,i) = sum(ctemp(1:3,i))/3d0
        end do
        ! --- return ---
        do i=1,nlist
            if(mod(irep,nsamp(i)) /= 0) cycle
            fcdata(i0fcdata(myfc(i))+1:i0fcdata(myfc(i))+ndata) = ctemp(1:4,i)
        end do
    end subroutine compute_temp
end module compute_temp_mod

