! compute temperature
module compute_compv_mod
    integer,allocatable :: myfc(:)
    integer,parameter :: ndata = 6
    integer :: nlist = 0
    integer,allocatable :: nsamp(:), atlist(:)
    double precision,allocatable :: comp(:,:), comv(:,:)
    !
contains
    subroutine init_compute_compv
        use update_mod
        use commn
        use file_mod
        use constant_mod
        use fix_compute_mod
        implicit none
        integer i, ios, iatype
        double precision,allocatable :: a(:)
        !
        nfc = nfc + 1
        ! === read conf file ===
        read(ifilecalc,*) iatype, i
        write(*,'(A,2I)') __FILE__, iatype, i
        !
        nsamp = [nsamp, i]
        atlist = [atlist, iatype]
        !
        comp = reshape( [0d0], [3,nlist+1], pad=[0d0] )
        comv = reshape( [0d0], [3,nlist+1], pad=[0d0] )
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
    end subroutine init_compute_compv
    !
    !
    ! --- compute temperature ---
    subroutine compute_compv
        use update_mod
        use commn
        use constant_mod
        use fix_compute_mod
        implicit none
        integer j,i,imol
        double precision random, work_now(1:nlist)
        !
        ! --- calculate temperature ---
        if(nlist==0) return
        comp(:,:) = 0d0
        comv(:,:) = 0d0
        do i=1,nlist
            if(mod(irep,nsamp(i)) /= 0) cycle
            do j=1,natseq(atlist(i))
                imol = imolatseq(i0atseq(atlist(i)) + j)
                ! calculate temperature
                comp(:,i) = comp(:,i) + p(:,imol)
                comv(:,i) = comv(:,i) + v(:,imol)
            end do
            comp(:,i) = comp(:,i)/natseq(atlist(i))
            comv(:,i) = comv(:,i)/natseq(atlist(i))
        end do
        ! --- return ---
        do i=1,nlist
            if(mod(irep,nsamp(i)) /= 0) cycle
            j = i0fcdata(myfc(i))
            fcdata(j+1:j+3) = comp(:,i)
            fcdata(j+4:j+ndata) = comv(:,i)
        end do
    end subroutine compute_compv
end module compute_compv_mod

