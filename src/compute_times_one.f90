! compute temperature
module compute_times_one_mod
    integer,allocatable :: myfc(:)
    integer :: nlist = 0
    integer :: ndata = 1
    integer,allocatable :: timesfc(:), timesfci(:), ntimesfc(:), itimes0fc(:)
    double precision,allocatable :: timesdata(:)
    !
contains
    subroutine init_compute_times_one
        use update_mod
        use commn
        use file_mod
        use constant_mod
        use fix_compute_mod
        use mpivar
        implicit none
        integer i, ios, iatype, ntimesfcnow
        double precision,allocatable :: a(:)
        character(64) form1
        character(128) char1
        !
        if(Iam/=master)then
            read(ifilecalc,*)
            read(ifilecalc,*)
            return
        end if
        nfc = nfc + 1
        ! === read conf file ===
        read(ifilecalc,*) ntimesfcnow
        ntimesfc = [ntimesfc, ntimesfcnow]
        if(nlist==0) itimes0fc = [0]
        itimes0fc = [itimes0fc, itimes0fc(nlist+1)+ntimesfcnow]
        do i=1,ntimesfcnow
            timesfc = [timesfc, -1]
            timesfci = [timesfci, -1]
        end do
        read(ifilecalc,*) (timesfc(itimes0fc(nlist+1)+i), timesfci(itimes0fc(nlist+1)+i),i=1,ntimesfcnow)
        write(form1,'(I0)') 2*ntimesfcnow
        form1 = "(A,I,"//trim(form1)//"I)"
        write(*,form1) __FILE__, ntimesfcnow, (timesfc(itimes0fc(nlist+1)+i),timesfci(itimes0fc(nlist+1)+i),i=1,ntimesfcnow)
        !
        timesdata = [timesdata, 0d0]
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
    end subroutine init_compute_times_one
    !
    !
    ! --- compute times ---
    subroutine compute_times_one
        use commn
        use constant_mod
        use fix_compute_mod
        implicit none
        integer j,i,imol,i0,n,k,l
        !
        if( nlist == 0 ) return
        ! --- calculate times ---
        timesdata(:) = 1d0
        do i=1,nlist
            do l=1,ntimesfc(i)
                j = fcnumi(timesfc(itimes0fc(i)+l))
                k = timesfci(itimes0fc(i)+l)
                i0 = i0fcdata(j)
                timesdata(i) = timesdata(i) * fcdata(i0+k)
            end do
            ! --- return ---
            fcdata(i0fcdata(myfc(i))+1:i0fcdata(myfc(i))+ndata) = timesdata(i)
            !
        end do
    end subroutine compute_times_one
end module compute_times_one_mod

