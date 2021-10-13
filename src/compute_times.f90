! compute temperature
module compute_times_mod
    integer,allocatable :: myfc(:)
    integer :: nlist = 0
    integer,allocatable :: ndata(:), i0timesdata(:)
    integer,allocatable :: timesfc(:), ntimesfc(:), itimes0fc(:)
    double precision,allocatable :: timesdata(:)
    !
contains
    subroutine init_compute_times
        use update_mod
        use commn
        use file_mod
        use constant_mod
        use fix_compute_mod
        use mpivar
        implicit none
        integer i, ios, iatype, ntimesfcnow, ndata_now
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
        end do
        read(ifilecalc,*) (timesfc(itimes0fc(nlist+1)+i),i=1,ntimesfcnow)
        write(form1,'(I0)') ntimesfcnow
        form1 = "(A,I,A,"//trim(form1)//"I)"
        write(*,form1) __FILE__, ntimesfcnow, (timesfc(itimes0fc(nlist+1)+i),i=1,ntimesfcnow)
        !
        ndata_now = 0
        do i=1,ntimesfcnow
            if(i==1)then
                ndata_now = nfcdata(fcnumi(timesfc(itimes0fc(nlist+1)+i)))
                cycle
            end if
            if(ndata_now /= nfcdata(fcnumi(timesfc(itimes0fc(nlist+1)+i)))) call error_msg("ndata_now mismatch @"//__FILE__)
            write(*,*) nfcdata(fcnumi(timesfc(itimes0fc(nlist+1)+i))), fcnumi(timesfc(itimes0fc(nlist+1)+i)), timesfc(itimes0fc(nlist+1)+i)
        end do
        ndata = [ndata, ndata_now]
        if(nlist==0) i0timesdata = [0]
        i0timesdata = [i0timesdata, i0timesdata(nlist+1)+ndata_now]
        do i=1,ndata_now
            timesdata = [timesdata, 0d0]
        end do
        nlist = nlist + 1
        ! --- return ---
        fcnum = [fcnum, ifc]
        myfc = [myfc, nfc]
        nfcdata = [nfcdata, ndata_now]
        ! --- calc fcnumi ---
        if(size(fcnumi) < ifc)then
            do i=1,ifc-size(fcnumi)
                fcnumi = [fcnumi, -1]
            end do
        end if
        fcnumi(ifc) = nfc
        ! ===
    end subroutine init_compute_times
    !
    !
    ! --- compute times ---
    subroutine compute_times
        use commn
        use constant_mod
        use fix_compute_mod
        implicit none
        integer j,i,imol,i0,n,k,l
        !
        if( nlist == 0 ) return
        ! --- calculate times ---
        k = 0
        timesdata(:) = 1d0
        do i=1,nlist
            n = ndata(i)
            do l=1,ntimesfc(i)
                j = fcnumi(timesfc(itimes0fc(i)+l))
                i0 = i0fcdata(j)
                timesdata(k+1:k+n) = timesdata(k+1:k+n) * fcdata(i0+1:i0+n)
            end do
            ! --- return ---
            fcdata(i0fcdata(myfc(i))+1:i0fcdata(myfc(i))+n) = timesdata(k+1:k+n)
            !
            k = k + n
        end do
    end subroutine compute_times
end module compute_times_mod

