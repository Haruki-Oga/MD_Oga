! compute temperature
module compute_sum_mod
    integer,allocatable :: myfc(:)
    integer :: nlist = 0
    integer,allocatable :: ndata(:), i0sumdata(:)
    integer,allocatable :: sumfc(:), nsumfc(:), isum0fc(:)
    double precision,allocatable :: sumdata(:)
    !
contains
    subroutine init_compute_sum
        use update_mod
        use commn
        use file_mod
        use constant_mod
        use fix_compute_mod
        implicit none
        integer i, ios, iatype, nsumfcnow, ndata_now
        double precision,allocatable :: a(:)
        character(64) form1
        character(128) char1
        !
        nfc = nfc + 1
        ! === read conf file ===
        read(ifilecalc,*) nsumfcnow
        nsumfc = [nsumfc, nsumfcnow]
        if(nlist==0) isum0fc = [0]
        isum0fc = [isum0fc, isum0fc(nlist+1)+nsumfcnow]
        do i=1,nsumfcnow
            sumfc = [sumfc, -1]
        end do
        read(ifilecalc,*) (sumfc(isum0fc(nlist+1)+i),i=1,nsumfcnow)
        write(form1,'(I0)') nsumfcnow
        form1 = "(A,I,A,"//trim(form1)//"I)"
        write(*,form1) __FILE__, nsumfcnow, (sumfc(isum0fc(nlist+1)+i),i=1,nsumfcnow)
        !
        ndata_now = 0
        do i=1,nsumfcnow
            if(i==1)then
                ndata_now = nfcdata(fcnumi(sumfc(isum0fc(nlist+1)+i)))
                cycle
            end if
            if(ndata_now /= nfcdata(fcnumi(sumfc(isum0fc(nlist+1)+i)))) call error_msg("ndata_now mismatch @"//__FILE__)
            write(*,*) nfcdata(fcnumi(sumfc(isum0fc(nlist+1)+i))), fcnumi(sumfc(isum0fc(nlist+1)+i)), sumfc(isum0fc(nlist+1)+i)
        end do
        ndata = [ndata, ndata_now]
        if(nlist==0) i0sumdata = [0]
        i0sumdata = [i0sumdata, i0sumdata(nlist+1)+ndata_now]
        do i=1,ndata_now
            sumdata = [sumdata, 0d0]
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
    end subroutine init_compute_sum
    !
    !
    ! --- compute sum ---
    subroutine compute_sum
        use commn
        use constant_mod
        use fix_compute_mod
        implicit none
        integer j,i,imol,i0,n,k,l
        !
        if( nlist == 0 ) return
        ! --- calculate sum ---
        k = 0
        sumdata(:) = 0d0
        do i=1,nlist
            n = ndata(i)
            do l=1,nsumfc(i)
                j = fcnumi(sumfc(isum0fc(i)+l))
                i0 = i0fcdata(j)
                sumdata(k+1:k+n) = sumdata(k+1:k+n) + fcdata(i0+1:i0+n)
            end do
            ! --- return ---
            fcdata(i0fcdata(myfc(i))+1:i0fcdata(myfc(i))+n) = sumdata(k+1:k+n)
            !
            k = k + n
        end do
    end subroutine compute_sum
end module compute_sum_mod

