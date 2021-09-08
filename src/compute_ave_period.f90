! compute temperature
module compute_ave_period_mod
    integer,allocatable :: myfc(:)
    integer :: nlist = 0
    integer,allocatable :: ndata(:), i0avedata(:)
    character(128),allocatable :: avefile(:)
    integer,allocatable :: avefc(:), navefc(:), iave0fc(:), period(:)
    double precision,allocatable :: avedata(:)
    !
contains
    subroutine init_compute_ave_period
        use update_mod
        use commn
        use file_mod
        use constant_mod
        use fix_compute_mod
        use mpivar
        implicit none
        integer i, j, ios, iatype, navefcnow, ndata_now, period_now
        double precision,allocatable :: a(:)
        character(64) form1
        character(128) char1
        !
        if(Iam/=master)then
            read(ifilecalc,*)
            read(ifilecalc,*)
            read(ifilecalc,*)
            return
        end if
        !
        nfc = nfc + 1
        ! === read conf file ===
        read(ifilecalc,*) navefcnow, period_now
        navefc = [navefc, navefcnow]
        period = [period, period_now]
        if(nlist==0) iave0fc = [0]
        iave0fc = [iave0fc, iave0fc(nlist+1)+navefcnow]
        read(ifilecalc,'(A)') char1
        char1 = trim(adjustl(char1))
        open(10,file=trim(char1))
        close(10)
        avefile = [avefile, char1]
        do i=1,navefcnow
            avefc = [avefc, -1]
        end do
        read(ifilecalc,*) (avefc(iave0fc(nlist+1)+i),i=1,navefcnow)
        write(form1,'(I0)') navefcnow
        form1 = "(A,I,A,"//trim(form1)//"I)"
        write(*,form1) __FILE__, navefcnow, trim(avefile(nlist+1)), (avefc(iave0fc(nlist+1)+i),i=1,navefcnow)
        !
        ndata_now = 0
        do i=1,navefcnow
            ndata_now = ndata_now + nfcdata(fcnumi(avefc(iave0fc(nlist+1)+i)))
            write(*,*) nfcdata(fcnumi(avefc(iave0fc(nlist+1)+i))), fcnumi(avefc(iave0fc(nlist+1)+i)), avefc(iave0fc(nlist+1)+i)
        end do
        ndata = [ndata, ndata_now]
        if(nlist==0) i0avedata = [0]
        i0avedata = [i0avedata, i0avedata(nlist+1)+ndata_now*period_now]
        do j=1,period_now
            do i=1,ndata_now
                avedata = [avedata, 0d0]
            end do
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
    end subroutine init_compute_ave_period
    !
    !
    ! --- compute average ---
    subroutine compute_ave_period
        use commn
        use constant_mod
        use fix_compute_mod
        implicit none
        integer j,i,imol,i0,n,k,l,iperiod
        !
        if( nlist == 0 ) return
        ! --- calculate average ---
        k = 0
        do i=1,nlist
            do l=1,navefc(i)
                j = fcnumi(avefc(iave0fc(i)+l))
                i0 = i0fcdata(j)
                n = nfcdata(j)
                iperiod = mod(irep,period(i))
                avedata(k+iperiod*n+1:k+iperiod*n+n) = avedata(k+iperiod*n+1:k+iperiod*n+n) &
                    + fcdata(i0+1:i0+n)
                ! --- return ---
                fcdata(i0fcdata(myfc(i))+1:i0fcdata(myfc(i))+n) = avedata(k+iperiod*n+1:k+iperiod*n+n)
                !
                k = k + n*period(i)
            end do
        end do
        ! --- output ---
        if(irep==nrep)then
            k = 0
            do i=1,nlist
                open(10,file=avefile(i))
                write(10,'(A15)',advance="no") "period"
                do l=1,navefc(i)
                    n = nfcdata(fcnumi(avefc(iave0fc(i)+l)))
                    do j=1,n
                        write(10,'(I10,I5)',advance="no") avefc(iave0fc(i)+l), j
                    end do
                end do
                write(10,*)
                !
                do iperiod=1,period(i)
                    write(10,'(I15)',advance="no") iperiod
                    do l=1,navefc(i)
                        n = nfcdata(fcnumi(avefc(iave0fc(i)+l)))
                        do j=1,n
                            write(10,'(E15.5e3)',advance="no") avedata(k+j)/((nrep-iperiod)/period(i)+1)
                        end do
                        k = k + n
                    end do
                    !
                    write(10,*)
                end do
                close(10)
            end do
        end if
    end subroutine compute_ave_period
end module compute_ave_period_mod

