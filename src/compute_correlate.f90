! compute temperature
module compute_correlate_mod
    integer myfc
    integer :: nlist = 0
    integer,allocatable :: nevery(:), nrepeat(:), i0cordata(:)
    character(128),allocatable :: corfile(:)
    integer,allocatable :: corfci(:), corfcj(:), corfcii(:), corfcji(:), ncorfc(:), i0corfc(:)
    double precision,allocatable :: cordata(:), cordata_Iam(:), timedatai(:), timedataj(:)
    integer,allocatable :: icortime0(:)
    !
contains
    subroutine init_compute_correlate
        use update_mod
        use commn
        use file_mod
        use constant_mod
        use fix_compute_mod
        use mpivar
        implicit none
        integer i, ios, iatype, ncorfcnow, ndata_now, neverynow, nrepeatnow
        double precision,allocatable :: a(:)
        character(64) form1
        character(128) char1
        !
        if(Iam==master) nfc = nfc + 1
        ! === read conf file ===
        read(ifilecalc,*) ncorfcnow, neverynow, nrepeatnow
        ncorfc = [ncorfc, ncorfcnow]
        nevery = [nevery, neverynow]
        nrepeat = [nrepeat, nrepeatnow]
        if(nlist==0) i0corfc = [0]
        i0corfc = [i0corfc, i0corfc(nlist+1)+ncorfcnow]
        read(ifilecalc,'(A)') char1
        char1 = trim(adjustl(char1))
        open(10,file=trim(char1))
        close(10)
        corfile = [corfile, char1]
        do i=1,ncorfcnow
            corfci = [corfci, -1]
            corfcj = [corfcj, -1]
            corfcii = [corfcii, -1]
            corfcji = [corfcji, -1]
            icortime0 = [icortime0, 0]
        end do
        read(ifilecalc,*) (corfci(i),corfcii(i),corfcj(i),corfcji(i),i=i0corfc(nlist+1)+1,i0corfc(nlist+1)+ncorfcnow)
        if(Iam==master)then
            write(form1,'(I0)') ncorfcnow*4
            form1 = "(A,I,A,"//trim(form1)//"I)"
            write(*,form1) __FILE__, ncorfcnow, trim(corfile(nlist+1)), &
                (corfci(i),corfcii(i),corfcj(i),corfcji(i),i=i0corfc(nlist+1)+1,i0corfc(nlist+1)+ncorfcnow)
        end if
        !
        ndata_now = ncorfcnow*(nrepeatnow+1)
        if(nlist==0) i0cordata = [0]
        i0cordata = [i0cordata, i0cordata(nlist+1)+ndata_now]
        do i=1,ndata_now
            cordata = [cordata, 0d0]
            cordata_Iam = [cordata_Iam, 0d0]
            timedatai = [timedatai, 0d0]
            timedataj = [timedataj, 0d0]
        end do
        nlist = nlist + 1
        ! --- return ---
        if(Iam==master)then
            fcnum = [fcnum, ifc]
            myfc = nfc
            nfcdata = [nfcdata, ndata_now]
            ! --- calc fcnumi ---
            if(size(fcnumi) < ifc)then
                do i=1,ifc-size(fcnumi)
                    fcnumi = [fcnumi, -1]
                end do
            end if
            fcnumi(ifc) = nfc
        end if
        ! ===
    end subroutine init_compute_correlate
    !
    !
    ! --- compute average ---
    subroutine compute_correlate
        use commn
        use constant_mod
        use fix_compute_mod
        use mpivar
        implicit none
        integer j,i,imol,i0,n,k,l,m,mt,ifcnumi, jfcnumi, iifcnumi, jifcnumi, idata,it0,it
        !
        if(nlist==0) return
        ! --- calculate correlation ---
        do i=1,nlist
            if(mod(irep,nevery(i))/=0) cycle
            ! --- calc icortime0 ---
            icortime0(i) = icortime0(i) + 1
            if(icortime0(i) > nrepeat(i)) icortime0(i) = 0
            !
            do l=1,ncorfc(i)
                if(Iam==master)then
                    ifcnumi = fcnumi(corfci(i0corfc(i)+l))
                    jfcnumi = fcnumi(corfcj(i0corfc(i)+l))
                    iifcnumi = corfcii(i0corfc(i)+l)
                    jifcnumi = corfcji(i0corfc(i)+l)
                end if
                idata = i0cordata(i) + (l-1)*(nrepeat(i)+1)
                ! --- read fcdata ---
                if(Iam==master)then
                    timedatai(idata+icortime0(i)+1) = fcdata(i0fcdata(ifcnumi)+iifcnumi)
                    timedataj(idata+icortime0(i)+1) = fcdata(i0fcdata(jfcnumi)+jifcnumi)
                end if
                call mpi_bcast(timedatai(idata+icortime0(i)+1), 1, mpi_real8, master, mpi_comm_world, ierr)
                call mpi_bcast(timedataj(idata+icortime0(i)+1), 1, mpi_real8, master, mpi_comm_world, ierr)
                ! --- calc correlation ---
                do m=Iam,nrepeat(i),Nproc
                    it0 = idata + icortime0(i) + 1
                    mt = icortime0(i) - m
                    if(mt < 0) mt = mt + nrepeat(i)+1
                    it = idata + mt + 1
                    cordata_Iam(idata+m+1) = cordata_Iam(idata+m+1) + timedatai(it0)*timedataj(it)
                end do
                !
            end do
        end do
        ! --- output ---
        if(irep/=nrep) return
        call mpi_reduce(cordata_Iam, cordata, size(cordata), mpi_real8,&
            mpi_sum, master, mpi_comm_world, ierr)
        if(Iam/=master) return
        !k = 0
        do i=1,nlist
            open(10,file=corfile(i))
            do l=1,ncorfc(i)
                idata = i0cordata(i) + (l-1)*(nrepeat(i)+1)
                do j=0,nrepeat(i)
                    cordata(idata+j+1) = cordata(idata+j+1)/(nrep/nevery(i)-j+1)
                end do
            end do
            ! --- write ---
            write(10,'(A15)',advance="no") "# time"
            do l=1,ncorfc(i)
                j = i0corfc(i)+l
                write(10,'(I5,I3,I4,I3)',advance="no") corfci(j),corfcii(j),corfcj(j),corfcji(j)
            end do
            write(10,*)
            !
            do j=0,nrepeat(i)
                write(10,'(E15.5e3)',advance="no") dt*nevery(i)*j
                do l=1,ncorfc(i)
                    idata = i0cordata(i) + (l-1)*(nrepeat(i)+1) + j + 1
                    write(10,'(E15.5e3)',advance="no") cordata(idata)
                end do
                write(10,*)
            end do
            close(10)
        end do
    end subroutine compute_correlate
end module compute_correlate_mod

