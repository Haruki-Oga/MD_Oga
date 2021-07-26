!output
module output_mod
    integer,allocatable :: outnum(:), nevery(:), fcout(:), nfcout(:), i0fcout(:)
    integer :: nout = 0
contains
    subroutine output
        use commn
        use mpivar
        implicit none
        !
        if(Iam==master)then
            if(mod(irep,nsamp_term)==0)then
                call output_term
                call output_posb
            end if
            call output_fix_compute
        end if
    end subroutine output

    subroutine init_output
        use mpivar
        use file_mod
        use fix_compute_mod
        use commn
        implicit none
        character(1) :: sharp = "#"
        character(64) outfile
        integer ios,i,j,k
        integer :: iopen = 100
        !
        if(Iam==master)then
            open(ifilecalc,file=trim(filecalc),action="read",status="old")
            do
                read(ifilecalc,*,iostat=ios) outfile
                if( trim(outfile)=="output" ) exit
                if(ios>0) call error_msg("read error @ init_output")
                if(ios<0) return
            end do
            nout = 0
            do
                read(ifilecalc,'(A)',iostat=ios) outfile
                outfile = trim(adjustl(outfile))
                !
                if(ios<0) exit
                if(outfile(1:1)==sharp .or. trim(outfile)=="") cycle
                if(ios>0) call error_msg("read error @init_output")
                !
                read(ifilecalc,*) i,j
                open(iopen, file=trim(outfile), action="write",status="replace")
                outnum = [outnum,iopen]
                nevery = [nevery, i]
                nfcout = [nfcout, j]
                if(nout==0) i0fcout = [0]
                i0fcout = [i0fcout, i0fcout(nout+1)+j]
                do i=1,j
                    fcout = [fcout, 0]
                end do
                read(ifilecalc,*) (fcout(i0fcout(nout+1)+i),i=1,j)
                nout = nout + 1
                iopen = iopen + 1
            end do
            close(ifilecalc)
        end if
    end subroutine init_output

    subroutine output_fix_compute
        use commn
        use mpivar
        use fix_compute_mod
        implicit none
        integer i,j,k,l,m
        character(32) form1, form0
        !
        do i=1,nout
            if(mod(irep,nevery(i))/=0) cycle
            !
            if(irep==0)then
                do k=1,nfcout(i)
                    l = i0fcout(i)+k
                    j = fcnumi(fcout(l))
                    !
                    write(form1,'(I0)') nfcdata(j)
                    !
                    if(k==1) write(outnum(i),'(A10)',advance="no") "# step"
                    form0 = "("//trim(form1)//"I15)"
                    write(outnum(i),form0,advance="no") (fcout(l),m=1,nfcdata(j))
                end do
                write(outnum(i),*)
            end if
            !
            do k=1,nfcout(i)
                l = i0fcout(i)+k
                j = fcnumi(fcout(l))
                !
                write(form1,'(I0)') nfcdata(j)
                !
                if(k==1) write(outnum(i),'(I10)',advance="no") irep
                form1 = "("//trim(form1)//"E15.5e3)"
                write(outnum(i),form1,advance="no") fcdata(i0fcdata(j)+1:i0fcdata(j)+nfcdata(j))
            end do
            write(outnum(i),*)
        end do
    end subroutine output_fix_compute

    subroutine output_term
        use commn
        use mpivar
        use update_mod
        use file_mod
        use constant_mod
        implicit none
        real(4) r43(3)
        integer i
        integer(4) i4nmol, i4step
        !
        if(Iam==master)then
            !--- terminal ---
            do i=1,nmol
                ek = ek + sum(v(:,i)**2d0)*wm(atype(i))*0.5
            end do
            temp = ek*2d0/(3d0*bk*nmol)
            if(irep==0)then
                write(*,'(A10,4A15)') "# step", "temp", "ep", "etot", "etot-eadd"
            end if
            write(*,'(I10,4E15.5e3)') irep, temp, ep, ek+ep, ek+ep-eadd
            !
        end if
    end subroutine output_term

    subroutine output_posb
        use commn
        use mpivar
        use update_mod
        use file_mod
        use constant_mod
        !use output_mod
        implicit none
        real(4) r43(3)
        integer i
        integer(4) i4nmol, i4step
        !
        if(Iam==master)then
            !--- posb.dat ---
            if(irep==0)then
                i4nmol = nmol
                i4step = nrep/nsamp_term+1
                write(fileposb) i4nmol, i4step
                write(fileposb) 0d0, dt*dble(nsamp_term)
                write(fileposb) lo(:),hi(:)
            end if
            do i=1,nmol
                r43(:) = p(:,i)
                write(fileposb) atype(i), r43(:)
            end do
            !
        end if
    end subroutine output_posb




    subroutine last_output
        use mpivar
        implicit none
        integer i
        !
        if(Iam==master)then
            do i=1,nout
                close(outnum(i))
            end do
        end if
    end subroutine last_output

end module output_mod
