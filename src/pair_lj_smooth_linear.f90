module pair_lj_smooth_linear_mod
    implicit none
    logical(1),allocatable ::  pairflag(:,:)
    double precision,allocatable :: sig(:,:), eps(:,:)
contains
    subroutine init
        use update_mod
        use file_mod
        use mpivar
        implicit none
        double precision r8a,r8b
        character(64) pairstyle
        integer ios,i,j
        allocate(sig(1:natype,1:natype),eps(1:natype,1:natype))
        allocate(pairflag(1:natype,1:natype))
        pairflag(1:natype,1:natype) = .false.
        eps(1:natype,1:natype) = -1d0
        sig(1:natype,1:natype) = -1d0
        open(ifilecalc,file=trim(filecalc),action="read",status="old")
        do
            read(ifilecalc,*) pairstyle
            if(trim(pairstyle)=="fix" .or. trim(pairstyle)=="compute")exit
            if(trim(pairstyle)=="smooth_linear" .or. trim(pairstyle)=="lj_smooth_linear"&
                .or. trim(pairstyle)=="pair_lj_smooth_linear")then
                do
                    read(ifilecalc,*,iostat=ios) i,j ,r8a,r8b
                    if(ios/=0)exit
                    if(Iam==master)then
                        write(*,'(A,2I,2E)') trim(pairstyle), i,j,r8a,r8b
                    end if
                    pairflag(i,j) = .true.
                    pairflag(j,i) = .true.
                    eps(i,j) = r8a
                    sig(i,j) = r8b
                    eps(j,i) = eps(i,j)
                    sig(j,i) = sig(i,j)
                end do
            end if
        end do
        close(ifilecalc)
    end subroutine init

    subroutine calc_force_one(atypei,atypej,rabs,rc,force_one,e)
        use update_mod
        implicit none
        integer(1) atypei,atypej
        double precision rabs,rc,force_one,e
        double precision rcsq,rabssq, sigsq, sbr6, sbrc6, dlja3, rrci
        if(.not. pairflag(atypei,atypej)) return
        rcsq = rc*rc
        rabssq = rabs*rabs
        sigsq = sig(atypei,atypej)**2.0
        sbr6 = (sigsq/rabssq)**3d0
        sbrc6 = (sigsq/rcsq)**3d0
        dlja3 = (2d0*sbrc6-1d0)*sbrc6
        rrci = rabs/rc
        force_one = 24d0*eps(atypei,atypej)*((2d0*sbr6-1d0)*sbr6 - dlja3*rrci)/rabssq
        e = e + 4d0*eps(atypei,atypej)*( (sbr6-1d0)*sbr6 - (sbrc6-1d0)*sbrc6 &
            - 6d0*(rrci-1d0)*dlja3 )
        return
    end subroutine calc_force_one
end module pair_lj_smooth_linear_mod
