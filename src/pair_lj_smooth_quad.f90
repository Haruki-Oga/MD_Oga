module pair_lj_smooth_quad_mod
    implicit none
    logical(1),allocatable ::  pairflag(:,:)
    double precision,allocatable :: sig(:,:), eps(:,:)
contains
    subroutine init
        use update_mod
        use file_mod
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
            if(trim(pairstyle)=="fix")exit
            if(trim(pairstyle)=="smooth_quad" .or. trim(pairstyle)=="lj_smooth_quad")then
                do
                    read(ifilecalc,*,iostat=ios) i,j ,r8a,r8b
                    if(ios/=0)exit
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

    function force_one(atypei,atypej,rabs,e,rc)
        use update_mod
        implicit none
        double precision rc ,rabs,e,force_one
        double precision rcsq,rabssq, sigsq, sbr6, sbrc6, dlja3
        integer(1) atypei,atypej
        if(.not. pairflag(atypei,atypej)) return
        rcsq = rc*rc
        rabssq = rabs*rabs
        sigsq = sig(atypei,atypej)**2.0
        sbr6 = (sigsq/rabssq)**3d0
        sbrc6 = (sigsq/rcsq)**3d0
        dlja3 = (2d0*sbrc6-1d0)*sbrc6*rabssq/rcsq
        !dljb = (-7.0d0*sbrc6 + 4.0d0)*sbrc6
        force_one = 24d0*eps(atypei,atypej)*((2d0*sbr6-1d0)*sbr6 - dlja3)/rabssq
        e = e + 4d0*eps(atypei,atypej)*( (sbr6-1d0)*sbr6 &
            + dlja3*3d0 + (-7d0*sbrc6+4d0)*sbrc6 )
        return
    end function force_one
end module pair_lj_smooth_quad_mod
