module force_mod
    double precision,allocatable :: f_Iam(:,:)
    double precision,allocatable :: f_pair(:,:), f_bonds(:,:)
    double precision ep_Iam, ep_pair, ep_bonds, dpssqmax
    double precision rc,rcsq
    logical(1),allocatable :: pairflag(:,:)
    logical(1) :: lj_smooth_quad = .false., smooth_linear = .false.
    integer,allocatable :: mode_pair(:,:)
    !
contains
    subroutine init_force
        use update_mod
        use file_mod
        use pair_lj_smooth_quad_mod&
            , only : init_pair_lj_smooth_quad => init, pair_lj_smooth_quad_pairflag => pairflag
        use pair_lj_smooth_linear_mod&
            , only : init_pair_lj_smooth_linear => init, pair_lj_smooth_linear_pairflag => pairflag
        implicit none
        double precision r8a, r8b
        integer i,j
        allocate(pairflag(1:natype,1:natype))
        pairflag(:,:) = .false.
        rc = 3.5d0
        rcsq = rc**2d0
        call init_pair_lj_smooth_quad
        call init_pair_lj_smooth_linear
        !
        do i=1,natype
            do j=1,natype
                if(pair_lj_smooth_quad_pairflag(i,j)) pairflag(i,j) = .true.
                if(pair_lj_smooth_linear_pairflag(i,j)) pairflag(i,j) = .true.
            end do
        end do
        !
    end subroutine init_force
    function fpair_one(atypei,atypej,rabs,e)
        use pair_lj_smooth_quad_mod, only : f_pair_lj_smooth_quad => calc_force_one 
        use pair_lj_smooth_linear_mod, only : f_pair_lj_smooth_linear => calc_force_one 
        implicit none
        integer(1) atypei,atypej
        double precision rabs,e
        double precision fpair_one
        !
        call f_pair_lj_smooth_quad(atypei,atypej,rabs,rc,fpair_one,e)
        call f_pair_lj_smooth_linear(atypei,atypej,rabs,rc,fpair_one,e)
    end function fpair_one
    !
    !!
    !function fij_lj_smooth_linear(eps,sig,rcsq,rabssq,e)
    !    implicit none
    !    double precision eps,sig,rcsq, rabssq, e
    !    double precision fij_lj_smooth_linear
    !    double precision sigsq, sbr6, sbrc6, dlja3, rrci
    !    sigsq = sig*sig
    !    sbr6 = (sigsq/rabssq)**3d0
    !    sbrc6 = (sigsq/rcsq)**3d0
    !    dlja3 = (2d0*sbrc6-1d0)*sbrc6
    !    rrci = sqrt(rabssq/rcsq)
    !    fij_lj_smooth_linear = 24d0*eps*((2d0*sbr6-1d0)*sbr6 - dlja3*rrci)/rabssq
    !    e = e + 4d0*eps*( (sbr6-1d0)*sbr6 - (sbrc6-1d0)*sbrc6 &
    !        - 6d0*(rrci-1d0)*dlja3 )
    !    return
    !end function fij_lj_smooth_linear
end module force_mod

