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
        use pair_lj_smooth_quad_mod, only : init_pair_lj_smooth_quad => init
        implicit none
        double precision r8a, r8b
        integer i,j
        rc = 3.5d0
        rcsq = rc**2d0
        call init_pair_lj_smooth_quad
        !
    end subroutine init_force
    function fpair_one(atypei,atypej,rabs,e)
        use pair_lj_smooth_quad_mod, only : f_pair_lj_smooth_quad => force_one 
        implicit none
        integer(1) atypei,atypej
        double precision rabs,e
        double precision fpair_one
        !
        fpair_one = f_pair_lj_smooth_quad(atypei,atypej,rabs,e,rc)
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

