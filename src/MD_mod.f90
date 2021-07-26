module commn
    integer nrep,irep,nsamp_term
    double precision dt
    integer nmol
    double precision vl(3), lo(3), hi(3), vlh(3), va(3)
    double precision, allocatable :: p(:,:), v(:,:), f(:,:)
    double precision ep, temp, ek, eadd
    logical(1) pflag, fflag, vflag 
end module commn

module mpivar
    use mpi
    integer &
        Iam, Nproc, dest, tag, master, ierr
    data master /0/
end module mpivar

module update_mod
    integer natype, nbonds, nbtype
    integer(1),allocatable :: atype(:)
    integer,allocatable :: btype(:), bondi(:), bondj(:)
    double precision,allocatable :: wm(:), k(:), r_eq(:)
    double precision,allocatable :: sig(:,:), eps(:,:)
    logical(1),allocatable :: ffixflag(:,:), vfixflag(:,:), pfixflag(:,:)
    double precision,allocatable :: dt_wmi2(:),dps(:,:),dp(:,:)
    double precision,allocatable :: fbuf(:,:)
    integer,allocatable :: natseq(:), i0atseq(:), imolatseq(:)
end module update_mod

module force_mod
    double precision,allocatable :: f_Iam(:,:)
    double precision,allocatable :: f_pair(:,:), f_bonds(:,:)
    double precision ep_Iam, ep_pair, ep_bonds, dpssqmax
    double precision rc,rcsq
    logical(1),allocatable :: pairflag(:,:)
end module force_mod

module book_keep_mod
    integer maxnlist, nlistall
    double precision rbk, rbksq
    integer,allocatable :: listi(:), listj(:)
    integer nlist
    integer,allocatable :: listi_gather(:), listj_gather(:)
    double precision rmar, rmarsq4
end module book_keep_mod

module file_mod
    integer fileposb
    integer,parameter :: ifilecalc = 11
    integer,parameter :: ifileinit = 12
    integer,parameter :: ifilelast = 13
    character(128) filecalc, fileinit, filelast, etcdir, resultdir
end module file_mod

module time_mod
    real(4) timeall, timeall0
    integer timealls,timeallm,timeallh
    real(4) timebk, timeforce
    !
contains
    function start_time
        use mpivar
        implicit none
        real(4) start_time
        call mpi_barrier(mpi_comm_world, ierr)
        start_time = mpi_wtime()
        return
    end function start_time
    function stop_time(time0)
        use mpivar
        implicit none
        real(4) time0, stop_time
        call mpi_barrier(mpi_comm_world, ierr)
        stop_time = mpi_wtime() - time0
        return
    end function stop_time
end module time_mod
