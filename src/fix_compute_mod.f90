module fix_compute_mod
    ! fcnum: imput=ifc, return=fcnum
    ! fcnumi: input=fcnum, return=ifc
    ! i0fcdata: input=ifc, return=first position of fcdata
    ! nfcdata: input=ifc, retrun=# of fcdata of fcnum
    integer ifc
    integer :: nfc = 0
    integer,allocatable :: fcnum(:), i0fcdata(:), nfcdata(:)
    integer,allocatable :: fcnumi(:)
    double precision,allocatable :: fcdata(:)
end module fix_compute_mod
