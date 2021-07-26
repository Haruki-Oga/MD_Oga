!compute
module compute_mpi_mod
contains
    subroutine init_compute_mpi
        use mpivar
        use file_mod
        use fix_compute_mod
        use compute_group_group_mod, only : init_compute_group_group
        use compute_group_group_hf_mod, only : init_compute_group_group_hf
        implicit none
        character(1) :: sharp = "#"
        character(64) comchar
        integer ios,i,j
        !
        open(ifilecalc,file=trim(filecalc),action="read",status="old")
        do
            read(ifilecalc,*,iostat=ios) comchar
            if( trim(comchar)=="compute" ) exit
            if(ios/=0) call error_msg("read error @ init_compute_mpi")
        end do
        do
            read(ifilecalc,*,iostat=ios) comchar, ifc
            !
            if(ios<0 .or. trim(comchar)=="output") exit
            if(comchar(1:1)==sharp) cycle
            if(ios>0) cycle !call error_msg("read error @init_compute_mpi")
            !
            if(comchar=="compute_group_group") call init_compute_group_group
            if(comchar=="compute_group_group_hf") call init_compute_group_group_hf
            !
        end do
        close(ifilecalc)
    end subroutine init_compute_mpi

    subroutine compute_mpi
        use mpivar
        use fix_compute_mod
        use compute_group_group_mod, only : compute_group_group
        use compute_group_group_hf_mod, only : compute_group_group_hf
        implicit none
        !
        call compute_group_group
        call compute_group_group_hf
        !
    end subroutine compute_mpi

end module compute_mpi_mod
