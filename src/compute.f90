!compute
module compute_mod
contains
    subroutine init_compute
        use mpivar
        use file_mod
        use fix_compute_mod
        use compute_temp_mod, only : init_compute_temp
        use compute_chunk_1d_mod, only : init_compute_chunk_1d
        use compute_compv_mod, only : init_compute_compv
        use compute_group_group_mod, only : init_compute_group_group
        use compute_group_group_hf_mod, only : init_compute_group_group_hf
        use compute_correlate_mod, only : init_compute_correlate
        use compute_ave_mod, only : init_compute_ave
        implicit none
        character(1) :: sharp = "#"
        character(64) comchar
        integer ios,i,j
        !
        open(ifilecalc,file=trim(filecalc),action="read",status="old")
        do
            read(ifilecalc,*,iostat=ios) comchar
            if( trim(comchar)=="compute" ) exit
            if(ios/=0) call error_msg("read error @ init_compute")
        end do
        do
            read(ifilecalc,*,iostat=ios) comchar, ifc
            !
            if(ios<0 .or. trim(comchar)=="output") exit
            if(comchar(1:1)==sharp) cycle
            if(ios>0) cycle ! call error_msg("read error @init_compute")
            ! --- init_compute ---
            if(Iam==master)then
                if(comchar=="compute_temp") call init_compute_temp
                if(comchar=="compute_chunk_1d") call init_compute_chunk_1d
                if(comchar=="compute_compv") call init_compute_compv
            end if
            if(comchar=="compute_group_group") call init_compute_group_group
            if(comchar=="compute_group_group_hf") call init_compute_group_group_hf
            !
            if(comchar=="compute_correlate") call init_compute_correlate
            if(Iam==master)then
                if(comchar=="compute_ave") call init_compute_ave
            end if
            ! ---
        end do
        close(ifilecalc)
        !
        if(Iam==master)then
            write(*,*) "all nfc=", nfc
            if(nfc == 0) return
            ! --- calc i0fcdata ---
            allocate(fcdata(1:sum(nfcdata(:))))
            fcdata(:) = 0d0
            allocate(i0fcdata(1:nfc))
            i0fcdata(1) = 0
            do i=1,nfc-1
                i0fcdata(i+1) = i0fcdata(i) + nfcdata(i)
            end do
            do i=1,nfc
                write(*,*) i,fcnum(i), nfcdata(i), i0fcdata(i)
            end do
        end if
    end subroutine init_compute

    subroutine compute
        use mpivar
        use compute_temp_mod, only : compute_temp
        use compute_chunk_1d_mod, only : compute_chunk_1d
        use compute_compv_mod, only : compute_compv
        use compute_group_group_mod, only : compute_group_group
        use compute_group_group_hf_mod, only : compute_group_group_hf
        use compute_correlate_mod, only : compute_correlate
        use compute_ave_mod, only : compute_ave
        implicit none
        !
        if(Iam==master)then
            call compute_temp
            call compute_chunk_1d
            call compute_compv
        end if
        !
        call compute_group_group
        call compute_group_group_hf
        !
        call compute_correlate
        if(Iam==master)then
            call compute_ave
        end if
    end subroutine compute

end module compute_mod
