!compute
module compute_mod
contains
    subroutine init_compute
        use mpivar
        use file_mod
        use fix_compute_mod
        use compute_temp_mod, only : init_compute_temp
        use compute_temp2_mod, only : init_compute_temp2
        use compute_chunk_1d_mod, only : init_compute_chunk_1d
        use compute_compv_mod, only : init_compute_compv
        use compute_gg_mod, only : init_compute_gg
        use compute_gg_hf_mod, only : init_compute_gg_hf
        use compute_gg_VAstress_mod, only : init_compute_gg_VAstress
        use compute_sum_mod, only : init_compute_sum
        use compute_times_mod, only : init_compute_times
        use compute_times_one_mod, only : init_compute_times_one
        use compute_correlate_mod, only : init_compute_correlate
        use compute_ave_mod, only : init_compute_ave
        use compute_ave_period_mod, only : init_compute_ave_period
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
            if(comchar=="compute_temp" .or. comchar=="ComputeTemp") call init_compute_temp
            if(comchar=="compute_temp2" .or. comchar=="ComputeTemp2") call init_compute_temp2
            if(comchar=="compute_chunk_1d" .or. comchar=="ComputeChunk1d") call init_compute_chunk_1d
            if(comchar=="compute_compv" .or. comchar=="ComputeCompv") call init_compute_compv
            if(comchar=="compute_group_group" .or. comchar=="compute_gg" .or. comchar=="ComputeGg") call init_compute_gg
            if(comchar=="compute_group_group_hf" .or. comchar=="compute_gg_hf" .or. comchar=="ComputeGgHf") call init_compute_gg_hf
            if(comchar=="compute_gg_VAstress" .or. comchar=="ComputeGg") call init_compute_gg_VAstress
            !
            if(comchar=="compute_sum" .or. comchar=="ComputeSum") call init_compute_sum
            if(comchar=="compute_times" .or. comchar=="Computetimes") call init_compute_times
            if(comchar=="compute_times_one" .or. comchar=="Computetimes") call init_compute_times_one
            if(comchar=="compute_correlate" .or. comchar=="ComputeCorrelate") call init_compute_correlate
            if(comchar=="compute_ave" .or. comchar=="ComputeAve") call init_compute_ave
            if(comchar=="compute_ave_period" .or. comchar=="ComputeAvePeriod") call init_compute_ave_period
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
                write(*,*) i, fcnumi(fcnum(i)), fcnum(i), nfcdata(i), i0fcdata(i)
            end do
        end if
    end subroutine init_compute

    subroutine compute
        use mpivar
        use commn
        use compute_temp_mod, only : compute_temp
        use compute_temp2_mod, only : compute_temp2
        use compute_chunk_1d_mod, only : compute_chunk_1d
        use compute_compv_mod, only : compute_compv
        use compute_gg_mod, only : compute_gg
        use compute_gg_hf_mod, only : compute_gg_hf
        use compute_gg_VAstress_mod, only : compute_gg_VAstress
        use compute_sum_mod, only : compute_sum
        use compute_times_mod, only : compute_times
        use compute_times_one_mod, only : compute_times_one
        use compute_correlate_mod, only : compute_correlate
        use compute_ave_mod, only : compute_ave
        use compute_ave_period_mod, only : compute_ave_period
        implicit none
        !
        if(Iam==master)then
            call compute_temp
            call compute_temp2
            call compute_chunk_1d
            call compute_compv
        end if
        !
        call compute_gg
        call compute_gg_hf
        call compute_gg_VAstress
        !
        if(Iam==master)then
            call compute_sum
            call compute_times
            call compute_times_one
        end if
        call compute_correlate
        if(Iam==master)then
            call compute_ave
            call compute_ave_period
        end if
    end subroutine compute
end module compute_mod
