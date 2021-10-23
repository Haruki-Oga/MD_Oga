! fix
module fix_mod
contains
    subroutine init_fix
        use mpivar
        use file_mod
        use fix_compute_mod
        use fix_freeze_mod, only : init_fix_freeze
        use fix_vel_scaling_mod, only : init_fix_vel_scaling
        use fix_langevin_mod, only : init_fix_langevin
        use fix_langevin2_mod, only : init_fix_langevin2
        use fix_langevin3_mod, only : init_fix_langevin3
        use fix_langevin_LAMMPS_mod, only : init_fix_langevin_LAMMPS
        use fix_NoseHoover_mod, only : init_fix_NoseHoover
        use fix_move_linear_mod, only : init_fix_move_linear
        use fix_move_wiggle_mod, only : init_fix_move_wiggle
        use fix_set_pos_mod, only : init_fix_set_pos
        use fix_addforce_mod, only : init_fix_addforce
        use fix_rigid_mod, only : init_fix_rigid
        use fix_rigid_rotate_mod, only : init_fix_rigid_rotate
        implicit none
        character(1) :: sharp = "#"
        character(64) fixchar
        integer ios,i,j
        !
        if(Iam==master)then
            open(ifilecalc,file=trim(filecalc),action="read",status="old")
            do
                read(ifilecalc,*,iostat=ios) fixchar
                if( trim(fixchar)=="fix" ) exit
                if(ios/=0) call error_msg("read error @ init_fix1")
            end do
            do
                read(ifilecalc,*,iostat=ios) fixchar, ifc
                !
                if(ios<0 .or. trim(fixchar)=="compute") exit
                if(fixchar(1:1)==sharp) cycle
                if(ios>0) call error_msg("read error @init_fix")
                !
                nfc = nfc + 1
                fcnum = [fcnum, ifc]
                if(fixchar=="fix_vel_scaling" .or. fixchar=="Fixvel_scaling") call init_fix_vel_scaling
                if(fixchar=="fix_langevin" .or. fixchar=="FixLangevin") call init_fix_langevin
                if(fixchar=="fix_langevin2" .or. fixchar=="FixLangevin2") call init_fix_langevin2
                if(fixchar=="fix_langevin3" .or. fixchar=="FixLangevin3") call init_fix_langevin3
                if(fixchar=="fix_langevin_LAMMPS" .or. fixchar=="FixLangevinLAMMPS") call init_fix_langevin_LAMMPS
                if(fixchar=="fix_NoseHoover" .or. fixchar=="FixNoseHoover") call init_fix_NoseHoover
                if(fixchar=="fix_addforce" .or. fixchar=="FixAddforce") call init_fix_addforce
                if(fixchar=="fix_move_linear" .or. fixchar=="FixMoveLinear") call init_fix_move_linear
                if(fixchar=="fix_move_wiggle" .or. fixchar=="FixMovewiggle") call init_fix_move_wiggle
                if(fixchar=="fix_set_pos" .or. fixchar=="FixSetPos") call init_fix_set_pos
                if(fixchar=="fix_rigid" .or. fixchar=="FixRigid") call init_fix_rigid
                if(fixchar=="fix_rigid_rotate" .or. fixchar=="FixRigidRotate") call init_fix_rigid_rotate
                if(fixchar=="fix_freeze" .or. fixchar=="FixFreeze") call init_fix_freeze
                ! --- calc fcnumi ---
                if(size(fcnumi) < ifc)then
                    do i=1,ifc-size(fcnumi)
                        fcnumi = [fcnumi, -1]
                    end do
                end if
                fcnumi(ifc) = nfc
            end do
            !write(*,*) "all nfc=", nfc
            close(ifilecalc)
        end if
    end subroutine init_fix

    subroutine fix
        use fix_freeze_mod, only : fix_freeze
        use fix_vel_scaling_mod, only : fix_vel_scaling
        use fix_langevin_mod, only : fix_langevin
        use fix_langevin2_mod, only : fix_langevin2
        use fix_langevin3_mod, only : fix_langevin3
        use fix_langevin_LAMMPS_mod, only : fix_langevin_LAMMPS
        use fix_NoseHoover_mod, only : fix_NoseHoover
        use fix_move_linear_mod, only : fix_move_linear
        use fix_move_wiggle_mod, only : fix_move_wiggle
        use fix_set_pos_mod, only : fix_set_pos
        use fix_addforce_mod, only : fix_addforce
        use fix_rigid_mod, only : fix_rigid
        use fix_rigid_rotate_mod, only : fix_rigid_rotate
        implicit none
        !
        call fix_freeze
        call fix_vel_scaling
        call fix_langevin
        call fix_langevin2
        call fix_langevin3
        call fix_langevin_LAMMPS
        call fix_NoseHoover
        call fix_move_linear
        call fix_move_wiggle
        call fix_set_pos
        call fix_addforce
        call fix_rigid
        call fix_rigid_rotate
    end subroutine fix

    subroutine fix_compute
        use commn
        !use fix_freeze_mod
        use fix_vel_scaling_mod, only : fix_vel_scaling_compute
        use fix_langevin_mod, only : fix_langevin_compute
        use fix_langevin2_mod, only : fix_langevin2_compute
        use fix_langevin3_mod, only : fix_langevin3_compute
        use fix_langevin_LAMMPS_mod, only : fix_langevin_LAMMPS_compute
        use fix_NoseHoover_mod, only : fix_NoseHoover_compute
        use fix_move_linear_mod, only : fix_move_linear_compute
        use fix_move_wiggle_mod, only : fix_move_wiggle_compute
        use fix_set_pos_mod, only : fix_set_pos_compute
        use fix_addforce_mod, only : fix_addforce_compute
        !use fix_rigid_mod
        implicit none

        call fix_vel_scaling_compute
        call fix_langevin_compute
        call fix_langevin2_compute
        call fix_langevin3_compute
        call fix_langevin_LAMMPS_compute
        call fix_NoseHoover_compute
        call fix_move_linear_compute
        call fix_move_wiggle_compute
        call fix_set_pos_compute
        call fix_addforce_compute

    end subroutine fix_compute
end module fix_mod
