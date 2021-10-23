! compute temperature
module compute_gg_VAstress_mod
    integer,allocatable :: myfc(:)
    integer,parameter :: ndata = 9
    integer :: nlist = 0
    integer ilist
    integer,allocatable :: nsamp(:)
    double precision,allocatable :: cforce(:,:), cforce_Iam(:,:), cforce_pair(:,:), cforce_bonds(:,:)
    double precision,allocatable :: stress_int_Iam(:,:,:)
    double precision,allocatable :: stress_kin(:,:,:),stress_int(:,:,:),stress(:,:,:)
    !
contains
    subroutine init_compute_gg_VAstress
        use update_mod
        use commn
        use file_mod
        use constant_mod
        use fix_compute_mod
        use mpivar
        implicit none
        integer i, ios, iatype, jatype
        double precision,allocatable :: a(:)
        !
        if(Iam==master) nfc = nfc + 1
        ! === read conf file ===
        read(ifilecalc,*) i
        if(Iam==master) write(*,'(A,I)') __FILE__, i
        !
        nsamp = [nsamp, i]
        !
        stress_int_Iam = reshape( [0d0], [3,3,nlist+1], pad=[0d0] )
        if(Iam==master)then
            stress_kin = reshape( [0d0], [3,3,nlist+1], pad=[0d0] )
            stress_int = reshape( [0d0], [3,3,nlist+1], pad=[0d0] )
            stress = reshape( [0d0], [3,3,nlist+1], pad=[0d0] )
        end if
        nlist = nlist + 1
        ! --- return ---
        if(Iam==master)then
            fcnum = [fcnum, ifc]
            myfc = [myfc, nfc]
            nfcdata = [nfcdata, ndata]
            ! --- calc fcnumi ---
            if(size(fcnumi) < ifc)then
                do i=1,ifc-size(fcnumi)
                    fcnumi = [fcnumi, -1]
                end do
            end if
            fcnumi(ifc) = nfc
        end if
        ! ===
    end subroutine init_compute_gg_VAstress
    !
    !
    ! --- compute temperature ---
    subroutine compute_gg_VAstress
        use update_mod
        use commn
        use constant_mod
        use fix_compute_mod
        use mpivar
        implicit none
        integer j,i,imol,l
        !
        ! --- calculate stress ---
        do ilist=1,nlist
            if(mod(irep,nsamp(ilist)) /= 0) cycle
            stress_int_Iam(:,:,ilist) = 0d0
            stress_int(:,:,ilist) = 0d0
            stress_kin(:,:,ilist) = 0d0
            call calc_force
            call mpi_reduce(stress_int_Iam(1,1,ilist), stress_int(1,1,ilist), 9, mpi_real8, &
                mpi_sum, master, mpi_comm_world, ierr)
            if(Iam==master) stress_int(:,:,ilist) = stress_int(:,:,ilist)/(vl(1)*vl(2)*vl(3))
            call calc_stress_kin
            if(Iam==master) stress(:,:,ilist) = stress_kin(:,:,ilist) + stress_int(:,:,ilist)
        end do
        !
        ! --- return ---
        if(Iam==master)then
            do i=1,nlist
                if(mod(irep,nsamp(i)) /= 0) cycle
                do j=1,3
                    do l=1,3
                        fcdata(i0fcdata(myfc(i))+(j-1)*3+l) = stress(j,l,i)
                    end do
                end do
            end do
        end if
    end subroutine compute_gg_VAstress

    subroutine calc_force
        use commn
        use mpivar
        use update_mod
        use force_mod
        implicit none
        !
        call calc_force_pair
        call calc_force_bonds
        !
    end subroutine calc_force

    subroutine calc_force_pair
        use commn
        use mpivar
        use update_mod
        use force_mod
        use book_keep_mod, only: nlistbk=>nlist, listbki=>listi, listbkj=>listj
        implicit none
        double precision r(3), rabssq, r8, rgomi
        integer ilistbk,li,lj
        integer :: typeflag = 0
        !
        cforce_Iam(:,ilist) = 0d0
        do ilistbk=1,nlistbk
            li = listbki(ilistbk)
            lj = listbkj(ilistbk)
            r(:) = p(:,li) - p(:,lj)
            if (r(1)>vlh(1)) r(1) = r(1) - vl(1)
            if (r(1)<-vlh(1)) r(1) = r(1) + vl(1)
            if (r(2)>vlh(2)) r(2) = r(2) - vl(2)
            if (r(2)<-vlh(2)) r(2) = r(2) + vl(2)
            if (r(3)>vlh(3)) r(3) = r(3) - vl(3)
            if (r(3)<-vlh(3)) r(3) = r(3) + vl(3)
            !
            rabssq = r(1)*r(1)+r(2)*r(2)+r(3)*r(3) 
            if(rabssq > rcsq) cycle
            !
            r8 = fpair_one(eps(atype(li),atype(lj)), sig(atype(li),atype(lj)), rcsq, rabssq, rgomi)
            call calc_stress_int(r,r8*r)
        end do
    end subroutine calc_force_pair

    subroutine calc_force_bonds
        use commn
        use mpivar
        use update_mod
        use force_mod
        implicit none
        double precision r(3), rabs, fij
        integer i,li,lj
        integer :: typeflag = 0
        !
        cforce_Iam(:,ilist) = 0d0
        do i=1+Iam,nbonds,Nproc
            li = bondi(i)
            lj = bondj(i)
            r(:) = p(:,li) - p(:,lj)
            if (r(1)>vlh(1)) r(1) = r(1) - vl(1)
            if (r(1)<-vlh(1)) r(1) = r(1) + vl(1)
            if (r(2)>vlh(2)) r(2) = r(2) - vl(2)
            if (r(2)<-vlh(2)) r(2) = r(2) + vl(2)
            if (r(3)>vlh(3)) r(3) = r(3) - vl(3)
            if (r(3)<-vlh(3)) r(3) = r(3) + vl(3)
            !
            rabs = sqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
            fij = - k(btype(i))*(rabs - r_eq(btype(i)))/rabs
            call calc_stress_int(r,fij*r)
        end do
    end subroutine calc_force_bonds

    subroutine calc_stress_int(r,f)
        implicit none
        real(8) r(3),f(3)
        integer i
        !
        do i=1,3
            stress_int_Iam(i,:,ilist) = stress_int_Iam(i,:,ilist) - r(i)*f(:)
        end do
        !
    end subroutine calc_stress_int

    subroutine calc_stress_kin
        use commn
        use mpivar
        use update_mod
        implicit none
        integer i,j,l
        !
        if(Iam==master)then
            do i=1,nmol
                do j=1,3
                    stress_kin(j,:,ilist) = stress_kin(j,:,ilist) - wm(atype(i))*v(j,i)*v(:,i)
                end do
            end do
            stress_kin(:,:,ilist) = stress_kin(:,:,ilist)/(vl(1)*vl(2)*vl(3))
        end if
    end subroutine calc_stress_kin

end module compute_gg_VAstress_mod

