! compute temperature
module compute_gg_mod
    integer,allocatable :: myfc(:)
    integer,parameter :: ndata = 4
    integer :: nlist = 0
    integer ilist
    integer,allocatable :: nsamp(:), atlisti(:), atlistj(:)
    double precision,allocatable :: cforce(:,:), cforce_Iam(:,:), cforce_pair(:,:), cforce_bonds(:,:)
    !
contains
    subroutine init_compute_gg
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
        read(ifilecalc,*) iatype, jatype, i
        if(Iam==master) write(*,'(A,3I)') __FILE__, iatype, jatype, i
        !
        nsamp = [nsamp, i]
        atlisti = [atlisti, iatype]
        atlistj = [atlistj, jatype]
        !
        cforce_Iam = reshape( [0d0], [4,nlist+1], pad=[0d0] )
        if(Iam==master)then
            cforce = reshape( [0d0], [4,nlist+1], pad=[0d0] )
            cforce_pair = reshape( [0d0], [4,nlist+1], pad=[0d0] )
            cforce_bonds = reshape( [0d0], [4,nlist+1], pad=[0d0] )
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
    end subroutine init_compute_gg
    !
    !
    ! --- compute temperature ---
    subroutine compute_gg
        use update_mod
        use commn
        use constant_mod
        use fix_compute_mod
        use mpivar
        implicit none
        integer j,i,imol
        !
        ! --- calculate force ---
        do ilist=1,nlist
            if(mod(irep,nsamp(ilist)) /= 0) cycle
            call calc_force
        end do
        !
        ! --- return ---
        if(Iam==master)then
            do i=1,nlist
                if(mod(irep,nsamp(i)) /= 0) cycle
                fcdata(i0fcdata(myfc(i))+1:i0fcdata(myfc(i))+ndata) = cforce(1:4,i)
            end do
        end if
    end subroutine compute_gg

    subroutine calc_force
        use commn
        use mpivar
        use update_mod
        use force_mod
        implicit none
        !
        call calc_force_pair
        call calc_force_bonds
        if(Iam==master) cforce(:,ilist) = cforce_pair(:,ilist) + cforce_bonds(:,ilist)
        !
    end subroutine calc_force

    subroutine calc_force_pair
        use commn
        use mpivar
        use update_mod
        use force_mod
        use book_keep_mod, only: nlistbk=>nlist, listbki=>listi, listbkj=>listj
        implicit none
        double precision r(3), rabssq, r8
        integer ilistbk,li,lj
        integer :: typeflag = 0
        !
        cforce_Iam(:,ilist) = 0d0
        do ilistbk=1,nlistbk
            li = listbki(ilistbk)
            lj = listbkj(ilistbk)
            !
            if(atlisti(ilist)==atype(li).and.atlistj(ilist)==atype(lj)) typeflag = 1
            if(atlisti(ilist)==atype(lj).and.atlistj(ilist)==atype(li)) typeflag = -1
            if(typeflag==0) cycle
            !
            r(:) = p(:,li) - p(:,lj)
            if (r(1)>vlh(1)) r(1) = r(1) - vl(1)
            if (r(1)<-vlh(1)) r(1) = r(1) + vl(1)
            if (r(2)>vlh(2)) r(2) = r(2) - vl(2)
            if (r(2)<-vlh(2)) r(2) = r(2) + vl(2)
            if (r(3)>vlh(3)) r(3) = r(3) - vl(3)
            if (r(3)<-vlh(3)) r(3) = r(3) + vl(3)
            !
            rabssq = r(1)*r(1)+r(2)*r(2)+r(3)*r(3) 
            if(rabssq > rcsq) then
                typeflag = 0
                cycle
            end if
            !
            r8 = fpair_one(eps(atype(li),atype(lj)), sig(atype(li),atype(lj)), rcsq, rabssq, cforce_Iam(4,ilist))
            cforce_Iam(1:3,ilist) = cforce_Iam(1:3,ilist) + typeflag*r8*r(:)
            typeflag = 0
        end do
        !
        call mpi_reduce(cforce_Iam(:,ilist), cforce_pair(:,ilist), 4, mpi_real8, &
            mpi_sum, master, mpi_comm_world, ierr)
        !
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
            !
            if(atlisti(ilist)==atype(li).and.atlistj(ilist)==atype(lj)) typeflag = 1
            if(atlisti(ilist)==atype(lj).and.atlistj(ilist)==atype(li)) typeflag = -1
            if(typeflag==0) cycle
            !
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
            cforce_Iam(4,ilist) = cforce_Iam(4,ilist) + 0.5d0*k(btype(i))*(rabs - r_eq(btype(i)))**2d0
            cforce_Iam(1:3,ilist) = cforce_Iam(1:3,ilist) + typeflag*fij*r(:)
            typeflag = 0
        end do
        !
        call mpi_reduce(cforce_Iam(:,ilist), cforce_bonds(:,ilist), 4, mpi_real8, &
            mpi_sum, master, mpi_comm_world, ierr)
        !
    end subroutine calc_force_bonds

end module compute_gg_mod

