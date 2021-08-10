program MD
    use commn
    use mpivar
    use force_mod, only : dpssqmax
    use book_keep_mod, only : rmarsq4
    use time_mod
    use fix_mod
    use output_mod, only : output, last_output
    implicit none
    !
    call mpi_init(ierr)                             !! starts mpi
    call mpi_comm_rank(mpi_comm_world, Iam, ierr)   !! get current proc id
    call mpi_comm_size(mpi_comm_world, Nproc, ierr) !! get number of procs
    !
    call random_seed_clock
    !
    call read_etc
    !
    call init_const
    call file_open
    !
    call init_book_keep
    !
    ! === main loop ===
    timeall0 = start_time()
    if(Iam==master) write(*,*) "--- main loop start ---"
    call book_keep
    call calc_force
    call calculation
    call output
    do irep=1,nrep
        call cp_buf
        ! state:  p(t-dt),f(t-dt),v(t-dt)
        !
        ! --- update atom position ---
        call update_pos
        ! state:  p(t),f(t-dt),v(t-dt)
        !
        ! --- book keeping ---
        if(dpssqmax > rmarsq4) call book_keep
        ! --- calculate atom force ---
        call calc_force
        ! state:  p(t),f(t),v(t-dt)
        !
        ! --- update atom velocity ---
        call update_vel
        ! state:  p(t),f(t),v(t)
        !
        ! --- calculation value ---
        call calculation
        !
        call output
    end do
    !
    if(Iam==master) write(*,*) "--- main loop end ---"
    timeall = stop_time(timeall0)
    ! =================
    !
    call output_last
    !
    call last_output
    call file_close
    !
    call mpi_finalize(ierr)   !!let mpi finish up ...
    !
    stop
end program MD


subroutine read_etc
    use commn
    use update_mod
    use book_keep_mod
    use mpivar
    use file_mod
    use force_mod
    use fix_mod, only : init_fix
    use compute_mod , only : init_compute
    use output_mod, only : init_output
    implicit none
    integer i,j,ios,imol
    double precision r8a, r8b
    character(128) char_cmt,chartest
    character(1) :: sharp = "#"
    !
    ! === read calc.dat (half way) ===
    if(command_argument_count()==0)then
        filecalc="./calc.dat"
    else
        call get_command_argument(1,filecalc)
    end if
    if(Iam==master) write(*,*) "calc: ",trim(filecalc)
    open(ifilecalc,file=trim(filecalc),action="read",status="old")
    !
    read(ifilecalc,'(A)') fileinit
    read(ifilecalc,'(A)') filelast
    read(ifilecalc,'(A)') resultdir
    fileinit = char_cmt(fileinit,sharp,len(sharp))
    filelast = char_cmt(filelast,sharp,len(sharp))
    resultdir = char_cmt(resultdir,sharp,len(sharp))
    read(ifilecalc,*) dt
    read(ifilecalc,*) nrep
    read(ifilecalc,*) nsamp_term
    if(nsamp_term > nrep) call error_msg("nsamp_term > nrep @read_etc")
    !
    ! === read init.dat ===
    open(ifileinit,file=trim(fileinit),action="read",status="old")
    !
    read(ifileinit,*)
    read(ifileinit,*)
    read(ifileinit,*) nmol
    read(ifileinit,*) natype
    write(*,*) Iam, "before nbonds"
    nbonds = 0
    read(ifileinit,*,iostat=ios) nbonds
    write(*,*) Iam, "after nbonds"
    if(ios/=0 .or. nbonds==0)then
        nbonds = 0
        nbtype = 0
        write(*,*) Iam, "aaa"
        rewind(ifileinit)
        read(ifileinit,*)
        read(ifileinit,*)
        read(ifileinit,*)
        read(ifileinit,*)
    else
        read(ifileinit,*) nbtype
        write(*,*) Iam, nbtype, nbonds
    end if
    !
    allocate(wm(1:natype), k(1:nbtype), r_eq(1:nbtype))
    allocate(sig(1:natype,1:natype),eps(1:natype,1:natype))
    allocate(pairflag(1:natype,1:natype))
    pairflag(1:natype,1:natype) = .false.
    allocate(atype(1:nmol), btype(1:nbonds), bondi(1:nbonds), bondj(1:nbonds))
    allocate(p(1:3,1:nmol))
    allocate(f_Iam(1:3,1:nmol))
    !
    if(Iam==master)then
        allocate(v(1:3,1:nmol),f(1:3,1:nmol),fbuf(1:3,1:nmol))
        allocate(f_pair(1:3,1:nmol), f_bonds(1:3,1:nmol))
        allocate(vfixflag(1:3,1:natype),ffixflag(1:3,1:natype),pfixflag(1:3,1:natype))
        vfixflag(:,:) = .false.
        ffixflag(:,:) = .false.
        pfixflag(:,:) = .false.
    end if
    !
    read(ifileinit,*)
    read(ifileinit,*) lo(1), hi(1)
    read(ifileinit,*) lo(2), hi(2)
    read(ifileinit,*) lo(3), hi(3)
    vl(1:3) = hi(1:3)-lo(1:3)
    vlh(1:3) = 0.5d0*vl(1:3)
    va(1) = vl(2)*vl(3)
    va(2) = vl(3)*vl(1)
    va(3) = vl(1)*vl(2)
    read(ifileinit,*)
    read(ifileinit,*)
    read(ifileinit,*)
    do j=1,natype
        read(ifileinit,*) i, wm(i)
    end do
    if(nbtype/=0)then
        read(ifileinit,*)
        read(ifileinit,*)
        read(ifileinit,*)
        do j=1,nbtype
            read(ifileinit,*) i,k(i), r_eq(i)
        end do
    end if
        read(ifileinit,*)
        read(ifileinit,*)
        read(ifileinit,*)
        do imol=1,nmol
            read(ifileinit,*) i, j, atype(i), p(:,i)
        end do
        read(ifileinit,*)
        read(ifileinit,*)
        read(ifileinit,*)
        do imol=1,nmol
            read(ifileinit,*) i, v(:,i)
        end do
    if(nbtype/=0)then
        read(ifileinit,*)
        read(ifileinit,*)
        read(ifileinit,*)
        do j=1,nbonds
            read(ifileinit,*) i,btype(i), bondi(i), bondj(i)
        end do
    end if
    close(ifileinit)
    !
    ! --- create atom type list ---
    ! imolatseq: sort imol by atype, impute=imol, return=imol
    ! i0atseq: initial position of imolatseq, imput=iatype, return=imol 
    ! natseq: number of atype molecules, imput=iatype, return=number_of_atype_molecule
    ! ex) atype(imolatseq(i0atseq(iatype))) = iatype ! initial postion of imolatseq
    ! ex) atype(imolatseq(i0atseq(iatype)+natseq(iatype)-1)) = iatype ! last postion of imolatseq
    if(Iam==master)then
        allocate(natseq(1:natype), i0atseq(1:natype), imolatseq(1:nmol))
        ! --- calculate natseq ---
        natseq(:) = 0
        do i=1,nmol
            natseq(atype(i)) = natseq(atype(i)) + 1
        end do
        ! --- calculate i0atseq ---
        i0atseq(:) = 0
        do i=2,natype
            i0atseq(i) = i0atseq(i-1) + natseq(i-1)
        end do
        ! --- calculate imolatseq ---
        imol = 1
        do j=1,natype
            do i=1,nmol
                if(atype(i) == j)then
                    imolatseq(imol) = i
                    imol = imol + 1
                end if
            end do
        end do
    end if
    !
    ! === calc.dat (second half) ===
    eps(1:natype,1:natype) = -1d0
    sig(1:natype,1:natype) = -1d0
    read(ifilecalc,*) pairstyle
    if(trim(pairstyle)=="smooth_quad") smooth_quad = .true.
    if(trim(pairstyle)=="smooth_linear") smooth_linear = .true.
    do
        read(ifilecalc,*,iostat=ios) i,j ,r8a,r8b
        if(ios/=0)exit
        pairflag(i,j) = .true.
        pairflag(j,i) = .true.
        eps(i,j) = r8a
        sig(i,j) = r8b
        eps(j,i) = eps(i,j)
        sig(j,i) = sig(i,j)
    end do
    close(ifilecalc)
    !
    call init_fix
    call init_compute
    call init_output
    !
end subroutine read_etc


subroutine init_const
    use commn
    use mpivar
    use update_mod
    use time_mod
    use force_mod
    implicit none
    !
    irep = 0
    allocate(dt_wmi2(1:natype))
    if(Iam==master)then
        allocate(dps(1:3,1:nmol))
        dps(1:3,1:nmol) = 0d0
        allocate(dp(1:3,1:nmol))
        dp(1:3,1:nmol) = 0d0
    end if
    dt_wmi2(1:natype) = 0.5d0*dt/wm(1:natype)
    timebk = 0d0
    timeforce = 0d0
    eadd = 0d0
    rc = 3.5d0
    rcsq = rc**2d0
end subroutine init_const


subroutine update_pos
    use commn
    use mpivar
    use update_mod
    use fix_mod
    use force_mod, only : dpssqmax
    implicit none
    integer i,j,l,m
    !
    if(Iam==master)then
        pflag = .true.
        fflag = .false.
        vflag = .false.
        ! --- calculate dp ---
        do i=1,nmol
            dp(:,i) = dt*(v(:,i) + f(:,i)*dt_wmi2(atype(i)))
        end do
        ! --- calculate dpssqmax ---
        dps(:,:) = dps(:,:) + dp(:,:)
        dpssqmax = 0d0
        do i=1,nmol
            if(dpssqmax < sum(dps(:,i)**2d0)) dpssqmax = sum(dps(:,i)**2d0)
        end do
        ! --- update position ---
        !
        call fix
        p(:,:) = p(:,:) + dp(:,:)
        !
        do i=1,nmol
            if (p(1,i)>hi(1)) p(1,i) = p(1,i) - vl(1)
            if (p(1,i)<lo(1)) p(1,i) = p(1,i) + vl(1)
            !
            if (p(2,i)>hi(2)) p(2,i) = p(2,i) - vl(2)
            if (p(2,i)<lo(2)) p(2,i) = p(2,i) + vl(2)
            !              
            if (p(3,i)>hi(3)) p(3,i) = p(3,i) - vl(3)
            if (p(3,i)<lo(3)) p(3,i) = p(3,i) + vl(3)
        end do
        !
    end if
    call mpi_bcast(p, nmol*3, mpi_real8, master, &
        mpi_comm_world, ierr)
    call mpi_bcast(dpssqmax, 1, mpi_real8, master, &
        mpi_comm_world, ierr)
    !
end subroutine update_pos

subroutine init_book_keep
    use book_keep_mod
    use force_mod
    implicit none
    integer i,j
    !
    maxnlist = 1
    allocate(listi(1:maxnlist))
    allocate(listj(1:maxnlist))
    allocate(listi_gather(1:maxnlist))
    allocate(listj_gather(1:maxnlist))
    rmar = 1.3d0
    rmarsq4 = rmar**2d0/4d0
    rbk = rc+rmar
    rbksq = rbk**2d0
end subroutine init_book_keep

subroutine book_keep
    use commn
    use mpivar
    use update_mod , only : dps, atype
    use book_keep_mod
    use time_mod
    use force_mod, only : pairflag
    implicit none
    double precision r(3),rabssq
    integer i,j,l,m,nlistmax , nlistmin, iatpairlist, iatype, jatype, m0
    integer nlist_gather(1:Nproc), displs(1:Nproc)
    real(4) t0
    !
    t0 = start_time()
    nlist = maxnlist+1
    do while(nlist>maxnlist)
        nlist = 0
        do i=1,nmol-1
            do j=i+1+Iam-master,nmol,Nproc
                if( .not. pairflag(atype(i),atype(j)))cycle
                r(3) = p(3,i) - p(3,j)
                if (r(3)>vlh(3)) r(3) = r(3) - vl(3)
                if (r(3)<-vlh(3)) r(3) = r(3) + vl(3)
                if (r(3)>rbk .or. r(3)<-rbk) cycle
                !
                r(1) = p(1,i) - p(1,j)
                if (r(1)>vlh(1)) r(1) = r(1) - vl(1)
                if (r(1)<-vlh(1)) r(1) = r(1) + vl(1)
                if (r(1)>rbk .or. r(1)<-rbk) cycle
                !
                r(2) = p(2,i) - p(2,j)
                if (r(2)>vlh(2)) r(2) = r(2) - vl(2)
                if (r(2)<-vlh(2)) r(2) = r(2) + vl(2)
                if (r(2)>rbk .or. r(2)<-rbk) cycle
                !
                rabssq = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
                if (rabssq<rbksq) then
                    ! save different book for defferent processors
                    nlist = nlist + 1
                    if(nlist<=maxnlist) then
                        listi(nlist) = i
                        listj(nlist) = j
                    end if
                    !
                end if
            end do
        end do
        call mpi_allreduce(nlist, nlistmax, 1, mpi_integer, &
            mpi_max, mpi_comm_world, ierr)
        if(nlistmax>maxnlist) then
            deallocate(listi)
            deallocate(listj)
            maxnlist = (nlistmax+1) * 1.1
            allocate(listi(1:maxnlist))
            allocate(listj(1:maxnlist))
            nlist = maxnlist + 1
        end if
    end do
    call mpi_allgather(nlist, 1, mpi_integer, nlist_gather, 1, mpi_integer&
        , mpi_comm_world, ierr)
    nlistall = sum(nlist_gather(1:Nproc))
    nlistmin = minval(nlist_gather(1:Nproc))
    nlist_gather(1:Nproc) = nlist_gather(1:Nproc) - nlistmin
    displs(1) = 0
    do i=2,Nproc
        displs(i) = displs(i-1) + nlist_gather(i-1)
    end do
    deallocate(listi_gather)
    deallocate(listj_gather)
    allocate(listi_gather(1:nlistall))
    allocate(listj_gather(1:nlistall))
    call mpi_gatherv(listi(nlistmin+1), nlist-nlistmin, mpi_integer, &
        listi_gather, nlist_gather, displs, mpi_integer, &
        master, mpi_comm_world, ierr)
    call mpi_gatherv(listj(nlistmin+1), nlist-nlistmin, mpi_integer, &
        listj_gather, nlist_gather, displs, mpi_integer, &
        master, mpi_comm_world, ierr)
    nlist = (nlistall+Iam-master)/Nproc
    do i=1,Nproc
        nlist_gather(i) = (nlistall - nlistmin*Nproc + i-1)/Nproc
    end do
    displs(1) = 0
    do i=2,Nproc
        displs(i) = displs(i-1) + nlist_gather(i-1)
    end do
    call mpi_scatterv(listi_gather, nlist_gather, displs, mpi_integer, &
        listi(nlistmin+1), nlist-nlistmin, mpi_integer, &
        master, mpi_comm_world, ierr)
    call mpi_scatterv(listj_gather, nlist_gather, displs, mpi_integer, &
        listj(nlistmin+1), nlist-nlistmin, mpi_integer, &
        master, mpi_comm_world, ierr)
    !
    if(Iam==master) dps(1:3,1:nmol) = 0d0
    timebk = timebk + stop_time(t0)
end subroutine book_keep


subroutine calc_force
    use commn
    use mpivar
    use update_mod
    use time_mod
    use fix_mod
    use force_mod
    implicit none
    real(4) t0
    integer i,j,l,m
    !
    t0 = start_time()
    call calc_force_pair
    call calc_force_bonds
    if(Iam==master)then
        pflag = .false.
        fflag = .true.
        vflag = .false.
        f(:,:) = f_pair(:,:) + f_bonds(:,:)
        ep = ep_pair + ep_bonds
        call fix
    end if
    timeforce = timeforce + stop_time(t0)
    !
end subroutine calc_force

subroutine calc_force_pair
    use commn
    use mpivar
    use update_mod
    use force_mod
    use book_keep_mod, only : listi, listj, nlist
    implicit none
    double precision r(3), rabssq, r8
    integer i,li,lj
    !
    ep_Iam = 0d0
    f_Iam(:,:) = 0d0
    do i=1,nlist
        li = listi(i)
        lj = listj(i)
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
        r8 = fij_lj(eps(atype(li),atype(lj)), sig(atype(li),atype(lj)), rcsq, rabssq, ep_Iam)
        f_Iam(:,li) = f_Iam(:,li) + r8*r(:)
        f_Iam(:,lj) = f_Iam(:,lj) - r8*r(:)
    end do
    !
    call mpi_reduce(f_Iam, f_pair, nmol*3, mpi_real8, &
        mpi_sum, master, mpi_comm_world, ierr)
    call mpi_reduce(ep_Iam, ep_pair, 1, mpi_real8, &
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
    !
    ep_Iam = 0d0
    f_Iam(1:3,1:nmol) = 0d0
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
        ep_Iam = ep_Iam + 0.5d0*k(btype(i))*(rabs - r_eq(btype(i)))**2d0
        f_Iam(:,li) = f_Iam(:,li) + fij*r(:)
        f_Iam(:,lj) = f_Iam(:,lj) - fij*r(:)
    end do
    !
    call mpi_reduce(f_Iam, f_bonds, nmol*3, mpi_real8, &
        mpi_sum, master, mpi_comm_world, ierr)
    call mpi_reduce(ep_Iam, ep_bonds, 1, mpi_real8, &
        mpi_sum, master, mpi_comm_world, ierr)
    !
end subroutine calc_force_bonds


subroutine update_vel
    use commn
    use mpivar
    use update_mod
    use fix_mod
    implicit none
    integer i,j,l,m
    !
    if(Iam==master)then
        pflag = .false.
        fflag = .false.
        vflag = .true.
        temp = 0d0
        ek = 0d0
        do i=1,nmol
            v(:,i) = v(:,i) + (fbuf(:,i)+f(:,i))*dt_wmi2(atype(i))
        end do
        call fix
    end if
    !
end subroutine update_vel

subroutine calculation
    use mpivar
    use fix_mod, only : fix_compute
    use compute_mod, only : compute
    implicit none
    !
    if(Iam==master) call fix_compute
    call compute
end subroutine calculation

subroutine cp_buf
    use commn
    use mpivar
    use update_mod
    implicit none
    !
    if(Iam==master)then
        fbuf(:,:) = f(:,:)
    end if
end subroutine cp_buf


subroutine output_last
    use commn
    use update_mod
    use mpivar
    use time_mod
    use file_mod
    implicit none
    integer i,j
    !
    if(Iam==master)then
        ! --- terminal ---
        timealls = int(timeall)
        write(*,'(A6,I5.4,":",I2.2,":",I2.2,":",I2.2)') "time: "&
            ,timealls/60/60, mod(timealls/60,60), mod(timealls,60),&
            int(100d0*(timeall-dble(timealls)))
        write(*,*) "--- time detail ---"
        write(*,'(4A10)') "all", "book_keep", "force", "other"
        write(*,'(4f10.3)')  timeall, timebk, timeforce, timeall-timebk-timeforce
        !
        ! --- last.dat ---
        open(ifilelast,file=trim(filelast),action="write",status="replace")
        write(ifilelast,*) "homemade init"
        write(ifilelast,*) 
        write(ifilelast,*) nmol, "atoms"
        write(ifilelast,*) natype, "atom types"
        write(ifilelast,*) nbonds, "bonds"
        write(ifilelast,*) nbtype, "bond types"
        write(ifilelast,*)
        write(ifilelast,*) lo(1), hi(1), "xlo xhi"
        write(ifilelast,*) lo(2), hi(2), "ylo yhi"
        write(ifilelast,*) lo(3), hi(3), "zlo zhi"
        write(ifilelast,*)
        write(ifilelast,*) "Masses"
        write(ifilelast,*)
        do i=1,natype
            write(ifilelast,*) i, wm(i)
        end do
        write(ifilelast,*)
        write(ifilelast,*) "Bond Coeffs # harmonic"
        write(ifilelast,*)
        do i=1,nbtype
            write(ifilelast,*) i,k(i), r_eq(i)
        end do
        write(ifilelast,*)
        write(ifilelast,*) "Atoms # bond"
        write(ifilelast,*)
        do i=1,nmol
            write(ifilelast,'(3I,3E)') i, 0, atype(i), p(:,i)!, 0,0,0
        end do
        write(ifilelast,*)
        write(ifilelast,*) "Velocities"
        write(ifilelast,*)
        do i=1,nmol
            write(ifilelast,'(I,3E)') i, v(:,i)
        end do
        write(ifilelast,*)
        write(ifilelast,*) "Bonds"
        write(ifilelast,*)
        do i=1,nbonds
            write(ifilelast,*) i,btype(i), bondi(i),bondj(i)
        end do
        close(ifilelast)
    end if
    !
end subroutine output_last


subroutine file_open
    use file_mod
    use mpivar
    implicit none
    real(4) r43(1:3)
    integer i
    !
    fileposb = 11
    if(Iam==master)then
        open(fileposb,file=trim(resultdir)//'posb.dat',access="stream",form="unformatted")
    end if
end subroutine file_open


subroutine file_close
    use file_mod
    use mpivar
    if(Iam==master)then
        close(fileposb)
    end if
end subroutine file_close

