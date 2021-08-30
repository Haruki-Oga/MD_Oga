! fix atom (dim, atype)
module fix_rigid_rotate_mod
    integer,allocatable :: myfc(:)
    integer,parameter :: ndata = 1
    integer :: nlist = 0
    integer,allocatable :: atlist(:)
    double precision,allocatable :: fsum(:,:), fsumbuf(:,:), dpcom(:,:), vcom(:,:), pcom(:,:)
    double precision,allocatable :: mass(:)
    double precision,allocatable :: pb(:,:),quat(:,:),omegab(:,:),omegabdot(:,:)
    double precision,allocatable :: torque(:,:),omegabbuf(:,:),omegabdotbuf(:,:)
    double precision,allocatable :: Ib(:,:)
    integer,parameter :: nrepomega = 20
    !
contains
    subroutine init_fix_rigid_rotate
        use update_mod
        use commn
        use file_mod
        use fix_compute_mod
        implicit none
        integer i,iatype,ndim,int3(1:3)
        integer,allocatable :: a(:)
        double precision,allocatable :: b(:)
        double precision vcomin(3),posr(3)
        integer imol,j
        double precision Im(3,3), Iminv(3,3), Lrot(1:3)
        logical(1) dim_true
        !
        ! === read conf file ===
        read(ifilecalc,*) iatype, ndim
        write(*,*) __FILE__, iatype, ndim
        int3(:) = 0
        do i=1,3
            if(dim_true(ndim,i)) int3(i) = 1
        end do
        !
        atlist = [atlist, iatype]
        !
        b = reshape([fsum(:,:)], [3*nlist])
        fsum = reshape( [b,0d0,0d0,0d0], [3,nlist+1] )
        fsumbuf = reshape( [b,0d0,0d0,0d0], [3,nlist+1] )
        dpcom = reshape( [b,0d0,0d0,0d0], [3,nlist+1] )
        vcom = reshape( [b,0d0,0d0,0d0], [3,nlist+1] )
        pcom = reshape( [b,0d0,0d0,0d0], [3,nlist+1] )
        pcom(:,nlist+1) = center_of_mass(iatype)
        mass = [mass, 0d0]
        quat = reshape( [0d0], [4,nlist+1], pad=[0d0] )
        quat(4,nlist+1) = 1d0
        omegab = reshape( [0d0], [4,nlist+1], pad=[0d0] )
        omegabdot = reshape( [0d0], [4,nlist+1], pad=[0d0] )
        omegabbuf = reshape( [0d0], [4,nlist+1], pad=[0d0] )
        omegabdotbuf = reshape( [0d0], [4,nlist+1], pad=[0d0] )
        ! calculate Ib
        Ib = reshape( [0d0], [3,nlist+1], pad=[0d0] )
        do i=1,natseq(iatype)
            imol = imolatseq(i0atseq(iatype) + i)
            posr = relative_pos(p(:,imol),pcom(:,nlist+1))
            Ib(1,nlist+1) = Ib(1,nlist+1) + wm(iatype)*(posr(2)**2d0+posr(3)**2d0)
            Ib(2,nlist+1) = Ib(2,nlist+1) + wm(iatype)*(posr(3)**2d0+posr(1)**2d0)
            Ib(3,nlist+1) = Ib(3,nlist+1) + wm(iatype)*(posr(1)**2d0+posr(2)**2d0)
        end do
        torque = reshape( [0d0], [3,nlist+1], pad=[0d0] )
        ! calculate pb : position of body coodinate
        if(nlist==0) allocate(pb(1:3,1:nmol))
        do j=1,natseq(iatype)
            imol = imolatseq(i0atseq(iatype) + j)
            pb(:,imol) = relative_pos(p(:,imol),center_of_mass(iatype))
        end do
        ! calculate omegab
        Im(:,:) = 0d0
        Im(1,1) = Ib(1,nlist+1)
        Im(2,2) = Ib(2,nlist+1)
        Im(3,3) = Ib(3,nlist+1)
        do i=1,natseq(iatype)
            imol = imolatseq(i0atseq(iatype) + i)
            posr = relative_pos(p(:,imol),pcom(:,nlist+1))
            Im(1,2) = Im(1,2) - wm(iatype)*posr(1)*posr(2)
            Im(2,3) = Im(2,3) - wm(iatype)*posr(2)*posr(3)
            Im(3,1) = Im(3,1) - wm(iatype)*posr(3)*posr(1)
        end do
        Im(2,1) = Im(1,2)
        Im(3,2) = Im(2,3)
        Im(1,3) = Im(3,1)
        call inverse(Im, Iminv, 3)
        do i=1,natseq(iatype)
            imol = imolatseq(i0atseq(iatype) + i)
            posr = relative_pos(p(:,imol),pcom(:,nlist+1))
            Lrot = Lrot + wm(iatype)*cross_product(posr, v(:,imol))
        end do
        omegab(1,nlist+1) = sum(Iminv(1,:)*Lrot(:))
        omegab(2,nlist+1) = sum(Iminv(2,:)*Lrot(:))
        omegab(3,nlist+1) = sum(Iminv(3,:)*Lrot(:))
        !
        nlist = nlist + 1
        do i=1,natseq(iatype)
            imol = imolatseq(i0atseq(iatype) + i)
            vcom(:,nlist) = vcom(:,nlist) + v(:,imol)
        end do
        vcom(:,nlist) = vcom(:,nlist)/natseq(iatype)
        mass(nlist) = natseq(iatype)*wm(iatype)
        !
        ! --- set fixflag ---
        do i=1,3
            if(pfixflag(i,atlist(nlist))) call error_msg("pfixflag true @"//__FILE__)
            if(vfixflag(i,atlist(nlist))) call error_msg("vfixflag true @"//__FILE__)
            pfixflag(i,atlist(nlist)) = .true.
            vfixflag(i,atlist(nlist)) = .true.
        end do
        ! --- return ---
        myfc = [myfc, nfc]
        nfcdata = [nfcdata, ndata]
        ! ===
    end subroutine init_fix_rigid_rotate
    !
    subroutine fix_rigid_rotate
        use update_mod
        use commn
        implicit none
        integer i,j,imol,jmol,irepomega
        double precision posi(3),r(3),M(4,4),s1,s2,s3,quatdot(4),quat2dot(4),Lambda
        !
        if( pflag )then
            do i=1,nlist
                ! ---translation---
                pcom(:,i) = pcom(:,i) + dt*(vcom(:,i) + fsum(:,i)*dt/mass(i)*0.5d0)
                ! ---rotation---
                ! calculate M(4,4)
                M = Mq(quat(:,i))
                ! dq/dt
                do j=1,4
                    quatdot(j) = 0.5d0*sum(M(j,:)*omegab(:,i))
                end do
                ! d^2q/dt^2
                do j=1,4
                    quat2dot(j) = 0.5d0*sum(M(j,:)*omegabdot(:,i)) - sum(quatdot(:)*quatdot(:))*quat(j,i)
                end do
                ! Lambda
                s1 = sum(quatdot(:)*quatdot(:))
                s2 = sum(quatdot(:)*quat2dot(:))
                s3 = sum(quat2dot(:)*quat2dot(:))
                Lambda = (1d0 - 0.5d0*dt*dt*s1 - sqrt(1d0 - dt*dt*(s1 + dt*(s2 + 0.25d0*dt*(s3-s1*s1)))))/(dt*dt)
                ! update q
                quat(:,i) = quat(:,i) + dt*quatdot(:) + 0.5d0*dt*dt*(quat2dot(:) - 2d0*Lambda*quat(:,i))
                ! write(*,*) sum(quat(:,i)**2d0)
                !
                do j=1,natseq(atlist(i))
                    imol = imolatseq(i0atseq(atlist(i)) + j)
                    dp(:,imol) = relative_pos(pcom(:,i)+body_to_space(pb(:,imol),quat(:,i)),p(:,imol))
                end do
            end do
        end if
        !
        if( fflag )then
            fsumbuf(:,:) = fsum(:,:)
            fsum(:,:) = 0d0
            torque(:,:) = 0d0
            do i=1,nlist
                do j=1,natseq(atlist(i))
                    imol = imolatseq(i0atseq(atlist(i)) + j)
                    fsum(:,i) = fsum(:,i) + f(:,imol)
                    torque(:,i) = torque(:,i) + cross_product(relative_pos(p(:,imol),pcom(:,i)), f(:,imol))
                end do
                torque(:,i) = space_to_body(torque(:,i), quat(:,i))
            end do
        end if
        !
        if( vflag )then
            do irepomega=1,nrepomega
                omegab(:,:) = omegabbuf(:,:) + 0.5d0*dt*(omegabdotbuf(:,:)+omegabdot(:,:))
                omegabdot(1,:) = (torque(1,:) + omegab(2,:)*omegab(3,:)*(Ib(2,:)-Ib(3,:)))/Ib(1,:)
                omegabdot(2,:) = (torque(2,:) + omegab(3,:)*omegab(1,:)*(Ib(3,:)-Ib(1,:)))/Ib(2,:)
                omegabdot(3,:) = (torque(3,:) + omegab(1,:)*omegab(2,:)*(Ib(1,:)-Ib(2,:)))/Ib(3,:)
                omegabdot(4,:) = 0d0
            end do
            omegabbuf(:,:) = omegab(:,:)
            omegabdotbuf(:,:) = omegabdot(:,:)
            do i=1,nlist
                vcom(:,i) = vcom(:,i) + (fsumbuf(:,i)+fsum(:,i))*dt/mass(i)*0.5d0
                !
                do j=1,natseq(atlist(i))
                    imol = imolatseq(i0atseq(atlist(i)) + j)
                    v(:,imol) = vcom(:,i) + body_to_space(cross_product(omegab(1:3,i), pb(:,imol)), quat(:,i))
                end do
            end do
        end if
    end subroutine fix_rigid_rotate

    function center_of_mass(atom_type)
        use update_mod
        use commn
        implicit none
        double precision center_of_mass(1:3)
        integer atom_type
        double precision posi(1:3), r(1:3)
        integer j,imol,jmol
        !
        center_of_mass(:) = 0d0
        do j=1,natseq(atom_type)
            imol = imolatseq(i0atseq(atom_type) + j)
            posi(:) = p(:,imol)
            if(j/=1)then
                jmol = imolatseq(i0atseq(atom_type) + j-1)
                r(:) = posi(:) - p(:,jmol)
                if(r(1) < -vlh(1)) posi(1) = posi(1) + vl(1)
                if(r(1) > +vlh(1)) posi(1) = posi(1) - vl(1)
                if(r(2) < -vlh(2)) posi(2) = posi(2) + vl(2)
                if(r(2) > +vlh(2)) posi(2) = posi(2) - vl(2)
                if(r(3) < -vlh(3)) posi(3) = posi(3) + vl(3)
                if(r(3) > +vlh(3)) posi(3) = posi(3) - vl(3)
            end if
            center_of_mass(:) = center_of_mass(:) + posi(:)
        end do
        center_of_mass(:) = center_of_mass(:)/natseq(atom_type)
        if(center_of_mass(1) < lo(1)) center_of_mass(1) = center_of_mass(1) + vl(1)
        if(center_of_mass(1) > hi(1)) center_of_mass(1) = center_of_mass(1) - vl(1)
        if(center_of_mass(2) < lo(2)) center_of_mass(2) = center_of_mass(2) + vl(2)
        if(center_of_mass(2) > hi(2)) center_of_mass(2) = center_of_mass(2) - vl(2)
        if(center_of_mass(3) < lo(3)) center_of_mass(3) = center_of_mass(3) + vl(3)
        if(center_of_mass(3) > hi(3)) center_of_mass(3) = center_of_mass(3) - vl(3)
        !
    end function center_of_mass

    function relative_pos(pj,pi)
        use commn
        implicit none
        double precision relative_pos(1:3)
        double precision pi(1:3),pj(1:3)
        !
        relative_pos(:) = pj(:) - pi(:)
        if(relative_pos(1) < -vlh(1)) relative_pos(1) = relative_pos(1) + vl(1)
        if(relative_pos(1) > +vlh(1)) relative_pos(1) = relative_pos(1) - vl(1)
        if(relative_pos(2) < -vlh(2)) relative_pos(2) = relative_pos(2) + vl(2)
        if(relative_pos(2) > +vlh(2)) relative_pos(2) = relative_pos(2) - vl(2)
        if(relative_pos(3) < -vlh(3)) relative_pos(3) = relative_pos(3) + vl(3)
        if(relative_pos(3) > +vlh(3)) relative_pos(3) = relative_pos(3) - vl(3)
    end function relative_pos

    function Mq(q)
        implicit none
        double precision Mq(4,4),q(4)
        !
        Mq(1,1) = -q(3)
        Mq(1,2) = -q(4)
        Mq(1,3) = +q(2)
        Mq(1,4) = +q(1)
        !
        Mq(2,1) = +q(4)
        Mq(2,2) = -q(3)
        Mq(2,3) = -q(1)
        Mq(2,4) = +q(2)
        !
        Mq(3,1) = +q(1)
        Mq(3,2) = +q(2)
        Mq(3,3) = +q(4)
        Mq(3,4) = +q(3)
        !
        Mq(4,1) = -q(2)
        Mq(4,2) = +q(1)
        Mq(4,3) = -q(3)
        Mq(4,4) = +q(4)
        !
    end function Mq

    function body_to_space(a,q)
        implicit none
        double precision a(3),q(4),body_to_space(3)
        !
        body_to_space(1) = &
            + (-q(1)*q(1)+q(2)*q(2)-q(3)*q(3)+q(4)*q(4)) *a(1)&
            - 2d0*(q(1)*q(2)+q(3)*q(4))                  *a(2)&
            + 2d0*(q(2)*q(3)-q(1)*q(4))                  *a(3)
        body_to_space(2) = &
            + 2d0*(q(3)*q(4)-q(1)*q(2))                  *a(1)&
            + (+q(1)*q(1)-q(2)*q(2)-q(3)*q(3)+q(4)*q(4)) *a(2)&
            - 2d0*(q(1)*q(3)+q(2)*q(4))                  *a(3)
        body_to_space(3) = &
            + 2d0*(q(2)*q(3)+q(1)*q(4))                  *a(1)&
            + 2d0*(q(2)*q(4)-q(1)*q(3))                  *a(2)&
            + (-q(1)*q(1)-q(2)*q(2)+q(3)*q(3)+q(4)*q(4)) *a(3)
    end function body_to_space

    function space_to_body(a,q)
        implicit none
        double precision a(3),q(4),space_to_body(3)
        !
        space_to_body(1) = &
            + (-q(1)*q(1)+q(2)*q(2)-q(3)*q(3)+q(4)*q(4)) *a(1)&
            + 2d0*(q(3)*q(4)-q(1)*q(2))                  *a(2)&
            + 2d0*(q(2)*q(3)+q(1)*q(4))                  *a(3)
        space_to_body(2) = &
            - 2d0*(q(1)*q(2)+q(3)*q(4))                  *a(1)&
            + (+q(1)*q(1)-q(2)*q(2)-q(3)*q(3)+q(4)*q(4)) *a(2)&
            + 2d0*(q(2)*q(4)-q(1)*q(3))                  *a(3)
        space_to_body(3) = &
            + 2d0*(q(2)*q(3)-q(1)*q(4))                  *a(1)&
            - 2d0*(q(1)*q(3)+q(2)*q(4))                  *a(2)&
            + (-q(1)*q(1)-q(2)*q(2)+q(3)*q(3)+q(4)*q(4)) *a(3)
    end function space_to_body

    function cross_product(a, b)
        implicit none
        double precision a(3),b(3), cross_product(3)
        !
        cross_product(1) = a(2)*b(3) - a(3)*b(2)
        cross_product(2) = a(3)*b(1) - a(1)*b(3)
        cross_product(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product
end module fix_rigid_rotate_mod

