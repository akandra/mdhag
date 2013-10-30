program mdhag
    ! Purpose:
    !       Do molecular dynamics calculations with the EMT potential.
    !
    ! Date          Author          History of Modification
    ! ====          ======          =======================
    ! 01.10.2013    Sascha&Svenja   Original
    !
    !
!    use atom_class
    use md_init
    use force
    use mdalgo
    use open_file

    implicit none

    type(atom), dimension(:), allocatable :: slab, teilchen

    real(8), dimension(:,:), allocatable :: rr,vv,aa,aao, aauo,ff,vp,vc, ac_temp
    real(8),dimension(:), allocatable    :: imass
    integer :: i,j, iter=3, q, bl
    real(8) :: Epot, rtemp, rdummy, Ekin(3)
    real(8) :: zeta, norm, delta


    call open_for_write(1,'testoutput2.dat')
    ! Call up procedure that gets us geometries
    call simbox_init(teilchen,slab)

    ! REMINDER:
    ! Print out what kind of potential has been chosen, the parameters and the
    ! initial conditions. DON'T FORGET!

 !   call pes(teilchen, slab, Epot)
    bl= spec_l%n+spec_p%n
    !teilchen%v(1) = 0.01
    !teilchen%v(2) = 0.02
    !teilchen%v(3) = -0.03


    allocate(rr(3,bl),vv(3,bl),vp(3,bl))
    allocate(aa(3,bl),aao(3,bl),ff(3,bl))
    allocate(vc(3,bl),aauo(3,bl),ac_temp(3,bl),imass(bl))
    rtemp=1.0d0/spec_p%mass
    do i=1,spec_p%n
        rr(:,i)     = teilchen(i)%r
        vv(:,i)     = teilchen(i)%v
        aa(:,i)     = teilchen(i)%a
        aao(:,i)    = teilchen(i)%aalt
        aauo(:,i)   = teilchen(i)%auralt
        ff(:,i)     = teilchen(i)%f
        imass(i)   = rtemp
    enddo
    rtemp=1.0d0/spec_l%mass
    do i=spec_p%n+1, bl
        j=i-spec_p%n
        rr(:,i)     = slab(j)%r
        vv(:,i)     = slab(j)%v
        aa(:,i)     = slab(j)%a
        aao(:,i)    = slab(j)%aalt
        aauo(:,i)   = slab(j)%auralt
        ff(:,i)     = slab(j)%f
        imass(i)    = rtemp
    enddo
    vv(3,1) = -0.3d0
    vv(3,2:bl-36) = 0.001d0
    !vv(1,2:bl-36) = 0.0003d0
    !vv(2,2:bl-36) = 0.0003d0
   ! vv(2,1) = 0.0042d0
    delta = 0.0010d0
    !vv(2,18) = 0.02d0
 !   vp = 0.0d0


do q = 1, nsteps

!    call beeman_1(rr,vv,aa,aao)
!    call  predict(vv,vp,aa,aao)

    do i=1,spec_p%n
        teilchen(i)%r = rr(:,i)
    enddo
    do i=spec_p%n+1, bl
        j=i-spec_p%n
        slab(j)%r = rr(:,i)
    enddo

    call pes(teilchen, slab, Epot)

    do j = 1, 3
    do i = 1, spec_l%n-36
        slab(i)%r(j) = rr(j,i+1)-delta
        call pes(teilchen, slab, rtemp)
        ff(j,i+1) = -(Epot-rtemp)/delta
        slab(i)%r(j) = rr(j,i+1)
    enddo
    do i = 1, spec_p%n
        teilchen(i)%r(j) = rr(j,i)-delta
        call pes(teilchen, slab, rtemp)
        ff(j,i) = -(Epot-rtemp)/delta
        teilchen(i)%r(j) = rr(j,i)
    enddo

    enddo

    rr(1,:) = rr(1,:) + vv(1,:)*step+step**2*imass*0.5d0*ff(1,:)
    rr(2,:) = rr(2,:) + vv(2,:)*step+step**2*imass*0.5d0*ff(2,:)
    rr(3,:) = rr(3,:) + vv(3,:)*step+step**2*imass*0.5d0*ff(3,:)

    vv(:,1) = vv(:,1) + ff(:,1)*step*imass(1)
    vv(:,2:bl-36) = vv(:,2:bl-36) + ff(:,2:bl-36)*step*imass(2)
    vv(:,bl-35:bl) =vv(:,bl-35:bl) + ff(:,bl-35:bl)*step*imass(2)
    !print *, rr(1,2), vv(1,2), ff(1,2)


!    call newton(ff,aa,imass)

!    if (spec_l%fric > 0) then
!        call fric_coef(rr,spec_l%fric,zeta)
!        do i = 1, iter
!            ac_temp = aa
!            call fric(ac_temp,vp,zeta)
!            call beeman_2(vv,vc,ac_temp,aao,aauo)
!            call norm_dist(vp,vc,norm)
!            if (norm < 1.0e-007) then
!                aa = ac_temp
!                exit
!            else
!                vp = vc
!            end if
!        end do
!    else
!        call beeman_2(vv,vc,aa,aao,aauo)
!    end if
!    vv=vc

    print *, q
    if (mod(q,1) == 0) then
     write(1,'(3f15.5)') rr
     Ekin(2)=(sum(vv(1,2:bl-36)**2)+sum(vv(2,2:bl-36)**2)+sum(vv(3,2:bl-36)**2))*0.5d0/imass(2)
     Ekin(3)=(sum(vv(1,bl-35:bl)**2)+sum(vv(2,bl-35:bl)**2)+sum(vv(3,bl-35:bl)**2))*0.5d0/imass(2)
     Ekin(1) = (vv(1,1)**2+vv(2,1)**2+vv(3,1)**2)*0.5d0/imass(1)
    write(1,*) Epot, Ekin
    end if
    if (rr(3,1) > 6.1 .or. rr(3,1) < -8.0) exit
end do

 !   write(*,'(3f10.5)') ( (vv(:,i)-vp(:,i)), i=1, spec_l%n)

 !   write(*,'(3e15.5)') (slab(i)%f, i=1, spec_l%n)
 !   write(*,'(3e15.5)') teilchen(1)%f

close(1)
end program mdhag
