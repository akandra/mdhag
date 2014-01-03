program mdhag
    ! Purpose:
    !       Do molecular dynamics calculations with the EMT potential.
    !
    ! Date          Author          History of Modification
    ! ====          ======          =======================
    ! 01.10.2013    Sascha&Svenja   Original
    !
    !
    use atom_class
    use md_init
    use force
!    use mdalgo
!    use open_file

    implicit none

    type(atoms) :: slab, teilchen   ! hold r, v and f for atoms in the box

    real(8) :: delta = 0.001d0      ! dr at numerical calculation of the forces
    real(8) :: E_ref, rtemp

!    real(8), dimension(:,:), allocatable :: rr,vv,aa,aao, aauo,ff,vp,vc, ac_temp
!    real(8),dimension(:), allocatable    :: imass
!    integer :: i,j, iter=3, q, bl, itraj, tw, s
!    real(8) :: Epot, rtemp, rdummy
!    real(8) :: zeta, norm
!    character(len=100)  :: filename
!    integer             :: randk = 13
!    real(8), dimension(2) :: xy_p
!    real(8),dimension(:), allocatable :: pe, ke_p, ke_l
!    real(8),dimension(:,:), allocatable :: rt, vt


    ! Call up procedure that gets us geometries
    call simbox_init(slab, teilchen)

    call pes(slab, teilchen, E_ref)

!    slab%r(3,1) = slab%r(3,1) + delta
    teilchen%r(3,1) = teilchen%r(3,1) + delta
    call pes(slab, teilchen, rtemp)
!    slab%r(3,1) = slab%r(3,1) - 2*delta
    teilchen%r(3,1) = teilchen%r(3,1) - 2*delta
    call pes(slab, teilchen, rtemp)
!    slab%r(3,1) = slab%r(3,1) + delta
    teilchen%r(3,1) = teilchen%r(3,1) + delta

!    print *, E_ref, rtemp, (rtemp - E_ref)/delta
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    !Call procedure that calcuates reference energy?
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    ! REMINDER:
!    ! Print out what kind of potential has been chosen, the parameters and the
!    ! initial conditions. DON'T FORGET!
!
!    bl= spec_l%n+spec_p%n
!
!    allocate(rr(3,bl),vv(3,bl),vp(3,bl))
!    allocate(aa(3,bl),aao(3,bl),ff(3,bl))
!    allocate(vc(3,bl),aauo(3,bl),ac_temp(3,bl),imass(bl))
!
!    tw = nsteps/wstep +1
!    allocate(pe(tw),ke_l(tw),ke_p(tw),rt(3,tw),vt(3,tw))
!
!    do itraj = start_tr, ntrajs+start_tr-1
!
!        print *, "Trajectory No. ", itraj
!        print *, "----------------------------------------------"
!
!        rtemp=1.0d0/spec_p%mass
!        do i=1,spec_p%n
!            rr(:,i)     = teilchen(i)%r
!            vv(:,i)     = teilchen(i)%v
!            aa(:,i)     = teilchen(i)%a
!            aao(:,i)    = teilchen(i)%aalt
!            aauo(:,i)   = teilchen(i)%auralt
!            ff(:,i)     = teilchen(i)%f
!            imass(i)   = rtemp
!        enddo
!
!        rtemp=1.0d0/spec_l%mass
!        do i=spec_p%n+1, bl
!            j=i-spec_p%n
!            rr(:,i)     = slab(j)%r
!            vv(:,i)     = slab(j)%v
!            aa(:,i)     = slab(j)%a
!            aao(:,i)    = slab(j)%aalt
!            aauo(:,i)   = slab(j)%auralt
!            ff(:,i)     = slab(j)%f
!            imass(i)    = rtemp
!        enddo
!
!        if (confname == 'POSCAR') then
!            !rr(1:2,1) = 0.0d0
!            !rr(3,1) = 6.0d0
!        else
!            ! Seed random number generator
!            do i=1,spec_p%n
!
!                call random_seed(size=randk)
!                call random_seed(put = itraj*randseed*i)
!
!                ! Initialize particle position
!                call random_number(xy_p)
!                xy_p =  matmul(celli(1:2,1:2),xy_p)
!                rr(1:2,i) = rr(1:2,i)+xy_p
!
!            enddo
!        end if
!
!        s=1
!
!        call pes(teilchen, slab, Epot)
!
!
!        ke_p(s) = (sum(vv(1,1:spec_p%n)**2+vv(2,1:spec_p%n)**2+vv(3,1:spec_p%n)**2))*0.5d0*spec_p%mass
!        ke_l(s) = (sum(vv(1,spec_p%n+1:bl)**2)+sum(vv(2,spec_p%n+1:bl)**2)+sum(vv(3,spec_p%n+1:bl)**2))*0.5d0*spec_l%mass
!        pe(s)   = Epot
!        rt(:,s) = rr(:,1)
!        vt(:,s) = vv(:,1)
!
!
!        do q = 1, nsteps
!            !print *, q
!        !    call beeman_1(rr,vv,aa,aao)
!        !    call  predict(vv,vp,aa,aao)
!            rr(1,:) = rr(1,:) + vv(1,:)*step + 0.5d0*step**2*aa(1,:)*imass
!            rr(2,:) = rr(2,:) + vv(2,:)*step + 0.5d0*step**2*aa(2,:)*imass
!            rr(3,:) = rr(3,:) + vv(3,:)*step + 0.5d0*step**2*aa(3,:)*imass
!
!
!            do i=1,spec_p%n
!                !write(*,'(3f15.5)') rr(:,i)
!                teilchen(i)%r = rr(:,i)
!            enddo
!
!            do i=spec_p%n+1, bl
!                j=i-spec_p%n
!                slab(j)%r = rr(:,i)
!            enddo
!
!            call pes(teilchen, slab, Epot)
!
!            do j = 1, 3
!            do i = 1, spec_l%n
!                    slab(i)%r(j) = rr(j,i+spec_p%n)-delta
!                    call pes(teilchen, slab, rtemp)
!
!                    ff(j,i+spec_p%n) = -(Epot-rtemp)/delta
!                    slab(i)%r(j) = rr(j,i+spec_p%n)
!            enddo
!            if (confname == 'POSCAR') then
!                ff(:,1:spec_p%n) = 0.0d0
!            else
!                do i = 1, spec_p%n
!                    teilchen(i)%r(j) = rr(j,i)-delta
!                    call pes(teilchen, slab, rtemp)
!                    ff(j,i) = -((Epot-rtemp)/delta)
!                    teilchen(i)%r(j) = rr(j,i)
!                enddo
!            end if
!
!            enddo
!
!
!            do j=1,3
!                vv(j,:) = vv(j,:) + 0.5d0*step*(ff(j,:) + aa(j,:))*imass
!            enddo
!
!            vv(:,bl-35:bl)= 0.0d0  ! Nail down the lower layer
!
!            aa = ff
!
!        !    call newton(ff,aa,imass)
!
!        !    if (spec_l%fric > 0) then
!        !        call fric_coef(rr,spec_l%fric,zeta)
!        !        do i = 1, iter
!        !            ac_temp = aa
!        !            call fric(ac_temp,vp,zeta)
!        !            call beeman_2(vv,vc,ac_temp,aao,aauo)
!        !            call norm_dist(vp,vc,norm)
!        !            if (norm < 1.0e-007) then
!        !                aa = ac_temp
!        !                exit
!        !            else
!        !                vp = vc
!        !            end if
!        !        end do
!        !    else
!        !        call beeman_2(vv,vc,aa,aao,aauo)
!        !    end if
!        !    vv=vc
!
!        !    if (mod(q,1) == 0) then
!        !    end if
!            if (mod(q,wstep) == 0) then
!                s = s+1
!                ke_p(s) = (vv(1,1)**2+vv(2,1)**2+vv(3,1)**2)*0.5d0*spec_p%mass
!                ke_l(s) = (sum(vv(1,2:bl)**2)+sum(vv(2,2:bl)**2)+sum(vv(3,2:bl)**2))*0.5d0*spec_l%mass
!                pe(s)   = Epot
!                rt(:,s) = rr(:,1)
!                vt(:,s) = vv(:,1)
!
!
!                if (confname == 'POSCAR') then
!                    write(*,'(4f12.3)') q*step, 2.0d0*ke_l(s)/(kB*3.0d0*(spec_l%n-36)),&
!                             ke_p(s)+ke_l(s)+pe(s)-(ke_p(1)+ke_l(1)+pe(1)), pe(s)
!                    if (q > 0/step) then
!                        write(filename,'(A,I7.7,A)') './traj191/conf',q,'.dat'
!                        call open_for_write(1,filename)
!                        write(1,'(3f18.8)') celli
!                        write(1,'(2I5)') spec_l%n,spec_l%n
!                        write(1,'(3f18.8)') rr(:,spec_p%n:bl)
!                        write(1,'(3f18.8)') vv(:,spec_p%n:bl)
!                        write(1,'(f18.8)') Epot
!                       print*,'print the configuration to', filename
!                    end if
!                end if
!
!
!            end if
!
!
!            if (confname .ne. 'POSCAR') then
!                if (rr(3,1) > 6.1 .or. rr(3,1) < -8.0) exit
!            end if
!
!        end do ! Cycle over time-steps
!!----------------------------------------------------------------------------------------------
!            if (confname == 'POSCAR') then
!                print *, 'fin'
!
!            else
!                write(filename,'(A,I7.7,A)') './trajs/traj',itraj,'.dat'
!                call open_for_write(1,filename)
!                write(1,'(10A15)') 'time(fs)','Ekin_p(eV)','Ekin_l(eV)','Epot(eV)',&
!                                    'x_p(A)','y_p(A)','z_p(A)',&
!                                    'vx_p(A)','vy_p(A)','vz_p(A)'
!                do q = 1,s
!                    write(1,'(14f15.5)') (q-1)*wstep*step, ke_p(q), ke_l(q), pe(q), rt(:,q), vt(:,q)
!                end do
!                close(1)
!            end if
!
!    end do ! Cycle over trajectories
!!--------------------------------------------------------------------------------------------
!
!deallocate(pe,ke_l,ke_p,rt,vt)

end program mdhag
