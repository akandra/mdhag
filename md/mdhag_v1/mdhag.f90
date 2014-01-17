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
use mdalgo

implicit none

integer :: i, j, itraj, q

real(8) :: imass_l, imass_p
real(8) :: Epot, Ekin_l, Ekin_p

!timing
!real(8) :: start, fin

type(atoms) :: slab, teil   ! hold r, v and f for atoms in the box

call open_for_write(10, 'dat/e_lan.dat')
call open_for_write(11, 'dat/pos_lan.dat')

! Construct simulation block and initialize everything
call simbox_init(slab, teil)

imass_l = 1.0d0/mass_l
imass_p = 1.0d0/mass_p

! This loop runs over the number of trajectories that are going to be calculated.
! We still need an option to read in different geometries
do itraj = start_tr, ntrajs+start_tr-1

    print *, "Trajectory No. ", itraj
    print *, "----------------------------------------------"

    if (teil%n_atoms > 0) then
        call emt(slab, teil, Epot)
    else
        call emt1(slab, Epot)
    end if

    call ldfa(slab)
    slab%a  = slab%f*imass_l
    slab%ao = slab%a
    slab%au = slab%ao
    if (teil%n_atoms > 0) then
       call ldfa(teil)
        teil%a  = teil%f*imass_p
        teil%ao = teil%a
        teil%au = teil%ao
    end if

    Ekin_l = 0.5d0*mass_l*sum(slab%v*slab%v)
    Ekin_p = 0.5d0*mass_p*sum(teil%v*teil%v)
    print '(i8,4f12.5)', 0, Ekin_l, Ekin_p, Ekin_l + Ekin_p + Epot, teil%r(3,1)
    write (10,'(i8,6f12.5)') 0, 2*Ekin_l/(3*slab%n_atoms*kB), 2*Ekin_p/(3*teil%n_atoms*kB), &
                  0.5d0*sum(teil%v(:,1)*teil%v(:,1))*mass_p, Ekin_l + Ekin_p + Epot, &
                  teil%r(3,1), teil%dens(1)*hbar


    write (11,'(i8)') 0
    write (11,'(3f12.5)') teil%r
    write (11,'(3f12.5)') slab%r

    !timing
    !call cpu_time(start)

    ! loop over timesteps
    do q = 1, nsteps

        call propagator_1(slab, md_algo_l, imass_l)         ! slab kick-drift

        if (teil%n_atoms > 0) then
            call propagator_1(teil, md_algo_p, imass_p)     ! projectile kick-drift
            call emt(slab, teil, Epot)                      ! slab-projectile forces
            call propagator_2(teil, md_algo_p, imass_p)     ! projectile kick
        else
            call emt1(slab, Epot)                           ! slab forces
        end if

        call propagator_2(slab, md_algo_l, imass_l)         ! slab kick


        if (mod(q,wstep)==0) then

            Ekin_l = 0.5d0*sum(slab%v*slab%v)*mass_l

            if (teil%n_atoms > 0) then
                Ekin_p = 0.5d0*sum(teil%v*teil%v)*mass_p
                write (10,'(i7,6f12.5)') q, 2*Ekin_l/(3*slab%n_atoms*kB),&
                                            2*Ekin_p/(3*teil%n_atoms*kB),&
                               0.5d0*sum(teil%v(:,1)*teil%v(:,1))*mass_p,&
                                                  Ekin_l + Ekin_p + Epot,&
                                                             teil%r(3,1)
                write (11,'(i8)') q
                write (11,'(3f12.5)') teil%r
                write (11,'(3f12.5)') slab%r
                print *, q

!                if (teil%r(3,1) > 6.1 .or. teil%r(3,1) < -8.0) exit
            else
                write (10,'(i8,2f12.5)') q, 2*Ekin_l/(3*slab%n_atoms*kB),&
                                                           Ekin_l + Epot
                write (11,'(i8)') q
                write (11,'(3f12.5)') slab%r
                print *, q
            end if
        end if



    end do ! steps

    !timing
    !call cpu_time(fin)
    !print *, fin - start, " seconds"

end do ! trajectories

close(10)
close(11)

end program mdhag
