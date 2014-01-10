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

real(8) :: imass_l, imass_p, norm_l, norm_p
real(8) :: Epot, Ekin_l, Ekin_p
real(8) :: delta = 0.001

type(atoms) :: slab, teil   ! hold r, v and f for atoms in the box

call open_for_write(10, 'e_lan.dat')
call open_for_write(11, 'pos_lan.dat')

! Construct simulation block and initialize everything
call simbox_init(slab, teil)

imass_l = 1.0d0/mass_l
imass_p = 1.0d0/mass_p

! This loop runs over the number of trajectories that are going to be calculated.
! We still need an option to read in different geometries
do itraj = start_tr, ntrajs+start_tr-1

    print *, "Trajectory No. ", itraj
    print *, "----------------------------------------------"

    call emt(slab, teil, Epot)

    slab%a  = slab%f*imass_l
    teil%a  = teil%f*imass_p
    slab%ao = slab%a
    teil%ao = teil%a
    slab%au = slab%ao
    teil%au = teil%ao

    Ekin_l = 0.5d0*mass_l*sum(slab%v*slab%v)
    Ekin_p = 0.5d0*mass_p*sum(teil%v*teil%v)
    print '(i8,4f12.5)', 0, Ekin_l, Ekin_p, Ekin_l + Ekin_p + Epot, teil%r(3,1)
    write (10,'(i8,4f12.5)') 0, 2*Ekin_l/(3*slab%n_atoms*kB), &
                  0.5d0*sum(teil%v(:,1)*teil%v(:,1))*mass_p, Ekin_l + Ekin_p + Epot, teil%r(3,1)


    write (11,'(i8)') 0
    write (11,'(3f12.5)') teil%r
    write (11,'(3f12.5)') slab%r

    ! loop over timesteps
    do q = 1, nsteps

        select case (md_algo)

            case (0) ! velocity Verlet Algorithm
                ! 1st and 2nd steps
                call verlet_1(slab)
                call verlet_1(teil)
                ! velocity-independent forces
                call emt(slab, teil, Epot)
                ! velocity-dependent forces

                ! 3rd step: new accelerations
                call newton(slab, imass_l)
                call newton(teil, imass_p)
                ! 4th step: update velocities
                call verlet_2(slab)
                call verlet_2(teil)

            case (1) ! Refson-Beeman Algorithm
                ! 1st and 2nd steps
                call beeman_1(slab)
                call beeman_1(teil)
                ! velocity-independent forces
                call emt(slab, teil, Epot)
                ! numerical forces
!                do j=1, 3
!                do i=1, slab%n_atoms
!                    slab%r(j,i) = slab%r(j,i) - delta
!                    call emt_e(slab, teil, Epot)
!                    slab%f(j,i) = Epot/(2*delta)
!                    slab%r(j,i) = slab%r(j,i) + 2*delta
!                    call emt_e(slab, teil, Epot)
!                    slab%r(j,i) = slab%r(j,i) - delta
!                    slab%f(j,i) =  slab%f(j,i) - Epot/(2*delta)
!                end do
!                do i=1, teil%n_atoms
!                    teil%r(j,i) = teil%r(j,i) - delta
!                    call emt_e(teil, teil, Epot)
!                    teil%f(j,i) = Epot/(2*delta)
!                    teil%r(j,i) = teil%r(j,i) + 2*delta
!                    call emt_e(teil, teil, Epot)
!                    teil%r(j,i) = teil%r(j,i) - delta
!                    teil%f(j,i) = teil%f(j,i) - Epot/(2*delta)
!                end do
!                end do
!
                ! predictor-corrector cycle
                do
                    ! velocity-dependent forces

                    ! 3rd step: new accelerations
                    call newton(slab, imass_l)
                    call newton(teil, imass_p)
                    ! Step 4: corrected velocities
                    call beeman_2(slab)
                    call beeman_2(teil)

                    call norm_dist(slab%vp, slab%vc, 3*slab%n_atoms, norm_l)
                    call norm_dist(teil%vp, teil%vc, 3*teil%n_atoms, norm_p)

                    if ((norm_l + norm_p) < 1.0e-007) then

                        slab%v  = slab%vc
                        teil%v  = teil%vc
                        slab%au = slab%ao
                        slab%ao = slab%a
                        teil%au = teil%ao
                        teil%ao = teil%a
                        exit

                    else

                        slab%vp = slab%vc
                        teil%vp = teil%vc

                    end if

                end do

            case (2) ! Langevin
                ! 1st and 2nd steps
                call verlet_1(slab)
                call verlet_1(teil)
                ! velocity-independent forces
                call emt(slab, teil, Epot)
                ! velocity-dependent forces

                ! 3rd step: new accelerations
                call newton(slab, imass_l)
                call newton(teil, imass_p)
                ! 4th step: update velocities
                call verlet_2(slab)
                call verlet_2(teil)

        end select
        if (mod(q,wstep)==0) then

            Ekin_l = 0.5d0*sum(slab%v*slab%v)*mass_l
            Ekin_p = 0.5d0*sum(teil%v*teil%v)*mass_p
            write (10,'(i8,4f12.5)') q, 2*Ekin_l/(3*slab%n_atoms*kB), &
                    0.5d0*sum(teil%v(:,1)*teil%v(:,1))*mass_p, Ekin_l + Ekin_p + Epot, teil%r(3,1)
            write (11,'(i8)') q
            write (11,'(3f12.5)') teil%r
            write (11,'(3f12.5)') slab%r
            print *, q

        end if

    if (teil%r(3,1) > 6.1 .or. teil%r(3,1) < -8.0) exit

    end do ! steps
end do ! trajectories

close(10)
close(11)

end program mdhag
