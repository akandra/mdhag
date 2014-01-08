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

integer :: i, itraj, q

real(8) :: imass_l, imass_p, norm_l, norm_p
real(8) :: Epot, Ekin_l, Ekin_p

type(atoms) :: slab, teil   ! hold r, v and f for atoms in the box

call open_for_write(10, 'e_beeman.dat')


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

    slab%ao = slab%f*imass_l
    teil%ao = teil%f*imass_p
    slab%au = slab%ao
    teil%au = teil%ao

    Ekin_l = 0.5d0*sum(slab%v*slab%v)*mass_l
    Ekin_p = 0.5d0*sum(teil%v*teil%v)*mass_p
    print '(i8,4f12.5)', 0, Ekin_l, Ekin_p, Ekin_l + Ekin_p + Epot, teil%r(3,1)

    ! loop over timesteps
    do q = 1, nsteps

        select case (md_algo)

            case (0)
                ! 1st and 2nd steps of velocity Verlet Algorithm
                call verlet_1(slab)
                call verlet_1(teil)

            case (1)
                ! 1st and 2nd steps of velocity Verlet Algorithm
                call beeman_1(slab)
                call beeman_1(teil)

        end select

        ! velocity-independent forces
        call emt(slab, teil, Epot)

        ! velocity-dependent forces

        ! New accelerations
        call newton(slab, imass_l)
        call newton(teil, imass_p)

        select case (md_algo)

            case (0)

                call verlet_2(slab)
                call verlet_2(teil)

            case (1)

                ! predictor-corrector cycle
                do

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

                        ! velocity-dependent forces

                        ! Step 3: new accelerations
                        call newton(slab, imass_l)
                        call newton(teil, imass_p)

                    end if

                end do

        end select

        if (mod(q,wstep)==0) then
            Ekin_l = 0.5d0*sum(slab%v*slab%v)*mass_l
            Ekin_p = 0.5d0*sum(teil%v*teil%v)*mass_p
            write (10,'(i8,4f12.5)'), q, 2*Ekin_l/(3*slab%n_atoms*kB), Ekin_p, Ekin_l + Ekin_p + Epot, teil%r(3,1)
        end if

    if (teil%r(3,1) > 6.1 .or. teil%r(3,1) < -8.0) exit

    end do ! steps
end do ! trajectories

close(10)

end program mdhag
