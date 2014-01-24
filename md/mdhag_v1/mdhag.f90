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
use output

implicit none

integer :: i, j, k, itraj, q, nwrites, ndata

real(8) :: imass_l, imass_p, rtemp
real(8), dimension(:,:), allocatable :: output_info
real(8), dimension(:,:), allocatable :: rmin_p              ! lowest particle position
integer, dimension(:), allocatable   :: col_start, col_end  ! collision time
logical :: exit_key

!timing
!real(8) :: start, fin

type(atoms) :: slab, teil   ! hold r, v and f for atoms in the box

! Construct simulation block and initialize everything
call simbox_init(slab, teil)
if (confname == 'mxt') call traj_init(slab, teil)

imass_l = 1.0d0/mass_l  ! inverse masses
imass_p = 1.0d0/mass_p

allocate(rmin_p(3,teil%n_atoms))
allocate(col_start(teil%n_atoms),col_end(teil%n_atoms))
nwrites=nsteps/wstep(2)
ndata= 4 + 7*teil%n_atoms ! Epot, Ekinl, Ekinp, density, r, v
allocate(output_info(ndata,nwrites))

! This loop runs over the number of trajectories that are going to be calculated.
! We still need an option to read in different geometries
do itraj = start_tr, ntrajs+start_tr-1


!------------------------------------------------------------------------------
!
!                       TRAJECTORY INITIALISATION
!
!------------------------------------------------------------------------------

!    call random_seed(put=randseed*itraj)
    do i = 1,itraj
        call random_number(rtemp)
    end do

    exit_key  = .false.
    overwrite = .true.
    rmin_p    = 6.1d0
    col_start = 0
    col_end   = 0
    nwrites   = 0
    ndata     = 0

    ! Continue run from given configuration
    if (confname == 'mxt') then
        call traj_init(slab, teil)
        if (teil%n_atoms > 0) call particle_init(teil)
    end if
    if (confname == 'POSCAR') then

        if (teil%n_atoms > 0) then
            call emt(slab, teil)
        else
            call emt1(slab)
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

    end if

    ! Assign projectile positions and velocities
!    if (n_p > 0) then
!
!        if (confname .ne. 'POSCAR') then
!            vinc = sqrt(2.0d0*einc/mass_p)
!            teil%v(1,:) =  vinc*sin(inclination)*cos(azimuth)
!            teil%v(2,:) =  vinc*sin(inclination)*sin(azimuth)
!            teil%v(3,:) = -vinc*cos(inclination)
!        end if
!    end if


    ! initial state
    if (wstep(1)==-1) call out_short(slab, teil, Epot, itraj, 0, rmin_p, col_end)

    print *, 'ntraj', itraj
    !timing
    !call cpu_time(start)

    ! loop over timesteps
    do q = 1, nsteps

!------------------------------------------------------------------------------
!
!                         PROPAGATION ROUTINE
!
!------------------------------------------------------------------------------
        call propagator_1(slab, md_algo_l, imass_l)         ! slab kick-drift

        if (teil%n_atoms > 0) then
           call propagator_1(teil, md_algo_p, imass_p)     ! projectile kick-drift
           call emt(slab, teil)                            ! slab-projectile forces
           call propagator_2(teil, md_algo_p, imass_p)     ! projectile kick
        else
            call emt1(slab)                           ! slab forces
        end if

        call propagator_2(slab, md_algo_l, imass_l)         ! slab kick

!------------------------------------------------------------------------------
!
!                      PROPERTIES CALCULATION ROUTINE
!
!------------------------------------------------------------------------------

        ! Lowest projectile positions
        if (teil%n_atoms > 0) then
            do i = 1,teil%n_atoms

                if (teil%r(3,i) < rmin_p(3,i)) rmin_p(:,i) = teil%r(:,i)

                if (teil%r(3,i) < 2.0d0) then
                    if (col_start(i) == 0) col_start(i) = q
                    col_end(i) = q
                end if

            end do
        end if


!------------------------------------------------------------------------------
!
!                            OUTPUT ROUTINE
!
!------------------------------------------------------------------------------

        if (mod(q,wstep(2))==0) then

            ! write out
            select case(wstep(1))

                case(0) ! save trajectory info
                        ! Epot, Ekinl, Ekinp, Etotal, density, r, v
                    nwrites = nwrites+1
                    output_info(1,nwrites) = Epot
                    output_info(2,nwrites) = E_kin(slab,mass_l)
                    output_info(3,nwrites) = E_kin(teil,mass_p)
                    output_info(4,nwrites) = Epot + output_info(2,nwrites)&
                                                  + output_info(3,nwrites)
                    j = teil%n_atoms
                    output_info(5    :4+  j,nwrites) = teil%dens
                    output_info(5+  j:4+4*j,nwrites) = reshape(teil%r,(/3*j/))
                    output_info(5+4*j:4+7*j,nwrites) = reshape(teil%v,(/3*j/))
                    ndata = ndata + 1

                case default ! full configuration of system
                    call full_conf(slab, teil,itraj)

            end select

        end if

        do i=1,teil%n_atoms
            if (teil%r(3,i) > 6.1 .or. teil%r(3,i) < -8.1) exit_key = .true.
        end do
        if (exit_key) exit

    end do ! steps

    col_end = col_end - col_start
    ! final state
    if (wstep(1)==-1) call out_short (slab, teil, Epot, itraj, q, rmin_p, col_end)
    if (wstep(1)== 0) call out_detail(output_info, ndata, itraj)


    !timing
    !call cpu_time(fin)
    !print *, fin - start, " seconds"
end do ! trajectories

deallocate(output_info)
deallocate(col_end, col_start)
deallocate(rmin_p, pars_l, pars_p)

end program mdhag
