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


type(atoms) :: slab, teil   ! hold r, v and f for atoms in the box

integer :: i, j, k, itraj, q, nwrites, ndata

real(8) :: imass_l, imass_p, rtemp, proj_downgone
real(8), dimension(:,:), allocatable :: output_info
real(8), dimension(:,:), allocatable :: rmin_p              ! lowest particle position
real(8), dimension(:,:), allocatable :: tock                ! collection array for momentum
integer, dimension(:), allocatable   :: col_start, col_end  ! collision time
integer, dimension(:), allocatable   :: imp, q_imp          ! collision number
logical :: exit_key, imp_switch
real(8), dimension(:), allocatable :: eed, eed_prec          ! embedded electron density for projectile

!timing
!real(8) :: start, fin


! Construct simulation block and initialize everything
call simbox_init(slab, teil)

! Inverse Mass
imass_l = 1.0d0/mass_l
imass_p = 1.0d0/mass_p

! Lowest permitted projectile position
proj_downgone = slab%r(3,slab%n_atoms) - 3.0d0

! Allocate useful arrays
allocate(rmin_p(3,teil%n_atoms))
allocate(col_start(teil%n_atoms),col_end(teil%n_atoms))
nwrites=nsteps/wstep(2)
ndata= 4 + 7*teil%n_atoms ! Epot, Ekinl, Ekinp, density, r, v
allocate(output_info(ndata,nwrites))
allocate(tock(3,teil%n_atoms))
allocate(imp(teil%n_atoms), q_imp(teil%n_atoms))
allocate(eed(teil%n_atoms),eed_prec(teil%n_atoms))

!------------------------------------------------------------------------------
!
!                      LOOP OVER TRAJECTORIES
!
!------------------------------------------------------------------------------
do itraj = start_tr, ntrajs+start_tr-1

!------------------- TRAJECTORY INITIALISATION ROUTINE ------------------------

    print *, 'Running trajectory number ... ', itraj

    call random_seed(put=randseed)  ! Seed random number generator
    do i = 1,100*(itraj - 1)
        call random_number(rtemp)   ! rotate it according to trajectory number
    end do

    exit_key    = .false.
    imp_switch  = .false.
    overwrite   = .true.
    rmin_p      = 6.10d0
    col_start   = 0
    col_end     = 0
    nwrites     = 0
    ndata       = 0
    imp         = 0
    q_imp       = 0
    eed_prec     = 0.0d0

    if (confname == 'mxt') call traj_init(slab, teil)

    if (teil%n_atoms > 0) then
        call particle_init(teil)
        call emt(slab, teil)
        if (md_algo_p == 3 .or. md_algo_p == 4 ) call ldfa(teil)
        teil%a  = teil%f*imass_p
        teil%ao = teil%a
        teil%au = teil%ao
    else
        call emt1(slab)
    end if

    if (md_algo_l == 3 .or. md_algo_l == 4 ) call ldfa(slab)
    slab%a  = slab%f*imass_l
    slab%ao = slab%a
    slab%au = slab%ao

    ! save initial state
    if (wstep(1)==-1) call out_short(slab, teil, Epot, itraj, 0, rmin_p, col_end, imp)
    tock = teil%v

    !print *, 'ntraj', itraj
    !timing
    !call cpu_time(start)

!------------------------------------------------------------------------------
!
!                         LOOP OVER TIME STEPS
!
!------------------------------------------------------------------------------
    do q = 1, nsteps

!----------------------- PROPAGATION ROUTINE ----------------------------------

        call propagator_1(slab, md_algo_l, imass_l)         ! slab kick-drift

        if (teil%n_atoms > 0) then
           call propagator_1(teil, md_algo_p, imass_p)     ! projectile kick-drift
           call emt(slab, teil)                            ! slab-projectile forces
           eed = teil%dens                                 ! keep eed values
           call propagator_2(teil, md_algo_p, imass_p)     ! projectile kick
        else
            call emt1(slab)                           ! slab forces
        end if


        call propagator_2(slab, md_algo_l, imass_l)         ! slab kick


!-------------------- PROPERTIES CALCULATION ROUTINE --------------------------

        if (teil%n_atoms > 0) then
            do i = 1,teil%n_atoms

                ! Lowest projectile position
                if (teil%r(3,i) < rmin_p(3,i)) rmin_p(:,i) = teil%r(:,i)

                ! Collision time
                if (teil%r(3,i) < 2.0d0) then
                    if (col_start(i) == 0) col_start(i) = q
                    col_end(i) = q
                end if
                ! Collision number. Criterium:
                ! embedded projectile electron density has maximum
                if (eed(i) > eed_prec(i))  then
                    imp_switch = .true.
                else if (imp_switch .and. q > q_imp(i) + step_imp) then
                    imp(i) = imp(i) + 1
                    imp_switch = .false.
                    q_imp(i) = q
                end if
                eed_prec = eed

                ! Collision number. Criterium:
                ! an angle between previous and current velocity larger than 5 degree
!                rtemp = sqrt(sum(teil%v(:,i)**2)*sum(tock(:,i)**2))
!                rtemp = dot_product(teil%v(:,i),tock(:,i))/rtemp
!                if (rtemp < crit_imp .or. rtemp < 0.0d0 ) then
!                    if  (q > q_imp(i) + step_imp) then
!                        q_imp(i) = q
!                        imp(i)   = imp(i) + 1
!                        tock(:,i) = teil%v(:,i)
!                    end if
!                end if

            end do
        end if

!---------------------------- OUTPUT ROUTINE ----------------------------------

        if (mod(q,wstep(2))==0) then

            ! write out
            select case(wstep(1))

                case(-1)

                case(0) ! save trajectory info
                        ! Epot, Ekinl, Ekinp, Etotal, density, r, v
                    nwrites = nwrites+1
                    output_info(1,nwrites) = Epot
                    output_info(2,nwrites) = E_kin(slab,mass_l)
                    output_info(3,nwrites) = E_kin(teil,mass_p)
                    output_info(4,nwrites) = Epot + output_info(2,nwrites)&
                                                  + output_info(3,nwrites)
                    j = teil%n_atoms
                    output_info(5    :4+  j,nwrites) = eed
                    output_info(5+  j:4+4*j,nwrites) = reshape(teil%r,(/3*j/))
                    output_info(5+4*j:4+7*j,nwrites) = reshape(teil%v,(/3*j/))
                    ndata = ndata + 1

                case default ! full configuration of system
                    if (q > wstep(1)) call full_conf(slab, teil,itraj)

            end select

        end if


        do i=1,teil%n_atoms
            if (teil%r(3,i) > proj_upgone .or. teil%r(3,i) < proj_downgone) then
                exit_key = .true.
            end if
        end do
        if (exit_key) exit

    end do ! steps

    col_end = col_end - col_start
    ! final state
    if (wstep(1)==-1) call out_short (slab, teil, Epot, itraj, q, rmin_p, col_end, imp)
    if (wstep(1)== 0) call out_detail(output_info, ndata, itraj)


    !timing
    !call cpu_time(fin)
    !print *, fin - start, " seconds"

end do ! trajectories

deallocate(eed_prec, eed, imp, tock)
deallocate(output_info)
deallocate(col_end, col_start)
deallocate(rmin_p, pars_l, pars_p)

end program mdhag
