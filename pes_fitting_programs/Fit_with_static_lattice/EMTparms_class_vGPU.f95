module EMTparms_class

    use cudafor

    implicit none

    real(8), private, parameter                 :: sqrt_2  = 1.41421356237
    real(8), private, parameter                 :: sqrt_3  = 1.73205080757
    real(8), private, parameter                 :: isqrt_2 = 0.70710678118
    real(8), private, parameter                 :: pi      = 3.14159265359
    real(8), private, parameter                 :: beta    = 1.8093997906
    real(8), private, parameter                 :: twelveth= 0.0833333333333333
    integer, private, dimension(3), parameter   :: b       = (/12, 6, 24/)

    ! storage for variables passed into emt_init that are needed by emt energy
    real(8), private                :: cell(3)      ! dimensions of cell in x,y and z
    integer, private                :: n_l          ! number of lattice atoms
    real(8), private,allocatable    :: r0_lat(:,:)  ! positions of lattice atoms
    real(8), private                :: Eref         ! reference energy

    type EMTparms
        character(2)::  name    = 'RV'

        real(8)       ::  eta2    = 0.d0   ! A^-1
        real(8)       ::  kappa   = 0.d0   ! A^-1
        real(8)       ::  lambda  = 0.d0   ! A^-1

        real(8)       ::  E0      = 0.d0   ! eV
        real(8)       ::  n0      = 0.d0   ! A^-3
        real(8)       ::  s0      = 0.d0   ! A
        real(8)       ::  V0      = 0.d0   ! eV

    end type EMTparms

contains

subroutine emt_init (cell_in, n_l_in, r0_lat_in, pars_p, pars_l, energy)

    use cudafor
!
! Purpose:
!           Here, the fitting procedure is just implemented and the reference
!           energy calculated.
!
implicit none

!------------------------------------------------------------------------------
!                                   PREAMBLE
!                                   ========
!------------------------------------------------------------------------------

! declare variables and parameters that are passed down from the program

    real(8), dimension(3), intent(in)   :: cell_in         ! dimensions of cell in x,y and z
    integer, intent(in)                 :: n_l_in          ! number of lattice atoms
    real(8), dimension(:,:), intent(in) :: r0_lat_in       ! positions of lattice atoms
    type(EMTparms), intent(inout)       :: pars_p       ! parameters of particle
    type(EMTparms), intent(inout)       :: pars_l       ! parameters of lattice atoms
    real(8), intent(out)                :: energy       ! calculated reference energy

! declare the variables that appear in the subroutine

    integer :: i,j, k                   ! running parameter
    real(8) :: r                        ! distance
    real(8) :: rcut, rr, acut           ! values to calculate cut-off
    real(8) :: igamma1l, igamma2l       ! inverse gamma for lattice atoms
    real(8) :: theta                    ! variable for cut-off calculation
    real(8) :: chilp                    ! mixing between lattice (l) and particle (p)
    real(8), dimension(n_l_in) :: sigma_ll ! contribution to ns, l only
    real(8), dimension(n_l_in) :: s_l      ! neutral sphere radius lattice atoms
    real(8) :: V_ll                     ! Pair potential contributions
    real(8) :: vref_l                   ! reference pair pot. contrib.
    real(8) :: Ecoh                     ! cohesive energy of part & lattice atoms
    real(8), dimension(3) :: xl         ! for cal. cut-off
    real(8), dimension(3) :: rnnl       ! next neighbour distance for cut-off
    real(8) :: betas0_l                 ! beta * s0= for l and p
    real(8) :: kappadbeta_l             ! beta * kappa for l and p
    real(8), dimension(3) :: r3temp     ! temporary array variable
    real(8) :: rtemp                    ! temporary variable


!-----------------------Save inputs in module for use by emt ------------------
    cell   = cell_in
    n_l    = n_l_in
    if (.not. allocated(r0_lat)) then
        allocate(r0_lat(3,n_l))
        print *, 'allocated r0_lat for', n_l, 'atoms'
    else
        print *, 'r0_lat is all ready allocated'
    end if
    r0_lat = r0_lat_in


!----------------------------VALUES OF FREQUENT USE----------------------------
! definition of a few values that appear frequently in calculation
    betas0_l = beta * pars_l%s0
    kappadbeta_l = pars_l%kappa / beta

! 'coupling' parameters between p and l
    chilp = pars_p%n0 / pars_l%n0


!------------------------------------------------------------------------------
!                                  CUT-OFF
!                                  =======
!------------------------------------------------------------------------------
! We use the distance to the next-next nearest neighbours as cut-off.
! FOR FUTURE REVISION:
!            cut-off should be defined via lattice constant _AND_ changeable.

    rcut = betas0_l * sqrt_3
    rr = 4 * rcut / (sqrt_3 + 2)
    acut = 9.21024/(rr -rcut) ! ln(10000)

! Distances to the considered neighbours
    rnnl(1) = betas0_l
    rnnl(2) = rnnl(1) * sqrt_2
    rnnl(3) = rnnl(1) * sqrt_3

    xl = b * twelveth / (1 + exp(acut*(rnnl-rcut)))

!-----------------------------------GAMMA--------------------------------------
! Gamma enforces the cut-off together with theta (see below)
! Gamma is defined as inverse.
    r3temp = rnnl-betas0_l
    igamma1l = 1.0 / sum(xl*exp(-pars_l%eta2 * r3temp))
    igamma2l = 1.0 /sum(xl*exp(-kappadbeta_l * r3temp))


!------------------------------------------------------------------------------
!                          INDIVIDUAL CONTRIBUTIONS
!                          ========================
!------------------------------------------------------------------------------
! The values for the sums are set to zero.

    sigma_ll = 0
    V_ll = 0

    do i = 1, n_l
        do j = i+1, n_l

        !-----------------PERIODIC BOUNDARY CONDITIONS LATTICE-----------------
        ! Because we want them.

            r3temp(1) = r0_lat(1,i)-r0_lat(1,j)
            r3temp(2) = r0_lat(2,i)-r0_lat(2,j)
            r3temp(1) = r3temp(1) - (cell(1)*ANINT(r3temp(1)/cell(1)))
            r3temp(2) = r3temp(2) - (cell(2)*ANINT(r3temp(2)/cell(2)))
            r3temp(3) = r0_lat(3,i)-r0_lat(3,j)
            r =  sqrt(sum(r3temp**2))


        !---------------------------THETA LATTICE------------------------------
        ! Theta enforces the cut-off together with gamma (see above). This
        ! function enacts cutoff by reducing contributions of atoms outside the
        ! cut-off to zero.

            theta = 1.0 / (1 + exp( acut * (r - rcut) ) )


        !----------------------------SIGMA LATTICE-----------------------------
        ! Sigma is a contribution to the neutral sphere radius.
        ! It is a list in which for each lattice atom, the contributions of the
        ! others are summed up. To enforce the cut-off, it will be later
        ! corrected by gamma.
        ! sigma_pp does not exist because there is only a single particle.

            rtemp = theta*exp(-pars_l%eta2 * (r - betas0_l) )
            sigma_ll(i) = sigma_ll(i) + rtemp
            sigma_ll(j) = sigma_ll(j) + rtemp



        !-----------------------PAIR POTENTIAL LATTICE-------------------------
        ! For the lattice only.
        ! Will later be subjected to gamma to complete the cut-off.
        ! The particle has no pair potential contribution since there is only
        ! one and thus does not have a partner to interact with.

            rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
            V_ll = V_ll + rtemp
        end do
    end do

!-------------------------------CUT-OFF ENACTION-------------------------------
! Don't forget the gamma!

    sigma_ll = sigma_ll * igamma1l
    V_ll = V_ll * pars_l%V0 * igamma2l

!-----------------------------NEUTRAL SPHERE RADIUS----------------------------
! The neutral sphere radius is the radius in which the entire density of the
! atom is included.

    s_l = -log( sigma_ll * twelveth ) &
            / ( beta * pars_l%eta2)

!----------------MIXED REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------
! These contributions have to be subtracted to account for the contributions
! that were included twice.

    vref_l = 12 * pars_l%V0 * sum( exp( -pars_l%kappa * s_l) )

!------------------------------------------------------------------------------
!                           CALCULATING THE ENERGY
!                           ======================
!------------------------------------------------------------------------------


!---------------------------COHESIVE ENERGY FUNCTION---------------------------
! Calculates and sums the contributions to the cohesive energy for both lattice
! and particle.

    Ecoh = sum( (1 + pars_l%lambda*s_l) * exp(-pars_l%lambda * s_l)-1 ) &
          * pars_l%E0

!-------------------------------OVERALL ENERGY---------------------------------
! Summation over all contributions.

   energy = Ecoh - V_ll + 0.5 * vref_l
   Eref   = energy

end subroutine emt_init



!subroutine emt (cell, n_l, r0_lat, r_part, pars_p, pars_l, energy)
subroutine emt (r_part, pars_p, pars_l, energy)
!
! Purpose:
!       emt calculates the energy according to the effective medium theory.
! Input variables are (in oder of appearance):
!           particle: incident particle
!           lattice : lattice atoms
!           r0      : equilibrium positions of lattice atoms
!           pars_p  : emt-parameters of particle
!           pars_l  : emt-parameters of lattice
! Output variables are:
!           energy  : emt-energy
implicit none

!------------------------------------------------------------------------------
!                                   PREAMBLE
!                                   ========
!------------------------------------------------------------------------------

! declare variables and parameters that are passed down from the program

!    real(8), dimension(3), intent(in)   :: cell     ! dimensions of cell in x,y and z
!    integer, intent(in)                 :: n_l      ! number of lattice atoms
    real(8), dimension(3), intent (in)  :: r_part  ! position of the particle (note:
                                                    ! only one position for reference)
!    real(8), dimension(:,:), intent(in) :: r0_lat   ! positions of lattice atoms
    type(EMTparms), intent(inout)       :: pars_p   ! parameters of particle
    type(EMTparms), intent(inout)       :: pars_l   ! parameters of lattice atoms
    real(8), intent(out)                :: energy ! calc. reference energy

! declare the variables that appear in the subroutine

    integer :: i,j, k                   ! running parameter
    real(8) :: r                        ! distance
    real(8) :: rcut, rr, acut           ! values to calculate cut-off
    real(8) :: igamma1p, igamma2p       ! inverse gamma for particle
    real(8) :: igamma1l, igamma2l       ! inverse gamma for lattice atoms
    real(8) :: theta                    ! variable for cut-off calculation
    real(8) :: chilp, chipl             ! mixing between lattice (l) and particle (p)
    real(8) :: sigma_pl                 ! mixed contribution to neutral sphere,
    real(8) :: s_p                      ! neutral sphere radius particle
    real(8), dimension(n_l) :: sigma_ll ! contribution to ns, l only
    real(8), dimension(n_l) :: sigma_lp, tmp ! mixed contribution to ns
    real(8), dimension(n_l) :: s_l      ! neutral sphere radius lattice atoms
    real(8) :: V_pl, V_lp, V_ll         ! Pair potential contributions
    real(8) :: vref_l, vref_p           ! reference pair pot. contrib.
    real(8) :: Ecoh                     ! cohesive energy of part & lattice atoms
    real(8), dimension(3) :: xl, xp     ! for cal. cut-off
    real(8), dimension(3) :: rnnl, rnnp ! next neighbour distance for cut-off
    real(8) :: betas0_l, betas0_p       ! beta * s0= for l and p
    real(8) :: kappadbeta_l, kappadbeta_p ! beta * kappa for l and p
    real(8), dimension(3) :: r3temp     ! temporary array variable
    real(8) :: rtemp                    ! temporary variable
                                        ! device array
    real(8), device, allocatable, dimension(:,:) ::  d_r0_lat
    real(8), device, allocatable, dimension(:)   ::  d_sigma_ll, d_cell
    real(8), device, allocatable, dimension(:)   ::  d_sum_V_ll, d_sigma_lp, d_r_part
    real(8), device, allocatable, dimension(:)   ::  d_sum_sigma_pl, d_sum_V_lp, d_sum_V_pl
    real(8), dimension(n_l)                      ::  h_sum_sigma_pl, h_sum_V_lp, h_sum_V_pl, h_sum_V_ll




!----------------------VALUES OF FREQUENT USE ---------------------------------
! definition of a few values that appear frequently in calculation
    betas0_l = beta * pars_l%s0
    betas0_p = beta * pars_p%s0
    kappadbeta_l = pars_l%kappa / beta
    kappadbeta_p = pars_p%kappa / beta

! 'coupling' parameters between p and l
    chilp = pars_p%n0 / pars_l%n0
    chipl = 1.0 / chilp

!------------------------------------------------------------------------------
!                                  CUT-OFF
!                                  =======
!------------------------------------------------------------------------------
! We use the distance to the next-next nearest neighbours as cut-off.
! FOR FUTURE REVISION:
!            cut-off should be defined via lattice constant _AND_ changeable.

    rcut = betas0_l * sqrt_3
    rr = 4 * rcut / (sqrt_3 + 2)
    acut = 9.21024/(rr -rcut) ! ln(10000)

! Distances to the considered neighbours
    rnnl(1) = betas0_l
    rnnl(2) = rnnl(1) * sqrt_2
    rnnl(3) = rnnl(1) * sqrt_3
    rnnp(1) = betas0_p
    rnnp(2) = rnnp(1) * sqrt_2
    rnnp(3) = rnnp(1) * sqrt_3

    xl = b * twelveth / (1 + exp(acut*(rnnl-rcut)))
    xp = b * twelveth/ (1 + exp(acut*(rnnp-rcut)))

!-----------------------------------GAMMA--------------------------------------
! Gamma enforces the cut-off together with theta (see below)
! Gamma is defined as inverse.
    r3temp = rnnl-betas0_l
    igamma1l = 1.0 / sum(xl*exp(-pars_l%eta2 * r3temp))
    igamma2l = 1.0 /sum(xl*exp(-kappadbeta_l * r3temp))

    r3temp = rnnp-betas0_p
    igamma1p = 1.0 / sum(xp*exp(-pars_p%eta2 * r3temp))
    igamma2p = 1.0 / sum(xp*exp(-kappadbeta_p * r3temp))


!------------------------------------------------------------------------------
!                          INDIVIDUAL CONTRIBUTIONS
!                          ========================
!------------------------------------------------------------------------------


!---------------------------BEGINNING OF GPU REGION----------------------------

    ! Allocate device memory
    allocate ( d_r0_lat(3, n_l), d_sigma_ll(n_l), d_cell(3), &
               d_sum_V_ll(n_l), d_sigma_lp(n_l), d_r_part(3), &
               d_sum_sigma_pl(n_l), d_sum_V_lp(n_l), d_sum_V_pl(n_l) )

    ! Copy to device
    d_r0_lat   = r0_lat(1:3, 1:n_l)
    d_cell     = cell(1:3)
    d_r_part   = r_part

    call GPU_double_sum<<<1, n_l>>>( d_r0_lat, d_sigma_ll, d_cell, rcut, acut, &
                                 kappadbeta_l, betas0_l, n_l, pars_l%eta2, &
                                 d_sum_V_ll, d_r_part, d_sigma_lp,         &
                                 pars_p%eta2, d_sum_sigma_pl, d_sum_V_lp,  &
                                 d_sum_V_pl, betas0_p, kappadbeta_p)

    ! Copy back to host, sum it up and deallocate
    tmp      = d_sum_sigma_pl
    sigma_pl = sum( tmp )

    tmp      = d_sum_V_ll
    V_ll     = sum( tmp )

    tmp      = d_sum_V_lp
    V_lp     = sum( tmp )

    tmp      = d_sum_V_pl
    V_pl     = sum( tmp )

    sigma_ll(1:n_l) = d_sigma_ll
    sigma_lp(1:n_l) = d_sigma_lp

    deallocate ( d_r0_lat, d_sigma_ll, d_cell, &
                 d_sum_V_ll, d_sigma_lp, d_r_part, &
                 d_sum_sigma_pl, d_sum_V_lp, d_sum_V_pl )


!-------------------------------END OF GPU REGION------------------------------

!-------------------------------CUT-OFF ENACTION-------------------------------
! Don't forget the gamma!

    sigma_ll = sigma_ll * igamma1l
    sigma_lp = sigma_lp * igamma1l
    sigma_pl = sigma_pl * igamma1p

    V_ll = V_ll * pars_l%V0 * igamma2l
    V_lp = V_lp *chilp * pars_l%V0 * igamma2l
    V_pl = V_pl * pars_p%V0 * igamma2p * chipl



!-----------------------------NEUTRAL SPHERE RADIUS----------------------------
! The neutral sphere radius is the radius in which the entire density of the
! atom is included.

    s_l = -log( (sigma_ll + chilp * sigma_lp) * twelveth ) &
            / ( beta * pars_l%eta2)
    s_p  = -log( sigma_pl * chipl * twelveth) / ( beta * pars_p%eta2)


!----------------MIXED REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------
! These contributions have to be substracted to account for the contributions
! that were included twice.

    vref_l = 12 * pars_l%V0 * sum( exp( -pars_l%kappa * s_l) )
    vref_p = 12 * pars_p%V0 * exp( -pars_p%kappa * s_p)



!------------------------------------------------------------------------------
!                           CALCULATING THE ENERGY
!                           ======================
!------------------------------------------------------------------------------


!---------------------------COHESIVE ENERGY FUNCTION---------------------------
! Calculates and sums the contributions to the cohesive energy for both lattice
! and particle.

    Ecoh = sum( (1 + pars_l%lambda*s_l) * exp(-pars_l%lambda * s_l)-1 ) &
          * pars_l%E0 &
          + (1 + pars_p%lambda*s_p) * exp(-pars_p%lambda * s_p)* pars_p%E0



!-------------------------------OVERALL ENERGY---------------------------------
! Summation over all contributions.

    energy = Ecoh - V_ll - 0.5 * ( V_lp + V_pl - vref_l - vref_p) - Eref

end subroutine emt

subroutine emt_parms2array (emt_parms, array)
    type(EMTparms)          :: emt_parms
    real(8), dimension(7)   :: array

    array(1) = emt_parms%eta2
    array(2) = emt_parms%kappa
    array(3) = emt_parms%lambda
    array(4) = emt_parms%E0
    array(5) = emt_parms%n0
    array(6) = emt_parms%s0
    array(7) = emt_parms%V0
end subroutine emt_parms2array

subroutine array2emt_parms (array, emt_parms)
    type(EMTparms)          :: emt_parms
    real(8), dimension(7)   :: array

    emt_parms%eta2   = array(1)
    emt_parms%kappa  = array(2)
    emt_parms%lambda = array(3)
    emt_parms%E0     = array(4)
    emt_parms%n0     = array(5)
    emt_parms%s0     = array(6)
    emt_parms%V0     = array(7)
end subroutine array2emt_parms


!-------------------------------GPU SUBROUTINE---------------------------------
attributes(device) subroutine GPU_double_sum(d_r0_lat, d_sigma_ll, d_cell, rcut, acut,  &
                                         kappadbeta_l, betas0_l, n_l, pars_l_eta2,      &
                                         d_sum_V_ll, d_r_part, d_sigma_lp, pars_p_eta2, &
                                         d_sum_sigma_pl, d_sum_V_lp, d_sum_V_pl,        &
                                         betas0_p, kappadbeta_p)
    ! 'device' means: implicitly private to the module
    ! cannot handle derived type variables
    ! INPUT
    integer, value                          ::  n_l
    real(8), device, dimension(:,:)         ::  d_r0_lat
    real(8), device, dimension(:)           ::  d_sigma_ll
    real(8), device, dimension(:)           ::  d_sum_V_ll, d_sum_V_lp, d_sum_V_pl
    real(8), device, dimension(3)           ::  d_cell
    real(8), device, dimension(:)           ::  d_sigma_lp, d_sum_sigma_pl
    real(8), value                          ::  rcut, acut, kappadbeta_p
    real(8), value                          ::  kappadbeta_l, betas0_l, betas0_p
    real(8), value                          ::  pars_l_eta2, pars_p_eta2
    real(8), device, dimension(3)           ::  d_r_part

    ! LOCAL VARIABLES
    real(8), device, dimension(3)           ::  r3temp
    real(8), device                         ::  r, theta, rtemp
    integer, value                          ::  thid, j
    real(8), shared, dimension(630)         ::  sh_sigma_ll, sh_sum_V_ll! n_l = 630 needs to be passed
                                                                        ! explicitly to make sure it fits
                                                                        ! into shared memory (max. 48 kB)

    ! RUN, RUN !
                          ! Only works if n_l <= 1024 lattice atoms.
    thid = threadIdx%x    ! Otherwise blockDim needs to be used additionally
                          ! e. g. (blockidx%x-1)*blockdim%x + threadidx%x

    ! initialize shared memory

    sh_sigma_ll(thid) = 0.0
    sh_sum_V_ll(thid) = 0.0

    call syncthreads()

    do j = thid+1, n_l

        !-----------------PERIODIC BOUNDARY CONDITIONS LATTICE-----------------
        ! Because we want them.

        r3temp(1) = d_r0_lat(1,thid)-d_r0_lat(1,j)
        r3temp(2) = d_r0_lat(2,thid)-d_r0_lat(2,j)
        r3temp(1) = r3temp(1) - (d_cell(1)*ANINT(r3temp(1)/d_cell(1)))
        r3temp(2) = r3temp(2) - (d_cell(2)*ANINT(r3temp(2)/d_cell(2)))
        r3temp(3) = d_r0_lat(3,thid)-d_r0_lat(3,j)
        r =  sqrt(sum(r3temp**2))


        !---------------------------THETA LATTICE------------------------------
        ! Theta enforces the cut-off together with gamma (see above). This
        ! function enacts cutoff by reducing contributions of atoms outside the
        ! cut-off to zero.

        theta = 1.0 / (1 + exp( acut * (r - rcut) ) )


        !----------------------------SIGMA LATTICE-----------------------------
        ! Sigma is a contribution to the neutral sphere radius.
        ! It is a list in which for each lattice atom, the contributions of the
        ! others are summed up. To enforce the cut-off, it will be later
        ! corrected by gamma.
        ! sigma_pp does not exist because there is only a single particle.

        rtemp = theta*exp(-pars_l_eta2 * (r - betas0_l) )

        sh_sigma_ll(thid) = sh_sigma_ll(thid) + rtemp
        call syncthreads()
        sh_sigma_ll(j)    = sh_sigma_ll(j) + rtemp
        call syncthreads()
        ! no problem, since j is f(thid)


        !-----------------------PAIR POTENTIAL LATTICE-------------------------
        ! For the lattice only.
        ! Will later be subjected to gamma to complete the cut-off.
        ! The particle has no pair potential contribution since there is only
        ! one and thus does not have a partner to interact with.

        rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
        sh_sum_V_ll(thid) = sh_sum_V_ll(thid) + rtemp
        call syncthreads()


    end do

    !-----------------PERIODIC BOUNDERY CONDITIONS PARTICLE--------------------

    r3temp(1) = d_r0_lat(1,thid)-d_r_part(1)
    r3temp(2) = d_r0_lat(2,thid)-d_r_part(2)
    r3temp(1) = r3temp(1) - (d_cell(1)*ANINT(r3temp(1)/d_cell(1)))
    r3temp(2) = r3temp(2) - (d_cell(2)*ANINT(r3temp(2)/d_cell(2)))
    r3temp(3) = d_r0_lat(3,thid)-d_r_part(3)
    r =  sqrt(sum(r3temp**2))


    !----------------------------THETA PARTICLE--------------------------------

    theta = 1.0 / (1 + exp( acut * (r - rcut) ) )


    !-------------------------------MIXED SIGMA--------------------------------
    ! Contributions of both particle and lattice to neutral sphere radius
    ! To fully include the cut-off, we correct them later by gamma.

    d_sigma_lp(thid) = theta*exp(-pars_p_eta2 * (r - betas0_p) )
    rtemp = theta*exp(-pars_l_eta2 * (r - betas0_l) )
    d_sum_sigma_pl(thid) = rtemp


    !--------------------MIXED PAIR POTENTIAL CONTRIUBUTION--------------------

    rtemp = theta*exp(-kappadbeta_p * (r - betas0_p))
    d_sum_V_lp(thid) = rtemp
    rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
    d_sum_V_pl(thid) = rtemp

    ! copy back to device memory which can be accessed by host routine

    call syncthreads()

    if (thid .eq. 1) d_sigma_ll = sh_sigma_ll
    if (thid .eq. 2) d_sum_V_ll = sh_sum_V_ll

end subroutine inner_loop

end module

