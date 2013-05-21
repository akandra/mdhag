module EMTparms_class

use atom_class
use math_functions

    real(8), parameter :: sqrt_2 = 1.41421356237
    real(8), parameter :: sqrt_3 = 1.73205080757
    real(8), parameter :: isqrt_2 = 0.70710678118
    real(8), parameter :: pi = 3.14159265359
    real(8), parameter :: beta = 1.8093997906
    real(8), parameter :: twelveth = 0.0833333333333333
    integer, dimension(3), parameter :: b = (/12, 6, 24/)

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

subroutine emt_init (cell, a_lat, n_l, r0_lat, pars_p, pars_l, energy)
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

    real(8), dimension(3), intent(in)   :: cell     ! dimensions of cell in x,y and z
    real(8),intent(in)                  :: a_lat    ! lattice constant of lattice
    integer, intent(in)                 :: n_l      ! number of lattice atoms
    real(8), dimension(:,:), intent(in) :: r0_lat   ! positions of lattice atoms
    type(EMTparms), intent(inout)       :: pars_p   ! parameters of particle
    type(EMTparms), intent(inout)       :: pars_l   ! parameters of lattice atoms
    real(8), intent(out)                :: energy    ! calc. reference energy

! declare the variables that appear in the subroutine

    integer :: i,j, k                   ! running parameter
    real(8) :: r                        ! distance
    real(8) :: rcut, rr, acut           ! values to calculate cut-off
    real(8) :: igamma1l, igamma2l       ! inverse gamma for lattice atoms
    real(8) :: theta                    ! variable for cut-off calculation
    real(8) :: chilp                    ! mixing between lattice (l) and particle (p)
    real(8), dimension(n_l) :: sigma_ll ! contribution to ns, l only
    real(8), dimension(n_l) :: s_l      ! neutral sphere radius lattice atoms
    real(8) :: V_ll                     ! Pair potential contributions
    real(8) :: vref_l                   ! reference pair pot. contrib.
    real(8) :: Ecoh                     ! cohesive energy of part & lattice atoms
    real(8), dimension(3) :: xl         ! for cal. cut-off
    real(8), dimension(3) :: rnnl       ! next neighbour distance for cut-off
    real(8) :: betas0_l                 ! beta * s0= for l and p
    real(8) :: kappadbeta_l             ! beta * kappa for l and p
    real(8), dimension(3) :: r3temp     ! temporary array variable
    real(8) :: rtemp                    ! temporary variable


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

    rcut = a_lat * sqrt_3 * isqrt_2
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
! These contributions have to be substracted to account for the contributions
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

! _____________________________________________________________________________

end subroutine emt_init



subroutine emt_fit (cell, a_lat, n_l, r_lat, r_part, pars_p, pars_l, energy)
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

    real(8), dimension(3), intent(in)   :: cell     ! dimensions of cell in x,y and z
    real(8), intent(in)                 :: a_lat   ! lattice constant of lattice
    integer, intent(in)                 :: n_l      ! number of lattice atoms
    real(8), dimension(3), intent (in)  :: r_part   ! positions of the particle
    real(8), dimension(:,:), intent(in) :: r_lat    ! positions of lattice atoms
    type(EMTparms), intent(inout)       :: pars_p   ! parameters of particle
    type(EMTparms), intent(inout)       :: pars_l   ! parameters of lattice atoms
    real(8), intent(out)                :: energy   ! calc. reference energy

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
    real(8), dimension(n_l) :: sigma_lp ! mixed contribution to ns
    real(8), dimension(n_l) :: s_l      ! neutral sphere radius lattice atoms
    real(8) :: V_pl, V_lp, V_ll         ! Pair potential contributions
    real(8) :: vref_l, vref_p           ! reference pair pot. contrib.
    real(8) :: Ecoh                     ! cohesive energy of part & lattice atoms
    real(8), dimension(3) :: xl, xp     ! for cal. cut-off
    real(8), dimension(3) :: rnnl, rnnp ! next neighbour distance for cut-off
    real(8), dimension(3) :: rnndbeta_l, rnndbeta_p ! divide rnn by beta
    real(8) :: betas0_l, betas0_p       ! beta * s0= for l and p
    real(8) :: kappadbeta_l, kappadbeta_p ! beta * kappa for l and p
    real(8), dimension(3) :: r3temp, r3temp1     ! temporary array variable
    real(8) :: rtemp, rtemp1                    ! temporary variable
    real(8), dimension(n_l) :: rn_ltemp

! Variables and Arrays for partial derivatives
    real(8), dimension(4) :: dtemp       ! array to handle fexp/fexplus
    real(8), dimension(2) :: dchilp, dchipl

    ! Order of parameters eta2, no, eo, lambda, vo, kappa, so.
    real(8), dimension(7) :: dgamma1l, dgamma2l, dgamma1p, dgamma2p
    real(8), dimension(7,n_l) :: dsigma_ll
    real(8), dimension(7,n_l) :: dsigma_lp_l, dsigma_lp_p
    real(8), dimension(7) :: dsigma_pl_l, dsigma_pl_p
    real(8), dimension(7) :: dV_ll
    real(8), dimension(7) :: dV_lp_l, dV_lp_p
    real(8), dimension(7) :: dV_pl_l, dV_pl_p

    real(8), dimension(7,n_l) :: ds_l_l, ds_l_p
    real(8), dimension(7) :: ds_p_l, ds_p_p
    real(8), dimension(7) :: dvref_l_l, dvref_l_p
    real(8), dimension(7) :: dvref_p_l, dvref_p_p
    real(8), dimension(7) :: dEcoh_l_l, dEcoh_l_p
    real(8), dimension(7) :: dEcoh_p_l, dEcoh_p_p




!----------------------VALUES OF FREQUENT USE ---------------------------------
! definition of a few values that appear frequently in calculation
    betas0_l = beta * pars_l%s0
    betas0_p = beta * pars_p%s0
    kappadbeta_l = pars_l%kappa / beta
    kappadbeta_p = pars_p%kappa / beta

! 'coupling' parameters between p and l
! derivatives: (1) derivative over nop, (2) over nol
    dchilp(1) = 1.0 / pars_l%n0         ! d chilp / d nop
    dchipl(2) = 1.0 / pars_p%n0         ! d chipl / d nol

    chilp = pars_p%n0 * dchilp(1)
    chipl = pars_l%n0 * dchipl(2)

    dchipl(1) = - chipl * dchipl(2)     ! d chipl / d nop
    dchilp(2) = - chilp * dchilp(1)     ! d chipl / d nol

!------------------------------------------------------------------------------
!                                  CUT-OFF
!                                  =======
!------------------------------------------------------------------------------
! We use the distance to the next-next nearest neighbours as cut-off.
! FOR FUTURE REVISION:
!            cut-off should be defined via lattice constant _AND_ changeable.

    rcut = a_lat * sqrt_3 * isqrt_2
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
    r3temp = rnnl - betas0_l

    r3temp1 = xl*exp(- pars_l%eta2*r3temp)
    igamma1l = 1.0 / sum(r3temp1)
    dgamma1l = 0.
    dgamma1l(1) = - sum(r3temp*r3temp1)
    dgamma1l(7) = sum(r3temp1)*pars_l%eta2*beta


    r3temp1 = xl*exp(-kappadbeta_l * r3temp)
    igamma2l = 1.0 / sum(r3temp1)
    dgamma2l = 0.
    dgamma2l(6) = - sum(r3temp*r3temp1) / beta
    dgamma2l(7) = sum(r3temp1)*pars_l%kappa


    r3temp = rnnp-betas0_p

    r3temp1 = xp*exp(- pars_p%eta2*r3temp)
    igamma1p = 1.0 / sum(r3temp1)
    dgamma1p = 0.
    dgamma1p(1) = - sum(r3temp*r3temp1)
    dgamma1p(7) = sum(r3temp1)*pars_p%eta2*beta


    r3temp1 = xp*exp(-kappadbeta_p * r3temp)
    igamma2p = 1.0 / sum(r3temp1)
    dgamma2p = 0.
    dgamma2p(6) = - sum(r3temp*r3temp1) / beta
    dgamma2p(7) = sum(r3temp1)*pars_p%kappa



!------------------------------------------------------------------------------
!                          INDIVIDUAL CONTRIBUTIONS
!                          ========================
!------------------------------------------------------------------------------
! The values for the sums are set to zero.

    sigma_ll = 0
    dsigma_ll = 0
    sigma_pl = 0
    dsigma_lp_l=0
    dsigma_lp_p=0
    dsigma_pl_l=0
    dsigma_pl_p=0
    V_ll = 0
    dV_ll = 0
    V_lp = 0
    V_pl = 0
    dV_lp_l = 0
    dV_lp_p = 0
    dV_pl_l = 0
    dV_pl_p = 0
    ds_l_l = 0
    ds_l_p = 0
    ds_p_l = 0
    ds_p_p = 0
    dvref_l_l=0
    dvref_l_p=0
    dvref_p_l=0
    dvref_p_p=0
    dEcoh_l_l = 0
    dEcoh_l_p = 0
    dEcoh_p_l = 0
    dEcoh_p_p = 0

    do i = 1, n_l
        do j = i+1, n_l

        !-----------------PERIODIC BOUNDARY CONDITIONS LATTICE-----------------
        ! Because we want them.

            r3temp(1) = r_lat(1,i)-r_lat(1,j)
            r3temp(2) = r_lat(2,i)-r_lat(2,j)
            r3temp(1) = r3temp(1) - (cell(1)*ANINT(r3temp(1)/cell(1)))
            r3temp(2) = r3temp(2) - (cell(2)*ANINT(r3temp(2)/cell(2)))
            r3temp(3) = r_lat(3,i)-r_lat(3,j)
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

            rtemp1 = rtemp*(r - betas0_l)
            dsigma_ll(1,i) = dsigma_ll(1,i) - rtemp1
            dsigma_ll(1,j) = dsigma_ll(1,j) - rtemp1

            rtemp1 = rtemp*beta*pars_l%eta2
            dsigma_ll(7,i) = dsigma_ll(7,i) + rtemp1
            dsigma_ll(7,j) = dsigma_ll(7,j) + rtemp1



        !-----------------------PAIR POTENTIAL LATTICE-------------------------
        ! For the lattice only.
        ! Will later be subjected to gamma to complete the cut-off.
        ! The particle has no pair potential contribution since there is only
        ! one and thus does not have a partner to interact with.

            rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
            V_ll = V_ll + rtemp

            rtemp1 = rtemp*(r - betas0_l)
            dV_ll(6) = dV_ll(6) + rtemp1

            rtemp1 = rtemp*pars_l%kappa
            dV_ll(7) = dV_ll(7) - rtemp1

        end do


    !-----------------PERIODIC BOUNDERY CONDITIONS PARTICLE--------------------

        r3temp(1) = r_lat(1,i)-r_part(1)
        r3temp(2) = r_lat(2,i)-r_part(2)
        r3temp(1) = r3temp(1) - (cell(1)*ANINT(r3temp(1)/cell(1)))
        r3temp(2) = r3temp(2) - (cell(2)*ANINT(r3temp(2)/cell(2)))
        r3temp(3) = r_lat(3,i)-r_part(3)
        r =  sqrt(sum(r3temp**2))

    !----------------------------THETA PARTICLE--------------------------------

        theta = 1.0 / (1 + exp( acut * (r - rcut) ) )


    !-------------------------------MIXED SIGMA--------------------------------
    ! Contributions of both particle and lattice to neutral sphere radius
    ! To fully include the cut-off, we correct them later by gamma.

        sigma_lp(i) = theta*exp(-pars_p%eta2 * (r - betas0_p) )
        rtemp = theta*exp(-pars_l%eta2 * (r - betas0_l) )
        sigma_pl = sigma_pl + rtemp

        dsigma_lp_p(1,i) = -(r-betas0_p)*sigma_lp(i)

        dsigma_pl_l(1) = dsigma_pl_l(1) - (r - betas0_l)*rtemp
        dsigma_pl_l(7) = pars_l%eta2*beta*sigma_pl


    !--------------------MIXED PAIR POTENTIAL CONTRIUBUTION--------------------

        rtemp = theta*exp(-kappadbeta_p * (r - betas0_p))
        V_lp= V_lp + rtemp
        dV_lp_p(6) = dV_lp_p(6) + rtemp*(r - betas0_p)
        dV_lp_p(7) = dV_lp_p(7) + rtemp*pars_p%kappa

        rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
        V_pl = V_pl + rtemp
        dV_pl_l(6) = dV_pl_l(6) + rtemp*(r - betas0_l)
        dV_pl_l(7) = dV_pl_l(7) + rtemp*pars_l%kappa




    end do


!-------------------------------CUT-OFF ENACTION-------------------------------
! Don't forget the gamma!

    sigma_ll = sigma_ll * igamma1l
    dsigma_ll(1,:) = (dsigma_ll(1,:) - sigma_ll*dgamma1l(1))*igamma1l
    dsigma_ll(7,:) = (dsigma_ll(7,:) - sigma_ll*dgamma1l(7))*igamma1l

    sigma_lp = sigma_lp * igamma1l
    dsigma_lp_l(1,:) = - sigma_lp*igamma1l*dgamma1l(1)
    dsigma_lp_l(7,:) = - sigma_lp*igamma1l*dgamma1l(7)
    dsigma_lp_p(1,:) = dsigma_lp_p(1,:)*igamma1l
    dsigma_lp_p(7,:) = sigma_lp(:)*beta*pars_p%eta2

    dsigma_pl_l(1) = dsigma_pl_l(1)*igamma1p
    dsigma_pl_l(7) = dsigma_pl_l(7)*igamma1p
    sigma_pl = sigma_pl * igamma1p
    dsigma_pl_p(1) = - sigma_pl*igamma1p*dgamma1p(1)
    dsigma_pl_p(7) = - sigma_pl*igamma1p*dgamma1p(7)


    V_ll = V_ll * pars_l%V0 * igamma2l
    dV_ll(5) = - V_ll/pars_l%V0
    dV_ll(6) = (dV_ll(6) * pars_l%V0/beta + V_ll*dgamma2l(6)) * igamma2l
    dV_ll(7) = (dV_ll(7) * pars_l%V0 + V_ll*dgamma2l(7)) * igamma2l

    V_lp = V_lp *chilp * pars_l%V0 * igamma2l
    V_pl = V_pl * pars_p%V0 * igamma2p * chipl

    dV_lp_l(2) = - V_lp/pars_l%n0
    dV_lp_l(5) = - V_lp/pars_l%V0
    dV_lp_l(7) = -V_lp*igamma2l
    dV_lp_l(6) = dV_lp_l(7)*dgamma2l(6)
    dV_lp_l(7) = dV_lp_l(7)*dgamma2l(7)
    dV_lp_p(2) = V_lp/pars_p%n0
    dV_lp_p(6) = dV_lp_p(6)*igamma2l*chilp * pars_l%V0 / beta
    dV_lp_p(7) = dV_lp_p(7)*igamma2l*chilp * pars_l%V0



    dV_pl_l(2) = - V_pl/pars_l%n0
    dV_pl_l(6) = dV_pl_l(6)*igamma2p*pars_p%V0*chipl/beta
    dV_pl_l(7) = dV_pl_l(7)*igamma2p*chipl * pars_p%V0
    dV_pl_p(2) = V_pl/pars_p%n0
    dV_pl_p(5) = - V_pl/pars_p%V0

    dV_pl_p(6) = V_pl*igamma2p*dgamma2p(6)
    dV_pl_p(7) = -V_pl*igamma2p*dgamma2p(7)





!-----------------------------NEUTRAL SPHERE RADIUS----------------------------
! The neutral sphere radius is the radius in which the entire density of the
! atom is included.

    rn_ltemp = sigma_ll + chilp*sigma_lp
    s_l = -log( rn_ltemp * twelveth ) &
            / ( beta * pars_l%eta2)

    rn_ltemp = 1.0/(rn_ltemp*pars_l%eta2*beta)
    ds_l_l(1,:) = -s_l/pars_l%eta2 &
                  - (dsigma_ll(1,:)+chilp*dsigma_lp_l(1,:))*rn_ltemp
    ds_l_l(2,:) = -sigma_lp*rn_ltemp*dchilp(2)
    ds_l_l(7,:) = -rn_ltemp*(dsigma_ll(7,:)+chilp*dsigma_lp_l(7,:))

    ds_l_p(1,:) = - rn_ltemp*chilp*dsigma_lp_p(1,:)
    ds_l_p(2,:) = -sigma_lp*rn_ltemp*dchilp(1)
    ds_l_p(7,:) = -rn_ltemp*chilp*dsigma_lp_p(7,:)


    rtemp=1.0/(beta*pars_p%eta2)
    s_p  = -log( sigma_pl * chipl * twelveth) *rtemp
    ds_p_l(1) = - rtemp /sigma_pl*dsigma_pl_l(1)
    ds_p_l(2) = -rtemp / pars_l%n0
    ds_p_l(7) = - rtemp / sigma_pl*dsigma_pl_l(7)

    ds_p_p(1) = -s_p/pars_p%eta2 - rtemp/(sigma_pl)*dsigma_pl_p(1)
    ds_p_p(2) = rtemp / pars_p%n0
    ds_p_p(7) = - rtemp/sigma_pl*dsigma_pl_p(7)



!----------------MIXED REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------
! These contributions have to be substracted to account for the contributions
! that were included twice.

    rn_ltemp = exp( -pars_l%kappa * s_l)
    rtemp = -12 * pars_l%V0 * pars_l%kappa
    vref_l = 12 * pars_l%V0 * sum(rn_ltemp)
    dvref_l_l(1) = rtemp*sum(rn_ltemp*ds_l_l(1,:))
    dvref_l_l(2) = rtemp*sum(rn_ltemp*ds_l_l(2,:))
    dvref_l_l(5) = vref_l/pars_l%V0
    dvref_l_l(6) = - 12 * pars_l%V0 * sum(rn_ltemp *s_l)
    dvref_l_l(7) = rtemp*sum(rn_ltemp*ds_l_l(7,:))
    dvref_l_p(1) = rtemp*sum(rn_ltemp*ds_l_p(1,:))
    dvref_l_p(2) = rtemp*sum(rn_ltemp*ds_l_p(2,:))
    dvref_l_p(7) = rtemp*sum(rn_ltemp*ds_l_p(7,:))

    rtemp = -pars_p%kappa
    vref_p = 12 * pars_p%V0 * exp( -pars_p%kappa * s_p)
    dvref_p_p(1) = rtemp*vref_p*ds_p_p(1)
    dvref_p_p(2) = rtemp*vref_p*ds_p_p(2)

    dvref_p_p(5) = vref_p/pars_p%V0
    dvref_p_p(6) = - vref_p*s_p

    dvref_p_p(7) = rtemp*vref_p*ds_p_p(7)

    dvref_p_l(1) = rtemp*vref_p*ds_p_l(1)
    dvref_p_l(2) = rtemp*vref_p*ds_p_l(2)
    dvref_p_l(7) = rtemp*vref_p*ds_p_l(7)



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

        ! Wir sind hier und morgen machen wir weiter :-D


    write(*,'(4f12.7)') vref_l/2, vref_p/2
    write(*,'(7f15.7)') dvref_l_l/2
    write(*,'(7f15.7)') dvref_l_p/2
    write(*,'(7f15.7)') dvref_p_p/2
    write(*,'(7f15.7)') dvref_p_l/2
!    write(*,'(7f15.7)') ds_p_p


    stop


!-------------------------------OVERALL ENERGY---------------------------------
! Summation over all contributions.

    energy = Ecoh - V_ll - 0.5 * ( V_lp + V_pl - vref_l - vref_p)

end subroutine emt_fit

end module
