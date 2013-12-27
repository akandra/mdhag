module EMTparms_class
    use emt_init_data

    implicit none
    save

    real(8), parameter                          :: sqrt_2  = 1.41421356237310d0
    real(8), private, parameter                 :: sqrt_3  = 1.73205080756888d0
    real(8), private, parameter                 :: isqrt_2 = 0.707106781186548d0
    real(8), private, parameter                 :: pi      = 3.14159265358979d0
    real(8),  parameter                 :: beta    = 1.8093997906d0
    real(8), private, parameter                 :: twelveth= 0.0833333333333333d0
    integer, private, dimension(3), parameter   :: b       = (/12, 6, 24/)

    ! Calculate the reference enrgy

    ! storage for variables passed into emt_init that are needed by emt energy
    real(8),allocatable             :: r0_lat(:,:)  ! positions of reference lattice atoms
    real(8), private                :: Eref         ! reference energy
    real(8)                         :: a_lat        ! global variable for the lattice constant
    real(8), dimension(14)          :: Ablei


    type EMTparms
        character(2)::  name    = 'RV'

        real(8)       ::  eta2    = 0.d0   ! A^-1
        real(8)       ::  n0      = 0.d0   ! A^-3
        real(8)       ::  E0      = 0.d0   ! eV

        real(8)       ::  lambda  = 0.d0   ! A^-1
        real(8)       ::  V0      = 0.d0   ! eV
        real(8)       ::  kappa   = 0.d0   ! A^-1
        real(8)       ::  s0      = 0.d0   ! A

    end type EMTparms

contains


subroutine emt(a_lat, cell, x_all, n_l, n_p, pars_p, pars_l, energy)
!
! Purpose
!       emt calculates the energy according to the effective medium theory.
! Input variables are (in oder of appearance):
!           a_lat: lattice constant
!           cell : conversion matrix from and to cartesian/direct coordinates
!           x_all : particle and lattice coordinates
!           n_l: number of lattice atoms
!           n_p: number of particle atoms
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

    real(8), intent(in)                 :: a_lat
    real(8), dimension(3,6), intent(in) :: cell   ! cell matrix and its inverse
    real(8), dimension(:,:), intent(in) :: x_all
    integer, intent(in)                 :: n_l, n_p ! number of lattice/particle atoms
    type(EMTparms), intent(inout)       :: pars_p   ! parameters of particle
    type(EMTparms), intent(inout)       :: pars_l   ! parameters of lattice atoms
    real(8), intent(out)                :: energy ! calc. reference energy

! declare the variables that appear in the subroutine

    integer :: i,j, k,q                   ! running parameter
    real(8) :: r                        ! distance
    real(8) :: rcut, rr, acut     ! values to calculate cut-off
    real(8) :: igamma1p, igamma2p       ! inverse gamma for particle
    real(8) :: igamma1l, igamma2l       ! inverse gamma for lattice atoms
    real(8) :: theta                    ! variable for cut-off calculation
    real(8) :: chilp, chipl             ! mixing between lattice (l) and particle (p)
    real(8), dimension(n_l) :: sigma_ll ! contribution to ns, l only
    real(8), dimension(n_p) :: sigma_pp ! contribution to ns, l only
    real(8), dimension(n_l) :: sigma_lp ! mixed contribution to ns
    real(8), dimension(n_p) :: sigma_pl ! mixed contribution to ns
    real(8), dimension(n_l) :: s_l      ! neutral sphere radius lattice atoms
    real(8), dimension(n_p) :: s_p      ! neutral sphere radius lattice atoms
    real(8) :: V_pl, V_lp, V_ll, V_pp   ! Pair potential contributions
    real(8) :: vref_l, vref_p           ! reference pair pot. contrib.
    real(8) :: Ecoh                     ! cohesive energy of part & lattice atoms
    real(8), dimension(3) :: xl, xp     ! for cal. cut-off
    real(8), dimension(3) :: rnnl, rnnp ! next neighbour distance for cut-off
    real(8) :: betas0_l, betas0_p       ! beta * s0= for l and p
    real(8) :: kappadbeta_l, kappadbeta_p ! beta * kappa for l and p
    real(8), dimension(3) :: r3temp     ! temporary array variable
    real(8) :: rtemp                    ! temporary variable

    real(8) :: E_ref                    ! reference energy
    real(8) :: Ecoh_ref                  !cohesive reference energy
    type(EMTparms)      :: particle_parms   ! parameters of particle
    type(EMTparms)      :: lattice_parms    ! parameters of lattice atoms

! For the reference energy
    real(8)                 :: rn_ltemp(n_l), r3temp1(3), rtemp1, teti

!----------------------VALUES OF FREQUENT USE ---------------------------------
! definition of a few values that appear frequently in calculation
    betas0_l = beta * pars_l%s0
    betas0_p = beta * pars_p%s0
    kappadbeta_l = pars_l%kappa / beta
    kappadbeta_p = pars_p%kappa / beta

! 'coupling' parameters between p and l
    chilp = pars_p%n0 / pars_l%n0
    chipl = 1.0d0 / chilp

!------------------------------------------------------------------------------
!                                  CUT-OFF
!                                  =======
!------------------------------------------------------------------------------
! We use the distance to the next-next nearest neighbours as cut-off.
! We only need one cut-off and we choose the one of the lattice atoms since s0
! is usually larger for them.

    rcut = betas0_l * sqrt_3
    rcut = a_lat * sqrt_3 * isqrt_2
    rr = 4 * rcut / (sqrt_3 + 2.0d0)
    acut = 9.210240d0/(rr -rcut) ! ln(10000)


! Distances to the considered neighbours
    rnnl(1) = betas0_l
    rnnl(2) = rnnl(1) * sqrt_2
    rnnl(3) = rnnl(1) * sqrt_3
    rnnp(1) = betas0_p
    rnnp(2) = rnnp(1) * sqrt_2
    rnnp(3) = rnnp(1) * sqrt_3

    xl = b * twelveth / (1.0d0 + exp(acut*(rnnl-rcut)))
    xp = b * twelveth/ (1.0d0 + exp(acut*(rnnp-rcut)))

!-----------------------------------GAMMA--------------------------------------
! Gamma enforces the cut-off together with theta (see below)
! Gamma is defined as inverse.
    r3temp = rnnl-betas0_l
    igamma1l = 1.0d0 / sum(xl*exp(-pars_l%eta2 * r3temp))
    igamma2l = 1.0d0 /sum(xl*exp(-kappadbeta_l * r3temp))

    r3temp = rnnp-betas0_p
    igamma1p = 1.0d0 / sum(xp*exp(-pars_p%eta2 * r3temp))
    igamma2p = 1.0d0 / sum(xp*exp(-kappadbeta_p * r3temp))



!------------------------------------------------------------------------------
!                          INDIVIDUAL CONTRIBUTIONS
!                          ========================
!------------------------------------------------------------------------------
! The values for the sums are set to zero.

    sigma_ll = 0.0d0
    sigma_pp = 0.0d0
    sigma_pl = 0.0d0
    sigma_lp = 0.0d0
    V_ll = 0.0d0
    V_pp = 0.0d0
    V_lp = 0.0d0
    V_pl = 0.0d0

    do i = 1, n_l+n_p
        do j = i+1, n_l+n_p


        !-----------------PERIODIC BOUNDARY CONDITIONS LATTICE-----------------
        ! Because we want them.

            r3temp(1) = x_all(1,i)-x_all(1,j)
            r3temp(2) = x_all(2,i)-x_all(2,j)
            r3temp(3) = x_all(3,i)-x_all(3,j)

            ! transform distances into direct coordinates
            r3temp= matmul(cell(1:3,4:6),r3temp)

            r3temp(1)=r3temp(1)-Anint(r3temp(1))
            r3temp(2)=r3temp(2)-Anint(r3temp(2))
            r3temp(3)=r3temp(3)-Anint(r3temp(3))

            r3temp=matmul(cell(1:3,1:3),r3temp)

            r =  sqrt(sum(r3temp**2))
            teti = teti+r

        !---------------------------THETA PARTICLE------------------------------
        ! Theta enforces the cut-off together with gamma (see above). This
        ! function enacts cutoff by reducing contributions of atoms outside the
        ! cut-off to zero.

            theta = 1.0d0 / (1.0d0 + exp( acut * (r - rcut) ) )


            if (i <= n_p .and. j <= n_p) then
                ! This loop deals with the contributions only due to the particle

                !----------------------------SIGMA PARTICLE-----------------------------
                ! Sigma is a contribution to the neutral sphere radius.
                ! It is a list in which for each lattice atom, the contributions of the
                ! others are summed up. To enforce the cut-off, it will be later
                ! corrected by gamma.

                rtemp = theta*exp(-pars_p%eta2 * (r - betas0_p) )
                sigma_pp(i) = sigma_pp(i) + rtemp
                sigma_pp(j) = sigma_pp(j) + rtemp



                !-----------------------PAIR POTENTIAL PARTICLE-------------------------
                ! Will later be subjected to gamma to complete the cut-off.

                rtemp = theta*exp(-kappadbeta_p * (r - betas0_p))
                V_pp = V_pp + rtemp


            else if (i <= n_p .and. j > n_p ) then
                ! This loop deals with the coupled terms
                k=j-n_p



                !-------------------------------MIXED SIGMA--------------------------------
                ! Contributions of both particle and lattice to neutral sphere radius
                ! To fully include the cut-off, we correct them later by gamma.

                sigma_lp(k) = sigma_lp(k)+theta*exp(-pars_p%eta2 * (r - betas0_p) )

                sigma_pl(i) = sigma_pl(i)+theta*exp(-pars_l%eta2 * (r - betas0_l) )


                !--------------------MIXED PAIR POTENTIAL CONTRIUBUTION--------------------

                rtemp = theta*exp(-kappadbeta_p * (r - betas0_p))
                V_lp= V_lp + rtemp
                rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
                V_pl = V_pl + rtemp


            else if (i > n_p .and. j > n_p) then
                ! This loop deals with the contributions only due to the lattice
                k=j-n_p
                q=i-n_p

                !----------------------------SIGMA LATTICE-----------------------------
                ! Sigma is a contribution to the neutral sphere radius.
                ! It is a list in which for each lattice atom, the contributions of the
                ! others are summed up. To enforce the cut-off, it will be later
                ! corrected by gamma.

                rtemp = theta*exp(-pars_l%eta2 * (r - betas0_l) )
                sigma_ll(q) = sigma_ll(q) + rtemp
                sigma_ll(k) = sigma_ll(k) + rtemp



                !-----------------------PAIR POTENTIAL LATTICE-------------------------
                ! Will later be subjected to gamma to complete the cut-off.

                rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
                V_ll = V_ll + rtemp

            else
                print *, 'Your selection between n_p and n_l is wrong'


            end if
        end do
    end do



!-------------------------------CUT-OFF ENACTION-------------------------------
! Don't forget the gamma!

    sigma_ll = sigma_ll * igamma1l
    sigma_pp = sigma_pp * igamma1p
    sigma_lp = sigma_lp * igamma1l
    sigma_pl = sigma_pl * igamma1p

    V_ll = V_ll * pars_l%V0 * igamma2l
    V_pp = V_pp * pars_p%V0 * igamma2p
    V_lp = V_lp *chilp * pars_l%V0 * igamma2l
    V_pl = V_pl * pars_p%V0 * igamma2p * chipl


!-----------------------------NEUTRAL SPHERE RADIUS----------------------------
! The neutral sphere radius is the radius in which the entire density of the
! atom is included.

    s_l = -log( (sigma_ll + chilp * sigma_lp ) * twelveth ) &
            / ( beta * pars_l%eta2)
    s_p = -log( (sigma_pp + chipl * sigma_pl ) * twelveth ) &
           / ( beta * pars_p%eta2)

!----------------MIXED REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------
! These contributions have to be substracted to account for the contributions
! that were included twice.

    vref_l = 12.0d0 * pars_l%V0 * sum( exp( -pars_l%kappa * s_l) )
    vref_p = 12.0d0 * pars_p%V0 * sum( exp( -pars_p%kappa * s_p) )


!------------------------------------------------------------------------------
!                           CALCULATING THE ENERGY
!                           ======================
!------------------------------------------------------------------------------


!---------------------------COHESIVE ENERGY FUNCTION---------------------------
! Calculates and sums the contributions to the cohesive energy for both lattice
! and particle.
! Pay attention: the -1.0d0 from Ecoh_l has been deleted, because we want the
!                VASP-reference energy.

    Ecoh =  sum( (1.0d0 + pars_l%lambda*s_l) * exp(-pars_l%lambda * s_l)-1.0d0) &
          * pars_l%E0 &
          + sum( (1.0d0 + pars_p%lambda*s_p) * exp(-pars_p%lambda * s_p)) &
          * pars_p%E0

!-------------------------------OVERALL ENERGY---------------------------------
! Summation over all contributions.

    energy = Ecoh - V_ll - V_pp - 0.50d0 * ( V_lp + V_pl - vref_l - vref_p)
    !energy=energy+E_pseudo


end subroutine emt

subroutine emt_fit(a_lat, cell, x_all, n_l, n_p, pars_p, pars_l, energy, denergy)
!
! Purpose
!       emt calculates the energy according to the effective medium theory.
! Input variables are (in oder of appearance):
!           a_lat: lattice constant
!           cell : conversion matrix from and to cartesian/direct coordinates
!           x_all : particle and lattice coordinates
!           n_l: number of lattice atoms
!           n_p: number of particle atoms
!           pars_p  : emt-parameters of particle
!           pars_l  : emt-parameters of lattice
! Output variables are:
!           energy  : emt-energy
!           denergy : emt-energy-derivatives
!
! General Remarks To Fit:
! =======================
! If you are using the potential for fitting, we recommand to refrain from
! fitting so. so_l is directly connected to the lattice constant. If you change
! it, you'll automatically change the lattice constant of your slab.
! so_l = a_lattice /(beta *sqrt_2)
! so_p is related to the DFT density, so, you should not change it, either, if
! you want to reproduce the DFT density.
!
implicit none

!------------------------------------------------------------------------------
!                                   PREAMBLE
!                                   ========
!------------------------------------------------------------------------------

! declare variables and parameters that are passed down from the program

    real(8), intent(in)                 :: a_lat        ! lattice constant of lattice
    real(8), dimension(3,6), intent(in) :: cell         ! cell matrix and its inverse
    integer, intent(in)                 :: n_l, n_p     ! number of lattice/particle atoms
    real(8), dimension(:,:), intent(in) :: x_all        ! positions of lattice atoms
    type(EMTparms), intent(inout)       :: pars_p       ! parameters of particle
    type(EMTparms), intent(inout)       :: pars_l           ! parameters of lattice atoms
    real(8), intent(out)                :: energy    ! calc. reference energy
    real(8), dimension(14), intent(out)  :: denergy  ! derivatives
                                                        ! with respect to to
                                                        ! eta2, followed by no,
                                                        ! eo, lambda, vo, kappa
                                                        ! and so.

! declare the variables that appear in the subroutine and initialise them.

    integer :: i,j, k,q                                 ! running parameter
    real(8) :: r                                        ! distance
    real(8) :: rcut= 0.0d0, rr= 0.0d0, acut= 0.0d0      ! values to calculate cut-off
    real(8) :: igamma1p= 0.0d0, igamma2p= 0.0d0         ! inverse gamma for particle
    real(8) :: igamma1l= 0.0d0, igamma2l= 0.0d0         ! inverse gamma for lattice atoms
    real(8) :: theta= 0.0d0                             ! variable for cut-off calculation
    real(8) :: chilp= 0.0d0, chipl= 0.0d0               ! mixing between lattice (l) and particle (p)
    real(8), dimension(n_p) :: sigma_pl          ! mixed contribution to neutral sphere,
    real(8), dimension(n_p) :: s_p               ! neutral sphere radius particle
    real(8), dimension(n_l) :: sigma_ll          ! contribution to ns, l only
    real(8), dimension(n_p) :: sigma_pp          ! contribution to ns, l only
    real(8), dimension(n_l) :: sigma_lp          ! mixed contribution to ns
    real(8), dimension(n_l) :: s_l               ! neutral sphere radius lattice atoms
    real(8) :: V_pl= 0.0d0, V_lp= 0.0d0, V_ll= 0.0d0, V_pp= 0.0d0   ! Pair potential contributions
    real(8) :: vref_l= 0.0d0, vref_p= 0.0d0             ! reference pair pot. contrib.
    real(8) :: Ecoh= 0.0d0                              ! cohesive energy of part & lattice atoms
    real(8), dimension(3) :: xl= 0.0d0, xp= 0.0d0       ! for cal. cut-off
    real(8), dimension(3) :: rnnl= 0.0d0, rnnp= 0.0d0   ! next neighbour distance for cut-off
    real(8), dimension(3) :: rnndbeta_l= 0.0d0, rnndbeta_p= 0.0d0 ! divide rnn by beta
    real(8) :: betas0_l= 0.0d0, betas0_p= 0.0d0         ! beta * s0= for l and p
    real(8) :: kappadbeta_l= 0.0d0, kappadbeta_p= 0.0d0 ! beta * kappa for l and p
    real(8), dimension(3) :: r3temp= 0.0d0, r3temp1= 0.0d0     ! temporary array variable
    real(8) :: rtemp= 0.0d0, rtemp1= 0.0d0              ! temporary variable
    real(8), dimension(n_l) :: rn_ltemp
    real(8), dimension(n_p) :: rn_ptemp

!-----------------------DECLARE VARIABLES FOR DERIVATIVES----------------------
! Variables and Arrays for partial derivatives
! Apart from chi, all derivatives are 7 long. The first place denotes the
! derivative with respect to to eta2, followed by no, eo, lambda, vo, kappa and so.
! The general notation is:
! e.g. dV_lp_l(1) : the derivative of V_lp with respect to eta2_l.
    real(8), dimension(2) :: dchilp= 0.0d0, dchipl= 0.0d0     ! First element: p, then l
    real(8), dimension(7) :: dgamma1l= 0.0d0, dgamma2l= 0.0d0
    real(8), dimension(7) :: dgamma1p= 0.0d0, dgamma2p= 0.0d0
    real(8), dimension(7,n_l) :: dsigma_ll
    real(8), dimension(7,n_p) :: dsigma_pp
    real(8), dimension(7,n_l) :: dsigma_lp_l
    real(8), dimension(7,n_l) :: dsigma_lp_p
    real(8), dimension(7, n_p) :: dsigma_pl_l, dsigma_pl_p
    real(8), dimension(7) :: dV_ll= 0.0d0, dV_pp= 0.0d0
    real(8), dimension(7) :: dV_lp_l= 0.0d0, dV_lp_p= 0.0d0
    real(8), dimension(7) :: dV_pl_l= 0.0d0, dV_pl_p= 0.0d0
    real(8), dimension(7,n_l) :: ds_l_l, ds_l_p
    real(8), dimension(7,n_p) :: ds_p_l, ds_p_p
    real(8), dimension(7) :: dvref_l_l= 0.0d0, dvref_l_p= 0.0d0
    real(8), dimension(7) :: dvref_p_l= 0.0d0, dvref_p_p= 0.0d0
    real(8), dimension(7) :: dEcoh_l= 0.0d0, dEcoh_p= 0.0d0

!______________________________________________________________________________

!----------------------VALUES OF FREQUENT USE ---------------------------------
! definition of a few values that appear frequently in calculation
    betas0_l = beta * pars_l%s0
    betas0_p = beta * pars_p%s0
    kappadbeta_l = pars_l%kappa / beta
    kappadbeta_p = pars_p%kappa / beta

! 'coupling' parameters between p and l
! In the original paper by Jacobsen et al (SS 366(1996),394), chi also contains
! exponential contributions. We conclude that those exponential contributions
! can be multiplied into the various no and are thus part of the fitting
! parameters.
! derivatives: (1) derivative over nop, (2) over nol
    dchilp(1) = 1.0d0 / pars_l%n0         ! d chilp / d nop
    dchipl(2) = 1.0d0 / pars_p%n0         ! d chipl / d nol

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
    rr = 4.0d0 * rcut / (sqrt_3 + 2.0d0)
    acut = 9.21024d0/(rr -rcut) ! ln(10000)


! Distances to the considered neighbours
    rnnl(1) = betas0_l
    rnnl(2) = rnnl(1) * sqrt_2
    rnnl(3) = rnnl(1) * sqrt_3
    rnnp(1) = betas0_p
    rnnp(2) = rnnp(1) * sqrt_2
    rnnp(3) = rnnp(1) * sqrt_3

    xl = b * twelveth / (1.0d0 + exp(acut*(rnnl-rcut)))
    xp = b * twelveth/ (1.0d0 + exp(acut*(rnnp-rcut)))


!-----------------------------------GAMMA--------------------------------------
! Gamma enforces the cut-off together with theta (see below)
! Gamma is defined as inverse.
! The derivative is not defined as the inverse and formed for each gamma
! individually.
    r3temp = rnnl - betas0_l

    r3temp1 = xl*exp(- pars_l%eta2*r3temp)
    igamma1l = 1.0d0 / sum(r3temp1)
    dgamma1l = 0.d0
    dgamma1l(1) = - sum(r3temp*r3temp1)
    dgamma1l(7) = sum(r3temp1)*pars_l%eta2*beta


    r3temp1 = xl*exp(-kappadbeta_l * r3temp)
    igamma2l = 1.0d0 / sum(r3temp1)
    dgamma2l = 0.
    dgamma2l(6) = - sum(r3temp*r3temp1) / beta
    dgamma2l(7) = sum(r3temp1)*pars_l%kappa


    r3temp = rnnp-betas0_p

    r3temp1 = xp*exp(- pars_p%eta2*r3temp)
    igamma1p = 1.0d0 / sum(r3temp1)
    dgamma1p = 0.
    dgamma1p(1) = - sum(r3temp*r3temp1)
    dgamma1p(7) = sum(r3temp1)*pars_p%eta2*beta


    r3temp1 = xp*exp(-kappadbeta_p * r3temp)
    igamma2p = 1.0d0 / sum(r3temp1)
    dgamma2p = 0.
    dgamma2p(6) = - sum(r3temp*r3temp1) / beta
    dgamma2p(7) = sum(r3temp1)*pars_p%kappa




!------------------------------------------------------------------------------
!                          INDIVIDUAL CONTRIBUTIONS
!                          ========================
!------------------------------------------------------------------------------
! The values for the sums and the derivatives are set to zero.

    sigma_ll = 0.0d0
    sigma_pp = 0.0d0
    dsigma_ll = 0.0d0
    dsigma_pp = 0.0d0
    sigma_pl = 0.0d0
    sigma_lp=0.0d0
    dsigma_lp_l=0.0d0
    dsigma_lp_p=0.0d0
    dsigma_pl_l=0.0d0
    dsigma_pl_p=0.0d0
    V_ll = 0.0d0
    V_pp = 0.0d0
    dV_ll = 0.0d0
    dV_pp = 0.0d0
    V_lp = 0.0d0
    V_pl = 0.0d0
    dV_lp_l = 0.0d0
    dV_lp_p = 0.0d0
    dV_pl_l = 0.0d0
    dV_pl_p = 0.0d0
    ds_l_l = 0.0d0
    ds_l_p = 0.0d0
    ds_p_l = 0.0d0
    ds_p_p = 0.0d0
    dvref_l_l=0.0d0
    dvref_l_p=0.0d0
    dvref_p_l=0.0d0
    dvref_p_p=0.0d0
    dEcoh_l = 0.0d0
    dEcoh_p = 0.0d0
    denergy = 0.0d0
    energy = 0.0d0
    rn_ltemp = 0.0d0
    rn_ptemp=0.0d0
    s_l = 0.0d0
    s_p = 0.0d0

    do i = 1,n_l+n_p
        do j = i+1, n_l+n_p


        !-----------------PERIODIC BOUNDARY CONDITIONS LATTICE-----------------
        ! Because we want them.

            r3temp(1) = x_all(1,i)-x_all(1,j)
            r3temp(2) = x_all(2,i)-x_all(2,j)
            r3temp(3) = x_all(3,i)-x_all(3,j)


        ! transform distances into direct coordinates
            r3temp= matmul(cell(1:3,4:6),r3temp)

            r3temp(1)=r3temp(1)-Anint(r3temp(1))
            r3temp(2)=r3temp(2)-Anint(r3temp(2))
            r3temp(3)=r3temp(3)-Anint(r3temp(3))

            r3temp=matmul(cell(1:3,1:3),r3temp)

            r =  sqrt(sum(r3temp**2))


        !---------------------------THETA LATTICE------------------------------
        ! Theta enforces the cut-off together with gamma (see above). This
        ! function enacts cutoff by reducing contributions of atoms outside the
        ! cut-off to zero.

            theta = 1.0d0 / (1 + exp( acut * (r - rcut) ) )




            if ( i <= n_p .and. j <= n_p) then
                !----------------------------SIGMA PARTICLE-----------------------------
                ! Sigma is a contribution to the neutral sphere radius.
                ! It is a list in which for each particle atom, the contributions of the
                ! others are summed up. To enforce the cut-off, it will be later
                ! corrected by gamma.

                rtemp = theta*exp(-pars_p%eta2 * (r - betas0_p) )
                sigma_pp(i) = sigma_pp(i) + rtemp
                sigma_pp(j) = sigma_pp(j) + rtemp

                rtemp1 = rtemp*(r - betas0_p)
                dsigma_pp(1,i) = dsigma_pp(1,i) - rtemp1
                dsigma_pp(1,j) = dsigma_pp(1,j) - rtemp1

                rtemp1 = rtemp*beta*pars_p%eta2
                dsigma_pp(7,i) = dsigma_pp(7,i) + rtemp1
                dsigma_pp(7,j) = dsigma_pp(7,j) + rtemp1



                !-----------------------PAIR POTENTIAL LATTICE-------------------------
                ! Will later be subjected to gamma to complete the cut-off.

                rtemp = theta*exp(-kappadbeta_p * (r - betas0_p))
                V_pp = V_pp + rtemp

                rtemp1 = rtemp*(r - betas0_p)  !müsste die ganze Geschichte nicht noch einmal durch beta geteilt werden?
                dV_pp(6) = dV_pp(6) + rtemp1

                rtemp1 = rtemp*pars_p%kappa
                dV_pp(7) = dV_pp(7) - rtemp1

            else if ( i <= n_p .and. j > n_p) then
                k = j-n_p   ! becauce we want all the arrays in here to start from 1


                !-------------------------------MIXED SIGMA--------------------------------
                ! Contributions of both particle and lattice to neutral sphere radius
                ! To fully include the cut-off, we correct them later by gamma.
                ! Each of the mixed sigmas depends on both l and p parameters. Not all
                ! contributions need to be under the loop.
                sigma_lp(k) = sigma_lp(k)+theta*exp(-pars_p%eta2 * (r - betas0_p) )

                rtemp = theta*exp(-pars_l%eta2 * (r - betas0_l) )
                sigma_pl(i) = sigma_pl(i)+rtemp

                dsigma_lp_p(1,k) = dsigma_lp_p(1,k)-(r-betas0_p)*sigma_lp(k)

                dsigma_pl_l(1,i) = dsigma_pl_l(1,i) - (r - betas0_l)*rtemp






                !--------------------MIXED PAIR POTENTIAL CONTRIUBUTION--------------------

                rtemp = theta*exp(-kappadbeta_p * (r - betas0_p))
                V_lp= V_lp + rtemp
                dV_lp_p(6) = dV_lp_p(6) + rtemp*(r - betas0_p)
                dV_lp_p(7) = dV_lp_p(7) + rtemp*pars_p%kappa

                rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
                V_pl = V_pl + rtemp
                dV_pl_l(6) = dV_pl_l(6) + rtemp*(r - betas0_l)
                dV_pl_l(7) = dV_pl_l(7) + rtemp*pars_l%kappa

            else if (i > n_p .and. j > n_p) then
                ! This loop deals with the contributions only due to the lattice
                k=j-n_p
                q=i-n_p

                !----------------------------SIGMA LATTICE-----------------------------
                ! Sigma is a contribution to the neutral sphere radius.
                ! It is a list in which for each lattice atom, the contributions of the
                ! others are summed up. To enforce the cut-off, it will be later
                ! corrected by gamma.

                rtemp = theta*exp(-pars_l%eta2 * (r - betas0_l) )
                sigma_ll(q) = sigma_ll(q) + rtemp
                sigma_ll(k) = sigma_ll(k) + rtemp

                rtemp1 = rtemp*(r - betas0_l)
                dsigma_ll(1,q) = dsigma_ll(1,q) - rtemp1
                dsigma_ll(1,k) = dsigma_ll(1,k) - rtemp1

                rtemp1 = rtemp*beta*pars_l%eta2
                dsigma_ll(7,q) = dsigma_ll(7,q) + rtemp1
                dsigma_ll(7,k) = dsigma_ll(7,k) + rtemp1



                !-----------------------PAIR POTENTIAL LATTICE-------------------------
                ! Will later be subjected to gamma to complete the cut-off.

                rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
                V_ll = V_ll + rtemp

                rtemp1 = rtemp*(r - betas0_l)  !müsste die ganze Geschichte nicht noch einmal durch beta geteilt werden?
                dV_ll(6) = dV_ll(6) + rtemp1

                rtemp1 = rtemp*pars_l%kappa
                dV_ll(7) = dV_ll(7) - rtemp1

            else
                print *, 'Your calculations in emt_fit are not running over the right loops'
            end if

        end do

    end do



!-------------------------------CUT-OFF ENACTION-------------------------------
! Don't forget the gamma!

! sigma and its derivatives.
    sigma_ll = sigma_ll * igamma1l
        dsigma_ll(1,:) = (dsigma_ll(1,:) - sigma_ll*dgamma1l(1))*igamma1l
        dsigma_ll(7,:) = (dsigma_ll(7,:) - sigma_ll*dgamma1l(7))*igamma1l

! sigma and its derivatives.
    sigma_pp = sigma_pp * igamma1p
        dsigma_pp(1,:) = (dsigma_pp(1,:) - sigma_pp*dgamma1p(1))*igamma1p
        dsigma_pp(7,:) = (dsigma_pp(7,:) - sigma_pp*dgamma1p(7))*igamma1p


    sigma_lp = sigma_lp * igamma1l
        ! Derivative with respect to l
        dsigma_lp_l(1,:) = - sigma_lp*igamma1l*dgamma1l(1)
        dsigma_lp_l(7,:) = - sigma_lp*igamma1l*dgamma1l(7)
        ! Derivative with respect to p
        dsigma_lp_p(1,:) = dsigma_lp_p(1,:)*igamma1l
        dsigma_lp_p(7,:) = sigma_lp(:)*beta*pars_p%eta2

    sigma_pl = sigma_pl * igamma1p
        ! Derivative with respect to l
        dsigma_pl_l(1,:) = dsigma_pl_l(1,:)*igamma1p
        dsigma_pl_l(7,:) = pars_l%eta2*beta*sigma_pl(:)
        ! Derivative with respect to p
        dsigma_pl_p(1,:) = - sigma_pl*igamma1p*dgamma1p(1)
        dsigma_pl_p(7,:) = - sigma_pl*igamma1p*dgamma1p(7)

! The pair potential and its derivatives
    V_ll = V_ll * pars_l%V0 * igamma2l
        dV_ll(5) = - V_ll/pars_l%V0
        dV_ll(6) = (dV_ll(6) * pars_l%V0/beta + V_ll*dgamma2l(6)) * igamma2l
        dV_ll(7) = -(dV_ll(7) * pars_l%V0 + V_ll*dgamma2l(7)) * igamma2l

    V_pp = V_pp * pars_p%V0 * igamma2p
        dV_pp(5) = - V_pp/pars_p%V0
        dV_pp(6) = (dV_pp(6) * pars_p%V0/beta + V_pp*dgamma2p(6)) * igamma2p
        dV_pp(7) = -(dV_pp(7) * pars_p%V0 + V_pp*dgamma2p(7)) * igamma2p


    V_lp = V_lp *chilp * pars_l%V0 * igamma2l
        ! Derivative with respect to l
        dV_lp_l(2) = V_lp/pars_l%n0
        dV_lp_l(5) = - V_lp/pars_l%V0
        dV_lp_l(7) = V_lp*igamma2l             ! this one is temporary
        dV_lp_l(6) = dV_lp_l(7)*dgamma2l(6)
        dV_lp_l(7) = dV_lp_l(7)*dgamma2l(7)
        ! Derivative with respect to p
        dV_lp_p(2) = -V_lp/pars_p%n0
        dV_lp_p(6) = dV_lp_p(6)*igamma2l*chilp * pars_l%V0 / beta
        dV_lp_p(7) = -dV_lp_p(7)*igamma2l*chilp * pars_l%V0

    V_pl = V_pl * pars_p%V0 * igamma2p * chipl
        ! Derivative with respect to l
        dV_pl_l(2) = - V_pl/pars_l%n0
        dV_pl_l(6) = dV_pl_l(6)*igamma2p*pars_p%V0*chipl/beta
        dV_pl_l(7) = -dV_pl_l(7)*igamma2p*chipl * pars_p%V0
        ! Derivative with respect to p
        dV_pl_p(2) = V_pl/pars_p%n0
        dV_pl_p(5) = - V_pl/pars_p%V0
        dV_pl_p(6) = V_pl*igamma2p*dgamma2p(6)
        dV_pl_p(7) = V_pl*igamma2p*dgamma2p(7)



!-----------------------------NEUTRAL SPHERE RADIUS----------------------------
! The neutral sphere radius is the radius in which the entire density of the
! atom is included.

    rn_ltemp = sigma_ll + chilp*sigma_lp
    s_l = -log( rn_ltemp * twelveth ) &
            / ( beta * pars_l%eta2)

    rn_ltemp = 1.0d0/(rn_ltemp*pars_l%eta2*beta)
        ! Derivative with respect to l
        ds_l_l(1,:) = -s_l/pars_l%eta2 &
                      - (dsigma_ll(1,:)+chilp*dsigma_lp_l(1,:))*rn_ltemp
        ds_l_l(2,:) = -sigma_lp*rn_ltemp*dchilp(2)
        ds_l_l(7,:) = -rn_ltemp*(dsigma_ll(7,:)+chilp*dsigma_lp_l(7,:))
        ! Derivative with respect to p
        ds_l_p(1,:) = - rn_ltemp*chilp*dsigma_lp_p(1,:)
        ds_l_p(2,:) = -sigma_lp*rn_ltemp*dchilp(1)
        ds_l_p(7,:) = -rn_ltemp*chilp*dsigma_lp_p(7,:)

    rn_ptemp= sigma_pp + chipl * sigma_pl
    s_p = -log( rn_ptemp * twelveth ) / (beta*pars_p%eta2)

    rn_ptemp = 1.0d0 / (rn_ptemp *pars_p%eta2*beta)
        ! Derivative with respect to p
        ds_p_p(1,:) = -s_p/pars_p%eta2 &
                      - (dsigma_pp(1,:)+chipl*dsigma_pl_p(1,:))*rn_ptemp
        ds_p_p(2,:) = - sigma_pl*rn_ptemp*dchipl(1)
        ds_p_p(7,:) = - rn_ptemp*(dsigma_pp(7,:)+chipl*dsigma_pl_p(7,:))
        ! Derivative with respect to p
        ds_p_l(1,:) = - rn_ptemp*chipl*dsigma_pl_l(1,:)
        ds_p_l(2,:) = - sigma_pl*rn_ptemp*dchipl(2)
        ds_p_l(7,:) = - rn_ptemp*chipl*dsigma_pl_l(7,:)



!----------------MIXED REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------
! These contributions have to be substracted to account for the contributions
! that were included twice.

    rn_ltemp = exp( -pars_l%kappa * s_l)
    rtemp = -12.0d0 * pars_l%V0 * pars_l%kappa
    vref_l = 12.0d0 * pars_l%V0 * sum(rn_ltemp)
    ! Derivative with respect to l
        dvref_l_l(1) = rtemp*sum(rn_ltemp*ds_l_l(1,:))
        dvref_l_l(2) = rtemp*sum(rn_ltemp*ds_l_l(2,:))
        dvref_l_l(5) = vref_l/pars_l%V0
        dvref_l_l(6) = - 12.0d0 * pars_l%V0 * sum(rn_ltemp *s_l)
        dvref_l_l(7) = rtemp*sum(rn_ltemp*ds_l_l(7,:))
        ! Derivative with respect to p
        dvref_l_p(1) = rtemp*sum(rn_ltemp*ds_l_p(1,:))
        dvref_l_p(2) = rtemp*sum(rn_ltemp*ds_l_p(2,:))
        dvref_l_p(7) = rtemp*sum(rn_ltemp*ds_l_p(7,:))


    rn_ptemp = exp( -pars_p%kappa * s_p)
    rtemp = -12.0d0 * pars_p%V0 * pars_p%kappa
    vref_p = 12.0d0 * pars_p%V0 * sum(rn_ptemp)
    ! Derivative with respect to l
        dvref_p_p(1) = rtemp*sum(rn_ptemp*ds_p_p(1,:))
        dvref_p_p(2) = rtemp*sum(rn_ptemp*ds_p_p(2,:))
        dvref_p_p(5) = vref_p/pars_p%V0
        dvref_p_p(6) = - 12.0d0 * pars_p%V0 * sum(rn_ptemp *s_p)
        dvref_p_p(7) = rtemp*sum(rn_ptemp*ds_p_p(7,:))
        ! Derivative with respect to p
        dvref_p_l(1) = rtemp*sum(rn_ptemp*ds_p_l(1,:))
        dvref_p_l(2) = rtemp*sum(rn_ptemp*ds_p_l(2,:))
        dvref_p_l(7) = rtemp*sum(rn_ptemp*ds_p_l(7,:))

!------------------------------------------------------------------------------
!                           CALCULATING THE ENERGY
!                           ======================
!------------------------------------------------------------------------------


!---------------------------COHESIVE ENERGY FUNCTION---------------------------
! Calculates and sums the contributions to the cohesive energy for both lattice
! and particle.

    Ecoh =  sum( (1.0d0 + pars_l%lambda*s_l) * exp(-pars_l%lambda * s_l)-1.0d0) &
          * pars_l%E0 &
          + sum(( 1.0d0 + pars_p%lambda*s_p) * exp(-pars_p%lambda * s_p)) &
          * pars_p%E0


    rn_ltemp = -pars_l%E0*pars_l%lambda*s_l*exp(-pars_l%lambda*s_l)
    rn_ptemp = -pars_p%lambda*s_p*pars_p%E0*exp(-pars_p%lambda*s_p)
        ! Derivative with respect to l
        dEcoh_l(1) =  sum(pars_l%lambda*rn_ltemp*ds_l_l(1,:))&
                    + sum(pars_p%lambda*rn_ptemp*ds_p_l(1,:))
        dEcoh_l(2) =  sum(pars_l%lambda*rn_ltemp*ds_l_l(2,:))&
                    + sum(pars_p%lambda*rn_ptemp*ds_p_l(2,:))
        dEcoh_l(3) =  sum((1.0d0 + pars_l%lambda*s_l) * exp(-pars_l%lambda * s_l)-1.0d0)
        dEcoh_l(4) =  sum(s_l*rn_ltemp)
        dEcoh_l(7) =  sum(pars_l%lambda*rn_ltemp*ds_l_l(7,:))&
                    + sum(pars_p%lambda*rn_ptemp*ds_p_l(7,:))

        ! Derivative with respect to p
        dEcoh_p(1) =  sum(pars_l%lambda*rn_ltemp*ds_l_p(1,:))&
                    + sum(pars_p%lambda*rn_ptemp*ds_p_p(1,:))
        dEcoh_p(2) =  sum(pars_l%lambda*rn_ltemp*ds_l_p(2,:))&
                    + sum(pars_p%lambda*rn_ptemp*ds_p_p(2,:))
        dEcoh_p(3) =  sum((1.0d0 + pars_p%lambda*s_p) * exp(-pars_p%lambda * s_p))
        dEcoh_p(4) =  sum(s_p*rn_ptemp)
        dEcoh_p(7) =  sum(pars_l%lambda*rn_ltemp*ds_l_p(7,:))&
                    + sum(rn_ptemp*pars_p%lambda*ds_p_p(7,:))


!-------------------------------OVERALL ENERGY---------------------------------
! Summation over all contributions.
! Overall energy
    energy = Ecoh - V_ll - V_pp - 0.5d0* ( V_lp + V_pl - vref_l - vref_p)
    !energy=energy+E_pseudo

    ! Derivative with respect to l
    denergy(8) = dEcoh_l(1) + 0.5d0*( dvref_l_l(1)+dvref_p_l(1))
    denergy(9) = dEcoh_l(2) &
                    +0.5d0*(dV_pl_l(2)+dvref_p_l(2)+dV_lp_l(2)+dvref_l_l(2))
    denergy(10) = dEcoh_l(3)
    denergy(11) = dEcoh_l(4)
    denergy(12) = dV_ll(5) + 0.5d0*(dvref_l_l(5)+dV_lp_l(5))
    denergy(13) = dV_ll(6) + 0.5d0*( dV_lp_l(6) + dV_pl_l(6) + dvref_l_l(6))
    denergy(14) = dEcoh_l(7) + dV_ll(7) &
                   + 0.5d0*(dV_lp_l(7)+dV_pl_l(7)+dvref_l_l(7)+dvref_p_l(7))

    ! Derivative with respect to p (no correction by dEref since those do not
    ! contain any p-contribution)
    denergy(1) = dEcoh_p(1) + 0.5d0*(dvref_l_p(1)+dvref_p_p(1))
    denergy(2) = dEcoh_p(2) &
                   + 0.5d0*(dV_pl_p(2)+dV_lp_p(2)+dvref_l_p(2)+dvref_p_p(2))
    denergy(3) = dEcoh_p(3)
    denergy(4) = dEcoh_p(4)
    denergy(5) = dV_pp(5)+0.5d0*(dV_pl_p(5) + dvref_p_p(5))
    denergy(6) = dV_pp(6)+0.5d0*(dV_lp_p(6)+dV_pl_p(6)+dvref_p_p(6))
    denergy(7) = dEcoh_p(7) + dV_pp(7)+ &
                   + 0.5d0*(dV_lp_p(7)+dV_pl_p(7)+dvref_l_p(7)+dvref_p_p(7))

    Ablei=denergy



end subroutine emt_fit




subroutine emt_parms2array (emt_parms, array)
    type(EMTparms)          :: emt_parms
    real(8), dimension(7)   :: array

    array(1) = emt_parms%eta2
    array(2) = emt_parms%n0
    array(3) = emt_parms%E0
    array(4) = emt_parms%lambda
    array(5) = emt_parms%V0
    array(6) = emt_parms%kappa
    array(7) = emt_parms%s0
end subroutine emt_parms2array

subroutine array2emt_parms (array, emt_parms)
    type(EMTparms)          :: emt_parms
    real(8), dimension(7)   :: array

    emt_parms%eta2   = array(1)
    emt_parms%n0     = array(2)
    emt_parms%E0     = array(3)
    emt_parms%lambda = array(4)
    emt_parms%V0     = array(5)
    emt_parms%kappa  = array(6)
    emt_parms%s0     = array(7)
end subroutine array2emt_parms


!subroutine emt_mixed (a_lat, cell, r_part, r_lat, n_l, pars_p, pars_l, energy)
!!
!! Purpose:
!!       emt calculates the energy according to the effective medium theory.
!! Input variables are (in oder of appearance):
!!           particle: incident particle
!!           lattice : lattice atoms
!!           r0      : equilibrium positions of lattice atoms
!!           pars_p  : emt-parameters of particle
!!           pars_l  : emt-parameters of lattice
!! Output variables are:
!!           energy  : emt-energy
!implicit none
!
!!------------------------------------------------------------------------------
!!                                   PREAMBLE
!!                                   ========
!!------------------------------------------------------------------------------
!
!! declare variables and parameters that are passed down from the program
!
!    real(8), intent(in)                 :: a_lat
!    real(8), dimension(3,6), intent(in) :: cell   ! cell matrix and its inverse
!    integer, intent(in)                 :: n_l      ! number of lattice atoms
!    real(8), dimension(:), intent (in)  :: r_part  ! position of the particle (note:
!                                                    ! only one position for reference)
!    real(8), dimension(:,:), intent(in) :: r_lat   ! positions of lattice atoms
!    type(EMTparms), intent(inout)       :: pars_p   ! parameters of particle
!    type(EMTparms), intent(inout)       :: pars_l   ! parameters of lattice atoms
!    real(8), intent(out)                :: energy ! calc. reference energy
!
!! declare the variables that appear in the subroutine
!
!    integer :: i,j, k                   ! running parameter
!    real(8) :: r                        ! distance
!    real(8) :: rcut, rr, acut           ! values to calculate cut-off
!    real(8) :: igamma1p, igamma2p       ! inverse gamma for particle
!    real(8) :: igamma1l, igamma2l       ! inverse gamma for lattice atoms
!    real(8) :: theta                    ! variable for cut-off calculation
!    real(8) :: chilp, chipl             ! mixing between lattice (l) and particle (p)
!    real(8) :: sigma_pl                 ! mixed contribution to neutral sphere,
!    real(8) :: s_p                      ! neutral sphere radius particle
!    real(8), dimension(n_l) :: sigma_ll ! contribution to ns, l only
!    real(8), dimension(n_l) :: sigma_lp ! mixed contribution to ns
!    real(8), dimension(n_l) :: s_l      ! neutral sphere radius lattice atoms
!    real(8) :: V_pl, V_lp, V_ll         ! Pair potential contributions
!    real(8) :: vref_l, vref_p           ! reference pair pot. contrib.
!    real(8) :: Ecoh                     ! cohesive energy of part & lattice atoms
!    real(8), dimension(3) :: xl, xp     ! for cal. cut-off
!    real(8), dimension(3) :: rnnl, rnnp ! next neighbour distance for cut-off
!    real(8) :: betas0_l, betas0_p       ! beta * s0= for l and p
!    real(8) :: kappadbeta_l, kappadbeta_p ! beta * kappa for l and p
!    real(8), dimension(3) :: r3temp     ! temporary array variable
!    real(8) :: rtemp                    ! temporary variable
!
!    real(8) :: E_ref                    ! reference energy
!    real(8) :: Ecoh_ref                  !cohesive reference energy
!!    type(EMTparms)      :: particle_parms   ! parameters of particle
!!    type(EMTparms)      :: lattice_parms    ! parameters of lattice atoms
!
!! For the reference energy
!    real(8)                 :: rn_ltemp(n_l), r3temp1(3), rtemp1
!    real(8), dimension(n_l) :: s_l_ref
!    real(8)                 :: vref_l_ref
!
!!----------------------VALUES OF FREQUENT USE ---------------------------------
!! definition of a few values that appear frequently in calculation
!    betas0_l = beta * pars_l%s0
!    betas0_p = beta * pars_p%s0
!    kappadbeta_l = pars_l%kappa / beta
!    kappadbeta_p = pars_p%kappa / beta
!
!! 'coupling' parameters between p and l
!    chilp = pars_p%n0 / pars_l%n0
!    chipl = 1.0d0 / chilp
!
!
!!------------------------------------------------------------------------------
!!                                  CUT-OFF
!!                                  =======
!!------------------------------------------------------------------------------
!! We use the distance to the next-next nearest neighbours as cut-off.
!! FOR FUTURE REVISION:
!!            cut-off should be defined via lattice constant _AND_ changeable.
!
!!    rcut = betas0_l * sqrt_3
!!    rr = 4 * rcut / (sqrt_3 + 2)
!!    acut = 9.21024/(rr -rcut) ! ln(10000)
!
!    rcut = a_lat * sqrt_3 * isqrt_2
!    rr = 4.0d0 * rcut / (sqrt_3 + 2.0d0)
!    acut = 9.21024d0/(rr -rcut) ! ln(10000)
!
!! Distances to the considered neighbours
!    rnnl(1) = betas0_l
!    rnnl(2) = rnnl(1) * sqrt_2
!    rnnl(3) = rnnl(1) * sqrt_3
!    rnnp(1) = betas0_p
!    rnnp(2) = rnnp(1) * sqrt_2
!    rnnp(3) = rnnp(1) * sqrt_3
!
!    xl = b * twelveth / (1.0d0 + exp(acut*(rnnl-rcut)))
!    xp = b * twelveth/ (1.0d0 + exp(acut*(rnnp-rcut)))
!
!!-----------------------------------GAMMA--------------------------------------
!! Gamma enforces the cut-off together with theta (see below)
!! Gamma is defined as inverse.
!    r3temp = rnnl-betas0_l
!    igamma1l = 1.0d0 / sum(xl*exp(-pars_l%eta2 * r3temp))
!    igamma2l = 1.0d0 /sum(xl*exp(-kappadbeta_l * r3temp))
!
!    r3temp = rnnp-betas0_p
!    igamma1p = 1.0d0 / sum(xp*exp(-pars_p%eta2 * r3temp))
!    igamma2p = 1.0d0 / sum(xp*exp(-kappadbeta_p * r3temp))
!
!!------------------------------------------------------------------------------
!!                          INDIVIDUAL CONTRIBUTIONS
!!                          ========================
!!------------------------------------------------------------------------------
!! The values for the sums are set to zero.
!
!    sigma_ll = 0
!    sigma_pl = 0
!    V_ll = 0
!    V_lp = 0
!    V_pl = 0
!
!    do i = 1, n_l
!        do j = i+1, n_l
!
!
!        !-----------------PERIODIC BOUNDARY CONDITIONS LATTICE-----------------
!        ! Because we want them.
!
!            r3temp(1) = r_lat(1,i)-r_lat(1,j)
!            r3temp(2) = r_lat(2,i)-r_lat(2,j)
!            r3temp(3) = r_lat(3,i)-r_lat(3,j)
!
!            ! transform distances into direct coordinates
!            r3temp= matmul(cell(1:3,4:6),r3temp)
!
!            r3temp(1)=r3temp(1)-Anint(r3temp(1))
!            r3temp(2)=r3temp(2)-Anint(r3temp(2))
!            r3temp(3)=r3temp(3)-Anint(r3temp(3))
!
!            r3temp=matmul(cell(1:3,1:3),r3temp)
!
!            r =  sqrt(sum(r3temp**2))
!
!
!        !---------------------------THETA LATTICE------------------------------
!        ! Theta enforces the cut-off together with gamma (see above). This
!        ! function enacts cutoff by reducing contributions of atoms outside the
!        ! cut-off to zero.
!
!            theta = 1.0d0 / (1.0d0 + exp( acut * (r - rcut) ) )
!
!
!        !----------------------------SIGMA LATTICE-----------------------------
!        ! Sigma is a contribution to the neutral sphere radius.
!        ! It is a list in which for each lattice atom, the contributions of the
!        ! others are summed up. To enforce the cut-off, it will be later
!        ! corrected by gamma.
!        ! sigma_pp does not exist because there is only a single particle.
!
!            rtemp = theta*exp(-pars_l%eta2 * (r - betas0_l) )
!            sigma_ll(i) = sigma_ll(i) + rtemp
!            sigma_ll(j) = sigma_ll(j) + rtemp
!
!
!
!        !-----------------------PAIR POTENTIAL LATTICE-------------------------
!        ! For the lattice only.
!        ! Will later be subjected to gamma to complete the cut-off.
!        ! The particle has no pair potential contribution since there is only
!        ! one and thus does not have a partner to interact with.
!
!            rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
!            V_ll = V_ll + rtemp
!
!
!        end do
!
!    !-----------------PERIODIC BOUNDERY CONDITIONS PARTICLE--------------------
!
!        r3temp(1) = r_lat(1,i)-r_part(1)
!        r3temp(2) = r_lat(2,i)-r_part(2)
!        r3temp(3) = r_lat(3,i)-r_part(3)
!    ! transform distances into direct coordinates
!        r3temp= matmul(cell(1:3,4:6),r3temp)
!
!        r3temp(1)=r3temp(1)-Anint(r3temp(1))
!        r3temp(2)=r3temp(2)-Anint(r3temp(2))
!        r3temp(3)=r3temp(3)-Anint(r3temp(3))
!
!        r3temp=matmul(cell(1:3,1:3),r3temp)
!        r =  sqrt(sum(r3temp**2))
!
!
!
!    !----------------------------THETA PARTICLE--------------------------------
!
!        theta = 1.0d0 / (1.0d0 + exp( acut * (r - rcut) ) )
!
!
!    !-------------------------------MIXED SIGMA--------------------------------
!    ! Contributions of both particle and lattice to neutral sphere radius
!    ! To fully include the cut-off, we correct them later by gamma.
!
!        sigma_lp(i) = theta*exp(-pars_p%eta2 * (r - betas0_p) )
!        rtemp = theta*exp(-pars_l%eta2 * (r - betas0_l) )
!        sigma_pl = sigma_pl + rtemp
!
!
!    !--------------------MIXED PAIR POTENTIAL CONTRIUBUTION--------------------
!
!        rtemp = theta*exp(-kappadbeta_p * (r - betas0_p))
!        V_lp= V_lp + rtemp
!        rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
!        V_pl = V_pl + rtemp
!
!    end do
!
!!-------------------------------CUT-OFF ENACTION-------------------------------
!! Don't forget the gamma!
!
!    sigma_ll = sigma_ll * igamma1l
!    sigma_lp = sigma_lp * igamma1l
!    sigma_pl = sigma_pl * igamma1p
!
!    V_ll = V_ll * pars_l%V0 * igamma2l
!    V_lp = V_lp *chilp * pars_l%V0 * igamma2l
!    V_pl = V_pl * pars_p%V0 * igamma2p * chipl
!
!
!
!!-----------------------------NEUTRAL SPHERE RADIUS----------------------------
!! The neutral sphere radius is the radius in which the entire density of the
!! atom is included.
!
!    s_l = -log( (sigma_ll + chilp * sigma_lp) * twelveth ) &
!            / ( beta * pars_l%eta2)
!    s_p  = -log( sigma_pl * chipl * twelveth) / ( beta * pars_p%eta2)
!
!!    s_l_ref = -log( sigma_ll*twelveth )/( beta*pars_l%eta2)
!
!
!!----------------MIXED REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------
!! These contributions have to be substracted to account for the contributions
!! that were included twice.
!
!    vref_l = 12.0d0 * pars_l%V0 * sum( exp( -pars_l%kappa * s_l) )
!    vref_p = 12.0d0 * pars_p%V0 * exp( -pars_p%kappa * s_p)
!
!    ! Reference energy
!!    vref_l_ref = 12.0d0 * pars_l%V0 * sum(exp( -pars_l%kappa * s_l_ref))
!
!!------------------------------------------------------------------------------
!!                           CALCULATING THE ENERGY
!!                           ======================
!!------------------------------------------------------------------------------
!
!
!!---------------------------COHESIVE ENERGY FUNCTION---------------------------
!! Calculates and sums the contributions to the cohesive energy for both lattice
!! and particle.
!
!    Ecoh = sum( (1.0d0 + pars_l%lambda*s_l) * exp(-pars_l%lambda * s_l)-1.d0 ) &
!          * pars_l%E0 &
!          + (1.0d0 + pars_p%lambda*s_p) * exp(-pars_p%lambda * s_p)* pars_p%E0
!
!    ! Reference energy
!! Ecoh_ref = sum((1.0d0 + pars_l%lambda*s_l_ref)*exp(-pars_l%lambda*s_l_ref)-1.0d0)&
!!          * pars_l%E0
!
!
!!-------------------------------OVERALL ENERGY---------------------------------
!! Summation over all contributions.
!
!
!    energy = Ecoh - V_ll - 0.50d0 * ( V_lp + V_pl - vref_l - vref_p)!-E_ref
!
!
!end subroutine emt_mixed




end module

!
!subroutine emt_fit_mixed (a_lat, cell, r_part, r_lat, n_l, pars_p, pars_l, energy, denergy)
!!subroutine emt_fit (cell, a_lat, n_l, r_lat, r_part, pars_p, pars_l, energy, denergy_l, denergy_p)
!!
!! Purpose:
!! ========
!!       emt calculates the energy according to the effective medium theory.
!!       This emt includes the derivitives for the Fit and the Forces
!! Input variables are (in oder of appearance):
!!           particle: incident particle
!!           lattice : lattice atoms
!!           r0      : equilibrium positions of lattice atoms
!!           pars_p  : emt-parameters of particle
!!           pars_l  : emt-parameters of lattice
!! Output variables are:
!!           energy  : emt-energy
!!
!! General Remarks To Fit:
!! =======================
!! If you are using the potential for fitting, we recommand to refrain from
!! fitting so. so_l is directly connected to the lattice constant. If you change
!! it, you'll automatically change the lattice constant of your slab.
!! so_l = a_lattice /(beta *sqrt_2)
!! so_p is related to the DFT density, so, you should not change it, either, if
!! you want to reproduce the DFT density.
!!
!implicit none
!
!!------------------------------------------------------------------------------
!!                                   PREAMBLE
!!                                   ========
!!------------------------------------------------------------------------------
!
!! declare variables and parameters that are passed down from the program
!
!    real(8), dimension(3,6), intent(in) :: cell     ! cell matrix and its inverse
!    real(8), intent(in)                 :: a_lat   ! lattice constant of lattice
!    integer, intent(in)                 :: n_l      ! number of lattice atoms
!    real(8), dimension(3), intent (in)  :: r_part   ! positions of the particle
!    real(8), dimension(:,:), intent(in) :: r_lat    ! positions of lattice atoms
!    type(EMTparms), intent(inout)       :: pars_p   ! parameters of particle
!    type(EMTparms), intent(inout)       :: pars_l   ! parameters of lattice atoms
!    real(8), intent(out)                :: energy   ! calc. reference energy
!    real(8), dimension(14), intent(out)  :: denergy ! derivatives
!                                                        ! with respect to to
!                                                        ! eta2, followed by no,
!                                                        ! eo, lambda, vo, kappa
!                                                        ! and so.
!
!! declare the variables that appear in the subroutine
!
!    integer :: i,j, k                   ! running parameter
!    real(8) :: r                        ! distance
!    real(8) :: rcut, rr, acut           ! values to calculate cut-off
!    real(8) :: igamma1p, igamma2p       ! inverse gamma for particle
!    real(8) :: igamma1l, igamma2l       ! inverse gamma for lattice atoms
!    real(8) :: theta                    ! variable for cut-off calculation
!    real(8) :: chilp, chipl             ! mixing between lattice (l) and particle (p)
!    real(8) :: sigma_pl                 ! mixed contribution to neutral sphere,
!    real(8) :: s_p                      ! neutral sphere radius particle
!    real(8), dimension(n_l) :: sigma_ll ! contribution to ns, l only
!    real(8), dimension(n_l) :: sigma_lp ! mixed contribution to ns
!    real(8), dimension(n_l) :: s_l      ! neutral sphere radius lattice atoms
!    real(8) :: V_pl, V_lp, V_ll         ! Pair potential contributions
!    real(8) :: vref_l, vref_p           ! reference pair pot. contrib.
!    real(8) :: Ecoh                     ! cohesive energy of part & lattice atoms
!    real(8), dimension(3) :: xl, xp     ! for cal. cut-off
!    real(8), dimension(3) :: rnnl, rnnp ! next neighbour distance for cut-off
!    real(8), dimension(3) :: rnndbeta_l, rnndbeta_p ! divide rnn by beta
!    real(8) :: betas0_l, betas0_p       ! beta * s0= for l and p
!    real(8) :: kappadbeta_l, kappadbeta_p ! beta * kappa for l and p
!    real(8), dimension(3) :: r3temp, r3temp1     ! temporary array variable
!    real(8) :: rtemp, rtemp1                    ! temporary variable
!    real(8), dimension(n_l) :: rn_ltemp
!
!!-----------------------DECLARE VARIABLES FOR DERIVATIVES----------------------
!! Variables and Arrays for partial derivatives
!! Apart from chi, all derivatives are 7 long. The first place denotes the
!! derivative with respect to to eta2, followed by no, eo, lambda, vo, kappa and so.
!! The general notation is:
!! e.g. dV_lp_l(1) : the derivative of V_lp with respect to eta2_l.
!    real(8), dimension(2) :: dchilp, dchipl     ! First element: p, then l
!    real(8), dimension(7) :: dgamma1l, dgamma2l
!    real(8), dimension(7) :: dgamma1p, dgamma2p
!    real(8), dimension(7,n_l) :: dsigma_ll
!    real(8), dimension(7,n_l) :: dsigma_lp_l
!    real(8), dimension(7,n_l) :: dsigma_lp_p
!    real(8), dimension(7) :: dsigma_pl_l, dsigma_pl_p
!    real(8), dimension(7) :: dV_ll
!    real(8), dimension(7) :: dV_lp_l, dV_lp_p
!    real(8), dimension(7) :: dV_pl_l, dV_pl_p
!    real(8), dimension(7,n_l) :: ds_l_l, ds_l_p
!    real(8), dimension(7) :: ds_p_l, ds_p_p
!    real(8), dimension(7) :: dvref_l_l, dvref_l_p
!    real(8), dimension(7) :: dvref_p_l, dvref_p_p
!    real(8), dimension(7) :: dEcoh_l, dEcoh_p
!
!!______________________________________________________________________________
!
!!----------------------VALUES OF FREQUENT USE ---------------------------------
!! definition of a few values that appear frequently in calculation
!    betas0_l = beta * pars_l%s0
!    betas0_p = beta * pars_p%s0
!    kappadbeta_l = pars_l%kappa / beta
!    kappadbeta_p = pars_p%kappa / beta
!
!! 'coupling' parameters between p and l
!! In the original paper by Jacobsen et al (SS 366(1996),394), chi also contains
!! exponential contributions. We conclude that those exponential contributions
!! can be multiplied into the various no and are thus part of the fitting
!! parameters.
!! derivatives: (1) derivative over nop, (2) over nol
!    dchilp(1) = 1.0d0 / pars_l%n0         ! d chilp / d nop
!    dchipl(2) = 1.0d0 / pars_p%n0         ! d chipl / d nol
!
!    chilp = pars_p%n0 * dchilp(1)
!    chipl = pars_l%n0 * dchipl(2)
!
!    dchipl(1) = - chipl * dchipl(2)     ! d chipl / d nop
!    dchilp(2) = - chilp * dchilp(1)     ! d chipl / d nol
!
!!------------------------------------------------------------------------------
!!                                  CUT-OFF
!!                                  =======
!!------------------------------------------------------------------------------
!! We use the distance to the next-next nearest neighbours as cut-off.
!! FOR FUTURE REVISION:
!!            cut-off should be defined via lattice constant _AND_ changeable.
!
!    rcut = a_lat * sqrt_3 * isqrt_2
!    rr = 4.0d0 * rcut / (sqrt_3 + 2.0d0)
!    acut = 9.21024d0/(rr -rcut) ! ln(10000)
!
!
!! Distances to the considered neighbours
!    rnnl(1) = betas0_l
!    rnnl(2) = rnnl(1) * sqrt_2
!    rnnl(3) = rnnl(1) * sqrt_3
!    rnnp(1) = betas0_p
!    rnnp(2) = rnnp(1) * sqrt_2
!    rnnp(3) = rnnp(1) * sqrt_3
!
!    xl = b * twelveth / (1.0d0 + exp(acut*(rnnl-rcut)))
!    xp = b * twelveth/ (1.0d0 + exp(acut*(rnnp-rcut)))
!
!
!!-----------------------------------GAMMA--------------------------------------
!! Gamma enforces the cut-off together with theta (see below)
!! Gamma is defined as inverse.
!! The derivative is not defined as the inverse and formed for each gamma
!! individually.
!    r3temp = rnnl - betas0_l
!
!    r3temp1 = xl*exp(- pars_l%eta2*r3temp)
!    igamma1l = 1.0d0 / sum(r3temp1)
!    dgamma1l = 0.d0
!    dgamma1l(1) = - sum(r3temp*r3temp1)
!    dgamma1l(7) = sum(r3temp1)*pars_l%eta2*beta
!
!
!    r3temp1 = xl*exp(-kappadbeta_l * r3temp)
!    igamma2l = 1.0d0 / sum(r3temp1)
!    dgamma2l = 0.
!    dgamma2l(6) = - sum(r3temp*r3temp1) / beta
!    dgamma2l(7) = sum(r3temp1)*pars_l%kappa
!
!
!    r3temp = rnnp-betas0_p
!
!    r3temp1 = xp*exp(- pars_p%eta2*r3temp)
!    igamma1p = 1.0d0 / sum(r3temp1)
!    dgamma1p = 0.
!    dgamma1p(1) = - sum(r3temp*r3temp1)
!    dgamma1p(7) = sum(r3temp1)*pars_p%eta2*beta
!
!
!    r3temp1 = xp*exp(-kappadbeta_p * r3temp)
!    igamma2p = 1.0d0 / sum(r3temp1)
!    dgamma2p = 0.
!    dgamma2p(6) = - sum(r3temp*r3temp1) / beta
!    dgamma2p(7) = sum(r3temp1)*pars_p%kappa
!
!
!
!
!!------------------------------------------------------------------------------
!!                          INDIVIDUAL CONTRIBUTIONS
!!                          ========================
!!------------------------------------------------------------------------------
!! The values for the sums and the derivatives are set to zero.
!
!    sigma_ll = 0
!    dsigma_ll = 0
!    sigma_pl = 0
!    dsigma_lp_l=0
!    dsigma_lp_p=0
!    dsigma_pl_l=0
!    dsigma_pl_p=0
!    V_ll = 0
!    dV_ll = 0
!    V_lp = 0
!    V_pl = 0
!    dV_lp_l = 0
!    dV_lp_p = 0
!    dV_pl_l = 0
!    dV_pl_p = 0
!    ds_l_l = 0
!    ds_l_p = 0
!    ds_p_l = 0
!    ds_p_p = 0
!    dvref_l_l=0
!    dvref_l_p=0
!    dvref_p_l=0
!    dvref_p_p=0
!    dEcoh_l = 0
!    dEcoh_p = 0
!    denergy = 0
!
!    do i = 1,n_l
!        do j = i+1, n_l
!
!
!        !-----------------PERIODIC BOUNDARY CONDITIONS LATTICE-----------------
!        ! Because we want them.
!
!            r3temp(1) = r_lat(1,i)-r_lat(1,j)
!            r3temp(2) = r_lat(2,i)-r_lat(2,j)
!            r3temp(3) = r_lat(3,i)-r_lat(3,j)
!
!
!        ! transform distances into direct coordinates
!            r3temp= matmul(cell(1:3,4:6),r3temp)
!
!            r3temp(1)=r3temp(1)-Anint(r3temp(1))
!            r3temp(2)=r3temp(2)-Anint(r3temp(2))
!            r3temp(3)=r3temp(3)-Anint(r3temp(3))
!
!            r3temp=matmul(cell(1:3,1:3),r3temp)
!
!            r =  sqrt(sum(r3temp**2))
!
!
!        !---------------------------THETA LATTICE------------------------------
!        ! Theta enforces the cut-off together with gamma (see above). This
!        ! function enacts cutoff by reducing contributions of atoms outside the
!        ! cut-off to zero.
!
!            theta = 1.0d0 / (1 + exp( acut * (r - rcut) ) )
!
!
!
!
!
!        !----------------------------SIGMA LATTICE-----------------------------
!        ! Sigma is a contribution to the neutral sphere radius.
!        ! It is a list in which for each lattice atom, the contributions of the
!        ! others are summed up. To enforce the cut-off, it will be later
!        ! corrected by gamma.
!        ! sigma_pp does not exist because there is only a single particle.
!
!            rtemp = theta*exp(-pars_l%eta2 * (r - betas0_l) )
!            sigma_ll(i) = sigma_ll(i) + rtemp
!            sigma_ll(j) = sigma_ll(j) + rtemp
!
!            rtemp1 = rtemp*(r - betas0_l)
!            dsigma_ll(1,i) = dsigma_ll(1,i) - rtemp1
!            dsigma_ll(1,j) = dsigma_ll(1,j) - rtemp1
!
!            rtemp1 = rtemp*beta*pars_l%eta2
!            dsigma_ll(7,i) = dsigma_ll(7,i) + rtemp1
!            dsigma_ll(7,j) = dsigma_ll(7,j) + rtemp1
!
!
!
!        !-----------------------PAIR POTENTIAL LATTICE-------------------------
!        ! For the lattice only.
!        ! Will later be subjected to gamma to complete the cut-off.
!        ! The particle has no pair potential contribution since there is only
!        ! one and thus does not have a partner to interact with.
!
!            rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
!            V_ll = V_ll + rtemp
!
!            rtemp1 = rtemp*(r - betas0_l)  ! müsste das hier nicht r/beta-sol sein?
!            dV_ll(6) = dV_ll(6) + rtemp1
!
!            rtemp1 = rtemp*pars_l%kappa
!            dV_ll(7) = dV_ll(7) - rtemp1
!
!        end do
!
!
!    !-----------------PERIODIC BOUNDARY CONDITIONS LATTICE-----------------
!        ! Because we want them.
!
!            r3temp(1) = r_lat(1,i)-r_part(1)
!            r3temp(2) = r_lat(2,i)-r_part(2)
!            r3temp(3) = r_lat(3,i)-r_part(3)
!
!
!        ! transform distances into direct coordinates
!            r3temp= matmul(cell(1:3,4:6),r3temp)
!
!            r3temp(1)=r3temp(1)-Anint(r3temp(1))
!            r3temp(2)=r3temp(2)-Anint(r3temp(2))
!            r3temp(3)=r3temp(3)-Anint(r3temp(3))
!
!            r3temp=matmul(cell(1:3,1:3),r3temp)
!
!            r =  sqrt(sum(r3temp**2))
!
!
!    !----------------------------THETA PARTICLE--------------------------------
!
!        theta = 1.0d0 / (1.0d0 + exp( acut * (r - rcut) ) )
!
!    !-------------------------------MIXED SIGMA--------------------------------
!    ! Contributions of both particle and lattice to neutral sphere radius
!    ! To fully include the cut-off, we correct them later by gamma.
!    ! Each of the mixed sigmas depends on both l and p parameters. Not all
!    ! contributions need to be under the loop.
!
!        sigma_lp(i) = theta*exp(-pars_p%eta2 * (r - betas0_p) )
!        rtemp = theta*exp(-pars_l%eta2 * (r - betas0_l) )
!        sigma_pl = sigma_pl + rtemp
!
!        dsigma_lp_p(1,i) = -(r-betas0_p)*sigma_lp(i)
!
!        dsigma_pl_l(1) = dsigma_pl_l(1) - (r - betas0_l)*rtemp
!        dsigma_pl_l(7) = pars_l%eta2*beta*sigma_pl
!
!
!    !--------------------MIXED PAIR POTENTIAL CONTRIUBUTION--------------------
!
!        rtemp = theta*exp(-kappadbeta_p * (r - betas0_p))
!        V_lp= V_lp + rtemp
!        dV_lp_p(6) = dV_lp_p(6) + rtemp*(r - betas0_p)
!        dV_lp_p(7) = dV_lp_p(7) + rtemp*pars_p%kappa
!
!        rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
!        V_pl = V_pl + rtemp
!        dV_pl_l(6) = dV_pl_l(6) + rtemp*(r - betas0_l)
!        dV_pl_l(7) = dV_pl_l(7) + rtemp*pars_l%kappa
!
!
!    end do
!
!
!!-------------------------------CUT-OFF ENACTION-------------------------------
!! Don't forget the gamma!
!
!! sigma and its derivatives.
!    sigma_ll = sigma_ll * igamma1l
!        dsigma_ll(1,:) = (dsigma_ll(1,:) - sigma_ll*dgamma1l(1))*igamma1l
!        dsigma_ll(7,:) = (dsigma_ll(7,:) - sigma_ll*dgamma1l(7))*igamma1l
!
!    sigma_lp = sigma_lp * igamma1l
!        ! Derivative with respect to l
!        dsigma_lp_l(1,:) = - sigma_lp*igamma1l*dgamma1l(1)
!        dsigma_lp_l(7,:) = - sigma_lp*igamma1l*dgamma1l(7)
!        ! Derivative with respect to p
!        dsigma_lp_p(1,:) = dsigma_lp_p(1,:)*igamma1l
!        dsigma_lp_p(7,:) = sigma_lp(:)*beta*pars_p%eta2
!
!    sigma_pl = sigma_pl * igamma1p
!        ! Derivative with respect to l
!        dsigma_pl_l(1) = dsigma_pl_l(1)*igamma1p
!        dsigma_pl_l(7) = dsigma_pl_l(7)*igamma1p
!        ! Derivative with respect to p
!        dsigma_pl_p(1) = - sigma_pl*igamma1p*dgamma1p(1)
!        dsigma_pl_p(7) = - sigma_pl*igamma1p*dgamma1p(7)
!
!! The pair potential and its derivatives
!    V_ll = V_ll * pars_l%V0 * igamma2l
!        dV_ll(5) = - V_ll/pars_l%V0
!        dV_ll(6) = (dV_ll(6) * pars_l%V0/beta + V_ll*dgamma2l(6)) * igamma2l
!        dV_ll(7) = -(dV_ll(7) * pars_l%V0 + V_ll*dgamma2l(7)) * igamma2l
!
!    V_lp = V_lp *chilp * pars_l%V0 * igamma2l
!        ! Derivative with respect to l
!        dV_lp_l(2) = V_lp/pars_l%n0
!        dV_lp_l(5) = - V_lp/pars_l%V0
!        dV_lp_l(7) = V_lp*igamma2l             ! this one is temporary
!        dV_lp_l(6) = dV_lp_l(7)*dgamma2l(6)
!        dV_lp_l(7) = dV_lp_l(7)*dgamma2l(7)
!        ! Derivative with respect to p
!        dV_lp_p(2) = -V_lp/pars_p%n0
!        dV_lp_p(6) = dV_lp_p(6)*igamma2l*chilp * pars_l%V0 / beta
!        dV_lp_p(7) = -dV_lp_p(7)*igamma2l*chilp * pars_l%V0
!
!    V_pl = V_pl * pars_p%V0 * igamma2p * chipl
!        ! Derivative with respect to l
!        dV_pl_l(2) = - V_pl/pars_l%n0
!        dV_pl_l(6) = dV_pl_l(6)*igamma2p*pars_p%V0*chipl/beta
!        dV_pl_l(7) = -dV_pl_l(7)*igamma2p*chipl * pars_p%V0
!        ! Derivative with respect to p
!        dV_pl_p(2) = V_pl/pars_p%n0
!        dV_pl_p(5) = - V_pl/pars_p%V0
!        dV_pl_p(6) = V_pl*igamma2p*dgamma2p(6)
!        dV_pl_p(7) = V_pl*igamma2p*dgamma2p(7)
!
!
!
!!-----------------------------NEUTRAL SPHERE RADIUS----------------------------
!! The neutral sphere radius is the radius in which the entire density of the
!! atom is included.
!
!    rn_ltemp = sigma_ll + chilp*sigma_lp
!    s_l = -log( rn_ltemp * twelveth ) &
!            / ( beta * pars_l%eta2)
!
!    rn_ltemp = 1.0d0/(rn_ltemp*pars_l%eta2*beta)
!        ! Derivative with respect to l
!        ds_l_l(1,:) = -s_l/pars_l%eta2 &
!                      - (dsigma_ll(1,:)+chilp*dsigma_lp_l(1,:))*rn_ltemp
!        ds_l_l(2,:) = -sigma_lp*rn_ltemp*dchilp(2)
!        ds_l_l(7,:) = -rn_ltemp*(dsigma_ll(7,:)+chilp*dsigma_lp_l(7,:))
!        ! Derivative with respect to p
!        ds_l_p(1,:) = - rn_ltemp*chilp*dsigma_lp_p(1,:)
!        ds_l_p(2,:) = -sigma_lp*rn_ltemp*dchilp(1)
!        ds_l_p(7,:) = -rn_ltemp*chilp*dsigma_lp_p(7,:)
!
!    rtemp=1.0d0/(beta*pars_p%eta2)
!    s_p  = -log( sigma_pl * chipl * twelveth) *rtemp
!        ! Derivative with respect to l
!        ds_p_l(1) = - rtemp /sigma_pl*dsigma_pl_l(1)
!        ds_p_l(2) = - rtemp / pars_l%n0
!        ds_p_l(7) = - rtemp / sigma_pl*dsigma_pl_l(7)
!        ! Derivative with respect to p
!        ds_p_p(1) = -s_p/pars_p%eta2 - rtemp/(sigma_pl)*dsigma_pl_p(1)
!        ds_p_p(2) = rtemp /pars_p%n0
!        ds_p_p(7) = - rtemp/sigma_pl*dsigma_pl_p(7)
!
!
!!----------------MIXED REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------
!! These contributions have to be substracted to account for the contributions
!! that were included twice.
!
!    rn_ltemp = exp( -pars_l%kappa * s_l)
!    rtemp = -12.0d0 * pars_l%V0 * pars_l%kappa
!    vref_l = 12.0d0 * pars_l%V0 * sum(rn_ltemp)
!    ! Derivative with respect to l
!        dvref_l_l(1) = rtemp*sum(rn_ltemp*ds_l_l(1,:))
!        dvref_l_l(2) = rtemp*sum(rn_ltemp*ds_l_l(2,:))
!        dvref_l_l(5) = vref_l/pars_l%V0
!        dvref_l_l(6) = - 12.0d0 * pars_l%V0 * sum(rn_ltemp *s_l)
!        dvref_l_l(7) = rtemp*sum(rn_ltemp*ds_l_l(7,:))
!        ! Derivative with respect to p
!        dvref_l_p(1) = rtemp*sum(rn_ltemp*ds_l_p(1,:))
!        dvref_l_p(2) = rtemp*sum(rn_ltemp*ds_l_p(2,:))
!        dvref_l_p(7) = rtemp*sum(rn_ltemp*ds_l_p(7,:))
!
!    rtemp = -pars_p%kappa
!    vref_p = 12.0d0 * pars_p%V0 * exp( -pars_p%kappa * s_p)
!        ! Derivative with respect to p
!        dvref_p_p(1) = rtemp*vref_p*ds_p_p(1)
!        dvref_p_p(2) = rtemp*vref_p*ds_p_p(2)
!        dvref_p_p(5) = vref_p/pars_p%V0
!        dvref_p_p(6) = - vref_p*s_p
!        dvref_p_p(7) = rtemp*vref_p*ds_p_p(7)
!        ! Derivative with respect to l
!        dvref_p_l(1) = rtemp*vref_p*ds_p_l(1)
!        dvref_p_l(2) = rtemp*vref_p*ds_p_l(2)
!        dvref_p_l(7) = rtemp*vref_p*ds_p_l(7)
!
!!------------------------------------------------------------------------------
!!                           CALCULATING THE ENERGY
!!                           ======================
!!------------------------------------------------------------------------------
!
!
!!---------------------------COHESIVE ENERGY FUNCTION---------------------------
!! Calculates and sums the contributions to the cohesive energy for both lattice
!! and particle.
!
!    Ecoh = sum( (1.0d0 + pars_l%lambda*s_l) * exp(-pars_l%lambda * s_l)-1.0d0 ) &
!          * pars_l%E0 &
!          + (1.0d0 + pars_p%lambda*s_p) * exp(-pars_p%lambda * s_p)* pars_p%E0
!
!
!    rn_ltemp=-pars_l%E0*pars_l%lambda*s_l*exp(-pars_l%lambda*s_l)
!    rtemp = -pars_p%lambda*s_p*pars_p%E0*exp(-pars_p%lambda*s_p)
!        ! Derivative with respect to l
!        dEcoh_l(1) = sum(pars_l%lambda*rn_ltemp*ds_l_l(1,:))&
!                    +pars_p%lambda*rtemp*ds_p_l(1)
!        dEcoh_l(2) = sum(pars_l%lambda*rn_ltemp*ds_l_l(2,:))&
!                    + pars_p%lambda*rtemp*ds_p_l(2)
!        dEcoh_l(3) = sum( (1.0d0 + pars_l%lambda*s_l) * exp(-pars_l%lambda * s_l)-1.0d0 )
!        dEcoh_l(4) = sum(s_l*rn_ltemp)
!        dEcoh_l(7) = sum(pars_l%lambda*rn_ltemp*ds_l_l(7,:))&
!                    +pars_p%lambda*rtemp*ds_p_l(7)
!
!    rtemp = -pars_p%lambda*s_p*pars_p%E0*exp(-pars_p%lambda*s_p)
!        ! Derivative with respect to p
!        dEcoh_p(1) = sum(pars_l%lambda*rn_ltemp*ds_l_p(1,:))&
!                    +pars_p%lambda*rtemp*ds_p_p(1)
!        dEcoh_p(2) = sum(pars_l%lambda*rn_ltemp*ds_l_p(2,:))&
!                    +pars_p%lambda*rtemp*ds_p_p(2)
!        dEcoh_p(3) = (1.0d0 + pars_p%lambda*s_p) * exp(-pars_p%lambda * s_p)
!        dEcoh_p(4) = s_p*rtemp
!        dEcoh_p(7) = sum(pars_l%lambda*rn_ltemp*ds_l_p(7,:))&
!                    +rtemp*pars_p%lambda*ds_p_p(7)
!
!
!!-------------------------------OVERALL ENERGY---------------------------------
!! Summation over all contributions.
!! Overall energy
!    energy = Ecoh - V_ll - 0.5d0* ( V_lp + V_pl - vref_l - vref_p)
!
!    ! Derivative with respect to l
!    denergy(8) = dEcoh_l(1) + 0.5d0*( dvref_l_l(1)+dvref_p_l(1))
!    denergy(9) = dEcoh_l(2) &
!                    +0.5d0*(dV_pl_l(2)+dvref_p_l(2)+dV_lp_l(2)+dvref_l_l(2))
!    denergy(10) = dEcoh_l(3)
!    denergy(11) = dEcoh_l(4)
!    denergy(12) = dV_ll(5) + 0.5d0*(dvref_l_l(5)+dV_lp_l(5))
!    denergy(13) = dV_ll(6) + 0.5d0*( dV_lp_l(6) + dV_pl_l(6) + dvref_l_l(6))
!    denergy(14) = dEcoh_l(7) + dV_ll(7) &
!                   + 0.5d0*(dV_lp_l(7)+dV_pl_l(7)+dvref_l_l(7)+dvref_p_l(7))
!
!    ! Derivative with respect to p (no correction by dEref since those do not
!    ! contain any p-contribution)
!    denergy(1) = dEcoh_p(1) + 0.5d0*(dvref_l_p(1)+dvref_p_p(1))
!    denergy(2) = dEcoh_p(2) &
!                   + 0.5d0*(dV_pl_p(2)+dV_lp_p(2)+dvref_l_p(2)+dvref_p_p(2))
!    denergy(3) = dEcoh_p(3)
!    denergy(4) = dEcoh_p(4)
!    denergy(5) = 0.5d0*(dV_pl_p(5) + dvref_p_p(5))
!    denergy(6) = 0.5d0*(dV_lp_p(6)+dV_pl_p(6)+dvref_p_p(6))
!    denergy(7) = dEcoh_p(7) &
!                   + 0.5d0*(dV_lp_p(7)+dV_pl_p(7)+dvref_l_p(7)+dvref_p_p(7))
!
!denergy(9) = ds_p_p(2)
!
!end subroutine emt_fit_mixed
!
!subroutine emt_l (a_lat, cell_0, n_l0, r0_lat_in, pars_p, pars_l, energy)
!!subroutine emt_init (cell_in, n_l_in, r0_lat_in, pars_p, pars_l, energy)
!!
!! Purpose:
!!           Here, the fitting procedure is just implemented and the reference
!!           energy calculated.
!!           Only the gold-parameters are calculated.
!!
!implicit none
!
!!------------------------------------------------------------------------------
!!                                   PREAMBLE
!!                                   ========
!!------------------------------------------------------------------------------
!
!! declare variables and parameters that are passed down from the program
!
!    real(8), intent(in)                 :: a_lat
!    real(8), dimension(3,6), intent(in)   :: cell_0         ! cell-matrix(1-3)
!                                                            !and inverse matrix(4-6)
!    integer, intent(in)                 :: n_l0          ! number of lattice atoms
!    real(8), dimension(:,:), intent(in) :: r0_lat_in       ! positions of lattice atoms
!    type(EMTparms), intent(inout)       :: pars_p       ! parameters of particle
!    type(EMTparms), intent(inout)       :: pars_l       ! parameters of lattice atoms
!    real(8), intent(out)                :: energy       ! calculated reference energy
!
!! declare the variables that appear in the subroutine
!
!    integer :: i,j, k                   ! running parameter
!    real(8) :: r                        ! distance
!    real(8) :: rcut, rr, acut           ! values to calculate cut-off
!    real(8) :: igamma1l, igamma2l       ! inverse gamma for lattice atoms
!    real(8) :: theta                    ! variable for cut-off calculation
!    real(8) :: chilp                    ! mixing between lattice (l) and particle (p)
!    real(8), dimension(n_l0) :: sigma_ll ! contribution to ns, l only
!    real(8), dimension(n_l0) :: s_l      ! neutral sphere radius lattice atoms
!    real(8) :: V_ll                     ! Pair potential contributions
!    real(8) :: vref_l                   ! reference pair pot. contrib.
!    real(8) :: Ecoh                     ! cohesive energy of part & lattice atoms
!    real(8), dimension(3) :: xl         ! for cal. cut-off
!    real(8), dimension(3) :: rnnl       ! next neighbour distance for cut-off
!    real(8) :: betas0_l                 ! beta * s0= for l and p
!    real(8) :: kappadbeta_l             ! beta * kappa for l and p
!    real(8), dimension(3) :: r3temp     ! temporary array variable
!    real(8) :: rtemp                    ! temporary variable
!
!    real(8) :: teti
!
!
!
!!-----------------------Save inputs in module for use by emt ------------------
!
!    if (.not. allocated(r0_lat)) then
!        allocate(r0_lat(3,n_l0))
!        write(*,*) 'allocated r0_lat for', n_l0, 'atoms'
!    end if
!    r0_lat = r0_lat_in
!
!!----------------------------VALUES OF FREQUENT USE----------------------------
!! definition of a few values that appear frequently in calculation
!    betas0_l = beta * pars_l%s0
!    kappadbeta_l = pars_l%kappa / beta
!
!! 'coupling' parameters between p and l
!    chilp = pars_p%n0 / pars_l%n0
!
!
!!------------------------------------------------------------------------------
!!                                  CUT-OFF
!!                                  =======
!!------------------------------------------------------------------------------
!! We use the distance to the next-next nearest neighbours as cut-off.
!! FOR FUTURE REVISION:
!!            cut-off should be defined via lattice constant _AND_ changeable.
!
!!    rcut = betas0_l * sqrt_3
!!    rr = 4 * rcut / (sqrt_3 + 2)
!!    acut = 9.21024/(rr -rcut) ! ln(10000)
!
!    rcut = a_lat * sqrt_3 * isqrt_2
!    rr = 4 * rcut / (sqrt_3 + 2.0d0)
!    acut = 9.21024d0/(rr -rcut) ! ln(10000)
!
!! Distances to the considered neighbours
!    rnnl(1) = betas0_l
!    rnnl(2) = rnnl(1) * sqrt_2
!    rnnl(3) = rnnl(1) * sqrt_3
!
!    xl = b * twelveth / (1.0d0 + exp(acut*(rnnl-rcut)))
!
!!-----------------------------------GAMMA--------------------------------------
!! Gamma enforces the cut-off together with theta (see below)
!! Gamma is defined as inverse.
!    r3temp = rnnl-betas0_l
!    igamma1l = 1.0d0 / sum(xl*exp(-pars_l%eta2 * r3temp))
!    igamma2l = 1.0d0 /sum(xl*exp(-kappadbeta_l * r3temp))
!
!
!!------------------------------------------------------------------------------
!!                          INDIVIDUAL CONTRIBUTIONS
!!                          ========================
!!------------------------------------------------------------------------------
!! The values for the sums are set to zero.
!
!    sigma_ll = 0.0d0
!    V_ll = 0.0d0
!
!    do i = 1, n_l0
!        do j = i+1, n_l0
!
!        !-----------------PERIODIC BOUNDARY CONDITIONS LATTICE-----------------
!        ! Because we want them.
!
!            r3temp(1) = r0_lat(1,i)-r0_lat(1,j)
!            r3temp(2) = r0_lat(2,i)-r0_lat(2,j)
!            r3temp(3) = r0_lat(3,i)-r0_lat(3,j)
!
!
!        ! transform distances into direct coordinates
!            r3temp= matmul(cell_0(1:3,4:6),r3temp)
!
!            r3temp(1)=r3temp(1)-Anint(r3temp(1))
!            r3temp(2)=r3temp(2)-Anint(r3temp(2))
!            r3temp(3)=r3temp(3)-Anint(r3temp(3))
!
!            r3temp=matmul(cell_0(1:3,1:3),r3temp)
!
!            r =  sqrt(sum(r3temp**2))
!
!
!        !---------------------------THETA LATTICE------------------------------
!        ! Theta enforces the cut-off together with gamma (see above). This
!        ! function enacts cutoff by reducing contributions of atoms outside the
!        ! cut-off to zero.
!
!            theta = 1.0 / (1 + exp( acut * (r - rcut) ) )
!
!
!        !----------------------------SIGMA LATTICE-----------------------------
!        ! Sigma is a contribution to the neutral sphere radius.
!        ! It is a list in which for each lattice atom, the contributions of the
!        ! others are summed up. To enforce the cut-off, it will be later
!        ! corrected by gamma.
!        ! sigma_pp does not exist because there is only a single particle.
!
!            rtemp = theta*exp(-pars_l%eta2 * (r - betas0_l) )
!            sigma_ll(i) = sigma_ll(i) + rtemp
!            sigma_ll(j) = sigma_ll(j) + rtemp
!
!
!
!        !-----------------------PAIR POTENTIAL LATTICE-------------------------
!        ! For the lattice only.
!        ! Will later be subjected to gamma to complete the cut-off.
!        ! The particle has no pair potential contribution since there is only
!        ! one and thus does not have a partner to interact with.
!
!            rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
!            V_ll = V_ll + rtemp
!
!        end do
!    end do
!
!
!!-------------------------------CUT-OFF ENACTION-------------------------------
!! Don't forget the gamma!
!
!    sigma_ll = sigma_ll * igamma1l
!    V_ll = V_ll * pars_l%V0 * igamma2l  ! Muss V_ll nicht auch nochmal durch 2 geteilt werden?
!
!!-----------------------------NEUTRAL SPHERE RADIUS----------------------------
!! The neutral sphere radius is the radius in which the entire density of the
!! atom is included.
!
!    s_l = -log( sigma_ll * twelveth ) &
!            / ( beta * pars_l%eta2)
!
!!----------------MIXED REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------
!! These contributions have to be subtracted to account for the contributions
!! that were included twice.
!
!    vref_l = 12.d0 * pars_l%V0 * sum( exp( -pars_l%kappa * s_l) )
!
!!------------------------------------------------------------------------------
!!                           CALCULATING THE ENERGY
!!                           ======================
!!------------------------------------------------------------------------------
!
!
!!---------------------------COHESIVE ENERGY FUNCTION---------------------------
!! Calculates and sums the contributions to the cohesive energy for both lattice
!! and particle.
!
!    Ecoh = sum( (1 + pars_l%lambda*s_l) * exp(-pars_l%lambda * s_l)-1 ) &
!          * pars_l%E0
!
!!-------------------------------OVERALL ENERGY---------------------------------
!! Summation over all contributions.
!
!   energy = Ecoh - V_ll + 0.5d0 * vref_l
!!   energy = acut
!   Eref   = energy
!
!end subroutine emt_l
!
!subroutine emt_fit_l (a_lat, cell_0, r0_lat, n_l0, pars_p, pars_l, E_dref, dE_ref)
!!subroutine emt_fit (cell, a_lat, n_l, r_lat, r_part, pars_p, pars_l, energy, denergy_l, denergy_p)
!!
!! Purpose:
!! ========
!!       emt_fit_init calculates the reference energy according to the effective
!!       medium theory including the derivatives for the reference energy.
!
!! Input variables are:
!!           lattice : lattice atoms
!!           r0      : equilibrium positions of lattice atoms
!!           pars_p  : emt-parameters of particle
!!           pars_l  : emt-parameters of lattice
!! Output variables are:
!!           energy  : emt-energy
!!
!! General Remarks To Fit:
!! =======================
!! If you are using the potential for fitting, we recommand to refrain from
!! fitting so. so_l is directly connected to the lattice constant. If you change
!! it, you'll automatically change the lattice constant of your slab.
!! so_l = a_lattice /(beta *sqrt_2)
!! so_p is related to the DFT density, so, you should not change it, either, if
!! you want to reproduce the DFT density.
!!
!implicit none
!
!!------------------------------------------------------------------------------
!!                                   PREAMBLE
!!                                   ========
!!------------------------------------------------------------------------------
!
!! declare variables and parameters that are passed down from the program
!
!    real(8), dimension(3,6), intent(in)   :: cell_0  ! cell matrix and its inverse
!    real(8), intent(in)                 :: a_lat   ! lattice constant of lattice
!    integer, intent(in)                 :: n_l0      ! number of lattice atoms
!    real(8), dimension(:,:), intent(in) :: r0_lat    ! positions of lattice atoms
!    type(EMTparms), intent(inout)       :: pars_p   ! parameters of particle
!    type(EMTparms), intent(inout)       :: pars_l   ! parameters of lattice atoms
!    real(8), intent(out)                :: E_dref   ! calc. reference energy
!    real(8), dimension(7), intent(out)  :: dE_ref ! derivatives
!                                                        ! with respect to to
!                                                        ! eta2, followed by no,
!                                                        ! eo, lambda, vo, kappa
!                                                        ! and so.
!
!! declare the variables that appear in the subroutine
!
!    integer :: i,j, k                   ! running parameter
!    real(8) :: r                        ! distance
!    real(8) :: rcut, rr, acut           ! values to calculate cut-off
!    real(8) :: igamma1l, igamma2l       ! inverse gamma for lattice atoms
!    real(8) :: theta                    ! variable for cut-off calculation
!    real(8), dimension(n_l0) :: sigma_ll ! contribution to ns, l only
!    real(8), dimension(n_l0) :: s_l_ref
!    real(8) ::  V_ll                    ! Pair potential contributions
!    real(8) :: vref_l_ref
!    real(8), dimension(3) :: xl         ! for cal. cut-off
!    real(8), dimension(3) :: rnnl       ! next neighbour distance for cut-off
!    real(8), dimension(3) :: rnndbeta_l ! divide rnn by beta
!    real(8) :: betas0_l                 ! beta * s0= for l and p
!    real(8) :: kappadbeta_l             ! beta * kappa for l and p
!    real(8), dimension(3) :: r3temp, r3temp1     ! temporary array variable
!    real(8) :: rtemp, rtemp1                    ! temporary variable
!    real(8), dimension(n_l0) :: rn_ltemp
!    real(8) :: Ecoh_ref, E_ref
!
!!-----------------------DECLARE VARIABLES FOR DERIVATIVES----------------------
!
!    real(8), dimension(7) :: dgamma1l, dgamma2l
!    real(8), dimension(7,n_l0) :: dsigma_ll
!    real(8), dimension(7) :: dV_ll
!    real(8), dimension(7,n_l0) :: ds_l_l_ref
!    real(8), dimension(7) :: dvref_l_l_ref
!    real(8), dimension(7) :: dEcoh_ref
!    real(8) :: chilp                    ! mixing between lattice (l) and particle (p)
!
!
!!______________________________________________________________________________
!
!!----------------------VALUES OF FREQUENT USE ---------------------------------
!! definition of a few values that appear frequently in calculation
!    betas0_l = beta * pars_l%s0
!    kappadbeta_l = pars_l%kappa / beta
!
!! 'coupling' parameters between p and l
!    chilp = pars_p%n0 / pars_l%n0
!
!!------------------------------------------------------------------------------
!!                                  CUT-OFF
!!                                  =======
!!------------------------------------------------------------------------------
!! We use the distance to the next-next nearest neighbours as cut-off.
!! FOR FUTURE REVISION:
!!            cut-off should be defined via lattice constant _AND_ changeable.
!
!    rcut = a_lat * sqrt_3 * isqrt_2
!    rr = 4.0d0 * rcut / (sqrt_3 + 2.0d0)
!    acut = 9.21024d0/(rr -rcut) ! ln(10000)
!
!
!! Distances to the considered neighbours
!    rnnl(1) = betas0_l
!    rnnl(2) = rnnl(1) * sqrt_2
!    rnnl(3) = rnnl(1) * sqrt_3
!
!    xl = b * twelveth / (1.0d0 + exp(acut*(rnnl-rcut)))
!
!!-----------------------------------GAMMA--------------------------------------
!! Gamma enforces the cut-off together with theta (see below)
!! Gamma is defined as inverse.
!! The derivative is not defined as the inverse and formed for each gamma
!! individually.
!    r3temp = rnnl - betas0_l
!    r3temp1 = xl*exp(- pars_l%eta2*r3temp)
!    igamma1l = 1.0d0 / sum(r3temp1)
!    dgamma1l = 0.d0
!    dgamma1l(1) = - sum(r3temp*r3temp1)
!    dgamma1l(7) = sum(r3temp1)*pars_l%eta2*beta
!
!
!    r3temp1 = xl*exp(-kappadbeta_l * r3temp)
!    igamma2l = 1.0d0 / sum(r3temp1)
!    dgamma2l = 0.
!    dgamma2l(6) = - sum(r3temp*r3temp1) / beta
!    dgamma2l(7) = sum(r3temp1)*pars_l%kappa
!
!
!!------------------------------------------------------------------------------
!!                          INDIVIDUAL CONTRIBUTIONS
!!                          ========================
!!------------------------------------------------------------------------------
!! The values for the sums and the derivatives are set to zero.
!
!    sigma_ll = 0
!    dsigma_ll = 0
!    V_ll = 0
!    dV_ll = 0
!    s_l_ref = 0
!    ds_l_l_ref = 0
!    dvref_l_l_ref=0
!    vref_l_ref=0
!    Ecoh_ref=0
!    dEcoh_ref=0
!    E_ref=0
!
!    dE_ref=0
!
!    do i = 1, n_l0
!        do j = i+1, n_l0
!
!
!        !-----------------PERIODIC BOUNDARY CONDITIONS LATTICE-----------------
!        ! Because we want them.
!
!            r3temp(1) = r0_lat(1,i)-r0_lat(1,j)
!            r3temp(2) = r0_lat(2,i)-r0_lat(2,j)
!            r3temp(3) = r0_lat(3,i)-r0_lat(3,j)
!
!
!        ! transform distances into direct coordinates
!            r3temp= matmul(cell_0(1:3,4:6),r3temp)
!
!            r3temp(1)=r3temp(1)-Anint(r3temp(1))
!            r3temp(2)=r3temp(2)-Anint(r3temp(2))
!            r3temp(3)=r3temp(3)-Anint(r3temp(3))
!
!            r3temp=matmul(cell_0(1:3,1:3),r3temp)
!
!            r =  sqrt(sum(r3temp**2))
!
!
!        !---------------------------THETA LATTICE------------------------------
!        ! Theta enforces the cut-off together with gamma (see above). This
!        ! function enacts cutoff by reducing contributions of atoms outside the
!        ! cut-off to zero.
!
!            theta = 1.0d0 / (1 + exp( acut * (r - rcut) ) )
!
!
!
!        !----------------------------SIGMA LATTICE-----------------------------
!        ! Sigma is a contribution to the neutral sphere radius.
!        ! It is a list in which for each lattice atom, the contributions of the
!        ! others are summed up. To enforce the cut-off, it will be later
!        ! corrected by gamma.
!        ! sigma_pp does not exist because there is only a single particle.
!
!            rtemp = theta*exp(-pars_l%eta2 * (r - betas0_l) )
!            sigma_ll(i) = sigma_ll(i) + rtemp
!            sigma_ll(j) = sigma_ll(j) + rtemp
!
!            rtemp1 = rtemp*(r - betas0_l)
!            dsigma_ll(1,i) = dsigma_ll(1,i) - rtemp1
!            dsigma_ll(1,j) = dsigma_ll(1,j) - rtemp1
!
!            rtemp1 = rtemp*beta*pars_l%eta2
!            dsigma_ll(7,i) = dsigma_ll(7,i) + rtemp1
!            dsigma_ll(7,j) = dsigma_ll(7,j) + rtemp1
!
!
!
!        !-----------------------PAIR POTENTIAL LATTICE-------------------------
!        ! For the lattice only.
!        ! Will later be subjected to gamma to complete the cut-off.
!        ! The particle has no pair potential contribution since there is only
!        ! one and thus does not have a partner to interact with.
!
!            rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
!            V_ll = V_ll + rtemp
!
!            rtemp1 = rtemp*(r - betas0_l)
!            dV_ll(6) = dV_ll(6) + rtemp1
!
!            rtemp1 = rtemp*pars_l%kappa
!            dV_ll(7) = dV_ll(7) - rtemp1
!
!        end do
!    end do
!
!!-------------------------------CUT-OFF ENACTION-------------------------------
!! Don't forget the gamma!
!
!! sigma and its derivatives.
!    sigma_ll = sigma_ll * igamma1l
!        dsigma_ll(1,:) = (dsigma_ll(1,:) - sigma_ll*dgamma1l(1))*igamma1l
!        dsigma_ll(7,:) = (dsigma_ll(7,:) - sigma_ll*dgamma1l(7))*igamma1l
!
! ! The pair potential and its derivatives
!    V_ll = V_ll * pars_l%V0 * igamma2l
!        dV_ll(5) = - V_ll/pars_l%V0
!        dV_ll(6) = (dV_ll(6) * pars_l%V0/beta + V_ll*dgamma2l(6)) * igamma2l
!        dV_ll(7) = -(dV_ll(7) * pars_l%V0 + V_ll*dgamma2l(7)) * igamma2l
!
!!-----------------------------NEUTRAL SPHERE RADIUS----------------------------
!! The neutral sphere radius is the radius in which the entire density of the
!! atom is included.
!
!    s_l_ref = -log( sigma_ll*twelveth )/( beta*pars_l%eta2)
!
!    ! Derivatives for reference energy:
!    rn_ltemp = 1.0d0/(sigma_ll*pars_l%eta2*beta)
!        ! Derivative with respect to l
!        ds_l_l_ref(1,:) = -s_l_ref/pars_l%eta2 &
!                      - dsigma_ll(1,:)*rn_ltemp
!        ds_l_l_ref(7,:) = -rn_ltemp*dsigma_ll(7,:)
!
!
!!----------------MIXED REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------
!! These contributions have to be substracted to account for the contributions
!! that were included twice.
!
!    rn_ltemp = exp( -pars_l%kappa * s_l_ref)
!    rtemp = -12.0d0 * pars_l%V0 * pars_l%kappa
!    vref_l_ref = 12.0d0 * pars_l%V0 * sum(rn_ltemp)
!    ! Derivative with respect to l
!        dvref_l_l_ref(1) = rtemp*sum(rn_ltemp*ds_l_l_ref(1,:))
!!        dvref_l_l(2) = rtemp*sum(rn_ltemp*ds_l_l_ref(2,:))
!        dvref_l_l_ref(5) = vref_l_ref/pars_l%V0
!        dvref_l_l_ref(6) = - 12.0d0 * pars_l%V0 * sum(rn_ltemp *s_l_ref)
!        dvref_l_l_ref(7) = -rtemp*sum(rn_ltemp*ds_l_l_ref(7,:))
!
!
!!------------------------------------------------------------------------------
!!                           CALCULATING THE ENERGY
!!                           ======================
!!------------------------------------------------------------------------------
!
!
!!---------------------------COHESIVE ENERGY FUNCTION---------------------------
!! Calculates and sums the contributions to the cohesive energy for both lattice
!! and particle.
!
!    Ecoh_ref = sum( (1.0d0 + pars_l%lambda*s_l_ref) * exp(-pars_l%lambda * s_l_ref)-1.0d0 ) &
!          * pars_l%E0
!
!
!
!    rn_ltemp=-pars_l%E0*pars_l%lambda*s_l_ref*exp(-pars_l%lambda*s_l_ref)
!        ! Derivative with respect to l
!        dEcoh_ref(1) = sum(pars_l%lambda*rn_ltemp*ds_l_l_ref(1,:))
!
!        !dEcoh_l(2) = sum(pars_l%lambda*rn_ltemp*ds_l_l_ref(2,:))
!
!        dEcoh_ref(3) = sum( (1.0d0 + pars_l%lambda*s_l_ref) * exp(-pars_l%lambda * s_l_ref)-1.0d0 )
!        dEcoh_ref(4) = sum(s_l_ref*rn_ltemp)
!        dEcoh_ref(7) = sum(pars_l%lambda*rn_ltemp*ds_l_l_ref(7,:))
!
!
!!-------------------------------OVERALL ENERGY---------------------------------
!! Summation over all contributions.
!! Reference energy
!     E_ref = Ecoh_ref - V_ll + 0.5d0 * vref_l_ref
!    ! Derivative with respect to l
!    dE_ref(1) = dEcoh_ref(1) + 0.5d0*dvref_l_l_ref(1)
!!    dE_ref(2) = dEcoh_ref(2) - 0.5*dvref_l_l_ref(2)
!    dE_ref(3) = dEcoh_ref(3)
!    dE_ref(4) = dEcoh_ref(4)
!    dE_ref(5) = dV_ll(5) + 0.5d0*dvref_l_l_ref(5)
!    dE_ref(6) = dV_ll(6) + 0.5d0*dvref_l_l_ref(6)
!    dE_ref(7) = dEcoh_ref(7) + dV_ll(7) + 0.5d0*dvref_l_l_ref(7)
!    E_dref = E_ref
!
!
!end subroutine emt_fit_l
!
!subroutine emt_fit_init (a_lat, cell, r0_lat, n_l, pars_p, pars_l, E_dref, dE_ref)
!!subroutine emt_fit (cell, a_lat, n_l, r_lat, r_part, pars_p, pars_l, energy, denergy_l, denergy_p)
!!
!! Purpose:
!! ========
!!       emt_fit_init calculates the reference energy according to the effective
!!       medium theory including the derivatives for the reference energy.
!
!! Input variables are:
!!           lattice : lattice atoms
!!           r0      : equilibrium positions of lattice atoms
!!           pars_p  : emt-parameters of particle
!!           pars_l  : emt-parameters of lattice
!! Output variables are:
!!           energy  : emt-energy
!!
!! General Remarks To Fit:
!! =======================
!! If you are using the potential for fitting, we recommand to refrain from
!! fitting so. so_l is directly connected to the lattice constant. If you change
!! it, you'll automatically change the lattice constant of your slab.
!! so_l = a_lattice /(beta *sqrt_2)
!! so_p is related to the DFT density, so, you should not change it, either, if
!! you want to reproduce the DFT density.
!!
!implicit none
!
!!------------------------------------------------------------------------------
!!                                   PREAMBLE
!!                                   ========
!!------------------------------------------------------------------------------
!
!! declare variables and parameters that are passed down from the program
!
!    real(8), dimension(3,6), intent(in)   :: cell  ! cell matrix and its inverse
!    real(8), intent(in)                 :: a_lat   ! lattice constant of lattice
!    integer, intent(in)                 :: n_l      ! number of lattice atoms
!    real(8), dimension(:,:), intent(in) :: r0_lat    ! positions of lattice atoms
!    type(EMTparms), intent(inout)       :: pars_p   ! parameters of particle
!    type(EMTparms), intent(inout)       :: pars_l   ! parameters of lattice atoms
!    real(8), intent(out)                :: E_dref   ! calc. reference energy
!    real(8), dimension(7), intent(out)  :: dE_ref ! derivatives
!                                                        ! with respect to to
!                                                        ! eta2, followed by no,
!                                                        ! eo, lambda, vo, kappa
!                                                        ! and so.
!
!! declare the variables that appear in the subroutine
!
!    integer :: i,j, k                   ! running parameter
!    real(8) :: r                        ! distance
!    real(8) :: rcut, rr, acut           ! values to calculate cut-off
!    real(8) :: igamma1l, igamma2l       ! inverse gamma for lattice atoms
!    real(8) :: theta                    ! variable for cut-off calculation
!    real(8), dimension(n_l) :: sigma_ll ! contribution to ns, l only
!    real(8), dimension(n_l) :: s_l
!    real(8) ::  V_ll                    ! Pair potential contributions
!    real(8) :: vref_l
!    real(8), dimension(3) :: xl         ! for cal. cut-off
!    real(8), dimension(3) :: rnnl       ! next neighbour distance for cut-off
!    real(8), dimension(3) :: rnndbeta_l ! divide rnn by beta
!    real(8) :: betas0_l                 ! beta * s0= for l and p
!    real(8) :: kappadbeta_l             ! beta * kappa for l and p
!    real(8), dimension(3) :: r3temp, r3temp1     ! temporary array variable
!    real(8) :: rtemp, rtemp1                    ! temporary variable
!    real(8), dimension(n_l) :: rn_ltemp
!    real(8) :: Ecoh, energy
!
!    real(8),parameter,dimension(3) :: r0_part=(/0.d0,0d0,6.0d0/)
!    real(8) :: igamma1p, igamma2p       ! inverse gamma for particle
!    real(8) :: chilp, chipl             ! mixing between lattice (l) and particle (p)
!    real(8) :: sigma_pl                 ! mixed contribution to neutral sphere,
!    real(8) :: s_p                      ! neutral sphere radius particle
!    real(8), dimension(n_l) :: sigma_lp ! mixed contribution to ns
!    real(8) :: V_pl, V_lp               ! Pair potential contributions
!    real(8) :: vref_p                   ! reference pair pot. contrib.
!    real(8), dimension(3) :: xp         ! for cal. cut-off
!    real(8), dimension(3) :: rnnp       ! next neighbour distance for cut-off
!    real(8), dimension(3) :: rnndbeta_p ! divide rnn by beta
!    real(8) :: betas0_p                 ! beta * s0= for l and p
!    real(8) :: kappadbeta_p             ! beta * kappa for l and p
!
!
!
!!-----------------------DECLARE VARIABLES FOR DERIVATIVES----------------------
!! Variables and Arrays for partial derivatives
!! Apart from chi, all derivatives are 7 long. The first place denotes the
!! derivative with respect to to eta2, followed by no, eo, lambda, vo, kappa and so.
!! The general notation is:
!! e.g. dV_lp_l(1) : the derivative of V_lp with respect to eta2_l.
!    real(8), dimension(2) :: dchilp, dchipl     ! First element: p, then l
!    real(8), dimension(7) :: dgamma1l, dgamma2l
!    real(8), dimension(7) :: dgamma1p, dgamma2p
!    real(8), dimension(7,n_l) :: dsigma_ll
!    real(8), dimension(7,n_l) :: dsigma_lp_l
!    real(8), dimension(7,n_l) :: dsigma_lp_p
!    real(8), dimension(7) :: dsigma_pl_l, dsigma_pl_p
!    real(8), dimension(7) :: dV_ll
!    real(8), dimension(7) :: dV_lp_l, dV_lp_p
!    real(8), dimension(7) :: dV_pl_l, dV_pl_p
!    real(8), dimension(7,n_l) :: ds_l_l, ds_l_p
!    real(8), dimension(7) :: ds_p_l, ds_p_p
!    real(8), dimension(7) :: dvref_l_l, dvref_l_p
!    real(8), dimension(7) :: dvref_p_l, dvref_p_p
!    real(8), dimension(7) :: dEcoh_l, dEcoh_p
!
!!______________________________________________________________________________
!
!!----------------------VALUES OF FREQUENT USE ---------------------------------
!! definition of a few values that appear frequently in calculation
!    betas0_l = beta * pars_l%s0
!    betas0_p = beta * pars_p%s0
!    kappadbeta_l = pars_l%kappa / beta
!    kappadbeta_p = pars_p%kappa / beta
!
!! 'coupling' parameters between p and l
!! In the original paper by Jacobsen et al (SS 366(1996),394), chi also contains
!! exponential contributions. We conclude that those exponential contributions
!! can be multiplied into the various no and are thus part of the fitting
!! parameters.
!! derivatives: (1) derivative over nop, (2) over nol
!    dchilp(1) = 1.0d0 / pars_l%n0         ! d chilp / d nop
!    dchipl(2) = 1.0d0 / pars_p%n0         ! d chipl / d nol
!
!    chilp = pars_p%n0 * dchilp(1)
!    chipl = pars_l%n0 * dchipl(2)
!
!    dchipl(1) = - chipl * dchipl(2)     ! d chipl / d nop
!    dchilp(2) = - chilp * dchilp(1)     ! d chipl / d nol
!
!!------------------------------------------------------------------------------
!!                                  CUT-OFF
!!                                  =======
!!------------------------------------------------------------------------------
!! We use the distance to the next-next nearest neighbours as cut-off.
!! FOR FUTURE REVISION:
!!            cut-off should be defined via lattice constant _AND_ changeable.
!
!    rcut = a_lat * sqrt_3 * isqrt_2
!    rr = 4.0d0 * rcut / (sqrt_3 + 2.0d0)
!    acut = 9.21024d0/(rr -rcut) ! ln(10000)
!
!
!! Distances to the considered neighbours
!    rnnl(1) = betas0_l
!    rnnl(2) = rnnl(1) * sqrt_2
!    rnnl(3) = rnnl(1) * sqrt_3
!    rnnp(1) = betas0_p
!    rnnp(2) = rnnp(1) * sqrt_2
!    rnnp(3) = rnnp(1) * sqrt_3
!
!    xl = b * twelveth / (1.0d0 + exp(acut*(rnnl-rcut)))
!    xp = b * twelveth/ (1.0d0 + exp(acut*(rnnp-rcut)))
!
!
!!-----------------------------------GAMMA--------------------------------------
!! Gamma enforces the cut-off together with theta (see below)
!! Gamma is defined as inverse.
!! The derivative is not defined as the inverse and formed for each gamma
!! individually.
!    r3temp = rnnl - betas0_l
!
!    r3temp1 = xl*exp(- pars_l%eta2*r3temp)
!    igamma1l = 1.0d0 / sum(r3temp1)
!    dgamma1l = 0.d0
!    dgamma1l(1) = - sum(r3temp*r3temp1)
!    dgamma1l(7) = sum(r3temp1)*pars_l%eta2*beta
!
!
!    r3temp1 = xl*exp(-kappadbeta_l * r3temp)
!    igamma2l = 1.0d0 / sum(r3temp1)
!    dgamma2l = 0.
!    dgamma2l(6) = - sum(r3temp*r3temp1) / beta
!    dgamma2l(7) = sum(r3temp1)*pars_l%kappa
!
!
!    r3temp = rnnp-betas0_p
!
!    r3temp1 = xp*exp(- pars_p%eta2*r3temp)
!    igamma1p = 1.0d0 / sum(r3temp1)
!    dgamma1p = 0.
!    dgamma1p(1) = - sum(r3temp*r3temp1)
!    dgamma1p(7) = sum(r3temp1)*pars_p%eta2*beta
!
!
!    r3temp1 = xp*exp(-kappadbeta_p * r3temp)
!    igamma2p = 1.0d0 / sum(r3temp1)
!    dgamma2p = 0.
!    dgamma2p(6) = - sum(r3temp*r3temp1) / beta
!    dgamma2p(7) = sum(r3temp1)*pars_p%kappa
!
!
!
!
!!------------------------------------------------------------------------------
!!                          INDIVIDUAL CONTRIBUTIONS
!!                          ========================
!!------------------------------------------------------------------------------
!! The values for the sums and the derivatives are set to zero.
!
!    sigma_ll = 0
!    dsigma_ll = 0
!    sigma_pl = 0
!    dsigma_lp_l=0
!    dsigma_lp_p=0
!    dsigma_pl_l=0
!    dsigma_pl_p=0
!    V_ll = 0
!    dV_ll = 0
!    V_lp = 0
!    V_pl = 0
!    dV_lp_l = 0
!    dV_lp_p = 0
!    dV_pl_l = 0
!    dV_pl_p = 0
!    ds_l_l = 0
!    ds_l_p = 0
!    ds_p_l = 0
!    ds_p_p = 0
!    dvref_l_l=0
!    dvref_l_p=0
!    dvref_p_l=0
!    dvref_p_p=0
!    dEcoh_l = 0
!    dEcoh_p = 0
!    dE_ref = 0
!
!    do i = 1,n_l
!        do j = i+1, n_l
!
!
!        !-----------------PERIODIC BOUNDARY CONDITIONS LATTICE-----------------
!        ! Because we want them.
!
!            r3temp(1) = r0_lat(1,i)-r0_lat(1,j)
!            r3temp(2) = r0_lat(2,i)-r0_lat(2,j)
!            r3temp(3) = r0_lat(3,i)-r0_lat(3,j)
!
!
!        ! transform distances into direct coordinates
!            r3temp= matmul(cell(1:3,4:6),r3temp)
!
!            r3temp(1)=r3temp(1)-Anint(r3temp(1))
!            r3temp(2)=r3temp(2)-Anint(r3temp(2))
!            r3temp(3)=r3temp(3)-Anint(r3temp(3))
!
!            r3temp=matmul(cell(1:3,1:3),r3temp)
!
!            r =  sqrt(sum(r3temp**2))
!
!
!        !---------------------------THETA LATTICE------------------------------
!        ! Theta enforces the cut-off together with gamma (see above). This
!        ! function enacts cutoff by reducing contributions of atoms outside the
!        ! cut-off to zero.
!
!            theta = 1.0d0 / (1 + exp( acut * (r - rcut) ) )
!
!
!
!
!
!        !----------------------------SIGMA LATTICE-----------------------------
!        ! Sigma is a contribution to the neutral sphere radius.
!        ! It is a list in which for each lattice atom, the contributions of the
!        ! others are summed up. To enforce the cut-off, it will be later
!        ! corrected by gamma.
!        ! sigma_pp does not exist because there is only a single particle.
!
!            rtemp = theta*exp(-pars_l%eta2 * (r - betas0_l) )
!            sigma_ll(i) = sigma_ll(i) + rtemp
!            sigma_ll(j) = sigma_ll(j) + rtemp
!
!            rtemp1 = rtemp*(r - betas0_l)
!            dsigma_ll(1,i) = dsigma_ll(1,i) - rtemp1
!            dsigma_ll(1,j) = dsigma_ll(1,j) - rtemp1
!
!            rtemp1 = rtemp*beta*pars_l%eta2
!            dsigma_ll(7,i) = dsigma_ll(7,i) + rtemp1
!            dsigma_ll(7,j) = dsigma_ll(7,j) + rtemp1
!
!
!
!        !-----------------------PAIR POTENTIAL LATTICE-------------------------
!        ! For the lattice only.
!        ! Will later be subjected to gamma to complete the cut-off.
!        ! The particle has no pair potential contribution since there is only
!        ! one and thus does not have a partner to interact with.
!
!            rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
!            V_ll = V_ll + rtemp
!
!            rtemp1 = rtemp*(r - betas0_l)  ! müsste das hier nicht r/beta-sol sein?
!            dV_ll(6) = dV_ll(6) + rtemp1
!
!            rtemp1 = rtemp*pars_l%kappa
!            dV_ll(7) = dV_ll(7) - rtemp1
!
!        end do
!
!
!    !-----------------PERIODIC BOUNDARY CONDITIONS LATTICE-----------------
!        ! Because we want them.
!
!            r3temp(1) = r0_lat(1,i)-r0_part(1)
!            r3temp(2) = r0_lat(2,i)-r0_part(2)
!            r3temp(3) = r0_lat(3,i)-r0_part(3)
!
!
!        ! transform distances into direct coordinates
!            r3temp= matmul(cell(1:3,4:6),r3temp)
!
!            r3temp(1)=r3temp(1)-Anint(r3temp(1))
!            r3temp(2)=r3temp(2)-Anint(r3temp(2))
!            r3temp(3)=r3temp(3)-Anint(r3temp(3))
!
!            r3temp=matmul(cell(1:3,1:3),r3temp)
!
!            r =  sqrt(sum(r3temp**2))
!
!
!    !----------------------------THETA PARTICLE--------------------------------
!
!        theta = 1.0d0 / (1.0d0 + exp( acut * (r - rcut) ) )
!
!    !-------------------------------MIXED SIGMA--------------------------------
!    ! Contributions of both particle and lattice to neutral sphere radius
!    ! To fully include the cut-off, we correct them later by gamma.
!    ! Each of the mixed sigmas depends on both l and p parameters. Not all
!    ! contributions need to be under the loop.
!
!        sigma_lp(i) = theta*exp(-pars_p%eta2 * (r - betas0_p) )
!        rtemp = theta*exp(-pars_l%eta2 * (r - betas0_l) )
!        sigma_pl = sigma_pl + rtemp
!
!        dsigma_lp_p(1,i) = -(r-betas0_p)*sigma_lp(i)
!
!        dsigma_pl_l(1) = dsigma_pl_l(1) - (r - betas0_l)*rtemp
!        dsigma_pl_l(7) = pars_l%eta2*beta*sigma_pl
!
!
!    !--------------------MIXED PAIR POTENTIAL CONTRIUBUTION--------------------
!
!        rtemp = theta*exp(-kappadbeta_p * (r - betas0_p))
!        V_lp= V_lp + rtemp
!        dV_lp_p(6) = dV_lp_p(6) + rtemp*(r - betas0_p)
!        dV_lp_p(7) = dV_lp_p(7) + rtemp*pars_p%kappa
!
!        rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
!        V_pl = V_pl + rtemp
!        dV_pl_l(6) = dV_pl_l(6) + rtemp*(r - betas0_l)
!        dV_pl_l(7) = dV_pl_l(7) + rtemp*pars_l%kappa
!
!
!    end do
!
!
!!-------------------------------CUT-OFF ENACTION-------------------------------
!! Don't forget the gamma!
!
!! sigma and its derivatives.
!    sigma_ll = sigma_ll * igamma1l
!        dsigma_ll(1,:) = (dsigma_ll(1,:) - sigma_ll*dgamma1l(1))*igamma1l
!        dsigma_ll(7,:) = (dsigma_ll(7,:) - sigma_ll*dgamma1l(7))*igamma1l
!
!    sigma_lp = sigma_lp * igamma1l
!        ! Derivative with respect to l
!        dsigma_lp_l(1,:) = - sigma_lp*igamma1l*dgamma1l(1)
!        dsigma_lp_l(7,:) = - sigma_lp*igamma1l*dgamma1l(7)
!        ! Derivative with respect to p
!        dsigma_lp_p(1,:) = dsigma_lp_p(1,:)*igamma1l
!        dsigma_lp_p(7,:) = sigma_lp(:)*beta*pars_p%eta2
!
!    sigma_pl = sigma_pl * igamma1p
!        ! Derivative with respect to l
!        dsigma_pl_l(1) = dsigma_pl_l(1)*igamma1p
!        dsigma_pl_l(7) = dsigma_pl_l(7)*igamma1p
!        ! Derivative with respect to p
!        dsigma_pl_p(1) = - sigma_pl*igamma1p*dgamma1p(1)
!        dsigma_pl_p(7) = - sigma_pl*igamma1p*dgamma1p(7)
!
!! The pair potential and its derivatives
!    V_ll = V_ll * pars_l%V0 * igamma2l
!        dV_ll(5) = - V_ll/pars_l%V0
!        dV_ll(6) = (dV_ll(6) * pars_l%V0/beta + V_ll*dgamma2l(6)) * igamma2l
!        dV_ll(7) = -(dV_ll(7) * pars_l%V0 + V_ll*dgamma2l(7)) * igamma2l
!
!    V_lp = V_lp *chilp * pars_l%V0 * igamma2l
!        ! Derivative with respect to l
!        dV_lp_l(2) = V_lp/pars_l%n0
!        dV_lp_l(5) = - V_lp/pars_l%V0
!        dV_lp_l(7) = V_lp*igamma2l             ! this one is temporary
!        dV_lp_l(6) = dV_lp_l(7)*dgamma2l(6)
!        dV_lp_l(7) = dV_lp_l(7)*dgamma2l(7)
!        ! Derivative with respect to p
!        dV_lp_p(2) = -V_lp/pars_p%n0
!        dV_lp_p(6) = dV_lp_p(6)*igamma2l*chilp * pars_l%V0 / beta
!        dV_lp_p(7) = -dV_lp_p(7)*igamma2l*chilp * pars_l%V0
!
!    V_pl = V_pl * pars_p%V0 * igamma2p * chipl
!        ! Derivative with respect to l
!        dV_pl_l(2) = - V_pl/pars_l%n0
!        dV_pl_l(6) = dV_pl_l(6)*igamma2p*pars_p%V0*chipl/beta
!        dV_pl_l(7) = -dV_pl_l(7)*igamma2p*chipl * pars_p%V0
!        ! Derivative with respect to p
!        dV_pl_p(2) = V_pl/pars_p%n0
!        dV_pl_p(5) = - V_pl/pars_p%V0
!        dV_pl_p(6) = V_pl*igamma2p*dgamma2p(6)
!        dV_pl_p(7) = V_pl*igamma2p*dgamma2p(7)
!
!
!
!!-----------------------------NEUTRAL SPHERE RADIUS----------------------------
!! The neutral sphere radius is the radius in which the entire density of the
!! atom is included.
!
!    rn_ltemp = sigma_ll + chilp*sigma_lp
!    s_l = -log( rn_ltemp * twelveth ) &
!            / ( beta * pars_l%eta2)
!
!    rn_ltemp = 1.0d0/(rn_ltemp*pars_l%eta2*beta)
!        ! Derivative with respect to l
!        ds_l_l(1,:) = -s_l/pars_l%eta2 &
!                      - (dsigma_ll(1,:)+chilp*dsigma_lp_l(1,:))*rn_ltemp
!        ds_l_l(2,:) = -sigma_lp*rn_ltemp*dchilp(2)
!        ds_l_l(7,:) = -rn_ltemp*(dsigma_ll(7,:)+chilp*dsigma_lp_l(7,:))
!        ! Derivative with respect to p
!        ds_l_p(1,:) = - rn_ltemp*chilp*dsigma_lp_p(1,:)
!        ds_l_p(2,:) = -sigma_lp*rn_ltemp*dchilp(1)
!        ds_l_p(7,:) = -rn_ltemp*chilp*dsigma_lp_p(7,:)
!
!    rtemp=1.0d0/(beta*pars_p%eta2)
!    s_p  = -log( sigma_pl * chipl * twelveth) *rtemp
!        ! Derivative with respect to l
!        ds_p_l(1) = - rtemp /sigma_pl*dsigma_pl_l(1)
!        ds_p_l(2) = -rtemp / pars_l%n0
!        ds_p_l(7) = - rtemp / sigma_pl*dsigma_pl_l(7)
!        ! Derivative with respect to p
!        ds_p_p(1) = -s_p/pars_p%eta2 - rtemp/(sigma_pl)*dsigma_pl_p(1)
!        ds_p_p(2) = rtemp / pars_p%n0
!        ds_p_p(7) = - rtemp/sigma_pl*dsigma_pl_p(7)
!
!
!
!!----------------MIXED REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------
!! These contributions have to be substracted to account for the contributions
!! that were included twice.
!
!    rn_ltemp = exp( -pars_l%kappa * s_l)
!    rtemp = -12.0d0 * pars_l%V0 * pars_l%kappa
!    vref_l = 12.0d0 * pars_l%V0 * sum(rn_ltemp)
!    ! Derivative with respect to l
!        dvref_l_l(1) = rtemp*sum(rn_ltemp*ds_l_l(1,:))
!        dvref_l_l(2) = rtemp*sum(rn_ltemp*ds_l_l(2,:))
!        dvref_l_l(5) = vref_l/pars_l%V0
!        dvref_l_l(6) = - 12.0d0 * pars_l%V0 * sum(rn_ltemp *s_l)
!        dvref_l_l(7) = rtemp*sum(rn_ltemp*ds_l_l(7,:))
!        ! Derivative with respect to p
!        dvref_l_p(1) = rtemp*sum(rn_ltemp*ds_l_p(1,:))
!        dvref_l_p(2) = rtemp*sum(rn_ltemp*ds_l_p(2,:))
!        dvref_l_p(7) = rtemp*sum(rn_ltemp*ds_l_p(7,:))
!
!    rtemp = -pars_p%kappa
!    vref_p = 12.0d0 * pars_p%V0 * exp( -pars_p%kappa * s_p)
!        ! Derivative with respect to p
!        dvref_p_p(1) = rtemp*vref_p*ds_p_p(1)
!        dvref_p_p(2) = rtemp*vref_p*ds_p_p(2)
!        dvref_p_p(5) = vref_p/pars_p%V0
!        dvref_p_p(6) = - vref_p*s_p
!        dvref_p_p(7) = rtemp*vref_p*ds_p_p(7)
!        ! Derivative with respect to l
!        dvref_p_l(1) = rtemp*vref_p*ds_p_l(1)
!        dvref_p_l(2) = rtemp*vref_p*ds_p_l(2)
!        dvref_p_l(7) = rtemp*vref_p*ds_p_l(7)
!
!!------------------------------------------------------------------------------
!!                           CALCULATING THE ENERGY
!!                           ======================
!!------------------------------------------------------------------------------
!
!
!!---------------------------COHESIVE ENERGY FUNCTION---------------------------
!! Calculates and sums the contributions to the cohesive energy for both lattice
!! and particle.
!
!    Ecoh = sum( (1.0d0 + pars_l%lambda*s_l) * exp(-pars_l%lambda * s_l)-1.0d0 ) &
!          * pars_l%E0 &
!          + (1.0d0 + pars_p%lambda*s_p) * exp(-pars_p%lambda * s_p)* pars_p%E0
!
!
!    rn_ltemp=-pars_l%E0*pars_l%lambda*s_l*exp(-pars_l%lambda*s_l)
!    rtemp = -pars_p%lambda*s_p*pars_p%E0*exp(-pars_p%lambda*s_p)
!        ! Derivative with respect to l
!        dEcoh_l(1) = sum(pars_l%lambda*rn_ltemp*ds_l_l(1,:))&
!                    +pars_p%lambda*rtemp*ds_p_l(1)
!        dEcoh_l(2) = sum(pars_l%lambda*rn_ltemp*ds_l_l(2,:))&
!                    +pars_p%lambda*rtemp*ds_p_l(2)
!        dEcoh_l(3) = sum( (1.0d0 + pars_l%lambda*s_l) * exp(-pars_l%lambda * s_l)-1.0d0 )
!        dEcoh_l(4) = sum(s_l*rn_ltemp)
!        dEcoh_l(7) = sum(pars_l%lambda*rn_ltemp*ds_l_l(7,:))&
!                    +pars_p%lambda*rtemp*ds_p_l(7)
!
!    rtemp = -pars_p%lambda*s_p*pars_p%E0*exp(-pars_p%lambda*s_p)
!        ! Derivative with respect to p
!        dEcoh_p(1) = sum(pars_l%lambda*rn_ltemp*ds_l_p(1,:))&
!                    +pars_p%lambda*rtemp*ds_p_p(1)
!        dEcoh_p(2) = sum(pars_l%lambda*rn_ltemp*ds_l_p(2,:))&
!                    +pars_p%lambda*rtemp*ds_p_p(2)
!        dEcoh_p(3) = (1.0d0 + pars_p%lambda*s_p) * exp(-pars_p%lambda * s_p)
!        dEcoh_p(4) = s_p*rtemp
!        dEcoh_p(7) = sum(pars_l%lambda*rn_ltemp*ds_l_p(7,:))&
!                    +rtemp*pars_p%lambda*ds_p_p(7)
!
!
!!-------------------------------OVERALL ENERGY---------------------------------
!! Summation over all contributions.
!! Overall energy
!    energy = Ecoh - V_ll - 0.5d0* ( V_lp + V_pl - vref_l - vref_p)
!    E_dref = energy
!
!    ! Derivative with respect to l
!    dE_ref(8) = dEcoh_l(1) + 0.5d0*( dvref_l_l(1)+dvref_p_l(1))
!    dE_ref(9) = dEcoh_l(2) &
!                    +0.5d0*(dV_pl_l(2)+dvref_p_l(2)+dV_lp_l(2)+dvref_l_l(2))
!    dE_ref(10) = dEcoh_l(3)
!    dE_ref(11) = dEcoh_l(4)
!    dE_ref(12) = dV_ll(5) + 0.5d0*(dvref_l_l(5)+dV_lp_l(5))
!    dE_ref(13) = dV_ll(6) + 0.5d0*( dV_lp_l(6) + dV_pl_l(6) + dvref_l_l(6))
!    dE_ref(14) = dEcoh_l(7) + dV_ll(7) &
!                   + 0.5d0*(dV_lp_l(7)+dV_pl_l(7)+dvref_l_l(7)+dvref_p_l(7))
!
!    ! Derivative with respect to p (no correction by dEref since those do not
!    ! contain any p-contribution)
!    dE_ref(1) = dEcoh_p(1) + 0.5d0*(dvref_l_p(1)+dvref_p_p(1))
!    dE_ref(2) = dEcoh_p(2) &
!                   + 0.5d0*(dV_pl_p(2)+dV_lp_p(2)+dvref_l_p(2)+dvref_p_p(2))
!    dE_ref(3) = dEcoh_p(3)
!    dE_ref(4) = dEcoh_p(4)
!    dE_ref(5) = 0.5d0*(dV_pl_p(5) + dvref_p_p(5))
!    dE_ref(6) = 0.5d0*(dV_lp_p(6)+dV_pl_p(6)+dvref_p_p(6))
!    dE_ref(7) = dEcoh_p(7) &
!                   + 0.5d0*(dV_lp_p(7)+dV_pl_p(7)+dvref_l_p(7)+dvref_p_p(7))
!
!
!
!
!end subroutine emt_fit_init
!
!subroutine emt_init (a_lat, cell, n_l, r0_lat, pars_p, pars_l, energy)
!!subroutine emt_init (cell_in, n_l_in, r0_lat_in, pars_p, pars_l, energy)
!!
!! Purpose:
!!           Here, the fitting procedure is just implemented and the reference
!!           energy calculated.
!!
!implicit none
!
!!------------------------------------------------------------------------------
!!                                   PREAMBLE
!!                                   ========
!!------------------------------------------------------------------------------
!
!! declare variables and parameters that are passed down from the program
!
!    real(8), intent(in)                 :: a_lat
!    real(8), dimension(3,6), intent(in)   :: cell         ! cell-matrix(1-3)
!                                                            !and inverse matrix(4-6)
!    integer, intent(in)                 :: n_l          ! number of lattice atoms
!    real(8), dimension(:,:), intent(in) :: r0_lat       ! positions of lattice atoms
!    type(EMTparms), intent(inout)       :: pars_p       ! parameters of particle
!    type(EMTparms), intent(inout)       :: pars_l       ! parameters of lattice atoms
!    real(8), intent(out)                :: energy       ! calculated reference energy
!
!! declare the variables that appear in the subroutine
!
!    integer :: i,j, k                   ! running parameter
!    real(8) :: r                        ! distance
!    real(8) :: rcut, rr, acut           ! values to calculate cut-off
!    real(8) :: igamma1l, igamma2l       ! inverse gamma for lattice atoms
!    real(8) :: theta                    ! variable for cut-off calculation
!    real(8) :: chilp                    ! mixing between lattice (l) and particle (p)
!    real(8), dimension(n_l) :: sigma_ll ! contribution to ns, l only
!    real(8), dimension(n_l) :: s_l      ! neutral sphere radius lattice atoms
!    real(8) :: V_ll                     ! Pair potential contributions
!    real(8) :: vref_l                   ! reference pair pot. contrib.
!    real(8) :: Ecoh                     ! cohesive energy of part & lattice atoms
!    real(8), dimension(3) :: xl         ! for cal. cut-off
!    real(8), dimension(3) :: rnnl       ! next neighbour distance for cut-off
!    real(8) :: betas0_l                 ! beta * s0= for l and p
!    real(8) :: kappadbeta_l             ! beta * kappa for l and p
!    real(8), dimension(3) :: r3temp     ! temporary array variable
!    real(8) :: rtemp                    ! temporary variable
!
!    real(8),parameter,dimension(3) :: r0_part=(/0.d0,0d0,6.0d0/)
!    real(8) :: igamma1p, igamma2p       ! inverse gamma for lattice atoms
!    real(8) :: chipl                    ! mixing between lattice (l) and particle (p)
!    real(8) :: sigma_pl                 ! mixed contribution to neutral sphere,
!    real(8) :: s_p                      ! neutral sphere radius particle
!    real(8), dimension(n_l) :: sigma_lp ! mixed contribution to ns
!    real(8) :: V_pl, V_lp               ! Pair potential contributions
!    real(8) :: vref_p                   ! reference pair pot. contrib.
!    real(8), dimension(3) :: xp         ! for cal. cut-off
!    real(8), dimension(3) :: rnnp       ! next neighbour distance for cut-off
!    real(8) :: betas0_p                 ! beta * s0= for l and p
!    real(8) :: kappadbeta_p             ! beta * kappa for l and p
!
!energy = 0.0d0
!
!!----------------------VALUES OF FREQUENT USE ---------------------------------
!! definition of a few values that appear frequently in calculation
!    betas0_l = beta * pars_l%s0
!    betas0_p = beta * pars_p%s0
!    kappadbeta_l = pars_l%kappa / beta
!    kappadbeta_p = pars_p%kappa / beta
!
!! 'coupling' parameters between p and l
!    chilp = pars_p%n0 / pars_l%n0
!    chipl = 1.0d0 / chilp
!
!
!!------------------------------------------------------------------------------
!!                                  CUT-OFF
!!                                  =======
!!------------------------------------------------------------------------------
!! We use the distance to the next-next nearest neighbours as cut-off.
!! FOR FUTURE REVISION:
!!            cut-off should be defined via lattice constant _AND_ changeable.
!
!!    rcut = betas0_l * sqrt_3
!!    rr = 4 * rcut / (sqrt_3 + 2)
!!    acut = 9.21024/(rr -rcut) ! ln(10000)
!
!    rcut = a_lat * sqrt_3 * isqrt_2
!    rr = 4.0d0 * rcut / (sqrt_3 + 2.0d0)
!    acut = 9.21024d0/(rr -rcut) ! ln(10000)
!
!! Distances to the considered neighbours
!    rnnl(1) = betas0_l
!    rnnl(2) = rnnl(1) * sqrt_2
!    rnnl(3) = rnnl(1) * sqrt_3
!    rnnp(1) = betas0_p
!    rnnp(2) = rnnp(1) * sqrt_2
!    rnnp(3) = rnnp(1) * sqrt_3
!
!    xl = b * twelveth / (1.0d0 + exp(acut*(rnnl-rcut)))
!    xp = b * twelveth/ (1.0d0 + exp(acut*(rnnp-rcut)))
!
!!-----------------------------------GAMMA--------------------------------------
!! Gamma enforces the cut-off together with theta (see below)
!! Gamma is defined as inverse.
!    r3temp = rnnl-betas0_l
!    igamma1l = 1.0d0 / sum(xl*exp(-pars_l%eta2 * r3temp))
!    igamma2l = 1.0d0 /sum(xl*exp(-kappadbeta_l * r3temp))
!
!    r3temp = rnnp-betas0_p
!    igamma1p = 1.0d0 / sum(xp*exp(-pars_p%eta2 * r3temp))
!    igamma2p = 1.0d0 / sum(xp*exp(-kappadbeta_p * r3temp))
!
!!------------------------------------------------------------------------------
!!                          INDIVIDUAL CONTRIBUTIONS
!!                          ========================
!!------------------------------------------------------------------------------
!! The values for the sums are set to zero.
!
!    sigma_ll = 0.0d0
!    sigma_pl = 0.0d0
!    V_ll = 0.0d0
!    V_lp = 0.0d0
!    V_pl = 0.0d0
!    sigma_lp= 0.0d0
!    vref_l= 0.0d0
!    vref_p=0.0d0
!    Ecoh=0.0d0
!
!    do i = 1, n_l
!        do j = i+1, n_l
!
!
!        !-----------------PERIODIC BOUNDARY CONDITIONS LATTICE-----------------
!        ! Because we want them.
!
!            r3temp(1) = r0_lat(1,i)-r0_lat(1,j)
!            r3temp(2) = r0_lat(2,i)-r0_lat(2,j)
!            r3temp(3) = r0_lat(3,i)-r0_lat(3,j)
!
!            ! transform distances into direct coordinates
!            r3temp= matmul(cell(1:3,4:6),r3temp)
!
!            r3temp(1)=r3temp(1)-Anint(r3temp(1))
!            r3temp(2)=r3temp(2)-Anint(r3temp(2))
!            r3temp(3)=r3temp(3)-Anint(r3temp(3))
!
!            r3temp=matmul(cell(1:3,1:3),r3temp)
!
!            r =  sqrt(sum(r3temp**2))
!
!
!        !---------------------------THETA LATTICE------------------------------
!        ! Theta enforces the cut-off together with gamma (see above). This
!        ! function enacts cutoff by reducing contributions of atoms outside the
!        ! cut-off to zero.
!
!            theta = 1.0d0 / (1.0d0 + exp( acut * (r - rcut) ) )
!
!
!        !----------------------------SIGMA LATTICE-----------------------------
!        ! Sigma is a contribution to the neutral sphere radius.
!        ! It is a list in which for each lattice atom, the contributions of the
!        ! others are summed up. To enforce the cut-off, it will be later
!        ! corrected by gamma.
!        ! sigma_pp does not exist because there is only a single particle.
!
!            rtemp = theta*exp(-pars_l%eta2 * (r - betas0_l) )
!            sigma_ll(i) = sigma_ll(i) + rtemp
!            sigma_ll(j) = sigma_ll(j) + rtemp
!
!
!
!        !-----------------------PAIR POTENTIAL LATTICE-------------------------
!        ! For the lattice only.
!        ! Will later be subjected to gamma to complete the cut-off.
!        ! The particle has no pair potential contribution since there is only
!        ! one and thus does not have a partner to interact with.
!
!            rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
!            V_ll = V_ll + rtemp
!
!
!        end do
!
!    !-----------------PERIODIC BOUNDERY CONDITIONS PARTICLE--------------------
!
!        r3temp(1) = r0_lat(1,i)-r0_part(1)
!        r3temp(2) = r0_lat(2,i)-r0_part(2)
!        r3temp(3) = r0_lat(3,i)-r0_part(3)
!    ! transform distances into direct coordinates
!        r3temp= matmul(cell(1:3,4:6),r3temp)
!
!        r3temp(1)=r3temp(1)-Anint(r3temp(1))
!        r3temp(2)=r3temp(2)-Anint(r3temp(2))
!        r3temp(3)=r3temp(3)-Anint(r3temp(3))
!
!        r3temp=matmul(cell(1:3,1:3),r3temp)
!        r =  sqrt(sum(r3temp**2))
!
!
!
!    !----------------------------THETA PARTICLE--------------------------------
!
!        theta = 1.0d0 / (1.0d0 + exp( acut * (r - rcut) ) )
!
!
!    !-------------------------------MIXED SIGMA--------------------------------
!    ! Contributions of both particle and lattice to neutral sphere radius
!    ! To fully include the cut-off, we correct them later by gamma.
!
!        sigma_lp(i) = theta*exp(-pars_p%eta2 * (r - betas0_p) )
!
!        rtemp = theta*exp(-pars_l%eta2 * (r - betas0_l) )
!        sigma_pl = sigma_pl + rtemp
!
!
!    !--------------------MIXED PAIR POTENTIAL CONTRIUBUTION--------------------
!
!        rtemp = theta*exp(-kappadbeta_p * (r - betas0_p))
!        V_lp= V_lp + rtemp
!
!        rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
!        V_pl = V_pl + rtemp
!
!    end do
!
!!-------------------------------CUT-OFF ENACTION-------------------------------
!! Don't forget the gamma!
!
!    sigma_ll = sigma_ll * igamma1l
!    sigma_lp = sigma_lp * igamma1l
!    sigma_pl = sigma_pl * igamma1p
!
!    V_ll = V_ll * pars_l%V0 * igamma2l
!    V_lp = V_lp *chilp * pars_l%V0 * igamma2l
!    V_pl = V_pl * pars_p%V0 * igamma2p * chipl
!
!
!
!!-----------------------------NEUTRAL SPHERE RADIUS----------------------------
!! The neutral sphere radius is the radius in which the entire density of the
!! atom is included.
!
!    s_l = -log( (sigma_ll + chilp * sigma_lp) * twelveth ) &
!            / ( beta * pars_l%eta2)
!    s_p  = -log( sigma_pl * chipl * twelveth) / ( beta * pars_p%eta2)
!
!
!!----------------MIXED REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------
!! These contributions have to be substracted to account for the contributions
!! that were included twice.
!
!    vref_l = 12.0d0 * pars_l%V0 * sum( exp( -pars_l%kappa * s_l) )
!    vref_p = 12.0d0 * pars_p%V0 * exp( -pars_p%kappa * s_p)
!
!
!!------------------------------------------------------------------------------
!!                           CALCULATING THE ENERGY
!!                           ======================
!!------------------------------------------------------------------------------
!
!
!!---------------------------COHESIVE ENERGY FUNCTION---------------------------
!! Calculates and sums the contributions to the cohesive energy for both lattice
!! and particle.
!
!    Ecoh = sum( (1.0d0 + pars_l%lambda*s_l) * exp(-pars_l%lambda * s_l)-1.d0 ) &
!          * pars_l%E0 &
!          + (1.0d0 + pars_p%lambda*s_p) * exp(-pars_p%lambda * s_p)* pars_p%E0
!
!
!!-------------------------------OVERALL ENERGY---------------------------------
!! Summation over all contributions.
!
!
!    energy = Ecoh - V_ll - 0.50d0 * ( V_lp + V_pl - vref_l - vref_p)
!    Eref   = energy
!
!end subroutine emt_init

