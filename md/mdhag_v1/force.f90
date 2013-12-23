module force
    !
    ! Purpose:
    !           Gather all potentials to calculate the forces for MD
    !           MAY THE FORCE BE WITH YOU!
    ! Definition of Zero of Energy:
    !                               The particle is at infinite distance to the
    !                               lattice.
    !
    ! Date          Author          History of Revison
    ! ====          ======          ==================
    ! 09.10.2013    Sascha&Svenja   Original
    !
    use atom_class
    use md_init

    implicit none
    save

    real(8), private, parameter                 :: beta    = 1.8093997906d0
    real(8), private, parameter                 :: twelfth= 0.0833333333333333d0
    integer, private, dimension(3), parameter   :: b       = (/12, 6, 24/)

contains


subroutine pes(slab, teilchen, Epot)

    implicit none

    type(atoms), intent(inout)    :: teilchen, slab
    real(8),     intent(out)      :: Epot

    integer :: i,j

    real(8) :: betas0_l, kappadbeta_l, chipl
    real(8) :: betas0_p, kappadbeta_p, chilp
    real(8) :: r, rcut, rr, acut, theta, rtemp
    real(8) :: igamma1p, igamma2p, igamma1l, igamma2l
    real(8) :: V_pl, V_lp, V_ll, V_pp, vref_l, vref_p, Ecoh

    real(8), dimension(3) :: rnnl, rnnp         ! nnn-distances
    real(8), dimension(3) :: xl, xp, r3temp
    real(8), dimension(3) :: dtheta

    real(8), dimension(:), allocatable :: sigma_ll, sigma_lp, s_l
    real(8), dimension(:), allocatable :: sigma_pp, sigma_pl, s_p

    real(8), dimension(:,:), allocatable :: dsigma_ll, dsigma_pp

!----------------------VALUES OF FREQUENT USE ---------------------------------

    ! beta * s0
    betas0_l = beta * pars_l(7)
    betas0_p = beta * pars_p(7)
    ! kappa / beta
    kappadbeta_l = pars_l(6) / beta
    kappadbeta_p = pars_p(6) / beta

    ! 'coupling' parameters between p and l
    chilp = pars_p(2) / pars_l(2)
    chipl = 1.0d0 / chilp

    ! Distances to the nearest, next-nearest and next-next-nearest neighbours
    rnnl(1) = betas0_l
    rnnl(2) = rnnl(1) * sqrt2
    rnnl(3) = rnnl(1) * sqrt3
    rnnp(1) = betas0_p
    rnnp(2) = rnnp(1) * sqrt2
    rnnp(3) = rnnp(1) * sqrt3

!------------------------------------------------------------------------------
!                                  CUT-OFF
!                                  =======
!------------------------------------------------------------------------------
! We use the distance to the next-next-nearest neighbours as cut-off.
! We only need one cut-off and we choose the one of the lattice atoms since s0
! is usually larger for them.

    rcut = betas0_l * sqrt3
    !rcut = a_lat * sqrt3 * isqrt2
    rr = 4 * rcut / (sqrt3 + 2.0d0)
    acut = 9.210240d0/(rr -rcut) ! ln(10000)

    xl = b * twelveth / (1.0d0 + exp(acut*(rnnl-rcut)))
    xp = b * twelveth / (1.0d0 + exp(acut*(rnnp-rcut)))

!-----------------------------------GAMMA--------------------------------------
! Gamma enforces the cut-off together with theta (see below)
! Gamma is defined as inverse.

    r3temp = rnnl - betas0_l
    igamma1l = 1.0d0 / sum(xl*exp(   -pars_l(1) * r3temp))
    igamma2l = 1.0d0 / sum(xl*exp(-kappadbeta_l * r3temp))

    r3temp = rnnp - betas0_p
    igamma1p = 1.0d0 / sum(xp*exp(   -pars_p(1) * r3temp))
    igamma2p = 1.0d0 / sum(xp*exp(-kappadbeta_p * r3temp))

!------------------------------------------------------------------------------
!                          Sigma and Pair-wise Contributions
!                          =================================
!------------------------------------------------------------------------------

    allocate(sigma_ll(slab%n_atoms), sigma_pp(teilchen%n_atoms))
    allocate(sigma_lp(slab%n_atoms), sigma_pl(teilchen%n_atoms))
    allocate(s_l(slab%n_atoms), s_p(teilchen%n_atoms))
    allocate(dsigma_ll(3,slab%n_atoms), dsigma_pp(3,teilchen%n_atoms))

    ! initialize some accumulators
    sigma_ll = 0.0d0
    sigma_pp = 0.0d0
    sigma_pl = 0.0d0
    sigma_lp = 0.0d0
    V_ll     = 0.0d0
    V_pp     = 0.0d0
    V_lp     = 0.0d0
    V_pl     = 0.0d0
    dsigma_ll= 0.0d0
    dsigma_pp= 0.0d0

    ! slab-slab
    do i = 1, slab%n_atoms
        do j = i+1, slab%n_atoms

            ! Applying PBCs
            r3temp = slab%r(:,i) - slab%r(:,j)   ! distance vector
            r3temp = matmul(cell_imat, r3temp)   ! transform to direct coordinates

            r3temp(1) = r3temp(1) - Anint(r3temp(1))! imaging
            r3temp(2) = r3temp(2) - Anint(r3temp(2))
            r3temp(3) = r3temp(3) - Anint(r3temp(3))
            r3temp    = matmul(cell_mat, r3temp)    ! back to cartesian coordinates

            r =  sqrt(sum(r3temp**2))               ! distance
            r3temp = r3temp/r                       ! unit vector j -> i

            ! cut-off function
            rtemp = exp(acut*(r - rcut))
            theta = 1.0d0 / (1.0d0 + rtemp)
            dtheta = (pars_l(1) + acut*rtemp*theta)*theta*r3temp

            rtemp = theta*exp(-pars_l(1)*(r - betas0_l))    ! sigma_ij*gamma1
            sigma_ll(i) = sigma_ll(i) + rtemp
            sigma_ll(j) = sigma_ll(j) + rtemp

            dtheta = dtheta*rtemp
            dsigma_ll(:,i) = dsigma_ll(:,i) - dtheta
            dsigma_ll(:,j) = dsigma_ll(:,j) + dtheta

            rtemp = theta*exp(-kappadbeta_l*(r - betas0_l)) ! V_ij*gamma2*V_0
            V_ll = V_ll + rtemp

        end do
    end do

    ! projectile-projectile
    do i = 1, teilchen%n_atoms
        do j = i+1, teilchen%n_atoms

            ! Applying PBCs
            r3temp = teilchen%r(:,i) - teilchen%r(:,j)   ! distance vector
            r3temp = matmul(cell_imat, r3temp)           ! transform to direct coordinates

            r3temp(1) = r3temp(1) - Anint(r3temp(1))! imaging
            r3temp(2) = r3temp(2) - Anint(r3temp(2))
            r3temp(3) = r3temp(3) - Anint(r3temp(3))
            r3temp    = matmul(cell_mat, r3temp)    ! back to cartesian coordinates

            r =  sqrt(sum(r3temp**2))               ! distance
            r3temp = r3temp/r                       ! unit vector j -> i

            ! cut-off function
            rtemp = exp(acut*(r - rcut))
            theta = 1.0d0 / (1.0d0 + rtemp)
            dtheta = (pars_p(1) + acut*rtemp*theta)*theta*r3temp

            rtemp = theta*exp(-pars_p(1) * (r - betas0_p) )     ! sigma_ij*gamma1
            sigma_pp(i) = sigma_pp(i) + rtemp
            sigma_pp(j) = sigma_pp(j) + rtemp

            dtheta = dtheta*rtemp
            dsigma_pp(:,i) = dsigma_pp(:,i) - dtheta
            dsigma_pp(:,j) = dsigma_pp(:,j) + dtheta

            rtemp = theta*exp(-kappadbeta_p * (r - betas0_p))   ! V_ij*gamma2*V_0
            V_pp = V_pp + rtemp

        end do
    end do

    ! projectile-slab
    do i = 1, teilchen%n_atoms
        do j = 1, slab%n_atoms

            ! Applying PBCs
            r3temp = teilchen%r(:,i) - slab%r(:,j)   ! distance vector
            r3temp = matmul(cell_imat, r3temp)       ! transform to direct coordinates

            r3temp(1) = r3temp(1) - Anint(r3temp(1))! imaging
            r3temp(2) = r3temp(2) - Anint(r3temp(2))
            r3temp(3) = r3temp(3) - Anint(r3temp(3))
            r3temp    = matmul(cell_mat, r3temp)    ! back to cartesian coordinates

            r =  sqrt(sum(r3temp**2))               ! distance
            r3temp = r3temp/r                       ! unit vector j -> i

            ! cut-off function
            rtemp = exp(acut*(r - rcut))
            theta = 1.0d0 / (1.0d0 + rtemp)
            dtheta = (pars_p(1) + acut*rtemp*theta)*theta*r3temp

            rtemp = theta*exp(-pars_p(1) * (r - betas0_p) )     ! sigma_ij*gamma1
            sigma_lp(j) = sigma_lp(j) + rtemp
            dsigma_lp(:,j) = dsigma_lp(:,j) - dtheta*rtemp

            dtheta = (pars_l(1) + acut*rtemp*theta)*theta*r3temp

            sigma_pl(i) = sigma_pl(i) + theta*exp(-pars_l(1)*(r - betas0_l))

            ! V_ij*gamma2*V_0
            V_lp = V_lp + theta*exp(-kappadbeta_p*(r - betas0_p))
            V_pl = V_pl + theta*exp(-kappadbeta_l*(r - betas0_l))

        end do
    end do

    ! divide by cut-off scaling factors
    sigma_ll =  sigma_ll*igamma1l
    V_ll     =     V_ll*igamma2l*pars_l(5)
    sigma_pp = sigma_pp*igamma1p
    V_pp     =     V_pp*igamma2p*pars_p(5)
    sigma_lp = sigma_lp*igamma1l
    sigma_pl = sigma_pl*igamma1p
    V_lp     =     V_lp*igamma2l*pars_l(5)*chilp
    V_pl     =     V_pl*igamma2p*pars_p(5)*chipl

    dsigma_ll = dsigma_ll*igamma1l
    dsigma_pp = dsigma_pp*igamma1p

    print *, dsigma_pp(:,1), sigma_pp(1)

!-----------------------------NEUTRAL SPHERE RADIUS----------------------------
! The neutral sphere radius is the radius in which the entire density of the
! atom is included.

    s_l = -log((sigma_ll + chilp*sigma_lp)*twelveth)/(beta*pars_l(1))
    s_p = -log((sigma_pp + chipl*sigma_pl)*twelveth)/(beta*pars_p(1))

!----------------REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------------

    vref_l = 12.0d0 * pars_l(5)*sum(exp(-pars_l(6)*s_l))*igamma2l
    vref_p = 12.0d0 * pars_p(5)*sum(exp(-pars_p(6)*s_p))*igamma2p

!---------------------------COHESIVE FUNCTION-----------------------------------

    Ecoh = sum((1.0d0 + pars_l(4)*s_l)*exp(-pars_l(4)*s_l) - 1.0d0)*pars_l(3) &
         + sum((1.0d0 + pars_p(4)*s_p)*exp(-pars_p(4)*s_p) - 1.0d0)*pars_p(3)

!-------------------------------TOTAL ENERGY---------------------------------

    Epot = Ecoh - V_ll - V_pp - 0.50d0*(V_lp + V_pl - vref_l - vref_p)

    deallocate(dsigma_pp, dsigma_ll)
    deallocate(s_p, s_l)
    deallocate(sigma_pl, sigma_lp, sigma_pp, sigma_ll)

end subroutine pes

end module force

