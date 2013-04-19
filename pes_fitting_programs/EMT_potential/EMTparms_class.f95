module EMTparms_class

use atom_class

    real(8), parameter :: sqrt_2 = 1.41421356237
    real(8), parameter :: sqrt_3 = 1.73205080757
    real(8), parameter :: isqrt_2 = 0.70710678118
    real(8), parameter :: pi = 3.14159265359
    real(8), parameter :: beta = 1.8093997906
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

subroutine emt_init (n_Au, r0_lat, r0_part, pars_p, pars_l, energy)
!
! Purpose:
!           Here, the fitting procedure is just implemented and the reference energy calculated.
!   Input variables:
!                   r0_lat : equilibrium positions of the lattice atoms
!                   r0_part: reference position of the particle (at 0,0, 6.0 A above the surface)
!                   pars_p : parameters for the particle
!                   pars_l : parameters for the lattice
!
implicit none

! declare variables and parameters
!    type(atom), intent (in)  :: particle, lattice(:)
!    real(8), dimension(:,:), intent (in) :: r0
!    type(EMTparms),   intent (in)  :: pars_p, pars_l
    integer, intent(in)                 :: n_Au
    real(8), dimension(3), intent (in)  :: r0_part
    real(8), dimension(:,:), intent(in) :: r0_lat
    type(EMTparms), intent(inout)       :: pars_p
    type(EMTparms), intent(inout)       :: pars_l
    real(8), intent(out)                :: energy

! declare the variables that appear in the program

    integer :: n_AuAu                        ! number of elements under diagonal
    integer :: i,j, k                        ! running parameter
    real(8), dimension(n_Au) :: r_HAu        ! distance between H and Au
    real(8), dimension(n_Au*(n_Au-1)/2) :: r_AuAu       ! distance between Au and Au
    real(8) :: rcut, rr, acut                ! values to calculate cut-off
    real(8) :: thetaH, gamma1H, gamma2H
    real(8) :: thetaAu, gamma1Au, gamma2Au
    real(8), dimension(3) :: xAu, xH, rnnAu, rnnH
    real(8) :: rtemp, betas0_l, betas0_p     ! temporary real variables
    real(8), dimension(3) :: r3temp



    n_AuAu = n_Au*(n_Au-1)/2
    ! calculate r_HAu and r_AuAu
    do i = 1, n_Au
        r_HAu(i) = sqrt((r0_lat(1,i)-r0_part(1))**2+&
                        (r0_lat(2,i)-r0_part(2))**2+&
                        (r0_lat(3,i)-r0_part(3))**2)
        k = (i - 1)*(n_Au - i/2) - i
        do j = i+1, n_Au
            r_AuAu(j+k) = sqrt((r0_lat(1,i)-r0_lat(1,j))**2+&
                               (r0_lat(2,i)-r0_lat(2,j))**2+&
                               (r0_lat(3,i)-r0_lat(3,j))**2)
        end do
    end do


! calculate cut-off
! FUTURE REVISION: cut-off should be defined via lattice constant _AND_ changeable.

    betas0_l = beta * pars_l%s0
    betas0_p = beta * pars_p%s0

    rcut = betas0_l * sqrt_3
    rr = 4 * rcut / (sqrt_3 + 2)
    acut = 9.21024/(rr -rcut) ! ln(10000)

    rnnAu(1) = betas0_l
    rnnAu(2) = rnnAu(1) * sqrt_2
    rnnAu(3) = rnnAu(1) * sqrt_3
    rnnH(1) = betas0_p
    rnnH(2) = rnnH(1) * sqrt_2
    rnnH(3) = rnnH(1) * sqrt_3

    xAu = b / (12 * (1 + exp(acut*(rnnAu-rcut))))
    xH = b / (12 * (1 + exp(acut*(rnnH-rcut))))

 ! Definition of gamma
    r3temp = rnnAu-betas0_l
    gamma1Au = sum(xAu*exp(-pars_l%eta2 * r3temp))
    gamma2Au = sum(xAu*exp(-pars_l%kappa/beta * r3temp))

    r3temp = rnnH-betas0_p
    gamma1H = sum(xH*exp(-pars_p%eta2 * r3temp))
    gamma2H = sum(xH*exp(-pars_p%kappa/beta * r3temp))

    print *, gamma1H, gamma2H

    Auloop: do i = 1, n_Au


! calculate theta
        thetaAu = 1 / (1 + exp( acut * (r_AuAu(i) - rcut) ) )!

        thetaH = 1 / (1 + exp( acut * (r_HAu(i) - rcut) ))

! Calculate Gamma (this is not finished and will turn out to be rather vexing)
!                    rAunn = pars%s0 * beta / 1
!                    xiAu(i,1,1) = 1 / ( 1 + exp( acut (rAunn) ) )
 !                   tempAu=
 !                   write(*,*) thetaAu, thetaH


    end do Auloop


end subroutine emt_init

subroutine emt (particle, lattice, r0, pars_p, pars_l, energy)
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

    type(atom), intent (in)  :: particle, lattice(:)
    real(8), dimension(:,:), intent (in) :: r0
    type(EMTparms),   intent (in)  :: pars_p, pars_l
    real(8),                 intent (out) :: energy

    energy = 1347

end subroutine emt

end module
