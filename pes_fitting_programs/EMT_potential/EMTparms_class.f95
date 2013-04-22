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

subroutine emt_init (cell, n_Au, r0_lat, r0_part, pars_p, pars_l, energy)
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
    real(8), dimension(3), intent(in)   :: cell
    integer, intent(in)                 :: n_Au
    real(8), dimension(3), intent (in)  :: r0_part
    real(8), dimension(:,:), intent(in) :: r0_lat
    type(EMTparms), intent(inout)       :: pars_p
    type(EMTparms), intent(inout)       :: pars_l
    real(8), intent(out)                :: energy

! declare the variables that appear in the program

    integer :: i,j, k                        ! running parameter
    real(8) :: r                             ! distance
    real(8) :: rcut, rr, acut                ! values to calculate cut-off
    real(8) :: gamma1H, gamma2H, Ecoh
    real(8) :: theta, gamma1Au, gamma2Au, chiAuH
    real(8) :: sigma_HAu, s_H
    real(8), dimension(n_Au) :: sigma_AuAu, sigma_AuH, s_Au
    real(8), dimension(3) :: xAu, xH, rnnAu, rnnH
    real(8) :: rtemp, betas0_l, betas0_p     ! temporary real variables
    real(8), dimension(3) :: r3temp

! calculate cut-off
! FUTURE REVISION: cut-off should be defined via lattice constant _AND_ changeable.

    betas0_l = beta * pars_l%s0
    betas0_p = beta * pars_p%s0

    chiAuH = pars_p%n0 / pars_l%n0

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


! Here, the main loop starts.

    sigma_AuAu = 0
    sigma_HAu = 0

    do i = 1, n_Au
    ! calculate r_HAu and r_AuAu. Now we have periodic boundery conditions, too. Muhahahaha!

        do j = i+1, n_Au

            r3temp(1) = r0_lat(1,i)-r0_lat(1,j)
            r3temp(2) = r0_lat(2,i)-r0_lat(2,j)
            r3temp(1) = r3temp(1) - (cell(1)*ANINT(r3temp(1)/cell(1)))
            r3temp(2) = r3temp(2) - (cell(2)*ANINT(r3temp(2)/cell(2)))
            r3temp(3) = r0_lat(3,i)-r0_lat(3,j)
            r =  sqrt(sum(r3temp**2))

            ! calculate theta for lattice
            theta = 1 / (1 + exp( acut * (r - rcut) ) )

            ! calculate sigma
            rtemp = theta*exp(-pars_l%eta2 * (r - betas0_l) )
            sigma_AuAu(i) = sigma_AuAu(i) + rtemp
            sigma_AuAu(j) = sigma_AuAu(j) + rtemp
        end do

        r3temp(1) = r0_lat(1,i)-r0_part(1)
        r3temp(2) = r0_lat(2,i)-r0_part(2)
        r3temp(1) = r3temp(1) - (cell(1)*ANINT(r3temp(1)/cell(1)))
        r3temp(2) = r3temp(2) - (cell(2)*ANINT(r3temp(2)/cell(2)))
        r3temp(3) = r0_lat(3,i)-r0_part(3)
        r =  sqrt(sum(r3temp**2))

        ! calculate theta
        theta = 1 / (1 + exp( acut * (r - rcut) ) )

        ! calculate sigma
        sigma_AuH(i) = theta*exp(-pars_p%eta2 * (r - betas0_p) )
        rtemp = theta*exp(-pars_l%eta2 * (r - betas0_l) )
        sigma_HAu = sigma_HAu + rtemp

    end do

    ! don't forget the gamma
    sigma_AuAu = sigma_AuAu / gamma1Au
    sigma_AuH = sigma_AuH / gamma1Au
    sigma_HAu = sigma_HAu / gamma1H

    ! calculation of s
    s_Au = -log( (sigma_AuAu + chiAuH * sigma_AuH) / 12 ) / ( beta * pars_l%eta2)
    s_H  = -log( sigma_HAu / (12*chiAuH) ) / ( beta * pars_p%eta2)

    ! calculation of cohesive function
    Ecoh = sum( (1 + pars_l%lambda*s_Au) * exp(-pars_l%lambda * s_Au)-1 ) * pars_l%E0&
            + (1 + pars_p%lambda*s_H) * exp(-pars_p%lambda * s_H)* pars_p%E0


    print *, Ecoh


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
