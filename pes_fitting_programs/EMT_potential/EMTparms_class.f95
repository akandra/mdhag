module EMTparms_class

use atom_class

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
    real(8), parameter :: beta = 1.8093997906
    real(8), parameter :: bohr = 0.529177


    integer :: i,j                           ! running parameter
    real(8), dimension(n_Au) :: r_HAu        ! distance between H and Au
    real(8), dimension(n_Au,1) :: r_AuAu       ! distance between Au and Au
    real(8) :: rcut, rr, acut                ! values to calculate cut-off
    real(8), dimension(n_Au,1) :: xAu              ! help-array to calculate theta
    real(8), dimension(n_Au) :: xH
    real(8), dimension(n_Au,1) :: thetaAu    ! cut-off theta
    real(8), dimension(n_Au) :: thetaH       ! cut-off theta
    real(8), dimension(n_Au,1,3) :: xiAu     ! Array for gamma-calculation for Au
    real(8), dimension(n_Au,3) :: xiH        ! Array for gamma-calculation for H

    real(8),dimension(n_Au,1) :: tempAu=0    ! dummy array for Gold
    real(8), dimension(n_Au) :: tempH=0      ! dummy array for hydrogen


! calculate cut-off
    rcut = pars_l%s0 * beta * sqrt(3.)
    rr = 4. * rcut / (sqrt(3.) + 2.)
    acut = log(9999.)/(rr -rcut)

    ! calculate r_HAu and r_AuAu
    rHauloop: do i = 1, n_Au

            r_HAu(i) = sqrt((r0_lat(1,i)+r0_part(1))**2+(r0_lat(2,i)+r0_part(2))**2+(r0_lat(3,i)+r0_part(3))**2)

    end do rHauloop

! Start here with first Au-loop. Which views everything from the perspective of the 45th gold atom since that one is pretty much in the middle of the
!   first layer and it's gonna be replaced by a j, anyway (or, by another layer-implementation).
    rAuAuloop1: do i = 1, n_Au
                    r_AuAu(i,1) = sqrt((r0_lat(1,i)-r0_lat(1,45))**2+(r0_lat(2,i)-r0_lat(2,45))**2+(r0_lat(3,i)-r0_lat(3,45))**2)


! calculate theta
                    xAu(i,1) =exp( acut * (r_AuAu(i,1) - rcut) )
                    thetaAu(i,1) = 1 / (1 + xAu(i,1))

                    xH(i) = exp( acut * (r_HAu(i) - rcut) )
                    thetaH(i) = 1 / (1 + xH(i))

! Calculate Gamma (this is not finished and will turn out to be rather vexing)
                    rAunn = pars%s0 * beta / 1
                    xiAu(i,1,1) = 1 / ( 1 + exp( acut (rAunn) ) )
                    tempAu=
                    write(*,*) thetaAu, thetaH


    end do rAuAuloop1


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
