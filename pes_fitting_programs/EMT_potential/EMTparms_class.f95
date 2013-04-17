module EMTparms_class

use atom_class

    type EMTparms
        character(2)::  name    = 'Au'

        real(8)       ::  eta2    = 3.84357     ! A^-1
        real(8)       ::  kappa   = 6.95507     ! A^-1
        real(8)       ::  lambda  = 4.1233825   ! A^-1

        real(8)       ::  E0      = ⁻3.8        ! eV
        real(8)       ::  n0      = 0.017325    ! A^-3
        real(8)       ::  s0      = 1.6417352   ! A
        real(8)       ::  V0      = 2.321       ! eV

    end type EMTparms

contains

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
