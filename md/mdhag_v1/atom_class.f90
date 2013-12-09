module atom_class
    !
    ! Purpose:
    !           This module containes the definitions of all the types and all
    !           constants that are used in the program.
    !           1. atoms class
    !           2. The EMT Parameter class.

    implicit none
    save

    real(8), parameter          :: sqrt2    = 1.41421356237d0
    real(8), parameter          :: isqrt2   = 0.70710678118d0
    real(8), parameter          :: sqrt3    = 1.73205080757d0
    real(8), parameter          :: pi       = 3.14159265359d0
    real(8), parameter          :: kB       = 0.0000861730d0
    real(8), parameter          :: twelveth = 0.0833333333333333d0
    integer, parameter          :: randseed(13) = (/8,6,7,5,3,11,9,1,17,2,9,6,4/)

    ! Unit conversion constants:
!
! The Basic Units are:
!           Length : Angström
!           Time   : fs
!           Energy : eV
!     Thus, the derived Units are:
!           Mass   : eV fs^2 / A^2 = 1/103.6382 amu
!           Angle  : radian = 180 deg
   real(8), parameter          :: amu2mass = 103.6382d0
   real(8), parameter          :: deg2rad  = pi/180.0d0



    type atom
        real(8), dimension(3)    :: r=0.0d0     ! position
        real(8), dimension(3)    :: v=0.0d0     ! velocity
        real(8), dimension(3)    :: vp=0.0d0     ! predicted velocity
        real(8), dimension(3)    :: a=0.0d0     ! acceleration
        real(8), dimension(3)    :: aalt=0.0d0  !  old acceleration
        real(8), dimension(3)    :: auralt=0.0d0  !  old acceleration

        real(8), dimension(3)    :: f=0.0d0     ! force
    end type atom

    type species
        character(len=10)           :: name     ! Element symbole
        real(8)                     :: mass     ! mass
        integer                     :: n        ! Number of atoms
        character(len=10)           :: pot      ! Name of Potential
        integer                     :: n_pars   ! Number of Parameters
        integer                     :: fric     ! Type of Friction
    end type species


end module
