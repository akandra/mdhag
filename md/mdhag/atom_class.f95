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


    type atom
        real(8), dimension(3)    :: r=0.0d0  ! position
        real(8), dimension(3)    :: v=0.0d0  ! velocity
        real(8), dimension(3)    :: f=0.0d0  ! force
    end type atom

    type species
        character(len=10)           :: name     ! Element symbole
        real(8)                     :: mass     ! mass
        integer                     :: n        ! Number of atoms
        character(len=10)           :: pot      ! Name of Potential
        integer                     :: n_pars   ! Number of Parameters
    end type species


end module
