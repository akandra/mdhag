module atom_class
    !
    ! Purpose:
    !           This module containes the definitions of all the types and all
    !           constants that are used in the program.
    !

    implicit none
    public         ! public for performance in accessing components
    save

    ! Various useful constants
    real(8), parameter          :: sqrt2    = 1.41421356237d0
    real(8), parameter          :: isqrt2   = 0.70710678118d0
    real(8), parameter          :: sqrt3    = 1.73205080757d0
    real(8), parameter          :: pi       = 3.14159265359d0
    real(8), parameter          :: kB       = 0.0000861730d0
    real(8), parameter          :: twelveth = 0.083d3
    integer, parameter          :: randseed(13) = (/8,6,7,5,3,11,9,1,17,2,9,6,4/)

    ! Conversion constants to program units
    !
    ! Program basic units
    !           Length : Ang
    !           Time   : fs
    !           Energy : eV
    ! Program derived units
    !           Mass   : eV fs^2 / A^2 = 1/103.6382 amu
    !           Angle  : radian = 180 deg
    real(8), parameter          :: amu2mass = 103.6382d0
    real(8), parameter          :: deg2rad  = pi/180.0d0

    !  Type atom (obsolete)
    !   structure to hold the position, velocity, force etc. for a single atom
    type atom
        real(8), dimension(3)    :: r       = 0.0d0 ! position
        real(8), dimension(3)    :: v       = 0.0d0 ! velocity
        real(8), dimension(3)    :: vp      = 0.0d0 ! predicted velocity
        real(8), dimension(3)    :: a       = 0.0d0 ! acceleration
        real(8), dimension(3)    :: aalt    = 0.0d0 ! old acceleration
        real(8), dimension(3)    :: auralt  = 0.0d0 ! old acceleration
        real(8), dimension(3)    :: f       = 0.0d0 ! force
    end type atom

    !  Type species (obsolete)
    !   structure to hold the species-specific parameters
    type species
        character(len=10)           :: name     ! Element symbole
        real(8)                     :: mass     ! mass
        integer                     :: n        ! Number of atoms
        character(len=10)           :: pot      ! Name of Potential
        integer                     :: n_pars   ! Number of Parameters
        integer                     :: fric     ! Type of Friction
    end type species

    !  Type atoms
    !   structure to hold the position, velocity, force etc. for multiple atoms
    !       use rank2 array so positions, velocities, forces etc. are stored
    !       in sequentional memory locations for efficient access
    !       each array should be allocated (3,n_atom)
    type atoms
        integer                                :: n_atoms  ! number of atoms
        real(8), dimension(:,:), allocatable   :: r        ! positions
        real(8), dimension(:,:), allocatable   :: v        ! velocities
        real(8), dimension(:,:), allocatable   :: f        ! forces
    end type atoms

    !-------------------------------------------------------------------------|
    !                  Constructor for type atoms                             |
    !-------------------------------------------------------------------------|
    !    input:  n_atoms                                                      |
    !    allocates arrays as (3,n_atom)                                       |
    !                                                                         |
    !-------------------------------------------------------------------------|

    interface atoms
        module procedure new_atoms
    end interface

contains

    function new_atoms(n_atoms)
        integer, intent(in) :: n_atoms
        type(atoms) new_atoms

        allocate(new_atoms%r(3,n_atoms))    !   allocate
        allocate(new_atoms%v(3,n_atoms))
        allocate(new_atoms%f(3,n_atoms))

        new_atoms%r = 0.0d0                 !   initialize
        new_atoms%v = 0.0d0
        new_atoms%f = 0.0d0

        new_atoms%n_atoms = n_atoms

    end function


end module
