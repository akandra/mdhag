module atom_class
    implicit none
    public                              ! sorry, for performance in accessing components

    type atom
        real, dimension(3)      :: r=1  ! position
        real, dimension(3)      :: v=2  ! velocity
        real, dimension(3)      :: f=3  ! force
    end type atom

    type atoms
        integer                             :: n_atoms=1!number of atoms
        real, dimension(:,:), allocatable   :: r        ! positions (3,natoms
        real, dimension(:,:), allocatable   :: v        ! velocities
        real, dimension(:,:), allocatable   :: f        ! forces
    end type atoms




    type atoms2
        integer                             :: n_atoms=1!number of atoms
        real, dimension(:,:), allocatable   :: r        ! positions (3,natoms
        real, dimension(:,:), allocatable   :: v        ! velocities
        real, dimension(:,:), allocatable   :: f        ! forces
    end type atoms2

    interface atoms2
        module procedure new_atoms2
    end interface

contains

    function new_atoms2(i)
        integer, intent(in) :: i
        type(atoms2) new_atoms2
!-------------------------------------------------
!   allocate r, v, and f
!-------------------------------------------------
        allocate(new_atoms2%r(3,i))
        allocate(new_atoms2%v(3,i))
        allocate(new_atoms2%f(3,i))

!-------------------------------------------------
!   initialize r, v, and f
!-------------------------------------------------
        new_atoms2%r = 1
        new_atoms2%v = 2
        new_atoms2%f = 3
    end function


end module
