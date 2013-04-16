module atom_class
    implicit none

    type atom
        real, dimension(3)    :: r=3  ! position
        real, dimension(3)    :: v=1  ! velocity
        real, dimension(3)    :: f=2  ! force
    contains

    end type atom


end module
