module atom_class
    implicit none

    type atom
        real(8), dimension(3)    :: r=3  ! position
        real(8), dimension(3)    :: v=1  ! velocity
        real(8), dimension(3)    :: f=2  ! force
    contains

    end type atom


end module
