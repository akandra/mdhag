module atom_class
    implicit none

    type atom
        real(8), dimension(3)    :: r=1  ! position
        real(8), dimension(3)    :: v=2  ! velocity
        real(8), dimension(3)    :: f=3  ! force
    end type atom


end module
