module atom_class
    implicit none
    public         ! public for performance in accessing components

    ! -------------------------------------------------------------------------
    !                  Type atom                                              |
    ! -------------------------------------------------------------------------
    !  structure to hold the position, velocity, and force on a single atom   |
    ! -------------------------------------------------------------------------
    type atom
        real(8), dimension(3)      :: r=(/1,2,3/)  ! position
        real(8), dimension(3)      :: v=(/4,5,6/)  ! velocity
        real(8), dimension(3)      :: f=(/7,8,9/)  ! force
    end type atom

    !-------------------------------------------------------------------------|
    !                  Type atoms                                             |
    !-------------------------------------------------------------------------|
    ! structure to hold the position, velocity, and force on a multiple atoms |
    !    use rank2 array so positions, velocities, and forces are stored      |
    !    in sequentional memory locations for efficient access                |
    !    each array should be allocated (3,n_atom)                            |
    !-------------------------------------------------------------------------|

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
    !                                                                         |
    !-------------------------------------------------------------------------|



    interface atoms
        module procedure new_atoms
    end interface

contains

    function new_atoms(n_atoms)
        integer, intent(in) :: n_atoms
        type(atoms) new_atoms
!-------------------------------------------------
!   allocate r, v, and f
!-------------------------------------------------
        allocate(new_atoms%r(3,n_atoms))
        allocate(new_atoms%v(3,n_atoms))
        allocate(new_atoms%f(3,n_atoms))

!-------------------------------------------------
!   initialize r, v, and f
!-------------------------------------------------
        new_atoms%r = 0
        new_atoms%v = 0
        new_atoms%f = 0
    end function


end module
