program mdhag
    ! Purpose:
    !       Do molecular dynamics calculations with the EMT potential.
    !
    ! Date          Author          History of Modification
    ! ====          ======          =======================
    ! 01.10.2013    Sascha&Svenja   Original
    !
    !
!    use atom_class
    use md_init
    use force

    implicit none

    real(8) :: Epot

    type(atom), dimension(:), allocatable :: slab, teilchen

    integer :: i

    ! Call up procedure that gets us geometries
    call simbox_init(teilchen,slab)

    ! REMINDER:
    ! Print out what kind of potential has been chosen, the parameters and the
    ! initial conditions. DON'T FORGET!

    call pes(teilchen, slab, Epot)

    write(*,'(10f10.5,2X,2A)') slab(1)
    write(*,'(10f10.5,2X,2A)') teilchen
    write(*,*) Epot

end program mdhag
