program testAtomClass

    use atom_class
    use EMTparms_class
    implicit none

    type (atom)              :: particle
    type (atom), allocatable, dimension(:) :: lattice
    ! Data directory: declare variables and types
    type (EMTparms) :: particle_pars, lattice_pars
    namelist / lattice_pars_list /  lattice_pars
    namelist / particle_pars_list / particle_pars

    ! Read in the lattice parameters
    open(8,file='parameters_Au_f119.nml')
    read(8,nml=lattice_pars_list)
    write(*,*) lattice_pars
    close(8)
    ! Read in the particle parameters
    open(8,file='parameters_H_f119.nml')
    read(8,nml=particle_pars_list)
    write(*,*) particle_pars
    close(8)



!    print *, "H positions with using defaults ", H%r
!    print *, "H velocities with using defaults ", H%v
!    print *, "H complete with using defaults ", H
!    H = atom(0)
!    print *, "H positions with using atom(0)  ", H%r
!    D = H
!    print *, "D positions with using D = H    ", H%r


! allocate the slab
!allocate(slab(100))

!slab(7)%r=(/0.0, 0.0, 6.0/)
!print *, 'This is he first atom coordinate', slab(1)%r
!print *, 'this is the second:', slab(7)%r

!deallocate(slab) ! not necessary but more beautiful

end program
