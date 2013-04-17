program testAtomClass

    use atom_class
    use EMTparms_class
    implicit none

    type (atom)              :: particle
    type (atom), allocatable, dimension(:) :: lattice
    real(8), allocatable, dimension(:,:) :: r0
    ! Data directory: declare variables and types
    type (EMTparms) :: particle_pars, lattice_pars
    namelist / lattice_pars_list /  lattice_pars
    namelist / particle_pars_list / particle_pars

    real(8) :: energy
    integer :: n_lat = 100            ! number of lattice atoms

    ! Read in the lattice parameters
    open(8,file='parameters_Au_f119.nml')
    read(8,nml=lattice_pars_list)
    close(8)
    ! Read in the particle parameters
    open(8,file='parameters_H_f119.nml')
    read(8,nml=particle_pars_list)
    close(8)

    ! allocate dynamic arrays
    allocate(lattice(n_lat), r0(3,n_lat))

    call emt (particle, lattice, r0, particle_pars, lattice_pars, energy)

    print *, energy

    deallocate(lattice, r0)

end program
