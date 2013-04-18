program testAtomClass

    use atom_class
    use EMTparms_class
    implicit none

    type (atom)                          :: particle
    type (atom), allocatable,dimension(:):: lattice
    real(8), allocatable, dimension(:,:) :: r0_lat
    real(8), dimension(3)                :: r0_part
    character, dimension(100)            :: reference_configuration_fname

    reference_configuration_fname='ref_conf_Au111.dat'

    ! Data directory: declare variables and types
    type (EMTparms) :: particle_pars, lattice_pars
    namelist / lattice_pars_list /  lattice_pars
    namelist / particle_pars_list / particle_pars

    real(8) :: energy
    integer :: n_lat = 2            ! number of lattice atoms



    ! Read in the lattice parameters
    open(8,file='parameters_Au_f119.nml')
    read(8,nml=lattice_pars_list)
    close(8)

    ! Read in the particle parameters
    open(8,file='parameters_H_f119.nml')
    read(8,nml=particle_pars_list)
    close(8)

    ! Read in the number of atoms in the lattice
    open(8,file=reference_configuration_fname)
    ! add spec that file exists and is read only.  error checking
    read(8, *)

    read

    ! allocate dynamic arrays
    allocate(lattice(n_lat))

!    call emt (particle, lattice, r0, particle_pars, lattice_pars, energy)

    print *, lattice
    deallocate(lattice)

end program
