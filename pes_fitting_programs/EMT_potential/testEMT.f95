program testAtomClass

    use atom_class
    use EMTparms_class
    implicit none
    ! Data directory: declare variables and types
    type (atom)                          :: particle
    type (atom), allocatable,dimension(:):: lattice
    real(8), allocatable, dimension(:,:) :: r0_lat
    real(8), dimension(3)                :: r0_part
    character(len=30)                      :: reference_configuration_fname


    type (EMTparms) :: particle_pars, lattice_pars
    namelist / lattice_pars_list /  lattice_pars
    namelist / particle_pars_list / particle_pars

    real(8) :: energy
    integer :: n_lat = 2            ! number of lattice atoms
    integer :: ierror               ! number to check if file opened correctly
    integer :: n_lat0_at            ! number of atoms in reference slab
    integer :: n_lay0               ! number of layers in reference slab
    integer :: i                    ! running integer to read in reference lattice
    real(8) :: nn0                  ! next neighbour distance in reference slab
    real(8) :: a0                   ! lattice constant


    reference_configuration_fname='ref_conf_Au111a.dat'

    ! Read in the lattice parameters
    open(8,file='parameters_Au_f119.nml')
    read(8,nml=lattice_pars_list)
    close(8)

    ! Read in the particle parameters
    open(8,file='parameters_H_f119.nml')
    read(8,nml=particle_pars_list)
    close(8)

    ! Read in the number of atoms in the lattice
    open(8, file=reference_configuration_fname, status='old', action='read', iostat=ierror)
    ! add spec that file exists and is read only.  error checking
        openif: if ( ierror == 0 ) then
                    read(8, *) n_lat0_at
                    read(8,*) n_lay0
                    read(8,*) nn0
                    ! read in the values for r0_lat
                    allocate(r0_lat(3,n_lat0_at))
!                    allocate(r0_lat(3,2))
                    readr0lat: do i = 1, n_lat0_at
                                read(8,*) r0_lat(1,i), r0_lat(2,i), r0_lat(3,i)

                    end do readr0lat
                    r0_part = (/0.0, 0.0, 6.0/)
!                    print *, r0_lat
!                    print *, n_lat0_at, n_lay0, nn0, r0_part
                    call emt_init(n_lat0_at,r0_lat, r0_part, particle_pars, lattice_pars, energy)
                    print *, energy

    end if openif

    ! allocate dynamic arrays
    allocate(lattice(n_lat))

 !   call emt_init (particle, lattice, r0, particle_pars, lattice_pars, energy)

!    print *, lattice
    deallocate(lattice)

end program
