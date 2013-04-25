program testaimd

    use atom_class
    use EMTparms_class
    implicit none

    ! Data directory: declare variables and types
    type (atom)                          :: particle
    type (atom), allocatable,dimension(:):: lattice
    type (EMTparms) :: particle_pars, lattice_pars
    namelist / lattice_pars_list /  lattice_pars
    namelist / particle_pars_list / particle_pars

! input for reference calculations
    real(8) :: E_ref, energy        ! reference energy and emt-energy
    integer :: n_lat = 2            ! number of lattice atoms
    integer :: ierror, ierror1      ! number to check if file opened correctly
    integer :: n_lat0_at            ! number of atoms in reference slab
    integer :: n_lay0               ! number of layers in reference slab
    integer :: i,k                  ! running integer to read in reference lattice
    real(8) :: nn0                  ! next neighbour distance in reference slab
    real(8) :: a0                   ! lattice constant
    real(8), dimension(3) :: cell   ! dimensions of the cell
    real(8), allocatable, dimension(:,:) :: r0_lat      ! lattice positions for reference calc.
    real(8), dimension(3)                :: r0_part     ! hydrogen position for reference calc.
    character(len=30)                    :: reference_configuration_fname


! input for emt-energy calculations
    integer :: n_lat_at             ! number of atoms slab
    integer :: n_lay                ! number of layers slab
    real(8) :: c11, c12, c22, c33   ! Matrix element for conversion into cartesian
    real(8), dimension(3)                :: r_H_d      ! H-position in direct coordinates
    real(8), dimension(3)                :: r_part     ! hydrogen position for reference calc.
    real(8), allocatable, dimension(:,:) :: r_lat_d    ! lattice positions in direct coordinates
    real(8), allocatable, dimension(:,:) :: r_lat      ! lattice positions
    character(len=30)                    :: position_fname
    integer :: ierror2, ierror3      ! number to check if file opened correctly

    position_fname = 'XDATCAR_ACC_fsv_005'
    reference_configuration_fname='ref_conf_Au111a.dat' ! File which contains Au coordinates
    ! Read in the lattice parameters
    open(8,file='parameters_Au_f119.nml')
    read(8,nml=lattice_pars_list)
    close(8)

    ! Read in the particle parameters
    open(8,file='parameters_H_f119.nml')
    read(8,nml=particle_pars_list)
    close(8)

    ! Read in the conversion to Cartesian coordinates and number of atoms in the lattice and the H-atom position
    open(8, file=reference_configuration_fname, status='old', action='read', iostat=ierror1)
    open(9, file=position_fname, status='old', action='read', iostat=ierror2)

    ! add spec that file exists and is read only.  error checking

    ierror = ierror1 + ierror2
    openif: if ( ierror == 0 ) then
                read(8, *) n_lat0_at
                read(8,*) n_lay0
                read(8,*) nn0
                read(8,*) cell(1)
                read(8,*) cell(2)

                ! read in the values for r0_lat
                allocate(r0_lat(3,n_lat0_at))
                readr0lat: do i = 1, n_lat0_at
                            read(8,*) r0_lat(1,i), r0_lat(2,i), r0_lat(3,i)

                end do readr0lat
                r0_part = (/0.0, 0.0, 6.0/)
                call emt_init(cell, n_lat0_at, r0_lat, particle_pars, lattice_pars, E_ref)
!---------------------------------------reference calc. finished --------------------------------
! Begin AIMD loop
                n_lat_at = 16
                n_lay= 4
                cell(1) = 2.9
                cell(2) = 2.5

                do i = 1, 1200

                    read(9,*) c11
                    read(9,*) c12, c22
                    read(9,*) c33
                    print *, c11, c12, c22, c33

!                        call emt(cell, n_lat0_at, r0_lat, r_part, particle_pars, lattice_pars, energy)

!                        write(*,'(1X, I2, 4F15.10)') loc(i), r_part(1), r_part(2), r_part(3), energy-E_ref
!                        write(7,'(1X, I2, 4F16.10)') loc(i), r_part(1), r_part(2), r_part(3), energy-E_ref


                end do


    end if openif
    close(8)
    close(9)


end program

