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
    integer :: i,k,j                ! running integer to read in reference lattice
    real(8) :: nn0                  ! next neighbour distance in reference slab
    real(8) :: a0                   ! lattice constant
    real(8), dimension(3) :: cell_ref   ! dimensions of the cell
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
    real(8), allocatable, dimension(:,:) :: temp      ! temporary array
    real(8), dimension(3) :: cell   ! dimensions of the cell
    character(len=30)                    :: position_fname
    integer :: ierror2, ierror3      ! number to check if file opened correctly
    real :: temp1

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
                read(8,*) cell_ref(1)
                read(8,*) cell_ref(2)

                ! read in the values for r0_lat
                allocate(r0_lat(3,n_lat0_at))
                readr0lat: do i = 1, n_lat0_at
                            read(8,*) r0_lat(1,i), r0_lat(2,i), r0_lat(3,i)

                end do readr0lat
                r0_part = (/0.0, 0.0, 6.0/)
                call emt_init(cell_ref, n_lat0_at, r0_lat, particle_pars, lattice_pars, E_ref)
!---------------------------------------reference calc. finished --------------------------------

! Begin AIMD loop
                n_lat_at = 16
                n_lay= 4
                cell(1) = 2*isqrt_2*4.205
                cell(2) = sqrt_3 * isqrt_2 * 4.205
                read(9,*)
                read(9,*)
                read(9,*) c11
                read(9,*) c12, c22
                read(9,*) temp1, temp1, c33
                do temp1= 1,2
                    read(9,*)
                end do
                allocate(r_lat_d(17,3))
                allocate(temp(17,3))
                allocate(r_lat(16,3))

! Loop over the steps in AIMD
                do i = 1, 1
                    read(9,*)
                    do j = 1, 17
                    ! Reads in direct coordinates and transfers them to cartesian
                        read(9,*) r_lat_d(j,1), r_lat_d(j,2), r_lat_d(j,3)
                        temp(j,1) = c11*r_lat_d(j,1) + c12 * r_lat_d(j,2)
                        temp(j,2) = c22* r_lat_d(j,2)
                        temp(j,3) = c33* r_lat_d(j,3)

                        if (temp(j,3) > 6.3) then
                            temp(j,3) = temp(j,3)-c33
                        end if

                    end do
                ! Asignes coordinates to gold atoms and hydrogen atom
                    do j = 1, 16
                        r_lat(j,1) = temp(j,1)
                        r_lat(j,2) = temp(j,2)
                        r_lat(j,3) = temp(j,3)
                    end do
                    r_part(1) = temp(17,1)
                    r_part(2) = temp(17,2)
                    r_part(3) = temp(17,3)
!                    print *, c11, c12, c22, c33
!                    print *, r_lat
                    print *, r_part
                    call emt(cell_ref, n_lat0_at, r0_lat, r_part, particle_pars, lattice_pars, energy)
                    print *, E_ref
                    print *, energy

                end do


    end if openif
    close(8)
    close(9)


end program

