module emt_init_data
    implicit none
    save

    integer                                 :: n_lat0_at         ! number of atoms in reference slab
    integer                                 :: n_lay0            ! number of layers in reference slab
    real(8)                                 :: nn0               ! next neighbour distance in reference slab
    real(8), dimension(3)                   :: cell              ! dimensions of the cell
    real(8), allocatable, dimension(:,:)    :: r0_lat            ! lattice positions for reference calc.

contains
subroutine variab_init
    integer :: i
    character(len=30) :: lattice_configuration_fname

    lattice_configuration_fname = 'ref_conf_Au111a.dat'

    call open_for_read (8, lattice_configuration_fname)
    read(8,*) n_lat0_at
    read(8,*) n_lay0
    read(8,*) nn0
    read(8,*) cell(1)
    read(8,*) cell(2)
    allocate(r0_lat(3,n_lat0_at))           ! allocate array to hold lattice coordinates
    readr0lat: do i = 1, n_lat0_at
        read(8,*) r0_lat(1,i), r0_lat(2,i), r0_lat(3,i)
        !write(7,*)r0_lat(1,i), r0_lat(2,i), r0_lat(3,i)
    end do readr0lat
end subroutine variab_init

end module emt_init_data
