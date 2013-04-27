program EMT_fit_1

   use EMTparms_class
   implicit none

   ! interface declarations of subroutines to open files
   ! if open fails, these subroutines abort with STOP 101
   interface
      subroutine open_for_read(lun, file_name)
         integer, intent(in)     :: lun
         character, intent(in)   :: file_name
      end subroutine

      subroutine open_for_write(lun, file_name)
         integer, intent(in)     :: lun
         character, intent(in)   :: file_name
      end subroutine
   endinterface

   ! Instantiate class for EMT parameters of lattice and particle
   type (EMTparms)               :: particle_pars, lattice_pars

   ! Namelists to handle I/O of parameters
   namelist / lattice_pars_list  / lattice_pars
   namelist / particle_pars_list / particle_pars

   ! lattice and particle variables
   integer                                :: n_lat0_at            ! number of atoms in reference slab
   integer                                :: n_lay0               ! number of layers in reference slab
   real(8)                                :: nn0                  ! next neighbour distance in reference slab
   real(8), dimension(3)                  :: cell                 ! dimensions of the cell
   real(8), allocatable, dimension(:,:)   :: r0_lat               ! lattice positions for reference calc.
   real(8), dimension(3)                  :: r0_part              ! hydrogen position for reference calc.

   integer                                :: i,k                  ! loop indicies


   ! variables for emt-energy calculations
   integer, dimension(560)                :: loc                  ! location number of site
                                                                  ! TODO define these in an enum
   real(8), dimension(3)                  :: r_part               ! hydrogen position for reference calc.
   real(8)                                :: energy               ! energy output from emt subroutines
   real(8)                                :: e_ref                ! reference energy with particle at infinity

   ! file names
   character(len=30)                      :: reference_configuration_fname
   character(len=30)                      :: H_position_fname
   character(len=30)                      :: H_dft_energy
   reference_configuration_fname = 'ref_conf_Au111a.dat'          ! file containing Au coordinates
   H_position_fname              = 'hau111_plot.E.dat'            ! file containing H coordinates
   H_dft_energy                  =  'hEMTfortran.dat'             ! file containing DFT energies

   !-----------------------------------------------------------------------------------------------------|
   !                     Variables that were declared in testEMT.f95                                      |
   !-----------------------------------------------------------------------------------------------------|
   !unused   integer                             :: n_lat_at             ! number of atoms slab          |
   !unused   integer                             :: n_lay                ! number of layers slab         |
   !unused   real(8)                             :: nn                   ! next neighbour distance slab  |
   !unused   real(8), allocatable, dimension(:,:):: r_lat                ! lattice positions             |
   !unused   real(8), allocatable, dimension(:,:):: r_part_v             ! hydrogen position             |
   !unused  real(8)                              :: a0                   ! lattice constant              |
   !unused   integer                             :: n_lat = 2            ! number of lattice atoms       |
   !-----------------------------------------------------------------------------------------------------|

   ! Read in the lattice EMT parameters
   call open_for_read(8,'parameters_Au_f119.nml')
   read(8,nml=lattice_pars_list)

   ! Read in the particle EMT parameters
   call open_for_read(8,'parameters_H_f119.nml')
   read(8,nml=particle_pars_list)


   ! Open files
   call open_for_write(7, H_dft_energy)

   ! Read in the lattice definition parameters
   call open_for_read (8, reference_configuration_fname)
   ! TO DO - replace these reads with namelist ?
   read(8,*) n_lat0_at
   read(8,*) n_lay0
   read(8,*) nn0
   read(8,*) cell(1)
   read(8,*) cell(2)

   ! allocate array to hold the lattice coordinates and read lattice coordinates
   allocate(r0_lat(3,n_lat0_at))
   r0_part = (/0.0, 0.0, 6.0/)  ! H atom at 6 Angstrom for reference energy calculation

   ! initialize the emt potential subroutine and calculate reference energy
   call emt_init(cell, n_lat0_at, r0_lat, particle_pars, lattice_pars, E_ref)

   k = 560

   call open_for_read (9, H_position_fname)
   readrH: do i = 1, k
      read(9,*) loc(i), r_part(1), r_part(2), r_part(3)

      call emt(cell, n_lat0_at, r0_lat, r_part, particle_pars, lattice_pars, energy)
      !
      !                        print *, r_part(1), r_part(2), r_part(3)
      write(*,'(1X, I2, 4F15.10)') loc(i), r_part(1), r_part(2), r_part(3), energy-E_ref
      write(7,'(1X, I2, 4F16.10)') loc(i), r_part(1), r_part(2), r_part(3), energy-E_ref
   end do readrH

print *, energy, energy-E_ref

close(7)
close(8)
close(9)

end program



