program EMT_fit_1

    use EMTparms_class
    use open_file
    implicit none
    integer, parameter          :: npts=560             ! number of points to be fitted
    integer, parameter          :: nparms=14            ! number of parameters in the fitting function
    integer, parameter          :: max_iterations=10    !maximum number of iterations

    real, dimension(npts)       :: Y                ! values of dependent variable at N data points.
    real, dimension(3,npts)     :: X                ! values of indepenent variable(s) at N data points.
    !       if M is number of independent variables
    !       M=1 => DIMENSION X(N)
    !       M>1 => DIMENSION X(N,M)

    real, dimension(nparms)     :: B                ! parameters for the fit
    !   At entry, B's are initial guesses for parameters
    !   On exit,  B's are the results of the fit.

    real, dimension(npts)       :: RRR
    !       Array used optionally for communication between program that calls
    !       NLLSQ and SUBROUTINE MODEL.
    !       used here to communicat the final values of the fit
    !
    !       the shape of RRR must be the same here and in SUBROUTINE MODEL
    !       RRR is not used directly by NLLSQ; only passed to model

    integer, dimension(8)          :: NARRAY = (/ npts, 3, nparms, 0, 1, 2, 10, max_iterations /)
    !                                              1    2    3     4  5  6   7    8
    ! Array of integer control parameters as follows:
    !       NARRAY(1) = N   = number of data points
    !       NARRAY(2) = M   = number of independent variables
    !       NARRAY(3) = K   = number of parameters
    !       NARRAY(4) = IP  = number of parameters to be held
    !                         constant during fit.
    !       NARRAY(5) = Intermediate Printout Option (0..3)
    !           0 =>  no intermediate printout
    !           1 =>  minimal printout
    !           2 =>  1 + correlation matrix
    !           3 =>  2 + PTP inverse
    !       NARRAY(6) = Final Printout Option (0..3)
    !          -1 =>  no printout except errors
    !           0 =>  parameters, statistical info.,correlation matrix,PTP inverse
    !           1 =>  0 + residuals and confidence regions
    !           2 =>  1 + PTP, PTP corr. coeff. nonlinear confidence regions.
    !       NARRAY(7) = Printout unit number (06)
    !       NARRAY(8) = KITER = Maximum number of iterations.

    real, dimension(8)          :: ARRAY = 0
    ! Array of statistical parameters. Use 0.0 to get default values.
    !       ARRAY(1) = AL     (.1)
    !       ARRAY(2) = DELTA  (.00001)
    !       ARRAY(3) = E      (.00005)
    !       ARRAY(4) = FF     (4.0)
    !       ARRAY(5) = GAMCR  (45 degrees)
    !       ARRAY(6) = T      (2.0)
    !       ARRAY(7) = TAU    (.001)
    !       ARRAY(8) = ZETA   (10**-31)

    integer, dimension(10)      :: IB = 0
    ! Integer Array containing the subscripts of parameters to be held constant.
    !       (number=NARRAY(4) ).
    !
    character(20)               :: TITLE    ! string of 20 characters used for title on intermediate and final printout.

    ! EMT parameters of lattice and particle
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

    ! not used real(8), dimension(3)       :: r0_part              ! hydrogen position for reference calc.

    integer                                :: i,k                  ! loop indicies


    ! variables for emt-energy calculations
    integer, dimension(560)                :: loc                  ! location number of site
    ! could define these in an enum

    real(8), dimension(3)                  :: r_part               ! hydrogen position for reference calc.
    real(8)                                :: energy               ! energy output from emt subroutines
    real(8)                                :: e_ref                ! reference energy with particle at infinity

    ! file names
    character(len=30)                      :: reference_configuration_fname
    character(len=30)                      :: H_position_fname
    character(len=30)                      :: H_dft_energy


    reference_configuration_fname = 'ref_conf_Au111a.dat'          ! file containing Au coordinates
    H_position_fname              = 'hau111_plot.E.dat'            ! file containing H coordinates
    H_dft_energy                  = 'hEMTfortran.dat'              ! file containing DFT energies
    TITLE = 'EMT NLLSQ Test'


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

    ! Maybe replace these reads with namelist ?

    read(8,*) n_lat0_at
    read(8,*) n_lay0
    read(8,*) nn0
    read(8,*) cell(1)
    read(8,*) cell(2)

    ! allocate array to hold the lattice coordinates and read lattice coordinates
    allocate(r0_lat(3,n_lat0_at))

    readr0lat: do i = 1, n_lat0_at
        read(8,*) r0_lat(1,i), r0_lat(2,i), r0_lat(3,i)
    end do readr0lat

    ! not used
    ! r0_part = (/0.0, 0.0, 6.0/)  ! H atom at 6 Angstrom for reference energy calculation

    ! initialize the emt potential subroutine and calculate reference energy

    call emt_init(cell, n_lat0_at, r0_lat, particle_pars, lattice_pars, E_ref)
    print *, 'reference energy=',E_ref

    k = 560

    call open_for_read (9, H_position_fname)
    readrH: do i = 1, k
        read(9,*) loc(i), r_part(1), r_part(2), r_part(3)

        call emt(r_part, particle_pars, lattice_pars, energy)
        !call emt(cell, n_lat0_at, r0_lat, r_part, particle_pars, lattice_pars, energy)
        !
        !                        print *, r_part(1), r_part(2), r_part(3)
        write(*,'(1X, I2, 4F15.10)') loc(i), r_part(1), r_part(2), r_part(3), energy
        write(7,'(1X, I2, 4F16.10)') loc(i), r_part(1), r_part(2), r_part(3), energy
    end do readrH

    close(7)
    close(8)
    close(9)

end program
    !-----------------------------------------------------------------------------------------------------|
    !                     Variables that were declared in testEMT.f95 but not usef                        |
    !-----------------------------------------------------------------------------------------------------|
    !unused   integer                             :: n_lat_at             ! number of atoms slab          |
    !unused   integer                             :: n_lay                ! number of layers slab         |
    !unused   real(8)                             :: nn                   ! next neighbour distance slab  |
    !unused   real(8), allocatable, dimension(:,:):: r_lat                ! lattice positions             |
    !unused   real(8), allocatable, dimension(:,:):: r_part_v             ! hydrogen position             |
    !unused  real(8)                              :: a0                   ! lattice constant              |
    !unused   integer                             :: n_lat = 2            ! number of lattice atoms       |
    !-----------------------------------------------------------------------------------------------------|
