program EMT_fit_1

!----------------------------------------------------------------------------------------------------------------------
!
!   File:       EMT_fit_1.f95
!
!   Purpose:    Nonlinear Least Squares fit of parameters of EMT potential to DFT points using static lattice
!
!   Author:     D.J. Auerbach
!   Date:       29 April, 2013
!   Copyright:  2013 Max Planck Institute for Biophysical Chemistry
!
!   Uses:       NLLSQ - a non linear least squres routine developed by D. Marquardt and modified at Bell Labs.
!               EMT_parms class developed by S. Janka and A. Kandrasenka.  Modified here for ease of use with NLLSQ
!
!   Description:
!
!   Additional Files:
!       NLLSQ.f                 the non linear least squares routine
!       model.f95               subroutine to implement the model being fitted
!       EMTparms_class_v2.f95   subroutines and structures for EMT energy calculations
!       open_file.f95           module with routines to open files for reading and writing with error checking
!

!----------------------------------------------------------------------------------------------------------------------
!   Revision Log
!
!   Date        Author      Change
!   30.04.2013  dja         Initial release
!
!
!
!----------------------------------------------------------------------------------------------------------------------

    use EMTparms_class
    use open_file

    implicit none



    !------------------------------------------------------------------------------------------------------------------
    !
    !                                   NLLSQ VARIABLES
    !
    !------------------------------------------------------------------------------------------------------------------
    ! todo make the arrays allocatable and automatically determine npts by reading to end of data
    !------------------------------------------------------------------------------------------------------------------
    integer                     :: npts             ! number of points to be fitted
    integer                     :: nparms           ! number of parameters in the fitting function
    integer                     :: IP               ! number of parameters to be held constant during fit
    integer                     :: max_iterations   ! maximum number of iterations
    real(8), allocatable        :: Y(:)             ! values of dependent variable at N data points.
    real(8), allocatable        :: X(:,:)           ! values of indepenent variable(s) at N data points.
    real(8), dimension(20)      :: B                ! parameters for the fit. 20 is the maximum hard coded in NLLSQ
                                                    !   At entry, B's are initial guesses for parameters
                                                    !   On exit,  B's are the results of the fit.
    real(8), allocatable        :: RRR(:)           ! Array used optionally for communication with MODEL
                                                    !   used here to communicat the final values of the fit
    integer, dimension(8)       :: NARRAY           ! Integer array of control parameters

    ! NARRAY(1) = N   = number of data points
    ! NARRAY(2) = M   = number of independent variables
    ! NARRAY(3) = K   = number of parameters
    ! NARRAY(4) = IP  = number of parameters to be held constant during fit.
    ! NARRAY(5) = Intermediate Printout Option (0..3)
    ! NARRAY(6) = Final Printout Option (0..3)
    ! NARRAY(7) = Printout unit number
    ! NARRAY(8) = Maximum number of iterations.

    real, dimension(8)          :: ARRAY = 0        ! Array of statistical parameters. Use 0.0 to get default values.

    !       ARRAY(1) = AL     (.1)
    !       ARRAY(2) = DELTA  (.00001)
    !       ARRAY(3) = E      (.00005)
    !       ARRAY(4) = FF     (4.0)
    !       ARRAY(5) = GAMCR  (45 degrees)
    !       ARRAY(6) = T      (2.0)
    !       ARRAY(7) = TAU    (.001)
    !       ARRAY(8) = ZETA   (10**-31)

    integer, dimension(14)      :: IB = 0           ! Integer Array containing the subscripts of parameters to be held constant.
                                                    ! number=NARRAY(4)
    character(20)               :: TITLE            ! string of 20 characters used for title
                                                    ! on intermediate and final printout.

    !------------------------------------------------------------------------------------------------------------------
    !
    !                                   EMT VARIABLES
    !
    !------------------------------------------------------------------------------------------------------------------
    type (EMTparms)                         :: particle_pars    ! EMT parameters of particle
    type (EMTparms)                         :: lattice_pars     ! EMT parameters of lattice

    type (EMTparms)                         :: test_pars        ! used for check conversion routines

    namelist / lattice_pars_list  / lattice_pars                ! Namelists to handle I/O of parameters
    namelist / particle_pars_list / particle_pars

    integer                                 :: n_lat0_at         ! number of atoms in reference slab
    integer                                 :: n_lay0            ! number of layers in reference slab
    real(8)                                 :: nn0               ! next neighbour distance in reference slab
    real(8), dimension(3)                   :: cell              ! dimensions of the cell
    real(8), allocatable, dimension(:,:)    :: r0_lat            ! lattice positions for reference calc.

    integer                                 :: i,j,k             ! loop indicies
    integer                                 :: ios               ! io status


    real(8), dimension(3)                   :: r_part            ! hydrogen position for reference calc.
    integer, allocatable, dimension(:)      :: site              ! impact site
    real(8)                                 :: energy            ! energy output from emt subroutines
    real(8)                                 :: e_ref             ! reference energy with particle at infinity

    ! file names
    character(len=30)                       :: reference_configuration_fname
    character(len=30)                       :: H_position_fname
    character(len=30)                       :: H_dft_energy


    reference_configuration_fname = 'ref_conf_Au111a.dat'          ! file containing Au coordinates
    H_position_fname              = 'hau111_plot.E.dat'            ! file containing H coordinates
    H_dft_energy                  = 'hEMTfortran.dat'              ! file containing DFT energies
    TITLE = 'EMT NLLSQ Test'

    !------------------------------------------------------------------------------------------------------------------
    !                       READ LATTICE AND PARTICLE EMT PARAMETERS
    !------------------------------------------------------------------------------------------------------------------
    call open_for_read(8,'parameters_Au_f119.nml')
    read(8,nml=lattice_pars_list)
    call open_for_read(8,'parameters_H_f119.nml')
    read(8,nml=particle_pars_list)

    !------------------------------------------------------------------------------------------------------------------
    !                       READ LATTICE DEFINITION PARAMETSR
    !------------------------------------------------------------------------------------------------------------------
    call open_for_read (8, reference_configuration_fname)
    read(8,*) n_lat0_at
    read(8,*) n_lay0
    read(8,*) nn0
    read(8,*) cell(1)
    read(8,*) cell(2)
    allocate(r0_lat(3,n_lat0_at))           ! allocate array to hold lattice coordinates
    readr0lat: do i = 1, n_lat0_at
        read(8,*) r0_lat(1,i), r0_lat(2,i), r0_lat(3,i)
    end do readr0lat

    !------------------------------------------------------------------------------------------------------------------
    !                       INITIALIZE EMT POTENTIAL SUBROUTINE AND CALCULATE REFERENCE ENERGY
    !------------------------------------------------------------------------------------------------------------------
    call emt_init(cell, n_lat0_at, r0_lat, particle_pars, lattice_pars, E_ref)
    print *, 'reference energy=',E_ref

    !------------------------------------------------------------------------------------------------------------------
    !                       READ THE PARICLE POSITIONS AND DFT ENERGIES
    !------------------------------------------------------------------------------------------------------------------
    ! first find the length of the file
    call open_for_read (8, H_position_fname)
    i=1
    do
        read(8,*,iostat=ios)
        if(ios <0) exit
        i=i+1
    end do

    rewind(8)
    npts = i-1
    print *,'the number of points in H position list=',npts

    ! now allocate arrays needed by NLLSQ
    allocate(X(3,npts))             ! particle coordinates - indepenet variables
    allocate(Y(npts))               ! DFT energies         - depenent variables
    allocate(RRR(npts))             ! Array for communication with subroutine model
    allocate(site(npts))            ! Impact site

    ! read in the site, particle coordinates, and DFT energiespositions
    j=1
    do i=1, npts
        read(8,*) site(j), X(1,i), X(2,j), X(3,j), Y(j)
        ! keep the point if the DFT energy is less than 10 eV
        if ( (Y(j)<=10) .and. (i<npts) ) j=j+1
    end do
    close(8)
    npts=j

    print *,'the number of points <= 10 eV=',npts


    !------------------------------------------------------------------------------------------------------------------
    !                       CHECK EMT POTENTIAL SUBROUTINE AND WRITE RESULTS
    !------------------------------------------------------------------------------------------------------------------
    k = 4
    call open_for_write(10, H_dft_energy)
    write(*,*)  'site       X              Y              Z             EMT           DFT'
    Write(10,*) 'site       X              Y              Z             EMT           DFT'
    readrH: do i = 1, k
        call emt(X(:,i), particle_pars, lattice_pars, energy)
        ! old form: call emt(cell, n_lat0_at, r0_lat, r_part, particle_pars, lattice_pars, energy)

        write(*,'(1X, I2, 5F15.10)')  site(i), X(1,i), X(2,i), X(3,i), energy, Y(i)
        write(10,'(1X, I2, 5F15.10)') site(i), X(1,i), X(2,i), X(3,i), energy, Y(i)
    end do readrH

    !------------------------------------------------------------------------------------------------------------------
    !                       SETUP FOR NLLSQ
    !------------------------------------------------------------------------------------------------------------------
    ! Usage: CALL NLLSQ ( Y , X , B , RRR , NARRAY , ARRAY , IB , TITLE)

    !------------------ SET UP B-----------------
    call emt_parms2array(particle_pars,B(1:7))
    call emt_parms2array(lattice_pars ,B(8:14))

    !-----------CHECK WE GOT IT RIGHT -----------
    !write(*,1000) 'particle_pars=',particle_pars
    !write(*,1000) 'lattice_pars =',lattice_pars
    !write(*,1010) 'B(1:14)      =',B(1:14)

    !---------------TEST ARRAY2EMT_PARMS--------
    !call array2emt_parms(B(1:7),test_pars)
    !write(*,1000) 'test_pars     =',test_pars

    !-----------SET UP NARRAY-------------------
    nparms         = 14
    IP             = 7
    max_iterations = 4

    !*********************************************************
    !npts = 3 !!!!!!!!   FOR TESTING  !!!!!!!!!!!!!!!!!!!!!!!!!!
    !*********************************************************

    NARRAY(1) = npts            ! number of data points
    NARRAY(2) = 3               ! number of independent variables
    NARRAY(3) = nparms          ! number of parameters
    NARRAY(4) = IP              ! number of parameters to be held constant during fit.
    NARRAY(5) = 0               ! intermediate printout option (0..3)
    NARRAY(6) = 0               ! final printout option (0..3)
    NARRAY(7) = 10              ! printout unit number
    NARRAY(8) = max_iterations  ! maximum number of iterations

    !-----------SET UP ARRAY-------------------
    ARRAY = 0                   ! use defaults for NLLSQ constrol parameters
    IB    = (/1,2,3,4,5,6,7,8,9,10,11,12,13,0/)   ! array of indicies of parameters held constant
    TITLE = 'Test of EMT_fit_1'


    CALL NLLSQ ( Y , X , B , RRR , NARRAY , ARRAY , IB , TITLE)


1000 format(a15,2x,a2,7F7.3)
1010 format(a15,4x,7F7.3,/,19x,7f7.3)

end program
    !-------------------------------------------------------------------------------------------------------|
    !                     Variables that were declared in testEMT.f95 but not usef                          |
    !-------------------------------------------------------------------------------------------------------|
    !unused   integer                               :: n_lat_at     ! number of atoms slab                  |
    !unused   integer                               :: n_lay        ! number of layers slab                 |
    !unused   real(8)                               :: nn           ! next neighbour distance slab          |
    !unused   real(8), allocatable, dimension(:,:)  :: r_lat        ! lattice positions                     |
    !unused   real(8), allocatable, dimension(:,:)  :: r_part_v     ! hydrogen position                     |
    !unused   real(8)                               :: a0           ! lattice constant                      |
    !unused   integer                               :: n_lat = 2    ! number of lattice atoms               |
    !unused   real(8), dimension(3)                 :: r0_part      ! particle position for reference calc. |
    !unused r0_part = (/0.0, 0.0, 6.0/)  ! H atom at 6 Angstrom for reference energy calculatio n           |
    !-------------------------------------------------------------------------------------------------------|
