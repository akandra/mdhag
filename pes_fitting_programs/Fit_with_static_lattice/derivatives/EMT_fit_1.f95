program EMT_fit_1

    !----------------------------------------------------------------------------------------------------------------------
    !
    ! File: EMT_fit_1.f95
    !
    ! Purpose: Nonlinear Least Squares fit of parameters of EMT potential to DFT points using static lattice
    !
    ! Author: D.J. Auerbach
    ! Date: 29 April, 2013
    ! Copyright: 2013 Max Planck Institute for Biophysical Chemistry
    !
    ! Uses: NLLSQ - a non linear least squres routine developed by D. Marquardt and modified at Bell Labs.
    ! EMT_parms class developed by S. Janka and A. Kandrasenka. Modified here for ease of use with NLLSQ
    !
    ! Description:
    !
    ! Additional Files:
    ! NLLSQ.f the non linear least squares routine
    ! model.f95 subroutine to implement the model being fitted
    ! EMTparms_class_v2.f95 subroutines and structures for EMT energy calculations
    ! open_file.f95 module with routines to open files for reading and writing with error checking
    !

    !----------------------------------------------------------------------------------------------------------------------
    ! Revision Log
    !
    ! Date Author Change
    ! 30.04.2013 dja Initial release
    ! 03.06.2013 smj derivatives into fit
    !
    !
    !
    !----------------------------------------------------------------------------------------------------------------------

    use EMTparms_class
    use open_file

    implicit none



    !------------------------------------------------------------------------------------------------------------------
    !
    ! NLLSQ VARIABLES
    !
    !------------------------------------------------------------------------------------------------------------------
    ! todo make the arrays allocatable and automatically determine npts by reading to end of data
    !------------------------------------------------------------------------------------------------------------------
    integer :: npts ! number of points to be fitted
    integer :: nparms ! number of parameters in the fitting function
    integer :: IP ! number of parameters to be held constant during fit
    integer :: max_iterations ! maximum number of iterations
    !real(8), allocatable :: Y(:) ! values of dependent variable at N data points.
    !real(8), allocatable :: X(:,:) ! values of indepenent variable(s) at N data points.
    real(8), dimension(20) :: B ! parameters for the fit. 20 is the maximum hard coded in NLLSQ
    ! At entry, B's are initial guesses for parameters
    ! On exit, B's are the results of the fit.
    !real(8), allocatable :: RRR(:) ! Array used optionally for communication with MODEL
    real(8) :: X(1000,3),Y(1000),RRR(1000)
    integer, dimension(8) :: NARRAY ! Integer array of control parameters
    real, dimension(8) :: ARRAY = 0 ! Array of statistical parameters. Use 0.0 to get default values.
    integer, dimension(14) :: IB = 0 ! Integer Array containing the subscripts of parameters to be held constant.
    ! number of parameters held constant = NARRAY(4)
    character(20) :: TITLE ! string of 20 characters used for title on printout



    !------------------------------------------------------------------------------------------------------------------
    !
    ! EMT VARIABLES
    !
    !------------------------------------------------------------------------------------------------------------------
    type (EMTparms) :: particle_pars ! EMT parameters of particle
    type (EMTparms) :: lattice_pars ! EMT parameters of lattice
    type (EMTparms) :: test_pars ! used for check conversion routines

    namelist / lattice_pars_list / lattice_pars ! Namelists to handle I/O of parameters
    namelist / particle_pars_list / particle_pars

    integer :: n_lat0_at ! number of atoms in reference slab
    integer :: n_lay0 ! number of layers in reference slab
    real(8) :: nn0 ! next neighbour distance in reference slab
! real(8) :: a_lat ! lattice constant
    real(8), dimension(3) :: cell ! dimensions of the cell
    real(8), allocatable, dimension(:,:) :: r0_lat ! lattice positions for reference calc.
    !integer, allocatable, dimension(:) :: site ! impact site
    integer site(1000)
    real(8) :: energy ! energy output from emt subroutines
    real(8) :: e_ref ! reference energy with particle at infinity
    real(8) :: e_max ! maximum DFT to use in the fit
    real(8) :: sumsq ! used to calculate rms error

    integer :: i,j,k ! loop indicies
    integer :: ios ! io status

    ! file names
    character(len=100) :: lattice_configuration_fname
    character(len=100) :: particle_position_and_DFT_energies_fname
    character(len=100) :: fit_results_fname
    character(len=100) :: particle_nml_in
    character(len=100) :: particle_nml_out
    character(len=100) :: lattice_nml_in
    character(len=100) :: lattice_nml_out

    logical debug
    common /debug/debug(5)
    debug=(/.false., .false., .false., .true., .false./)
    call open_for_overwrite(7,'fit_debug.out')


    lattice_configuration_fname = 'ref_conf_Au111a.dat'
    particle_position_and_DFT_energies_fname= 'hau111_plot.E.dat'

! These paths are for windows, they don't work under linux
! fit_results_fname = 'parameters_and_fit_results\f119.02.NLLSQ.out'

! particle_nml_in = 'parameters_and_fit_results\f119.00.H.nml'
! particle_nml_out = 'parameters_and_fit_results\f119.02.H.nml'

! lattice_nml_in = 'parameters_and_fit_results\f119.00.Au.nml'
! lattice_nml_out = 'parameters_and_fit_results\f119.02.Au.nml'

! These paths are for linux, they probably won't work under windows.
! fit_results_fname = 'parameters_and_fit_results/f119.02.NLLSQ.out'

! particle_nml_in = 'parameters_and_fit_results/f119.00.H.nml'
! particle_nml_out = 'parameters_and_fit_results/f119.02.H.nml'

! lattice_nml_in = 'parameters_and_fit_results/f119.00.Au.nml'
! lattice_nml_out = 'parameters_and_fit_results/f119.02.Au.nml'

! Strömqvist
! Strömqvist parameters modified in so, so they'll give a good fit.
    fit_results_fname = 'parameters_and_fit_results/stroem.02.NLLSQ.out'
    particle_nml_in = 'parameters_and_fit_results/stroem.00.H.nml'
    particle_nml_out = 'parameters_and_fit_results/stroem.02.H.nml'
    lattice_nml_in = 'parameters_and_fit_results/stroem.00.Au.nml'
    lattice_nml_out = 'parameters_and_fit_results/stroem.02.Au.nml'


    TITLE = 'EMT NLLSQ Test'

    !------------------------------------------------------------------------------------------------------------------
    ! READ LATTICE AND PARTICLE EMT PARAMETERS
    !------------------------------------------------------------------------------------------------------------------
    call open_for_read(8,lattice_nml_in)
    read(8,nml=lattice_pars_list)
    write(7,*) 'lattice_pars_list'
    write(7,nml=lattice_pars_list)

    call open_for_read(8,particle_nml_in)
    read(8,nml=particle_pars_list)
    write(7,*) 'particle_pars_list'
    write(7,nml=particle_pars_list)



    !------------------------------------------------------------------------------------------------------------------
    ! READ LATTICE DEFINITION PARAMETERS
    !------------------------------------------------------------------------------------------------------------------
    call open_for_read (8, lattice_configuration_fname)
    read(8,*) n_lat0_at
    read(8,*) n_lay0
    read(8,*) nn0
    read(8,*) cell(1)
    read(8,*) cell(2)
    allocate(r0_lat(3,n_lat0_at)) ! allocate array to hold lattice coordinates
    readr0lat: do i = 1, n_lat0_at
        read(8,*) r0_lat(1,i), r0_lat(2,i), r0_lat(3,i)
        !write(7,*)r0_lat(1,i), r0_lat(2,i), r0_lat(3,i)
    end do readr0lat
    a_lat = nn0*sqrt_2

    !------------------------------------------------------------------------------------
    ! INITIALIZE EMT POTENTIAL SUBROUTINE AND CALCULATE REFERENCE ENERGY
    !------------------------------------------------------------------------------------
    call emt_init(a_lat, cell, n_lat0_at, r0_lat, particle_pars, lattice_pars, E_ref)
    write(*,'(//(a),F9.5,(a)//)') 'the reference energy = ',E_ref, ' eV'
    write(7,'(//(a),F9.5,(a)//)') 'the reference energy = ',E_ref, ' eV'
    !------------------------------------------------------------------------------------
    ! READ THE PARTICLE POSITIONS AND DFT ENERGIES
    !------------------------------------------------------------------------------------
    ! FIRST FIND THE NUMBER OF POINTS BY READING TO EOF
    !------------------------------------------------------------------------------------
    call open_for_read (8, particle_position_and_DFT_energies_fname)
    i=1
    do
read(8,*,iostat=ios)
        if(ios <0) exit
i=i+1
    end do

rewind(8)
    npts = i-1
    print '((a),i4)','the number of particle positions and energies =',npts


    !------------------------------------------------------------------------------------
    ! NOW ALLOCATE ARRAYS TO STORE THE POINTS
    !------------------------------------------------------------------------------------
    ! now allocate arrays needed by NLLSQ
    !allocate(X(npts,3)) ! particle coordinates - indepenet variables
    !allocate(Y(npts)) ! DFT energies - depenent variables
    !allocate(RRR(npts)) ! Array for communication with subroutine model
    !allocate(site(npts)) ! Impact site

    !------------------------------------------------------------------------------------
    ! READ IN POINTS, KEEPING ONLY POINTS WITH ENERGY < E_MAX
    !------------------------------------------------------------------------------------
    j=1
    e_max = 10 ! maximum energy to fit in eV

    do i=1, npts
        read(8,*) site(j), X(j,1), X(j,2), X(j,3), Y(j)
     ! write(7,*)i,j,site(j), X(j,1), X(j,2), X(j,3), Y(j)
        if ( (abs(Y(j))<=10) .and. (i<npts) ) j=j+1
    end do
close(8)

    npts=j

    !write(7,'(//(a)/)') 'filtered data list'
    !do i=1,npts
    ! write(7,*)i,site(j), X(i,1), X(i,2), X(i,3), Y(i)
    !end do

    call open_for_write(10, fit_results_fname)


    !*********************************************************************************
    !npts = 50 !************************** FOR TESTING **************************
    !*********************************************************************************

    !------------------------------------------------------------------------------------------------------------------
    ! CHECK EMT POTENTIAL SUBROUTINE AND WRITE RESULTS
    !------------------------------------------------------------------------------------------------------------------
    k = npts

    if(k>0) then
write(*,'(/(a))')'CHECK EMT ENERGY CALCULATION IS WORKING'
        write(*,*) 'site X Y Z EMT DFT'
        Write(10,*) 'site X Y Z EMT DFT'
        sumsq = 0
        do i = 1, k
            call emt(a_lat, X(i,:), particle_pars, lattice_pars, energy)
            if(i<10) write( *,'(1X, I2, 5F15.10)') site(i), X(i,1), X(i,2), X(i,3), energy, Y(i)
            write(10,'(1X, I2, 5F15.10)') site(i), X(i,1), X(i,2), X(i,3), energy, Y(i)
            write(7, '(1X, I2, 4F16.10)') site(i), X(i,1), X(i,2), X(i,3), energy

            sumsq=sumsq+(energy-Y(i))**2
        end do
write(*,*)
        write(10,*)
        write(*,*) 'rms error using starting parameters =',sqrt(sumsq/npts), ' Eref=', E_ref
        write(10,*) 'rms error using starting parameters =',sqrt(sumsq/npts), ' Eref=', E_ref
        write(*,*)
        write(10,*)
    end if
! Here, the fitting procedure starts. So, for debugging, you might want to comment in the 'stop' .
!stop
    !------------------------------------------------------------------------------------------------------------------
    ! SETUP FOR NLLSQ
    !------------------------------------------------------------------------------------------------------------------
    ! Usage: CALL NLLSQ ( Y , X , B , RRR , NARRAY , ARRAY , IB , TITLE)

    !------------------ SET UP B-----------------
    call emt_parms2array(particle_pars,B(1:7))
    call emt_parms2array(lattice_pars ,B(8:14))

    !-----------CHECK WE GOT IT RIGHT -----------
    if (debug(1)) then
write(*,'(/(a))')'CHECK CONVERSION SUBROUTINE emt_parms2array'
        write(*,1000) 'particle_pars=',particle_pars
        write(*,1010) 'B(1:7) =',B(1:7)
        write(*,1000) 'lattice_pars =',lattice_pars
        write(*,1010) 'B(8:14) =',B(8:14)

        !---------------TEST ARRAY2EMT_PARMS--------
        call array2emt_parms(B(1:7),test_pars)
        write(*,'(/(a))') 'CHECK CONVERSION SUBROUTINE array2emt_parms'
        write(*,'( (a))') ' convert B(1:7) to emt_parms structure and print'
        write(*,1000) 'test_pars =',test_pars
    end if

    !--------------------------------------------------------------------------
    ! SET UP PARAMETERS TO HOLD CONSTANT
    !--------------------------------------------------------------------------
    ! Name part lat H Au
    !
    ! eta2 1 8
    ! kappa 2 9
    ! lambda 3 10
    ! E0 4 11 x x destabilizes fit
    ! n0 5 12 ? ? may destabilize fit
    ! s0 6 13 x x shouldn't change
    ! V0 7 14 x shouldn't be <0

    IB = (/6,13,14,4,11,0,0,0,0,0,0,0,0,0/) ! indicies of parameters held constant
    IP = 5 ! number of parameters held constant

    !--------------------------------------------------------------------------
    ! SET UP NARRAY
    !--------------------------------------------------------------------------
    nparms = 14
    max_iterations = 30
    NARRAY(1) = npts ! number of data points
    NARRAY(2) = 3 ! number of independent variables (cartesian coordinates)
    NARRAY(3) = nparms ! number of parameters
    NARRAY(4) = IP ! number of parameters to be held constant during fit.
    NARRAY(5) = 1 ! intermediate printout option (0..3)
    NARRAY(6) = 1 ! final printout option (0..3)
    NARRAY(7) = 10 ! printout unit number
    NARRAY(8) = max_iterations ! maximum number of iterations
    !*********************** ARRAY ***************************
    ! variable default
    ARRAY(1) = 0.1 ! AL .1
    ARRAY(2) = 0.00001 ! DELTA - FOR DERIVATIVES .00001
    ARRAY(3) = 0.00005 ! E - CONVERGENCE CRITERION .00005
    ARRAY(4) = 4.0 ! FF 4.0
    ARRAY(5) =45.0 ! GAMCR - CRITICAL ANGLE 45.
    ARRAY(6) = 2.0 ! T 2.
    ARRAY(7) = 0.001 ! TAU .001
    ARRAY(8) = 1E-31 ! ZETA - CRITERION FOR 1E-31
                        ! SINGULAR MATRIX
    ARRAY=0 ! use default values


    TITLE = 'Test of EMT_fit_1'

    if(debug(2)) then
print '(//(a))', 'BEFORE CALL TO NLLSQ'
        print '((a),4f10.5)', ' x1,y1= ', X(1,1), X(1,2), X(1,3), Y(1)
        print '((a),4f10.5)', ' x2,y2= ', X(2,1), X(2,2), X(2,3), Y(2)
        print '((a),3I10)', ' loc(x1))=', loc(x(1,1)), loc(x(1,2)), loc(x(1,3))
        print '((a),3I10)', ' loc(x2))=', loc(x(2,1)), loc(x(2,2)), loc(x(2,3))
        print '((a),8i5)', ' NARRAY=',NARRAY
        print '((a),7f7.3/8x,7f7.3/)',' B=',B(1:14)
    end if

CALL NLLSQ ( Y , X , B , RRR , NARRAY , ARRAY , IB , TITLE)

    call array2emt_parms(B(1:07),particle_pars)
    call array2emt_parms(B(8:14),lattice_pars)
    call open_for_write(11,particle_nml_out)
    call open_for_write(12,lattice_nml_out)
    write (11,nml=particle_pars_list)
    write (12,nml=lattice_pars_list)


    1000 format(2x,a15,2x,a2,7F7.3)
    1010 format(2x,a15,4x,7F7.3,/,19x,7f7.3)

end program
!-------------------------------------------------------------------------------------------------------|
! Variables that were declared in testEMT.f95 but not usef |
!-------------------------------------------------------------------------------------------------------|
!unused integer :: n_lat_at ! number of atoms slab |
!unused integer :: n_lay ! number of layers slab |
!unused real(8) :: nn ! next neighbour distance slab |
!unused real(8), allocatable, dimension(:,:) :: r_lat ! lattice positions |
!unused real(8), allocatable, dimension(:,:) :: r_part_v ! hydrogen position |
!unused real(8) :: a0 ! lattice constant |
!unused integer :: n_lat = 2 ! number of lattice atoms |
!unused real(8), dimension(3) :: r0_part ! particle position for reference calc. |
!unused r0_part = (/0.0, 0.0, 6.0/) ! H atom at 6 Angstrom for reference energy calculatio n |
!-------------------------------------------------------------------------------------------------------|
