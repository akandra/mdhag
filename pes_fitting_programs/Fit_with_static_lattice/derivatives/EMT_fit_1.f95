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
    use emt_init_data

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
    real(8) :: X(1000,3,1000),Y(1000),RRR(1000)
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
    type (EMTparms) :: part_new ! EMT parameters of particle
    type (EMTparms) :: latt_new ! EMT parameters of lattice

    namelist / lattice_pars_list_new / latt_new ! Namelists to handle I/O of parameters
    namelist / particle_pars_list_new / part_new

    namelist / lattice_pars_list / lattice_pars ! Namelists to handle I/O of parameters
    namelist / particle_pars_list / particle_pars

    integer :: n_l_in ! number of atoms in reference slab
    integer :: n_lay0 ! number of layers in reference slab
    real(8) :: nn0 ! next neighbour distance in reference slab
    real(8), dimension(3) :: cell_in

    real(8), allocatable, dimension(:,:) :: r1 ! lattice positions for reference calc.

    integer site(1000)
    real(8) :: energy ! energy output from emt subroutines
    real(8), dimension(14) :: denergy, dE_ref
    real(8) :: e_ref ! reference energy with particle at infinity
    real(8) :: e_max ! maximum DFT to use in the fit
    real(8) :: sumsq,se ! used to calculate rms error

    integer :: i,j,k ! loop indicies
    integer :: ios ! io status

! For AIMD readin
!    integer                  :: time            ! number of different configurations
    real(8), dimension(:,:,:), allocatable:: x_all,x_ball      ! Position of all the atoms. Is the same as x(:)
    real(8),allocatable, dimension(:)     :: E_all    ! energy from AIMD and DFT calc
    !integer                 :: rep      ! Repetitions of lattice (1,2,3) 2 should be sufficient.
    integer, dimension(3):: cell_b  ! contains cell geometry (2x2x4)
!    real(8), dimension(:,:), allocatable :: celli       ! extension of cell once it's expanded
!    integer :: n_l ! number of gold atoms in extended slab
    integer :: l_aimd ! number of aimd configurations in  all the configurations
    integer :: control  ! controls which fraction of aimd-atoms is used (1-100)
                        ! control = 200 : uses only frozen lattice DFT-points
                        ! control = 201 : uses only aimd-points
    real(8), dimension(:,:,:), allocatable :: r_l
    real(8), dimension(:,:), allocatable :: r_p
    integer :: q
    real(8) :: e_aimd_max
    real(8), dimension(:,:), allocatable :: check


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
    call open_for_overwrite(7,'data/fit_debug.out')


    lattice_configuration_fname = 'data/ref_conf_Au111a.dat'
    particle_position_and_DFT_energies_fname= 'data/hau111_plot.E.dat'


! Strömqvist parameters modified in so, so they'll give a good fit.
    fit_results_fname = 'data/parameters_and_fit_results/stroem_der.237.NLLSQ.out'
    particle_nml_out  = 'data/parameters_and_fit_results/stroem_der.237.H.nml'
    lattice_nml_out   = 'data/parameters_and_fit_results/stroem_der.237.Au.nml'

    particle_nml_in = 'data/parameters_and_fit_results/stroem_der.233.H.nml' !stroem.00.H.nml'
    lattice_nml_in  = 'data/parameters_and_fit_results/stroem_der.233.Au.nml' !stroem.00.Au.nml'



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
    read(8,*) n_l_in
    read(8,*) n_lay0
    read(8,*) nn0
    a_lat = nn0*sqrt_2

!    call gold_pos(r1, n_l_in, r0_lat, cell_0)

    ! routine gets the gold positions set in case they differ from positions of reference
    ! system

    !------------------------------------------------------------------------------------
    !                   GET INFORMATION OF ATOMS THAT SHALL BE FITTED
    !                   =============================================
    !------------------------------------------------------------------------------------
    !
    ! Usage of input parameters:
    ! =========================
    ! rep = 1, 2 or 3   How many times shall the cell from DFT be reproduced.
    !                   Recommendation: rep=2
    ! cell_b (/2,2,4/)  DFT-cell-geometry
    ! control =         1-100 : Fraction of AIMD points that go into fit
    !                   200   : only DFT points
    !                   201   : only AIMD points

    rep = 1
    cell_b=(/2,2,4/)
    control=2
    e_aimd_max=0.00
    just_l = .false.
    one_p = .false.


    call l_p_position(a_lat, rep, cell_b, control, e_aimd_max, just_l, one_p, time, &
                                                l_aimd, n_l,n_p,celli, x_all, E_all)

    print *, 'l_aimd', l_aimd

    !time = time-1
    X(1:time,:,1:n_l+n_p)=x_all(1:time,:,1:n_l+n_p)
    Y(1:time)=E_all(1:time)!+E_pseudo
    !write(*,'(3f15.5)') x(47,:,1:n_l+1)
    !stop

!    allocate(x_ball(time,3,n_l+1))


    call open_for_write(10, fit_results_fname)

!    write(*,*) 'q   ', 'EDFT-E_multiple_H   ', 'EDFT-E_single_H'
!    do q=1,30
!        x_ball(q,:,1) = x_all(q,:,5)
!        X(q,:,1)= x_all(q,:,5)
!        x_ball(q,:,2:n_l+1) = x_all(q,:,n_p+1:n_p+n_l)
!        X(q,:,2:n_l+1) = x_all(q,:,n_p+1:n_p+n_l)
!        call emt(a_lat, celli, x_all(q,:,:), n_l, n_p, particle_pars, lattice_pars, energy)
!        call emt_mixed(a_lat, celli, x_all(q,:,5), x_all(q,:,n_p+1:n_p+n_l), n_l, particle_pars, lattice_pars, e_ref)
!        write(*,'(I5,5f15.5)') q, x_all(q,3,5), energy-18.14163, e_ref-18.14163, Y(q)+24.99569, Y(q)* (2*rep+1)**2+224.9612
!    end do

    !------------------------------------------------------------------------------------------------------------------
    ! CHECK EMT POTENTIAL SUBROUTINE AND WRITE RESULTS
    !------------------------------------------------------------------------------------------------------------------

    k = time
    npts=time

    if(k>0) then
        write(*,'(/(a))')'CHECK EMT ENERGY CALCULATION IS WORKING'
        write(*,*) 'site X Y Z EMT DFT'
        Write(10,*) 'site X Y Z EMT DFT'
        sumsq = 0
        call emt_fit(a_lat, celli, x_ref, n_l, n_p, particle_pars, lattice_pars, e_ref, Ablei)
        write(*,'(7f15.5)'), Ablei




        do q=1,time
            call emt(a_lat, celli, x_all(q,:,:), n_l, n_p, particle_pars, lattice_pars, energy)


            if(q<10) write( *,'(1X, 5F15.8)') X(q,1,1), X(q,2,1), X(q,3,1), (energy-e_ref)/(2*rep+1)**2, Y(q)
            write(10,'(1X, 5F15.8)')  X(q,1,1), X(q,2,1), X(q,3,1), (energy-e_ref)/(2*rep+1)**2, Y(q)
            write(7, '(1X, 4F16.8)')  X(q,1,1), X(q,2,1), X(q,3,1), (energy-e_ref)/(2*rep+1)**2
            if (q > 483) then
            !write( *,'(1X, 5F15.10)') X(q,1,1), X(q,2,1), X(q,3,1), energy-E_ref, Y(q)

            end if
            sumsq=sumsq+((energy-e_ref)/(2*rep+1)**2-Y(q))**2
            !sumsq=sumsq+(energy-e_ref-Y(q))**2
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
    ! Name part lat H Au (be careful! Changed with regard to old procedure!)
    !
    ! eta2      1 8
    ! n0        2 9
    ! E0        3 10
    ! lambda    4 11
    ! V0        5 12    x x shouldn't be <0
    ! kappa     6 13
    ! s0        7 14    x x shouldn't change
    IB = (/3,7,10,14,0,0,0,0,0,0,0,0,0,0/) ! indicies of parameters held constant
    IP = 4 ! number of parameters held constant


    !--------------------------------------------------------------------------
    ! SET UP NARRAY
    !--------------------------------------------------------------------------
    nparms = 14
    max_iterations = 200
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
    ARRAY(1) = 0.1d0 ! AL .1
    ARRAY(2) = 0.00001d0 ! DELTA - FOR DERIVATIVES .00001
    ARRAY(3) = 0.00005d0 ! E - CONVERGENCE CRITERION .00005
    ARRAY(4) = 4.0d0 ! FF 4.0
    ARRAY(5) =45.0d0 ! GAMCR - CRITICAL ANGLE 45.
    ARRAY(6) = 2.0d0 ! T 2.
    ARRAY(7) = 0.001d0 ! TAU .001
    ARRAY(8) = 1d-31 ! ZETA - CRITERION FOR 1E-31
                        ! SINGULAR MATRIX
    ARRAY=0 ! use default values


    TITLE = 'Test of EMT_fit_1'

    if(debug(2)) then
print '(//(a))', 'BEFORE CALL TO NLLSQ'
        print '((a),4f10.5)', ' x1,y1= ', X(1,1,1), X(1,2,1), X(1,3,1), Y(1)
        print '((a),4f10.5)', ' x2,y2= ', X(2,1,1), X(2,2,1), X(2,3,1), Y(2)
        print '((a),3I20)', ' loc(x1))=', loc(x(1,1,1)), loc(x(1,2,1)), loc(x(1,3,1))
        print '((a),3I20)', ' loc(x2))=', loc(x(2,1,1)), loc(x(2,2,1)), loc(x(2,3,1))
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

    !deallocate(r1)

    1000 format(2x,a15,2x,a2,7F7.3)
    1010 format(2x,a15,4x,7F7.3,/,19x,7f7.3)


    call array2emt_parms( B(1:7 ), part_new)
    call array2emt_parms( B(8:14), latt_new )

    sumsq=0.0d0
    se=0.0d0
    call emt(a_lat, celli, x_ref, n_l, n_p, part_new, latt_new, e_ref)
    call open_for_write(1234,'deviation_199.dat')
    do q=1,time
            call emt(a_lat, celli, x_all(q,:,:), n_l, n_p, part_new, latt_new, energy)
            write(1234,'(1X, I5, 5F15.8)') q, energy-e_ref-Y(q)
            sumsq=sumsq+((energy-e_ref)/(2*rep+1)**2-Y(q))**2
            se = se+Sqrt(((energy-e_ref)/(2*rep+1)**2-Y(q))**2)
    end do
    print *, 'rms =',  sqrt(sumsq/time)*1000, 'meV'
    print*,  'SE =', (se/time)*1000, 'meV'
    write(*,'(7f15.5)'), Ablei

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
