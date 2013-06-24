!----------------------------------------------------------------------------------------------------------------------
!
!   File:       demoNLLSQ.f95
!
!   Purpose:    Demonstrate use of NLLSQ,the Bell Labs modification of SHARE #3094 NLIN by D. W. Marquardt
!
!   Author:     D.J. Auerbach
!   Date:       22 April, 2013
!   Copyright:  2013 Max Planck Institute for Biophysical Chemistry
!
!   Description:
!       1) Allocate arrays needed by the NLLSQ routine.
!
!       2) Populate these arrays with some test data.  In the initial version the data is calculated from a
!          2nd order polynomial with small arbitrary errors.  Also populate the arrays that control the NLLSQ program
!
!       3) Call NLLSQ
!
!       4) Print results
!
!   The project needs 2 additional files:
!       NLLSQ.f     the non linear least squares routine
!       model.f95   subroutine to implement the model being fitted
!
!       Communication betwen NLLSQ and model is via common blocks and based arguments
!
!   Documentation:
!       Brief documentation is given in NLLSQDOC.txt
!----------------------------------------------------------------------------------------------------------------------
!   Revision Log
!
!   Date        Author      Change
!   22.04.2013  dja         Initial release
!
!
!
!----------------------------------------------------------------------------------------------------------------------
program Demo_NLLSQ
    implicit none

    logical             :: debug, first

    integer, parameter  :: npts=100     ! number of points to be fitted
    integer, parameter  :: n_ind=2      ! number of independent variables
    integer, parameter  :: nparms=4     ! number of parameters in the fitting function
    integer             :: max_iter     ! maximum number of iteration
    integer             :: IP           ! number of parameters to be held constant


    !--------------------- VARIABLES FOR SUBROUTINE NLLSQ -------------------------------------------------------------
    ! Usage: CALL NLLSQ ( Y , X , B , RRR , NARRAY , ARRAY , IB , TITLE)
    !------------------------------------------------------------------------------------------------------------------
    real(8)         :: Y(npts)=0        ! Values of dependent variable at npts data points.
    real(8)         :: X(npts,n_ind)=0  ! Values of n_ind indepenent variables at npts data points.
                                        !   for N points and M independent variables
                                        !   M=1 -> dimension B(N)
                                        !   M>N    Dimnesion B(N.M), for M>1
    real(8)         :: B=0              ! Array of parameters, dimmension defined in common/BLK1/
                                        !   at entry, B gives the initial guesses for parameters
                                        !   at exit, B gives the value of the fit parameters
    real(8)         :: RRR(npts)        ! Array for final result from model
    integer         :: NARRAY(8)        ! Integer array of control parameters, defined in detail below
    real(8)         :: ARRAY(8)         ! Array of statistical control parameters, defined in detail below
    integer         :: IB(20)           ! Integer Array containing the subscripts of parameters to be held constant.

    character(20)   :: TITLE            ! string used for title on intermediate and final printout


    !--------------------- COMMON BLOCKS FOR COMMUNICATION WITH MODEL -------------------------------------------------
    real(8)  P, RES
    integer N, M                        ! number of data points, independed variables, and iteration returned by model
    common/BLK1/ B(20),P(20),RES,N,M
    common/DJA1/ debug(5), first


    !--------------------- VARIABLES USED TO DEFINE FIT DATA ----------------------------------------------------------
    real(8)         :: noise(npts)      ! array hold random numbers to introduce noise into the "data"
    integer         :: i
    real(8)         :: b0(4), b1(4)     ! start parameters, guess parameters
    real(8)         :: F


    debug=(/.false.,.false.,.false.,.false.,.false./)
    first=.true.
    !--------------------------------------------------------------------------
    !                  SET UP NARRAY
    !--------------------------------------------------------------------------
    IP       = 0                ! number of parameters to be held constant during fit.
    max_iter = 20

    NARRAY(1) = npts            ! number of data points
    NARRAY(2) = 3               ! number of independent variables (cartesian coordinates)
    NARRAY(3) = nparms          ! number of parameters
    NARRAY(4) = IP              ! number of parameters to be held constant during fit.
    NARRAY(5) = 1               ! intermediate printout option (0..3)
    NARRAY(6) = 1               ! final printout option (0..3)
    NARRAY(7) = 10              ! printout unit number
    NARRAY(8) = max_iter        ! maximum number of iterations
    !***********************  ARRAY ***************************
    !                     variable                      default
    ARRAY(1)  = 0.1     ! AL                             .1
    ARRAY(2)  = 0.00005 ! DELTA - FOR DERIVATIVES        .00001
    ARRAY(3)  = 0.00005 ! E     - CONVERGENCE CRITERION  .00005
    ARRAY(4)  = 4.0     ! FF                             4.0
    ARRAY(5)  = 45.0    ! GAMCR - CRITICAL ANGLE         45.
    ARRAY(6)  = 2.0     ! T                              2.
    ARRAY(7)  = 0.001   ! TAU                            .001
    ARRAY(8)  = 1E-31   ! ZETA  - CRITERION FOR          1E-31
                        !         SINGULAR MATRIX


    TITLE = 'Demonstrate NLLSQ'

    !--------------------- SET UP THE DATA ARRAYS -----------------------------
    call random_seed()
    call random_number(noise)

    !--------------------- PARAMETERS USED TO DEFINE DATA ---------------------
    b0(1)=1.0   ! A
    b0(2)=0.5   ! E
    b0(3)=0.1   ! W
    b0(4)=0.0   ! Y0


    !---------------- COMMON/BLK1 VARIABLE FOR COMMUNICATION TO MODEL ---------
    B(1:4)=b0   ! pass parameters to model in blk1
    N=npts      ! pass number of points to model in blk1
    M=n_ind     ! pass number of independent variable to model

    !---------------- DEFINE FIT DATA USING MODEL -----------------------------
    do i=1,npts
        x(i,1)=real(i)/100.
        x(i,2)=real(i)/1000.
    end do
    ! 2 do loops used so that x is completely initialized before first call to model
    do i=1,npts
        call model(F, Y, X, RRR, I, 1)
        y(i)=F !+ noise(i)*0.1
    enddo

    if(debug(1)) then
        print '(//(a))', '-------------MAIN AFTER DATA SETUP-------------'
        print '((a),4f10.5)', '  x1,y1=   ', X(1,1), X(1,2), Y(1)
        print '((a),4f10.5)', '  x2,y2=   ', X(2,1), X(2,2), Y(2)
        print '((a),3I10)',   '  loc(x1))=', loc(x(1,1)), loc(x(1,2))
        print '((a),3I10)',   '  loc(x2))=', loc(x(2,1)), loc(x(2,2))
        print '((a),8i5)',    '  NARRAY=',NARRAY
        print '((a),7f7.3/4x,7f7.3/)','  B=',B(1:14)
        write(*,'((a)/(i5,3f10.5))'), '    i     x1        x2        y', (i, x(i,:), y(i),i=1,npts)
        print '((a)//)', '-----------END MAIN AFTER DATA SETUP-----------'
    end if

! PAUSE 'before call to NLLSQ'





    !--------------- INITIAL GUESS OF VALUES OF PARAMETRS --------------------
    b1(1) = .5
    b1(2) = .3
    b1(3) = .5
    b1(4) = .1
    B(1:4)= b1

    call nllsq (Y, X, B, RRR, NARRAY, ARRAY, IB, TITLE)

    write( *, 1000) (i, x(i,:), y(i), RRR(i), i=1,npts)

    write(*,1010) 'start parameters=', b0
    write(*,1010) 'guess parameters=', b1
    write(*,1010) 'fit   parameters=', B

1000 format(i3,4f10.3)
1010 format((a), 4f10.4)

end program
