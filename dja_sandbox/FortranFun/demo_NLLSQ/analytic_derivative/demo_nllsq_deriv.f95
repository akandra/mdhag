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

    integer, parameter          :: npts=101             ! number of points to be fitted
    integer, parameter          :: nparms=4             ! number of parameters in the fitting function
    integer, parameter          :: max_iterations=30    !maximum number of iterations

    ! Usage: CALL NLLSQ ( Y , X , B , RRR , NARRAY , ARRAY , IB , TITLE)


    real(8), dimension(npts)       :: Y                ! values of dependent variable at N data points.

    real(8), dimension(npts)       :: X                ! values of indepenent variable(s) at N data points.
    !       if M is number of independent variables
    !       M=1 => DIMENSION X(N)
    !       M>1 => DIMENSION X(N,M)

    real(8), dimension(nparms)     :: B                ! parameters for the fit
    !   At entry, B's are initial guesses for parameters
    !   On exit,  B's are the results of the fit.

    real(8), dimension(nparms)     :: P
    !   Array used by NLLSQ and Model for partial derivatives wrt each parameter
    !   P is defined here because we need to have it set up in COMMON/BLK1/ to use
    !     model to calculated the values of the function

    real(8), dimension(npts)       :: RRR
    !       Array used optionally for communication between program that calls
    !       NLLSQ and SUBROUTINE MODEL.
    !       used here to communicat the final values of the fit
    !
    !       the shape of RRR must be the same here and in SUBROUTINE MODEL
    !       RRR is not used directly by NLLSQ; only passed to model

    integer, dimension(8)          :: NARRAY = (/ npts, 1, nparms, 0, 1, 2, 10, max_iterations /)
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

    real(8), dimension(8)          :: ARRAY = 0
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
    character(20)               :: TITLE    ! string of 20 characters used for title
                                            ! on intermediate and final printout.

    real(8), dimension(npts)    :: noise    ! array hold random numbers to introduce noise into
                                            ! the "data"
    integer                     :: i
    real(8), dimension(4)       :: b0,b1

    real(8)                     :: RES      ! for common/BLK1/
    integer                     :: N,M,K    ! for common/BLK1/

    common/BLK1/B,P,RES,N,M,K               ! for communication with model
                                            ! note B and P are defined above with dim(20)

    TITLE = 'Demonstrate NLLSQ'

    ! ---------------------------------------------------------------
    ! Set up an array of randeom numbers to use as "noise"
    ! ---------------------------------------------------------------
    call random_seed()
    call random_number(noise)

    ! ---------------------------------------------------------------
    ! Initialize B for subroutine model
    ! ---------------------------------------------------------------
    B(1)=1      ! A
    B(2)=.5     ! E0
    B(3)=.1     ! W
    B(4)=0      ! Y0
    b0=B        !save initial values for printing later

    ! ---------------------------------------------------------------
    ! Set x, call model and save model value + noise in array y
    ! ---------------------------------------------------------------
    N = npts
    do i=1,npts
        x(i)=real(i-1)/100
        call model(Y(i),Y,X,RRR,i,1)
        y(i)=y(i)+noise(i)*0.1
    enddo
    ! ---------------------------------------------------------------
    ! Initial guess of parameters
    ! ---------------------------------------------------------------
    B(1) = .5
    B(2) = .5
    B(3) = .5
    B(4) = .5
    b1=B       !save guess values for printing late

    ! ---------------------------------------------------------------
    ! hold E0 fixed for the fit
    ! ---------------------------------------------------------------
    NARRAY(4)= 1
    IB(1)    = 2

    ! ---------------------------------------------------------------
    ! Call subroutine NLLSQ to do the fit and print the results
    !   Unit 10: output of NLLSQ
    !   Unit 11: data for plotting
    ! ---------------------------------------------------------------
    call nllsq (Y, X, B, RRR, NARRAY, ARRAY, IB, TITLE)

    print '(/(a))', '-------------FIT RESULTS-----------------'
    print '(/(a))', '  i       x    fit data       y'
    write(* , '(i3,3f10.3)') (i, x(i), y(i), RRR(i), i=1,npts)
    write(11, '(i3,3f10.3)') (i, x(i), y(i), RRR(i), i=1,npts)

    print '(/(a),7F7.3 / 3x,f7.3)', 'start parameters=', b0
    print '(/(a),7F7.3 / 3x,f7.3)', 'guess parameters=', b1
    print '(/(a),7F7.3 / 3x,f7.3)', 'fit t parameters=', B

end program
