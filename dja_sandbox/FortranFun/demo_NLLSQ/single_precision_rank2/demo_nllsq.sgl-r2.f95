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

    real, dimension(npts)       :: Y        ! values of dependent variable at N data points.
    real, dimension(npts,2)     :: X        ! values of indepenent variable(s) at N data points.
    !real, dimension(npts)       :: X        ! values of indepenent variable(s) at N data points.

    real                        :: B, P, RES
    integer                     :: N, M, KK
    common /BLK1/ B(20),P(20),RES,N,M,KK

    real, dimension(npts)       :: RRR

   !integer, dimension(8)          :: NARRAY = (/ npts, 1, nparms, 0, 1, 2, 10, max_iterations /)
    integer, dimension(8)          :: NARRAY = (/ npts, 2, nparms, 0, 1, 2, 10, max_iterations /)
    !                                              1    2    3     4  5  6   7    8
    ! Array of integer control parameters as follows:

    real, dimension(8)          :: ARRAY = 0
    ! Array of statistical parameters. Use 0.0 to get default values.

    integer, dimension(nparms)  :: IB = 0
    ! Integer Array containing the subscripts of parameters to be held constant.
    !       (number=NARRAY(4) ).
    !
    character(20)               :: TITLE    ! string of 20 characters used for title on intermediate and final printout.
    real, dimension(npts)       :: noise    ! array hold random numbers to introduce noise into the "data"
    integer                     :: i
    real, dimension(4)          :: b0, b1   ! start parameters, guess parameters
    real                        :: F

    TITLE = 'Demonstrate NLLSQ'

    ! Set up the data arrays
    call random_seed()
    call random_number(noise)

    b0(1)=1.0
    b0(2)=0.5
    b0(3)=0.1
    b0(4)=0.0
    B(1:4)=b0   ! pass parameters to model in blk1
    N=npts      ! pass number of points to model in blk1

    do i=1,npts
       !X(i)  =real(i-1)/100.
        X(i,1)=real(i-1)/100.
        X(i,2)=real(i-1)/1000.

        call model(F, Y, X, RRR, I, 1)
        y(i)=F + noise(i)*0.1
        !write(*,'(i5,3f10.5)'),i,xdat(i,1),xdat(i,2) !,ydat(i)
        write(*,'((a)i5,3f10.5)'), 'i, x1, x2, y=', i, x(i,:), y(i)
    enddo

    ! Initial guess of parameters
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
