subroutine model( F, YDAT, XDAT, RRR, I, JP)
    use EMTparms_class
    use emt_init_data

    implicit none

    logical debug

    type(EMTparms)      :: particle_parms   ! parameters of particle
    type(EMTparms)      :: lattice_parms    ! parameters of lattice atoms
    real(8)             :: F
    real(8)             :: YDAT(1000)   !YDAT(N)
    real(8)             :: XDAT(1000,3,1000) !XDAT(N,M)
    real(8)             :: RRR(1000)
    integer             :: I, JP,ij

    real(8)             :: B, P, RES
    integer             :: IP,IB,J
    integer             :: N, M, KK

    integer             :: iteration
    real(8)             :: energy,  E_dref
    real(8), dimension(14):: denergy
    real(8), dimension(14) ::dE_ref
    real(8)             :: Enew,dEnew(14)
    real(8)             :: c44, para, diff

    logical, save       :: first_run=.true.

    real(8), dimension(:,:,:), allocatable :: r_l
    real(8), dimension(:,:), allocatable   :: r_p
    real(8), dimension(:,:), allocatable   :: denergy1
    real(8) :: Eref
    integer :: q


    COMMON /BLK1/ B(20),P(20),RES,N,M,KK
    COMMON/ DJA1/ iteration
    COMMON/BLK5/IB(20),IP
    COMMON/debug/ debug(5)

    if (debug(3) .and. first_run) then
        print *
        print '((a))', 'FIRST RUN OF MODEL'
        print '((a),4i5)', '  i jp n m=', i, jp, n, m
        print '((a),4f20.5)', '  xi,yi=   ', xdat(i,1,1), xdat(i,2,1), xdat(i,3,1), ydat(i)
        print '((a),3I20)',   '  loc(xi))=', loc(xdat(i,1,1)), loc(xdat(i,2,1)), loc(xdat(i,3,1))
        print '((a),3I20)',   '  loc(x2))=', loc(xdat(2,1,1)), loc(xdat(2,2,1)), loc(xdat(2,3,1))
        print *
        write(7,*) 'it i jp    x       y      z       DFT     EMT    RES'
        !pause 'first run of model';
    end if

    call array2emt_parms( B(1:7 ), particle_parms)
    call array2emt_parms( B(8:14), lattice_parms )
!    r_part=XDAT(i,:)

    allocate(r_l(time,3,n_l))
    allocate(r_p(time,3))
    allocate(denergy1(time,14))

    diff=1.0d0/(2.0d0*rep+1.0d0)**2.0d0


!    ! This loop calculates the parameters that are dependend on the fitting parameter. It should only
!    ! be enacted when all but one parameter are held fixed.
!    ! This is for fitting the gold surface only.
!    c44=2.6210d0*10**9      !C44-constant of Au
!    para=2.35890d0       ! c44 times parameter=9*10^-10
!
!    if (IP == 13) then
!
!        if (105-sum(IB) == 8) then
!            ! correct kappa
!            B(13) = B(8)*beta/1.07
!            ! correct lambda
!            B(11) = B(13)/1.446
!            ! correct vo
!            B(12) = para*B(14)/(B(13)*B(8)*beta-B(13)**2)
!        end if
!        if (105-sum(IB) == 11) then
!            ! correct kappa
!            B(13) = B(11)*1.446
!            ! correct eta
!            B(8) = B(13)*1.07/beta
!            ! correct vo
!            B(12) = para*B(14)/(B(13)*B(8)*beta-B(13)**2)
!
!        end if
!        if (105-sum(IB) == 12) then
!            ! correct eta
!            B(8) = Sqrt( para*B(14)*1.07/(B(12)*beta*beta*(1-1/1.07)) )
!            ! correct kappa
!            B(13) = B(8)*beta/1.07
!            ! correct lambda
!            B(11) = B(13)/1.446
!
!        end if
!        if (105-sum(IB) == 13) then
!            ! correct eta
!            B(8) = B(13)*1.07/beta
!            ! correct lambda
!            B(11) = B(13)/1.446
!            ! correct vo
!            B(12) = para*B(14)/(B(13)*B(8)*beta-B(13)**2)
!        end if
!
!    end if

    ! Select if derivatives shall be called or not.
    select case(jp)
    case(1)
        call emt(a_lat, celli, x_ref, n_l, n_p, particle_parms,lattice_parms, Eref)
        call emt(a_lat, celli, XDAT(I,:,:), n_l, n_p, particle_parms,lattice_parms, energy)
        energy=(energy-Eref)*diff
        F   = energy
        RES = YDAT(I) - F


    case(2)
        call emt_fit(a_lat, celli, x_ref, n_l, n_p, particle_parms, lattice_parms, Eref, dE_ref)
        call emt_fit(a_lat, celli, XDAT(I,:,:), n_l, n_p, particle_parms, lattice_parms, energy, denergy)

        energy=(energy-Eref)*diff
        denergy=(denergy-dE_ref)*diff
        !write(*,'(7f15.5)') denergy
        !stop
        F   = energy
        RES = YDAT(I) - F
        P(1:14)=denergy

        ! Set derivatives of those parameters zero that aren't supposed to contribute to fit.
        do ij=1,IP
            P(IB(ij)) = 0.0d0
        end do


        !--------WRITE ITERATION AND POINT TO SHOW STATUS ------------

        if ( ((mod(i,10)==0) .or. (i==N))) then
            !write(*,1000) iteration, i, YDAT(i), F, B(1:14)
            write(*,1000) iteration, i, YDAT(i), YDAT(i)-F, B(1:14)
        end if

        if (debug(4)) then
            write(7,1010) iteration,i, xdat(i,:,1),YDAT(i), F, RES, B(1:14)
        end if


    case(4)
        call emt(a_lat, celli, XDAT(I,:,:), n_l, n_p, particle_parms,lattice_parms, energy)
        F   = energy
        RRR(I) = F

    end select

       ! if ( ((mod(i,10)==0) .or. (i==N))) then
       !     write(*,1000) iteration, i, YDAT(i), F, B(1:14)
       ! end if

        !if (debug(4)) then
        !    write(7,1010) iteration,i, xdat(i,:,1),YDAT(i), F, RES, B(1:14)
        !end if

!    F   = energy
!    RES = YDAT(I) - F


    !write(*,1000) 'iter, i, jp, DFT EMT RES',iteration,i, JP, YDAT(i), F, RES

    ! filetered out above IF (JP.EQ.4) RRR(I) = F

    first_run=.false.
    RETURN

    1000 format(2i4, 2f12.4, 7f6.2 / 32x,7f6.2)
    1010 format(2i4,3f6.2,3E12.3,14f6.2)

end subroutine MODEL

