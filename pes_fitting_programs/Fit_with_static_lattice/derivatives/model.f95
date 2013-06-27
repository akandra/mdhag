subroutine model( F, YDAT, XDAT, RRR, I, JP)
    use EMTparms_class

    implicit none

    logical debug

    type(EMTparms)      :: particle_parms   ! parameters of particle
    type(EMTparms)      :: lattice_parms    ! parameters of lattice atoms
    real(8)             :: F
    real(8)             :: YDAT(1000)   !YDAT(N)
    real(8)             :: XDAT(1000,3) !XDAT(N,M)
    real(8)             :: RRR(1000)
    integer             :: I, JP,ij

    real(8)             :: B, P, RES
    integer             :: IP,IB
    integer             :: N, M, KK

    integer             :: iteration
    real(8)             :: energy
    real(8), dimension(14):: denergy
    real(8)             :: r_part(3)

    logical, save       :: first_run=.true.


    COMMON /BLK1/ B(20),P(20),RES,N,M,KK
    COMMON/ DJA1/ iteration
    COMMON/BLK5/IB(20),IP
    COMMON/debug/ debug(5)

    if (debug(3) .and. first_run) then
        print *
        print '((a))', 'FIRST RUN OF MODEL'
        print '((a),4i5)', '  i jp n m=', i, jp, n, m
        print '((a),4f20.5)', '  xi,yi=   ', xdat(i,1), xdat(i,2), xdat(i,3), ydat(i)
        print '((a),3I20)',   '  loc(xi))=', loc(xdat(i,1)), loc(xdat(i,2)), loc(xdat(i,3))
        print '((a),3I20)',   '  loc(x2))=', loc(xdat(2,1)), loc(xdat(2,2)), loc(xdat(2,3))
        print *
        write(7,*) 'it i jp    x       y      z       DFT     EMT    RES'
        !pause 'first run of model';
    end if

    call array2emt_parms( B(1:7 ), particle_parms)
    call array2emt_parms( B(8:14), lattice_parms )
    r_part=XDAT(i,:)

    select case(jp)
    case(1)
        call emt (a_lat, r_part, particle_parms, lattice_parms, energy)
        F   = energy
        RES = YDAT(I) - F

    case(2)
        call emt_fit (a_lat, r_part, particle_parms, lattice_parms, energy, denergy)
! don't forget to hold parameters fast
        F   = energy
        RES = YDAT(I) - F
        P(1:14) = denergy


        do ij=1,IP
            P(IB(ij)) = 0.0
        end do

        !--------WRITE ITERATION AND POINT TO SHOW STATUS ------------

        if ( ((mod(i,10)==0) .or. (i==N))) then
            write(*,1000) iteration, i, YDAT(i), F, B(1:14)
        end if

        if (debug(4)) then
            write(7,1010) iteration,i, xdat(i,:),YDAT(i), F, RES, B(1:14)
        end if


    case(4)
        call emt (a_lat, r_part, particle_parms, lattice_parms, energy)
        F   = energy
        RRR(I) = F

    end select

    F   = energy
    RES = YDAT(I) - F


    !write(*,1000) 'iter, i, jp, DFT EMT RES',iteration,i, JP, YDAT(i), F, RES

    ! filetered out above IF (JP.EQ.4) RRR(I) = F

    first_run=.false.
    RETURN

    1000 format(2i4, 2f6.2, 7f6.2 / 20x,7f6.2)
    1010 format(2i4,3f6.2,3E12.3,14f6.2)

end subroutine MODEL

subroutine MODEL2606old( F, YDAT, XDAT, RRR, I, JP, denergy )
    use EMTparms_class

    implicit none

    logical debug

    type(EMTparms)      :: particle_parms   ! parameters of particle
    type(EMTparms)      :: lattice_parms    ! parameters of lattice atoms
    real(8)             :: F
    real(8)             :: YDAT(1000)   !YDAT(N)
    real(8)             :: XDAT(1000,3) !XDAT(N,M)
    real(8)             :: RRR(1000)
    integer             :: I, JP

    real(8)             :: B, P, RES
    integer             :: N, M, KK

    integer             :: iteration
    real(8)             :: energy
    real(8)             :: r_part(3)
    real(8), dimension(14) :: denergy


    logical, save       :: first_run=.true.


    COMMON /BLK1/ B(20),P(20),RES,N,M,KK
    COMMON/ DJA1/ iteration
    COMMON/debug/ debug(5)

    if (debug(3) .and. first_run) then
        print *
        print '((a))', 'FIRST RUN OF MODEL'
        print '((a),4i5)', '  i jp n m=', i, jp, n, m
        print '((a),4f20.5)', '  xi,yi=   ', xdat(i,1), xdat(i,2), xdat(i,3), ydat(i)
        print '((a),3I20)',   '  loc(xi))=', loc(xdat(i,1)), loc(xdat(i,2)), loc(xdat(i,3))
        print '((a),3I20)',   '  loc(x2))=', loc(xdat(2,1)), loc(xdat(2,2)), loc(xdat(2,3))
        print *
        write(7,*) 'it i jp    x       y      z       DFT     EMT    RES'
        !pause 'first run of model';
    end if

    call array2emt_parms( B(1:7 ), particle_parms)
    call array2emt_parms( B(8:14), lattice_parms )
    r_part=XDAT(i,:)

    select case(jp)
    case(2)
        call emt_fit (a_lat, r_part, particle_parms, lattice_parms, energy, denergy)
    case(1)
        call emt (a_lat, r_part, particle_parms, lattice_parms, energy)
    end select

    F   = energy
    RES = YDAT(I) - F


    !--------WRITE ITERATION AND POINT TO SHOW STATUS ------------
    if ( jp.eq.2 .and. ((mod(i,10)==0) .or. (i==N))) then
        write(*,1000) iteration, i, YDAT(i), F, B(1:14)
    end if

    if (debug(4).and.JP==2) then
        write(7,1010) iteration,i, xdat(i,:),YDAT(i), F, RES, B(1:14)
    end if
    !write(*,1000) 'iter, i, jp, DFT EMT RES',iteration,i, JP, YDAT(i), F, RES

    ! filetered out above IF (JP.EQ.4) RRR(I) = F

    if(jp .eq. 2) jp=3
    if(jp .eq. 4) rrr(i)=F
    first_run=.false.
    RETURN

    1000 format(2i4, 2f6.2, 7f6.2 / 20x,7f6.2)
    1010 format(2i4,3f6.2,3E12.3,14f6.2)

end subroutine MODEL2606old
