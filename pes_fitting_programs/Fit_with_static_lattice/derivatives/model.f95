subroutine model( F, YDAT, XDAT, RRR, I, JP, PX )
!subroutine model( F, YDAT, XDAT, RRR, I, JP)
    use EMTparms_class
    use emt_init_data

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

!    real(8)             :: a_lat       ! doesn't need to be declared because common
    real(8), dimension(7) :: denergy_l, denergy_p
    real(8)             :: PX(20)
    real(8)             :: E_ref

    logical, save       :: first_run=.true.


    COMMON /BLK1/ B(20),P(20),RES,N,M,KK
    COMMON/ DJA1/ iteration
    COMMON/debug/ debug(5)

    ! Just put the lattice constant here... perhaps we will find a better place to put it
    ! later.


    if ((debug(3)) .and. first_run) then
        print *
        print '((a))', 'FIRST RUN OF MODEL'
        print '((a),4i5)', '  i jp n m=', i, jp, n, m
        print '((a),4f10.5)', '  xi,yi=   ', xdat(i,1), xdat(i,2), xdat(i,3), ydat(i)
        print '((a),3I10)',   '  loc(xi))=', loc(xdat(i,1)), loc(xdat(i,2)), loc(xdat(i,3))
        print '((a),3I10)',   '  loc(x2))=', loc(xdat(2,1)), loc(xdat(2,2)), loc(xdat(2,3))
        print *
        write(7,*) 'it i jp    x       y      z       DFT     EMT    RES'
        !pause 'first run of model';
    end if

    call array2emt_parms( B(1:7 ), particle_parms)
    call array2emt_parms( B(8:14), lattice_parms )
    r_part=XDAT(i,:)
!    print *, lattice_parms
!    print *, particle_parms

    ! Cycle which decides via JP if only energy is calculated or also derivatives
    if ( JP == 1 ) then
        call emt (a_lat, r_part, particle_parms, lattice_parms, energy)
    else if ( JP == 2 ) then
        call emt_fit(a_lat, r_part, particle_parms, lattice_parms, energy, denergy_l, denergy_p)
        PX(1:7) = denergy_p
        PX(8:14) = denergy_l
    end if

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
    if(jp .eq. 4) rrr(i)=f
    first_run=.false.
    RETURN

    1000 format(2i4, 2f6.2, 7f6.2 / 20x,7f6.2)
    1010 format(2i4,3f6.2,3E12.3,14f6.2)

end subroutine MODEL
