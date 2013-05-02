subroutine model( F, YDAT, XDAT, RRR, I, JP )
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


    COMMON /BLK1/ B(20),P(20),RES,N,M,KK
    COMMON/ DJA1/ iteration
    COMMON/debug/ debug(5)

    if ((debug(3))) then
        debug(3)=.false.
        print *
        print '((a))', 'FIRST RUN OF MODEL'
        print '((a),4i5)', '  i jp n m=', i, jp, n, m
        print '((a),4f10.5)', '  xi,yi=   ', xdat(i,1), xdat(i,2), xdat(i,3), ydat(i)
        print '((a),3I10)', '  loc(xi))=', loc(xdat(i,1)), loc(xdat(i,2)), loc(xdat(i,3))
        print *
        pause 'first run of model';
    end if

    call array2emt_parms( B(1:7 ), particle_parms)
    call array2emt_parms( B(8:14), lattice_parms )
    r_part=XDAT(i,:)

    call emt (r_part, particle_parms, lattice_parms, energy)

    F   = energy
    RES = YDAT(I) - F
    if(jp==4) RRR(i)=RES

    !--------WRITE ITERATION AND POINT TO SHOW STATUS ------------
    if ( (mod(i,10)==0) .or. (i==N)) then
        write(*,1000) 'i, jp, DFT EMT RES',iteration,i, JP, YDAT(i), F, RES
        write(7,1000) 'i, jp, DFT EMT RES',iteration,i, JP, YDAT(i), F, RES
    end if

    !write(*,1000) 'iter, i, jp, DFT EMT RES',iteration,i, JP, YDAT(i), F, RES

    ! filetered out above IF (JP.EQ.4) RRR(I) = F

    if(jp .eq. 2) jp=3
    if(jp .eq. 4) rrr(i)=f

    RETURN

    1000 format((a) 3i5, 3E12.3)

end subroutine MODEL
