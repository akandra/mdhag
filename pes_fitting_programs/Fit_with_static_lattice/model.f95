subroutine model( F, YDAT, XDAT, RRR, I, JP )
    use EMTparms_class

    implicit none

    type(EMTparms)      :: particle_parms   ! parameters of particle
    type(EMTparms)      :: lattice_parms    ! parameters of lattice atoms
    real(8)             :: F, YDAT(1), XDAT(3,1), RRR(1)
    integer             :: I, JP

    real(8)             :: B, P, RES
    integer             :: N, M, KK

    real(8)             :: energy
    real(8)             :: b_old(20)

    COMMON /BLK1/ B(20),P(20),RES,N,M,KK

    if(jp .eq. 2) jp=3
    if(jp .gt. 2) return

    call array2emt_parms( B(1:7 ), particle_parms)
    call array2emt_parms( B(8:14), lattice_parms )

    !write(*,1000) 'particle_pars=',particle_parms



    if (mod(i,10)==0) write(*,1000) i,jp,B(8:14)
                      write(7,1000) i,jp,B(8:14)
    b_old=B

    call emt (XDAT(:,i), particle_parms, lattice_parms, energy)

    F   = energy
    RES = YDAT(I) - F
    if (YDAT(I) > 5) RES=0       ! ignore points that have energy > 5 eV


    IF (JP.EQ.4) RRR(I) = F


    RETURN

1000 format(i3,i2,2x,7F10.6/7x7F10.7)

1010 format(a15,4x,7F7.3,/,19x,7f7.3)
1020 format((a), 4F7.3)
end subroutine MODEL
