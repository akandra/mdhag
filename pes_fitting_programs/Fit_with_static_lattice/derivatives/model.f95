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
    integer             :: IP,IB
    integer             :: N, M, KK

    integer             :: iteration
    real(8)             :: energy,  E_dref
    real(8), dimension(14):: denergy
    real(8), dimension(7) ::dE_ref
    real(8)             :: Enew,dEnew(14)

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

! Select if derivatives shall be called or not.
    select case(jp)
    case(1)

        !call emt_init(a_lat,cell_0, n_l0, r0_lat, particle_parms,lattice_parms, Eref)

        r_l(I,:,:)=XDAT(I,:,2:n_l+1)
        r_p(I,:)=XDAT(I,:,1)
        call emt_init(a_lat, celli(1,:), n_l, XDAT(1,:,2:n_l+1), particle_parms,lattice_parms, Eref)
        call emt(a_lat, celli(I,:), r_p(I,:), r_l(I,:,:), n_l, particle_parms, lattice_parms, energy)
!        call emt_old(a_lat, celli(I,:), r_p(I,:), r_l(I,:,:), n_l, particle_parms, lattice_parms, energy)
! Get emt init implemented here and get it out of energy.
        energy= energy-Eref

        F   = energy
        RES = YDAT(I) - F



    case(2)
!        call emt_fit_init(a_lat, cell_0, r0_lat,n_l0, particle_parms, lattice_parms, E_dref, dE_ref)
        r_l(I,:,:)=XDAT(I,:,2:n_l+1)
        r_p(I,:)=XDAT(I,:,1)
        call emt_fit_init(a_lat, celli(1,:), XDAT(1,:,2:n_l+1), n_l, particle_parms, lattice_parms, E_dref, dE_ref)
        call emt_fit(a_lat, celli(I,:), r_p(I,:), r_l(I,:,:), n_l, particle_parms, lattice_parms, energy, denergy)
        energy=energy-E_dref
!        call emt_fit_old(a_lat, celli(I,:), r_p(I,:), r_l(I,:,:), n_l, particle_parms, lattice_parms, energy, denergy)

        F   = energy
        RES = YDAT(I) - F
        denergy(8:14) = denergy(8:14)-dE_ref
        P(1:14)=denergy


        do ij=1,IP
            P(IB(ij)) = 0.0d0
        end do
        !print *, 'der', P(1:7)

        !--------WRITE ITERATION AND POINT TO SHOW STATUS ------------

        if ( ((mod(i,10)==0) .or. (i==N))) then
            write(*,1000) iteration, i, YDAT(i), F, B(1:14)
        end if

        if (debug(4)) then
            write(7,1010) iteration,i, xdat(i,:,1),YDAT(i), F, RES, B(1:14)
        end if


    case(4)
!        call emt (a_lat, r_part, particle_parms, lattice_parms, energy)
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
!        call emt_fit (a_lat, r_part, particle_parms, lattice_parms, energy, denergy)
        F   = energy
        RES = YDAT(I) - F
    case(1)
!        call emt (a_lat, r_part, particle_parms, lattice_parms, energy)
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
