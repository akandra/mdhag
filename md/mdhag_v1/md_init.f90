module md_init
!
! Purpose:
!          This module prepares the system for the md simulations
!           It should be able to do the following things:
!           1. Read in configuration
!           2. Construct simulations cell (that is, multiply 1. as often as necessary), depending
!                on the form of the read-in conf.
!           3. Assign velocities to atoms (i.e.: get temperature)

    use atom_class
    use open_file
    implicit none
    save

    real(8) :: a_lat
    real(8) :: step
    integer :: nsteps, wstep
    integer :: start_tr
    integer :: ntrajs
    real(8),dimension(3,6) :: celli
    type(species) :: spec_l, spec_p
    real(8), dimension(:), allocatable :: pars_l, pars_p
    real(8), dimension(:,:), allocatable :: x_ref
    character(len=7)        :: confname



contains



! COMMENT:  To improve the performance, it is better to redefine arrays with coordinates
!           in the way (coordinate, point) to escape the non-contiguous array problem
!           by passing a deferred-array to a subroutine


subroutine simbox_init(teilchen,slab)
!
! Purpose:
!           Initialises the entire system fully:
!               1. Geometry
!               2. Interaction Potential
!               3. Velocities
!
! Date          Author          History of Revision
! ====          ======          ===================
!01.10.2013     Sascha&Svenja   Implementation of 1.
!09.10.2013     Sascha&Svenja   Implementation of 2.
!
    implicit none

! Declare in and output

    type(atom), dimension(:), allocatable,intent(out) :: slab, teilchen

! control parameters default values
    character(len=35)       :: pos_init_file
    integer, dimension(3)   :: celldim=(/2,2,4/)  ! contains cell geometry (2x2x4)
    integer                 :: rep=2      ! Repetitions of lattice (1,2,3)
    integer                 :: fric_l, fric_p ! 0: no friction, 1: fixed coefficent
! other variables
    character(len=100)      :: buffer, label
    character(len=10):: name_p, name_l
    character(len=10):: pot_p, pot_l
    character(len=100) key_p, key_l
    real(8) :: mass_p, mass_l, Tsurf
    integer :: pos1, ios = 0, line = 0
    real(8) :: einc, inclination, azimuth, temp, E_pdof, v_pdof
    real(8), dimension(3,3) :: c_matrix, d_matrix
    integer :: i, j, k,l, s,r
    integer :: n_l0, itemp
    character(len=1) coord_sys
    real(8), dimension(:,:), allocatable :: start_l, start_p, d_l, pos_l, v_l, pos_p
    real(8), dimension(:,:), allocatable :: d_ref, d_p,  v_p
    integer :: n_l, n_p=1
    integer :: npars_p, npars_l
    logical :: exists
    character(len=100) :: confname_file
    integer             :: randk = 13
!______________________________________________________________________________


    call getarg(1, pos_init_file)
!    pos_init_file = 'mdhag.inp'

!------------------------------------------------------------------------------
!                       READ IN GEOMETRIES
!                       ===============================
!------------------------------------------------------------------------------

    call open_for_read(38,pos_init_file )

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.

    do while (ios == 0)
        read(38, '(A)', iostat=ios) buffer
        if (ios == 0) then
            line = line + 1

        ! Find the first instance of whitespace.  Split label and data.
            pos1 = scan(buffer, ' ')
            label = buffer(1:pos1)
            buffer = buffer(pos1+1:)

            select case (label)
            case('wstep')
                read(buffer,*,iostat=ios) wstep
            case('Tsurf')
                read(buffer,*,iostat=ios) Tsurf
            case('Einc')
                read(buffer,*,iostat=ios) einc
            case('start')
                read(buffer,*,iostat=ios) start_tr
            case('ntrajs')
                read(buffer,*,iostat=ios) ntrajs
            case('inclination')
                read(buffer,*,iostat=ios) inclination
                inclination = inclination*deg2rad
            case('azimuth')
                read(buffer,*,iostat=ios) azimuth
                azimuth = azimuth*deg2rad
            case('step')
                read(buffer,*,iostat=ios) step
            case('nsteps')
                read(buffer,*,iostat=ios) nsteps
            case ('particle')
                read(buffer, *, iostat=ios) name_p, mass_p, pot_p, npars_p, key_p, fric_p
                mass_p=mass_p*amu2mass
            case ('lattice')
                read(buffer, *, iostat=ios) name_l, mass_l, pot_l, npars_l, key_l, fric_l
                mass_l=mass_l*amu2mass
            case ('celldim')
                read(buffer, *, iostat=ios) celldim
            case ('rep')
                read(buffer, *, iostat=ios) rep
            case ('conf')
                read(buffer, *, iostat=ios) confname, confname_file

           case default
!                print *, 'Skipping invalid label at line', line, label
            end select
        end if
    end do
    close(38)

    call open_for_read(17, confname_file)

    if (confname == 'POSCAR') then

        read(17,*) buffer
        read(17,*) temp
        read(17,*) c_matrix
        read(17,*) n_l0, n_p
        read(17,*) coord_sys

        ! Make sure that zero-arguments in c_matrix are indeed zeros
        c_matrix=Nint(c_matrix*10000)/10000.0d0

        ! Construct the Inverse matrix of c_matrix
        ! CHANGE THIS TO INTEL-PROCEDURE
        d_matrix = 0.0d0
        d_matrix(1,1) = 1.0d0/c_matrix(1,1)
        d_matrix(2,2) = 1.0d0/c_matrix(2,2)
        d_matrix(3,3) = 1.0d0/c_matrix(3,3)
        d_matrix(1,2) = -d_matrix(2,2)*c_matrix(1,2)*d_matrix(1,1)

        ! Celli contains both c_ and d_matrix for pbc-procedure
        celli=0.0d0
        celli(1:2,1:2)=c_matrix(1:2,1:2)*(2*rep+1)
        celli(3,3)=c_matrix(3,3)

        celli(1:2,4:5)=d_matrix(1:2,1:2)/(2*rep+1)
        celli(3,6)=d_matrix(3,3)

        allocate(start_l(3,n_l0))
        read(17,*) start_l
        allocate(start_p(3,n_p))
        read(17,*) start_p

        ! Transform the read in coordinates into direct if they are cartesians:
        if (coord_sys == 'C' .or. coord_sys == 'c') then
            start_l=matmul(d_matrix,start_l)
            start_p=matmul(d_matrix,start_p)
        else
            !start_p=matmul(c_matrix,start_p)
            !start_l=matmul(c_matrix,start_l)
        end if

        ! Define geometry for reference energy
        allocate(d_ref(3,n_l0+1))
        d_ref(:,1) = (/0.0d0,0.0d0,0.29618952180d0/)
        d_ref(:,2:n_l0+1) = start_l

    !------------------------------------------------------------------------------
    !                    TRANSLATE LATTICE INTO BIGGER CELLS
    !                    ===================================
    !------------------------------------------------------------------------------
    ! 1. Allocate d_l-array:
    ! The array size is determined by the amounts of repetition of the cell image:
    ! the number of au-atoms in the cell has to be multiplied by the number of
    ! permutations and then, the number of gold atoms has to be added again, since
    ! one also wants to keep the original image in the new array.
    ! temp is the number of gold atoms
    ! if the cell is translated by rep =1, then 8 new images are formed
    !                              rep =2, then 8 + 16
    !                              rep =3, then 8 + 16 + 24
    ! We have decided to replicate quadratically around the original lattice,
    ! otherwise, rep = 1 would not make much sense.
    !
        itemp=celldim(1)*celldim(2)
        n_l=itemp*celldim(3)*(2*rep+1)**2
        n_p = (2*rep+1)**2

    ! allocate arrays
        allocate(d_l(3,n_l))
        allocate(d_p(3,n_p))
        allocate(x_ref(3,n_p+n_l))
        d_l=0.d0
        d_p = 0.0d0
        x_ref=0.0d0

    ! Translation of the entire story


        i = 1
        s = 1
        do l = 1, celldim(3)
            do j =-rep, rep
                do k=-rep, rep
                    d_l(1,i:i+itemp-1) = start_l(1,(l-1)*itemp+1:l*itemp)+j
                    d_l(2,i:i+itemp-1) = start_l(2,(l-1)*itemp+1:l*itemp)+k
                    d_l(3,i:i+itemp-1) = start_l(3,(l-1)*itemp+1:l*itemp)

                    i = i+itemp

                end do
            end do
        end do

        itemp=itemp*celldim(3)
        i=1
        j=1
        do r =-rep, rep
            do s=-rep, rep
                ! Set the running parameters
                x_ref(1,n_p+i:n_p+n_l) = d_ref(1,2:n_l0+1)+r
                x_ref(2,n_p+i:n_p+n_l) = d_ref(2,2:n_l0+1)+s
                x_ref(3,n_p+i:n_p+n_l) = d_ref(3,2:n_l0+1)
                i = i+itemp

                x_ref(1,j) = d_ref(1,1)+r
                x_ref(2,j) = d_ref(2,1)+s
                x_ref(3,j) = d_ref(3,1)

                d_p(1,j) = start_p(1,1)+r
                d_p(2,j) = start_p(2,1)+s
                d_p(3,j) = start_p(3,1)

                j=j+1

            end do
        end do

        allocate(pos_l(3,n_l),v_l(3,n_l))
        allocate(pos_p(3,n_p),v_p(3,n_p))

        pos_l = matmul(c_matrix,d_l)
        pos_p = matmul(c_matrix,d_p)
        x_ref = matmul(c_matrix,x_ref)

        !write(*,'(3f15.5)') x_ref1
        !stop

    ! Now, we implement the velocities of the Au atoms
    ! This is the doubled energy per atom
        E_pdof = Tsurf*kB
        v_pdof = sqrt(2.0d0*E_pdof/(mass_l))
        n_l0 = n_l/celldim(3)*(celldim(3)-1)
        v_l = 0.0d0
        call random_seed(size=randk)
        call random_seed(put=randseed)

        do i=1,n_l0
            v_l(1,i) = normal(0.0d0,v_pdof)
            v_l(2,i) = normal(0.0d0,v_pdof)
            v_l(3,i) = normal(0.0d0,v_pdof)
        enddo




        v_p = 0.0d0
        ! Set c.-of-m. velocity to zero
        v_l(1,1:n_l0) = v_l(1,1:n_l0) - sum(v_l(1,1:n_l0))/n_l0
        v_l(2,1:n_l0) = v_l(2,1:n_l0) - sum(v_l(2,1:n_l0))/n_l0
        v_l(3,1:n_l0) = v_l(3,1:n_l0) - sum(v_l(3,1:n_l0))/n_l0



    else

        read(17,*) celli
        read(17,*) n_l, n_p

        allocate(pos_l(3,n_l),v_l(3,n_l),start_p(3,n_p))
        read(17,*) start_p
        read(17,*) pos_l
        read(17,*) v_l

    endif

    close(17)

    allocate(pars_l(npars_l),pars_p(npars_p))

    call open_for_read(23,key_l)
    read(23,'(/)')
    do i = 1, npars_l
        read(23,*) buffer, pars_l(i)
    end do
    close(23)
!    print *, 'The parameters for the lattice are:'
!    print *, pars_l

    call open_for_read(23,key_p)
    read(23,'(/)')
    do i = 1, npars_p
        read(23,*) buffer, pars_p(i)
    end do
    close(23)
!    print *, pars_p

    ! Species definition
    spec_l%name = name_l
    spec_l%mass = mass_l
    spec_l%n    = n_l
    spec_l%pot  = pot_l
    spec_l%n_pars= npars_l
    spec_l%fric = fric_l

    spec_p%name = name_p
    spec_p%mass = mass_p
    spec_p%n    = n_p
    spec_p%pot  = pot_p
    spec_p%n_pars= npars_p
    spec_p%fric = fric_p


    allocate(teilchen(n_p), slab(n_l))
    do i = 1, n_p
        teilchen(i)%r=pos_p(:,i)
    end do
    if (confname .ne. 'POSCAR') then
        einc = sqrt(2.0d0*einc/spec_p%mass) ! projectile speed
        teilchen(:)%v(1) = einc*sin(inclination)*cos(azimuth)
        teilchen(:)%v(2) = einc*sin(inclination)*sin(azimuth)
        teilchen(:)%v(3) = - einc*cos(inclination)
    else
        teilchen(1)%v = 0.0d0
    end if

    do i = 1, n_l
        slab(i)%r=pos_l(:,i)
        slab(i)%v=v_l(:,i)
    end do

! Create a directory for trajectory data
    inquire(file='trajs',exist=exists)
    if (.not. exists) then
        call system('mkdir trajs')
    end if


    deallocate(start_p)
    deallocate(pos_l,d_l,pos_p, d_p, d_ref,start_l,v_l)

end subroutine simbox_init

function ran1()  !returns random number between 0 - 1
 implicit none
        real(8) ran1,x
        call random_number(x) ! built in fortran 90 random number function
        ran1=x
end function ran1

function normal(mean,sigma) !returns a normal distribution
 implicit none
        real(8) normal,tmp
        real(8) mean,sigma   ! Sigma is the velocity we want to achieve
        integer flag
        real(8) fac,gsave,rsq,r1,r2
        save flag,gsave
        data flag /0/
        if (flag.eq.0) then
        rsq=2.0d0
            do while(rsq.ge.1.0d0.or.rsq.eq.0.0d0) ! new from for do
                r1=2.0d0*ran1()-1.0d0
                r2=2.0d0*ran1()-1.0d0
                rsq=r1*r1+r2*r2
            enddo
            fac=sqrt(-2.0d0*log(rsq)/rsq)
            gsave=r1*fac        ! shouldn't those two values be below zero?
            tmp=r2*fac          !
            flag=1
        else
            tmp=gsave
            flag=0
        endif
        normal=tmp*sigma+mean
        return
end function normal

end module md_init

