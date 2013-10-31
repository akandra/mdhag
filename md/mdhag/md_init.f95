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
    integer :: nsteps
    real(8),dimension(3,6) :: celli
    type(species) :: spec_l, spec_p
    real(8), dimension(:), allocatable :: pars_l, pars_p


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
    character(len=7)        :: confname='POSCAR'
    integer                 :: fric_l, fric_p ! 0: no friction, 1: fixed coefficent
! other variables
    character(len=100)      :: buffer, label
    character(len=10):: name_p, name_l
    character(len=10):: pot_p, pot_l
    character(len=100) key_p, key_l
    real(8) :: mass_p, mass_l
    integer :: pos1, ios = 0, line = 0
    real(8) :: einc, inclination, azimuth, temp
    real(8), dimension(3,3) :: c_matrix, d_matrix
    integer :: i, j, k,l
    integer :: n_l0, itemp
    character(len=1) coord_sys
    real(8), dimension(:,:), allocatable :: start_l, start_p, d_l, pos_l, pos_p
    integer :: n_l, n_p=1
    integer :: npars_p, npars_l




!______________________________________________________________________________


    pos_init_file = 'au111_2x2x4_emt.in'

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
            case('Einc')
                read(buffer,*,iostat=ios) einc
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
            case ('lattice')
                read(buffer, *, iostat=ios) name_l, mass_l, pot_l, npars_l, key_l, fric_l
            case ('celldim')
                read(buffer, *, iostat=ios) celldim
            case ('rep')
                read(buffer, *, iostat=ios) rep
            case ('conf')
                read(buffer, *, iostat=ios) confname
                read(38,*) buffer
                read(38,*) temp
                read(38,*) c_matrix
                read(38,*) n_l0, n_p
                read(38,*) coord_sys

                allocate(start_l(3,n_l0))
                read(38,*) start_l
                allocate(start_p(3,n_p))
                read(38,*) start_p

           case default
!                print *, 'Skipping invalid label at line', line, label
            end select
        end if
    end do
    close(38)

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

    ! Transform the read in coordinates into direct if they are cartesians:
    if (coord_sys == 'C' .or. coord_sys == 'c') then
        start_l=matmul(d_matrix,start_l)
!        start_l=Nint(start_l*10000)/10000.0d0
    else
        start_p=matmul(c_matrix,start_p)
    end if

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

! allocate arrays
    allocate(d_l(3,n_l))
    d_l=0.d0

! Translation of the entire story


    i = 1
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



    allocate(pos_l(3,n_l))
    pos_l = matmul(c_matrix,d_l)
!    write(*,'(3f15.5)') pos_l
!    stop

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
    spec_l%mass = mass_l*amu2mass
    spec_l%n    = n_l
    spec_l%pot  = pot_l
    spec_l%n_pars= npars_l
    spec_l%fric = fric_l

    spec_p%name = name_p
    spec_p%mass = mass_p*amu2mass
    spec_p%n    = n_p
    spec_p%pot  = pot_p
    spec_p%n_pars= npars_p
    spec_p%fric = fric_p


    allocate(teilchen(n_p), slab(n_l))
    do i = 1, n_p
        teilchen(i)%r=start_p(:,i)
    end do
    einc = sqrt(2.0d0*einc/spec_p%mass) ! projectile speed
    teilchen(1)%v(1) = einc*sin(inclination)*cos(azimuth)
    teilchen(1)%v(2) = einc*sin(inclination)*sin(azimuth)
    teilchen(1)%v(3) = - einc*cos(inclination)

   do i = 1, n_l
        slab(i)%r=pos_l(:,i)
    end do

    deallocate(start_p)
    deallocate(pos_l,d_l,start_l)

end subroutine simbox_init


end module md_init

