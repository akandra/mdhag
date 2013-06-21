module emt_init_data
    use EMTparms_class
    implicit none
    save

!    integer                                 :: n_lat0_at         ! number of atoms in reference slab
!    integer                                 :: n_lay0            ! number of layers in reference slab
!    real(8)                                 :: nn0               ! next neighbour distance in reference slab
!    real(8), dimension(3)                   :: cell              ! dimensions of the cell
!    real(8), allocatable, dimension(:,:)    :: r0_lat            ! lattice positions for reference calc.

contains

subroutine l_p_position(time, r_l,r_p, cell_in, rep)
!
! Purpose:
!           Reads in the gold and hydrogen positions from AIMD
!
! Date          Author          History of Revision
! ====          ======          ===================
! 20.06.2013    S. M. Janke     original
!
    implicit none

! Declare in and output
    integer, intent(out)                  :: time        ! timestep of AIMD
    real(8), dimension(:,:), allocatable, intent(out):: r_l      ! Position of lattice atoms
    real(8), dimension(3), intent(out)    :: r_p        ! Position of particle
    real(8),allocatable                   :: E_dft(:)    ! DFT energy
    integer, intent(in)                 :: rep      ! Repetitions of lattice (1,2,3)
    integer, dimension(3), intent(in)   :: cell_in  ! contains cell geometry (2x2x4)
    integer                             :: run      ! which geometry of file, should be replaced by t

! other variables
    character(len=35)   :: position_of_l_and_p
    character(len=35)    :: energy_l_and_p
    integer             :: i, j,k, l, m,n, ios, start, ende, ende2
    real(8)             :: c11,c12,c22, c33
    real(8)             :: temp
    character(len=35)   :: empty
    real(8), dimension(:,:), allocatable :: read_l ! array of au after read in
    real(8), dimension(:),allocatable    :: d_p     ! array of h
    real(8), dimension(:,:), allocatable :: d_l     ! array after multiplying lattice image
    real(8), dimension(3)                :: cell, cell_min, cell_max

    position_of_l_and_p = 'traj005/XDATCAR_ACC_fsv_005.dat'
    energy_l_and_p = 'traj005/analyse_005.out'


! read in energies and timesteps
    call open_for_read(17,energy_l_and_p)
    i = 0
    do
        read(17,*,iostat=ios)
        if(ios <0) exit
        i=i+1
    end do
    rewind(17)
    time = i

    allocate(E_dft(time))
    read(17,'(A)') empty
    do i=2,time,1
        read(17,*) temp, temp, E_dft(i)
!        read(7,'(3f15.14)') E_dft(i), E_dft(i), E_dft(i)
    end do
    close(17)
    E_dft=E_dft+25.024789d0

! read the transformation matrix in
    call open_for_read(18, position_of_l_and_p)
    read(18,'(A,/)') empty
    read(18,'(f13.10)') c11
    read(18,'(2f13.10)') c12, c22
    read(18,'(28X,f13.10,//)') c33

!    print*, c11,c12, c22, c33
! reads geometries in.
    ende = (time-2)*18
    allocate(read_l(3,ende))
    allocate(d_p(time-2))
!    allocate(read_l(3,16))
!    allocate(d_p(1))

    do j=1,10!time-2
        start = 1 + j*18
        ende2 = 16 + j*18
        do k =start, ende2
            read(18,*)  read_l(:,k)
            print *, read_l(:,k)
        end do
        read(18,*) d_p(j)
!        read(8,*) empty
    end do
    close(18)


stop
! translate lattice into bigger super cell
! 1. Allocate d_l-array:
! The array size is determined by the amounts of repetition of the cell image:
! the number of au-atoms in the cell has to be multiplied by the number of
! permutations and then, the number of gold atoms has to be added again, since
! one also wants to keep the original image in the new array.
    temp=cell_in(1)*cell_in(2)*cell_in(3)
    if (rep==1) then
        n = temp*8 + temp
    elseif (rep==2) then
        n = temp*(8 + 16) + temp
    else
        n = temp*(8 + 16 + 24) + temp
    end if

    allocate(d_l(3,n))

! Translation of lattice in x and y
    ! set the running parameters
    i=1
    k=temp
    l=1
    j=temp
    m = temp


    ! Translation for rep = 1 (and also all the others)
    ! identity
    d_l(:,i:k)=read_l(:,l:j)

    ! in x
    i=i+m
    k=k+m
    d_l(1,i:k)=read_l(1,l:j)+1
    d_l(2,i:k)= read_l(2,l:j)
    d_l(3,i:k)= read_l(3,l:j)

    i=i+m
    k=k+m
    d_l(1,i:k)=read_l(1,l:j)-1
    d_l(2,i:k)= read_l(2,l:j)
    d_l(3,i:k)= read_l(3,l:j)


    ! in y
    i=i+m
    k=k+m
    d_l(:,i:k)=read_l(:,l:j)
    d_l(2,i:k)=read_l(2,l:j)+1

    i=i+m
    k=k+m
    d_l(:,i:k)=read_l(:,l:j)
    d_l(2,i:k)=read_l(2,l:j)-1


    ! in x and y
    i=i+m
    k=k+m
    d_l(:,i:k)=read_l(:,l:j)
    d_l(1,i:k)=read_l(1,l:j)+1
    d_l(2,i:k)=read_l(2,l:j)-1

    i=i+m
    k=k+m
    d_l(:,i:k)=read_l(:,l:j)
    d_l(1,i:k)=read_l(1,l:j)-1
    d_l(2,i:k)=read_l(2,l:j)+1

    i=i+m
    k=k+m
    d_l(:,i:k)=read_l(:,l:j)
    d_l(1,i:k)=read_l(1,l:j)+1
    d_l(2,i:k)=read_l(2,l:j)+1

    i=i+m
    k=k+m
    d_l(:,i:k)=read_l(:,l:j)
    d_l(1,i:k)=read_l(1,l:j)-1
    d_l(2,i:k)=read_l(2,l:j)-1


! Tranlation for rep=2
    if (rep .gt. 1) then
    ! in x
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)+2

        i=i+m
        k=k+m
            d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)-2

        ! in y
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(2,i:k)=read_l(2,l:j)+2

        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(2,i:k)=read_l(2,l:j)-2

        ! in x and y
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)+2
        d_l(2,i:k)=read_l(2,l:j)-2
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)+1
        d_l(2,i:k)=read_l(2,l:j)-2
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)+2
        d_l(2,i:k)=read_l(2,l:j)-1


        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)-2
        d_l(2,i:k)=read_l(2,l:j)+2
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)-1
        d_l(2,i:k)=read_l(2,l:j)+2
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)-2
        d_l(2,i:k)=read_l(2,l:j)+1

        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)+2
        d_l(2,i:k)=read_l(2,l:j)+2
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)+2
        d_l(2,i:k)=read_l(2,l:j)+1
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)+1
        d_l(2,i:k)=read_l(2,l:j)+2

        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)-2
        d_l(2,i:k)=read_l(2,l:j)-2
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)-2
        d_l(2,i:k)=read_l(2,l:j)-1
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)-1
        d_l(2,i:k)=read_l(2,l:j)-2
    end if

    if (rep .gt. 2) then
    !a =3
    ! in x
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)+3
        d_l(1,i:k)=read_l(1,l:j)+3

        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)-3

        ! in y
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(2,i:k)=read_l(2,l:j)+3

        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(2,i:k)=read_l(2,l:j)-3


        ! in x and y
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)+3
        d_l(2,i:k)=read_l(2,l:j)-3
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)+3
        d_l(2,i:k)=read_l(2,l:j)-2
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)+3
        d_l(2,i:k)=read_l(2,l:j)-1
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)+2
        d_l(2,i:k)=read_l(2,l:j)-3
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)+1
        d_l(2,i:k)=read_l(2,l:j)-3



        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)-3
        d_l(2,i:k)=read_l(2,l:j)+3
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)-3
        d_l(2,i:k)=read_l(2,l:j)+2
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)-3
        d_l(2,i:k)=read_l(2,l:j)+1
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)-1
        d_l(2,i:k)=read_l(2,l:j)+3
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)-2
        d_l(2,i:k)=read_l(2,l:j)+3



        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)+3
        d_l(2,i:k)=read_l(2,l:j)+3
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)+3
        d_l(2,i:k)=read_l(2,l:j)+2
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)+3
        d_l(2,i:k)=read_l(2,l:j)+1
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)+1
        d_l(2,i:k)=read_l(2,l:j)+3
        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)+2
        d_l(2,i:k)=read_l(2,l:j)+3



        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)-3
        d_l(2,i:k)=read_l(2,l:j)-3

        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)-3
        d_l(2,i:k)=read_l(2,l:j)-2

        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)-3
        d_l(2,i:k)=read_l(2,l:j)-1

        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)-1
        d_l(2,i:k)=read_l(2,l:j)-3

        i=i+m
        k=k+m
        d_l(:,i:k)=read_l(:,l:j)
        d_l(1,i:k)=read_l(1,l:j)-2
        d_l(2,i:k)=read_l(2,l:j)-3

    end if


    do m= 1,k
 !       print *, d_l(:,m)
    end do

    ! allocate r_l
    allocate(r_l(3,k))

! transform from direct into cartesian coordinates and determine extension of
! cell
    cell_max=0
    cell_min=0

    do m= 1, k
        r_l(1,m) = c11*d_l(1,m)+c12*d_l(2,m)
        if ( r_l(1,m) .lt. cell_min(1))then
            cell_min(1)=r_l(1,m)
        else if ( r_l(1,m) .gt. cell_max(1))then
            cell_max(1)=r_l(1,m)
        end if

        r_l(2,m) = c22*d_l(2,m)
        if ( r_l(2,m) .lt. cell_min(2))then
            cell_min(2)=r_l(2,m)
        else if ( r_l(2,m) .gt. cell_max(2))then
            cell_max(2)=r_l(2,m)
        end if

        r_l(3,m) = c33*d_l(3,m)
        if ( r_l(3,m) .lt. cell_min(3))then
            cell_min(3)=r_l(3,m)
        else if ( r_l(3,m) .gt. cell_max(3))then
            cell_max(3)=r_l(3,m)
        end if


    end do
    cell = cell_max-cell_min




    r_p(1) = c11*d_p(1)+c12*d_p(2)
    r_p(2) = c22*d_p(2)
    r_p(3) = c33*d_p(3)

    do m= 1,k
!        print *, r_l(:,m)
    end do
    print *, cell

end subroutine


end module emt_init_data
