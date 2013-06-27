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
    real(8), dimension(:,:,:), allocatable, intent(out):: r_l      ! Position of lattice atoms
    real(8), dimension(:,:),allocatable, intent(out)    :: r_p        ! Position of particle
    real(8),allocatable                   :: E_dft(:)    ! DFT energy
    integer, intent(in)                 :: rep      ! Repetitions of lattice (1,2,3)
    integer, dimension(3), intent(in)   :: cell_in  ! contains cell geometry (2x2x4)
    integer                             :: run      ! which geometry of file, should be replaced by t

! other variables
    character(len=35)   :: position_of_l_and_p
    character(len=35)    :: energy_l_and_p
    integer             :: i, j,k, l, m,n,q, ios, start, ende, ende2,o
    !integer, dimension(:),allocatable :: o
    real(8)             :: c11,c12,c22, c33
    real(8)             :: temp
    character(len=35)   :: empty
    ! Arrays. Structure: timestep, x,y,z of atoms
    real(8), dimension(:,:,:), allocatable :: read_l ! array of au after read in
    real(8), dimension(:,:),allocatable    :: d_p     ! array of h
    real(8), dimension(:,:,:), allocatable :: d_l     ! array after multiplying lattice image
    real(8), dimension(:,:), allocatable :: celli, cell_min, cell_max

    position_of_l_and_p = 'traj005/XDATCAR_ACC_fsv_005.dat'
    energy_l_and_p = 'traj005/analyse_005.out'


! read in energies and timesteps
! This is an other file than that for the coordinates
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

! Note to self: Do everything with multidimensional array (3D)... will be easier to distinguish between things.
! read the transformation matrix in
    call open_for_read(18, position_of_l_and_p)
    read(18,'(A,/)') empty
    read(18,'(f13.10)') c11
    read(18,'(2f13.10)') c12, c22
    read(18,'(28X,f13.10,//)') c33

!    print*, c11,c12, c22, c33
! reads geometries in.

    ende = (time-2)*18
    allocate(read_l(time,3,ende))
    allocate(d_p(time-2,3))
!    allocate(read_l(1,3,16))
!    allocate(d_p(1,3))


    do j=1,10!time-2
        start = 1 + j*18
        ende2 = 16 + j*18
        do k =start, ende2
            read(18,*)  read_l(j,:,k)
!            print *, read_l(:,k)
        end do
        read(18,*) d_p(j,:)
!        read(8,*) empty
    end do
    close(18)
time =3
! translate lattice into bigger super cell
! 1. Allocate d_l-array:
! The array size is determined by the amounts of repetition of the cell image:
! the number of au-atoms in the cell has to be multiplied by the number of
! permutations and then, the number of gold atoms has to be added again, since
! one also wants to keep the original image in the new array.

    temp=cell_in(1)*cell_in(2)*cell_in(3)
    if (rep==1) then
        n = (temp*8 + temp)
    elseif (rep==2) then
        n = (temp*(8 + 16) + temp)
    else
        n = (temp*(8 + 16 + 24) + temp)
    end if

    allocate(d_l(time,3,n))
    allocate(r_l(time,3,n))
    allocate(celli(time,3))
    allocate(cell_max(time,3))
    allocate(cell_min(time,3))
    allocate(r_p(time,3))

    print *, read_l(1,3,:)

! loop so entire positions are converted
do q= 1, time-2
! Translation of lattice in x and y
    ! set the running parameters
    i=1!+(q-1)*temp
    k=temp!+(q-1)*temp
    l=1!+(q-1)*temp
    j=temp!+(q-1)*temp
    m = temp!+(q-1)*temp
!    print *, 'q', q

    ! Translation for rep = 1 (and also all the others)
    ! identity
    d_l(q,:,i:k)=read_l(q,:,l:j)
    print*, d_l

    ! in x
    i=i+m
    k=k+m
    d_l(q,1,i:k)=read_l(q,1,l:j)+1
    d_l(q,2,i:k)= read_l(q,2,l:j)
    d_l(q,3,i:k)= read_l(q,3,l:j)

    i=i+m
    k=k+m
    d_l(q,1,i:k)=read_l(q,1,l:j)-1
    d_l(q,2,i:k)= read_l(q,2,l:j)
    d_l(q,3,i:k)= read_l(q,3,l:j)


    ! in y
    i=i+m
    k=k+m
    d_l(q,:,i:k)=read_l(q,:,l:j)
    d_l(q,2,i:k)=read_l(q,2,l:j)+1

    i=i+m
    k=k+m
    d_l(q,:,i:k)=read_l(q,:,l:j)
    d_l(q,2,i:k)=read_l(q,2,l:j)-1


    ! in x and y
    i=i+m
    k=k+m
    d_l(q,:,i:k)=read_l(q,:,l:j)
    d_l(q,1,i:k)=read_l(q,1,l:j)+1
    d_l(q,2,i:k)=read_l(q,2,l:j)-1

    i=i+m
    k=k+m
    d_l(q,:,i:k)=read_l(q,:,l:j)
    d_l(q,1,i:k)=read_l(q,1,l:j)-1
    d_l(q,2,i:k)=read_l(q,2,l:j)+1

    i=i+m
    k=k+m
    d_l(q,:,i:k)=read_l(q,:,l:j)
    d_l(q,1,i:k)=read_l(q,1,l:j)+1
    d_l(q,2,i:k)=read_l(q,2,l:j)+1

    i=i+m
    k=k+m
    d_l(q,:,i:k)=read_l(q,:,l:j)
    d_l(q,1,i:k)=read_l(q,1,l:j)-1
    d_l(q,2,i:k)=read_l(q,2,l:j)-1


! Tranlation for rep=2
    if (rep .gt. 1) then
    ! in x
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)+2

        i=i+m
        k=k+m
            d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)-2

        ! in y
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,2,i:k)=read_l(q,2,l:j)+2

        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,2,i:k)=read_l(q,2,l:j)-2

        ! in x and y
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)+2
        d_l(q,2,i:k)=read_l(q,2,l:j)-2
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)+1
        d_l(q,2,i:k)=read_l(q,2,l:j)-2
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)+2
        d_l(q,2,i:k)=read_l(q,2,l:j)-1


        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)-2
        d_l(q,2,i:k)=read_l(q,2,l:j)+2
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)-1
        d_l(q,2,i:k)=read_l(q,2,l:j)+2
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)-2
        d_l(q,2,i:k)=read_l(q,2,l:j)+1

        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)+2
        d_l(q,2,i:k)=read_l(q,2,l:j)+2
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)+2
        d_l(q,2,i:k)=read_l(q,2,l:j)+1
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)+1
        d_l(q,2,i:k)=read_l(q,2,l:j)+2

        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)-2
        d_l(q,2,i:k)=read_l(q,2,l:j)-2
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)-2
        d_l(q,2,i:k)=read_l(q,2,l:j)-1
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)-1
        d_l(q,2,i:k)=read_l(q,2,l:j)-2
    end if

! Translation for rep=3
    if (rep .gt. 2) then
    !a =3
    ! in x
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)+3
        d_l(q,1,i:k)=read_l(q,1,l:j)+3

        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)-3

        ! in y
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,2,i:k)=read_l(q,2,l:j)+3

        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,2,i:k)=read_l(q,2,l:j)-3


        ! in x and y
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)+3
        d_l(q,2,i:k)=read_l(q,2,l:j)-3
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)+3
        d_l(q,2,i:k)=read_l(q,2,l:j)-2
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)+3
        d_l(q,2,i:k)=read_l(q,2,l:j)-1
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)+2
        d_l(q,2,i:k)=read_l(q,2,l:j)-3
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)+1
        d_l(q,2,i:k)=read_l(q,2,l:j)-3



        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)-3
        d_l(q,2,i:k)=read_l(q,2,l:j)+3
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)-3
        d_l(q,2,i:k)=read_l(q,2,l:j)+2
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)-3
        d_l(q,2,i:k)=read_l(q,2,l:j)+1
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)-1
        d_l(q,2,i:k)=read_l(q,2,l:j)+3
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)-2
        d_l(q,2,i:k)=read_l(q,2,l:j)+3



        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)+3
        d_l(q,2,i:k)=read_l(q,2,l:j)+3
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)+3
        d_l(q,2,i:k)=read_l(q,2,l:j)+2
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)+3
        d_l(q,2,i:k)=read_l(q,2,l:j)+1
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)+1
        d_l(q,2,i:k)=read_l(q,2,l:j)+3
        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)+2
        d_l(q,2,i:k)=read_l(q,2,l:j)+3



        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)-3
        d_l(q,2,i:k)=read_l(q,2,l:j)-3

        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)-3
        d_l(q,2,i:k)=read_l(q,2,l:j)-2

        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)-3
        d_l(q,2,i:k)=read_l(q,2,l:j)-1

        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)-1
        d_l(q,2,i:k)=read_l(q,2,l:j)-3

        i=i+m
        k=k+m
        d_l(q,:,i:k)=read_l(q,:,l:j)
        d_l(q,1,i:k)=read_l(q,1,l:j)-2
        d_l(q,2,i:k)=read_l(q,2,l:j)-3

    end if

    o=k              ! Saves the length of the translated array for the qth timestep


end do

! transform from direct into cartesian coordinates and determine extension of
! cell
    cell_max=0
    cell_min=0

    do q = 1, time
        do m= 1, o
            r_l(q,1,m) = c11*d_l(q,1,m)+c12*d_l(q,2,m)
            if ( r_l(q,1,m) .lt. cell_min(q,1))then
                cell_min(q,1)=r_l(q,1,m)
            else if ( r_l(q,1,m) .gt. cell_max(q,1))then
                cell_max(q,1)=r_l(q,1,m)
            end if

            r_l(q,2,m) = c22*d_l(q,2,m)
            if ( r_l(q,2,m) .lt. cell_min(q,2))then
                cell_min(q,2)=r_l(q,2,m)
            else if ( r_l(q,2,m) .gt. cell_max(q,2))then
                cell_max(q,2)=r_l(q,2,m)
            end if

            r_l(q,3,m) = c33*d_l(q,3,m)
            if ( r_l(q,3,m) .lt. cell_min(q,3))then
                cell_min(q,3)=r_l(q,3,m)
            else if ( r_l(q,3,m) .gt. cell_max(q,3))then
                cell_max(q,3)=r_l(q,3,m)
            end if
        end do

    celli = cell_max-cell_min
    r_p(1,q) = c11*d_p(1,q)+c12*d_p(2,q)
    r_p(2,q) = c22*d_p(2,q)
    r_p(3,q) = c33*d_p(3,q)

    end do

    print *, celli


    ! DON'T FORGET TO DEALLOCATE EVERYTHING!!!!!

end subroutine


end module emt_init_data
