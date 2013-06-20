module emt_init_data
    implicit none
    save

!    integer                                 :: n_lat0_at         ! number of atoms in reference slab
!    integer                                 :: n_lay0            ! number of layers in reference slab
!    real(8)                                 :: nn0               ! next neighbour distance in reference slab
!    real(8), dimension(3)                   :: cell              ! dimensions of the cell
!    real(8), allocatable, dimension(:,:)    :: r0_lat            ! lattice positions for reference calc.

contains

subroutine l_p_position(time, r_l,r_p)
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
    real(8), dimension(3,16), intent(out):: r_l      ! Position of lattice atoms
    real(8), dimension(3), intent(out)    :: r_p        ! Position of particle
    real(8),allocatable                   :: E_dft(:)    ! DFT energy
    integer                             :: rep      ! Repetitions of lattice
    integer                             :: run      ! which geometry of file, should be replaced by t

! other variables
    character(len=35)   :: position_of_l_and_p
    character(len=35)    :: energy_l_and_p
    integer             :: i, j, ios, k, start, ende, ende2
    real(8)             :: c11,c12,c22, c33
    real(8)             :: temp
    character(len=35)   :: empty
    real(8), dimension(:,:), allocatable :: d_l
    real(8), dimension(:),allocatable    :: d_p
    real(8), allocatable, dimension(:,:) :: rep_l   ! array after replication into all directions

    position_of_l_and_p = 'traj005/XDATCAR_ACC_fsv_005.dat'
    energy_l_and_p = 'traj005/analyse_005.out'


! read in energies and timesteps
    call open_for_read(7,energy_l_and_p)
    i = 0
    do
        read(7,*,iostat=ios)
        if(ios <0) exit
        i=i+1
    end do
    rewind(7)
    time = i
    print*, time

    allocate(E_dft(time))
    read(7,'(A)') empty
    do i=2,time,1
        read(7,*) temp, temp, E_dft(i)
!        read(7,'(3f15.14)') E_dft(i), E_dft(i), E_dft(i)
    end do
    E_dft=E_dft+25.024789d0

! read the transformation matrix in
    call open_for_read(8, position_of_l_and_p)
    read(8,'(A,/)') empty
    read(8,'(f13.10)') c11
    read(8,'(2f13.10)') c12, c22
    read(8,'(28X,f13.10,//)') c33
!    print*, c12, c11, c22, c33
! reads geometries in.
    ende = (time-2)*18
    allocate(d_l(3,ende))
    allocate(d_p(time-2))
    do j=1,time-2
        start = 1 + j*18
        ende2 = 16 + j*18
        do k =start, ende2
            read(8,*)  d_l(:,k)
!            print *, d_l(:,k)
        end do
        read(8,*) d_p(j)
!        read(8,*) empty
    end do
    print*, j,ende, d_p
!    read(8,*) d_p

    print *, 'temp', temp

stop



! transform from direct into cartesian coordinates

    do k= 1, 16
        r_l(1,k) = c11*d_l(1,k)+c12*d_l(2,k)
        r_l(2,k) = c22*d_l(2,k)
        r_l(3,k) = c22*d_l(3,k)

    end do
    r_p(1) = c11*d_p(1)+c12*d_p(2)
    r_p(2) = c22*d_p(2)
    r_p(3) = c22*d_p(3)
!    print *, r_l

end subroutine


end module emt_init_data
