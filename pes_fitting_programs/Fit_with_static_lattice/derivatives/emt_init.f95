module emt_init_data
    use EMTparms_class
    implicit none
    save

    integer :: time
    integer :: n_l
    real(8), allocatable, dimension(:,:) :: celli

contains

! COMMENT:  To improve the performance, it is better to redefine arrays with coordinates in the way (coordinate, point)
!           to escape the non-contiguous array problem by passing a deferred-array to a subroutine


subroutine l_p_position(a_lat, rep, cell_in, control,e_aimd_max, time, l_aimd, n_l, celli, x_all, E_all)
!
! Purpose:
!           Reads in the gold and hydrogen positions from AIMD
!
! Date          Author          History of Revision
! ====          ======          ===================
! 01.07.2013    S. M. Janke     original
!
    implicit none

! Declare in and output
    real(8), intent(in)                 :: a_lat    ! lattice constant
    integer, intent(in)                 :: rep      ! Repetitions of lattice (1,2,3)
    integer, dimension(3), intent(in)   :: cell_in  ! contains cell geometry (2x2x4)
    integer, intent(out)                  :: time        ! different configuration and energies
    real(8),allocatable,dimension(:), intent(out)                 :: E_all    ! DFT energy
    real(8),allocatable,dimension(:,:), intent(out)               :: celli
    integer, intent(out)                :: n_l, l_aimd  ! number of l atoms, number of aimd contributions
    real(8), dimension(:,:,:), allocatable, intent(out) :: x_all
    integer , intent(in)                        :: control
    real(8), intent(in) :: e_aimd_max

! other variables
    character(len=35)   :: position_of_l_and_p, fix_position
    character(len=35)    :: energy_l_and_p, fix_energy
    integer             :: i, j,k,u,l, m,n,q, ios, start, ende, ende2,o,npts
    real(8)             :: c11,c12,c22, c33
    real(8)             :: temp, e_max
    real(8)             :: isqrt_2
    character(len=35)   :: empty
    ! Arrays. Structure: timestep, x,y,z of atoms
    real(8), dimension(:,:,:), allocatable :: aimd_l, aimd_cor_l, prae_read_l, read_l ! array of au after read in
    real(8), dimension(:,:),allocatable    :: aimd_p, aimd_cor_p, prae_d_p, d_p     ! array of h
    real(8), dimension(:,:,:), allocatable :: d_l     ! array after multiplying lattice image
    real(8), dimension(:,:,:), allocatable :: fix_l ! array of au after read in
    real(8), dimension(:,:),allocatable    :: fix_p     ! array of h
    real(8), dimension(:,:,:), allocatable:: r_l      ! Position of lattice atoms
    real(8), dimension(:,:),allocatable   :: r_p        ! Position of particle

    real(8), dimension(:,:), allocatable :: cell_min, cell_max
    real(8), dimension(483)             :: E_fix    ! energy of fixed lattice geometires
    integer, dimension(:), allocatable  :: chaos
    real(8),allocatable,dimension(:)                 :: E_dft1, E_dft, prae_E_dft    ! read-in-dft-energy


    position_of_l_and_p = 'data/traj833/XDATCAR_833.dat'
    energy_l_and_p =      'data/traj833/analyse_833.out'

    fix_position = 'data/DFT_aimd_form2x2x4.dat'
    fix_energy = 'data/hau111_plot.E.dat'

!------------------------------------------------------------------------------
!                       READ IN GEOMETRIES AND ENERGIES
!                       ===============================
!------------------------------------------------------------------------------

! -------------------------READ IN FIXED LATTICE-------------------------------
! read in energies
    call open_for_read(39,fix_energy)
        i=1
    do
        read(39,*,iostat=ios)
        if(ios <0) exit
        i=i+1
    end do
    rewind(39)
    npts = i-1

!    npts=10
    e_max=10.0
    E_fix=0.0

    j=1
    do i=1, npts
        read(39,*) empty, empty, empty, empty, E_fix(j)
        if (abs(E_fix(j))<=e_max) j=j+1
    end do
    close(39)
    npts=j-1

! Read in the geometry of the fixed lattice
    allocate(fix_l(npts,3,16))
    allocate(fix_p(npts,3))
    call open_for_read(38,fix_position)
    fix_p = 0.0
    fix_l = 0.0

    do j = 1, npts
        do k=1,16
            read(38,*)  fix_l(j,1,k), fix_l(j,2,k), fix_l(j,3,k)
!            print *, fix_l(:,k)
        end do
        read(38,*) fix_p(j,1), fix_p(j,2), fix_p(j,3)
    end do
!        read(8,*) empty
    close(38)
    !write(*,'(3f13.10)') fix_l(483,:,:)
!    write(*,'(3f13.10)') transpose(fix_p)
!    print *, shape(fix_p)


!---------------------------READ IN AIMD GEOMETRIES----------------------------

! read in energies and timesteps of the AIMD trajectory
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

    allocate(E_dft1(time))
    allocate(chaos(time))
    E_dft1=0.0
    chaos= 0

    read(17,'(A)') empty
    read(17,*) temp, temp, E_dft1(1)
    j=2
    chaos(1) = 1


    do i=2,time-2
        read(17,*) temp, temp, E_dft1(j)
        temp = E_dft1(j-1)-E_dft1(j)
        if (abs(temp) >= e_aimd_max) then   ! Throws out energies where energies changes too little
            j = j+1
            chaos(i) = 1
        end if
    end do
    close(17)
    ende2 = j-1
    E_dft1=E_dft1+25.019988

! Read in AIMD Energy
    allocate(E_dft(j))
    do i=1,j
        E_dft(i)=E_dft1(i)
    end do

! The coefficients of the transformation matrix.
    call open_for_read(18, position_of_l_and_p)
    read(18,'(A,/)') empty
    read(18,'(f13.10)') c11
    read(18,'(2f13.10)') c12, c22
    read(18,'(28X,f13.10,//)') c33

! correction because aimd a=4.205, not 4.201
!    isqrt_2=0.7071067812d0
!    c12=-a_lat*isqrt_2*cell_in(1)/2
!    c11=a_lat*isqrt_2*cell_in(1)
!    c22=a_lat*sqrt(3.)*isqrt_2*cell_in(1)/2
!    c33=(cell_in(3)-1)*a_lat/sqrt(3.)+13



! read AIMD geometries in
    allocate(aimd_l(time,3,16))
    allocate(aimd_p(time,3))
    aimd_p=0

    ende=time-2
    do j = 1, ende
        do k=1,16
            read(18,*)  aimd_l(j,:,k)
!            print *, read_l(:,k)
        end do
        read(18,*) aimd_p(j,:)
    end do
    close(18)

! Trouble with periodic boundery conditions might be solved here
!    do j = 1, ende
!        if (aimd_l(j,1,2)>0.5 .and. aimd_l(j,2,2)>0.99 .and. aimd_l(j,3,2)>0.99) then
!            aimd_l(j,2,2)=aimd_l(j,2,2) -1.0
!        end if
!    end do


! Throw those geometries out between which the energy changes little
    allocate(aimd_cor_l(ende2,3,16))
    allocate(aimd_cor_p(ende2,3))
    aimd_cor_l=0.0d0
    aimd_cor_p=0.0d0
    j = 0
    do i = 1, ende
        if (chaos(i) == 1) then
            j=j+1
            aimd_cor_l(j,:,:)=aimd_l(i,:,:)
            aimd_cor_p(j,:)=aimd_p(i,:)
        end if
    end do

!------------------------------------------------------------------------------
!                   HOW MUCH AIMD CONTRIBUTION DO WE WANT?
!                   ======================================
!------------------------------------------------------------------------------

    if (control < 100) then
        o=1 + ende2/control
        allocate(prae_read_l(o,3,16))
        allocate(prae_d_p(o,3))
        allocate(prae_E_dft(o))

        u = 1
        do q=1,ende2
            prae_read_l(q,:,:)=aimd_cor_l(u,:,:)
            prae_d_p(q,:) = aimd_cor_p(u,:)
            prae_E_dft(q)=E_dft(u)
            u=u+control
            if (u .gt. ende2) exit
        end do
        l_aimd=o


    ! combine the DFT and the AIMD-array
        time=o+npts
        allocate(read_l(time,3,16))
        allocate(d_p(time,3))
        allocate(E_all(time))
        read_l=0
        d_p=0

        read_l(1:npts,:,:) = fix_l(1:npts,:,:)
        read_l(npts+1:time,:,:) = prae_read_l(1:o,:,:)
        d_p(1:npts,:) = fix_p(1:npts,:)
        d_p(npts+1:time,:) = prae_d_p(1:o,:)
        E_all(1:npts) = E_fix(1:npts)
        E_all(npts+1:time)=prae_E_dft(1:o)

    ! only original DFT points for fit
    else if (control == 200) then
        time = npts
        allocate(read_l(time,3,16))
        allocate(d_p(time,3))
        allocate(E_all(time))
        read_l=0
        d_p=0
        read_l = fix_l
        d_p = fix_p
        l_aimd=0
        E_all=E_fix

    ! only AIMD points for fit
    else if (control == 201) then
        time=ende2
        allocate(read_l(ende2,3,16))
        allocate(d_p(ende2,3))
        allocate(E_all(time))
        read_l=0.0
        d_p=0
        read_l = aimd_cor_l
        d_p = aimd_cor_p
        l_aimd = time
        E_all=E_dft
    end if
 !   write(*,'(3f13.10)') transpose(d_p)
    ende = time


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
    temp=cell_in(1)*cell_in(2)*cell_in(3)
    if (rep==1) then
        n = temp*8 + temp
    elseif (rep==2) then
        n = temp*(8 + 16) + temp
    elseif (rep==3) then
        n = temp*(8 + 16 + 24) + temp
    elseif (rep==0) then
        n= temp
    end if

    n_l=n
! allocate arrays
    allocate(d_l(ende,3,n))
    d_l=0.d0

if (rep>0) then
! Translation of the entire story
    do q=1, ende
        ! Set the running parameters
        i = 1
        k = temp
        l = 1
        j = temp
        m = temp

        ! keep identiy
        d_l(q,:,i:k) = read_l(q,:,l:k)

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

    end do
end if
    if (rep==0) then
        d_l = read_l

    end if

! allocate important arrays
    allocate(r_l(ende,3,n))
    allocate(r_p(ende,3))
    allocate(cell_min(ende,3))
    allocate(cell_max(ende,3))
    allocate(celli(ende,3))
! transform the cell into cartesian coordinates
    cell_min=1000
    cell_max=-1000
    do q=1,ende
        do m=1, n
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
!       The selection of the minimum and maximum is not necessary for z since
!       we already know that it will be within 0 and c33 (due to the nature of
!       direct coordinates)
!            if ( r_l(q,3,m) .lt. cell_min(q,3))then
!                cell_min(q,3)=r_l(q,3,m)
!            else if ( r_l(q,3,m) .gt. cell_max(q,3))then
!                cell_max(q,3)=r_l(q,3,m)
!            end if
        end do
        ! calculate the coordinates for the hydrogen atom
        r_p(q,1) = c11*d_p(q,1)+c12*d_p(q,2)
        r_p(q,2) = c22*d_p(q,2)
        r_p(q,3) = c33*d_p(q,3)
    end do
    celli=cell_max-cell_min
    celli(:,3)=c33
    ! Write the overall array that contains both H and Au positions

    k=n+1
    allocate(x_all(ende,3,k))

    x_all(:,:,1)=r_p
    x_all(:,:,2:k)=r_l



    ! DON'T FORGET TO DEALLOCATE EVERYTHING!!!!!
    deallocate(read_l)
    deallocate(d_p)
    deallocate(d_l)
    deallocate(aimd_cor_l)
    deallocate(aimd_l)
    deallocate(aimd_p)

end subroutine

subroutine make_dft_formate(p_dft, l_dft)
!
! Purpose:
!       get the dft-data into same shape as aimd-data
!
    implicit none

    real(8), allocatable, dimension(:,:), intent(out)    :: p_dft
    real(8), allocatable, dimension(:,:), intent(out)    :: l_dft
    real(8), dimension(3,3) :: c_matrix, d_matrix
    real(8), allocatable, dimension(:,:)                 :: p_dir, aimd_form
    real(8), allocatable, dimension(:)                   :: E_dft

    real(8) :: c11, c12, c22, c33
    integer :: q, i, k, npts, ios, j, l
    real(8) :: empty, temp, e_max
    real(8), dimension(3,16) :: au_pos1


    character(len=100) :: pos_Au
    character(len=100) :: pos_H

    ! This procedure will read out two different files.
    ! the XDATA-file will contain the positions of the H and Au atoms
    ! in direct coordinates.
    ! The other file will contain the energies. Really, one does not need to make
    ! this file. We can just have a seperate section in the l_p_position program.


    pos_Au='file.dat'
    pos_H='hau111_plot.E.dat'

    ! Read in Au-positions
    c_matrix=0.d0
    d_matrix=0.0d0
    call open_for_read(28, pos_Au)
    read(28,'(A)') empty
    read(28,'(f15.10)') c_matrix(1,1)
    read(28,'(2f15.10)') c_matrix(1,2), c_matrix(2,2)
    read(28,'(28X,f13.10,/)') c_matrix(3,3)

    d_matrix(1,1) = 1/c_matrix(1,1)
    d_matrix(2,2) = 1/c_matrix(2,2)
    d_matrix(3,3) = 1/c_matrix(3,3)
    d_matrix(1,2) = -d_matrix(2,2)*c_matrix(1,2)*d_matrix(1,1)



    do i=1,16
        read(28,*) au_pos1(:,i)
    end do

    call open_for_read(29,pos_H)
        i=1
    do
        read(29,*,iostat=ios)
        if(ios <0) exit
        i=i+1
    end do
    rewind(29)
    npts = i-1

    allocate(p_dft(3,483))
    allocate(E_dft(483))
    e_max=10

    j=1
    do i=1, npts
        read(29,*) empty, p_dft(1,j), p_dft(2,j), p_dft(3,j), E_dft(j)
        if (abs(E_dft(j))<=e_max) j=j+1
    end do
    close(29)

    allocate(p_dir(3,483))

    do i=1,483
        p_dir(1,i) = p_dft(1,i)*d_matrix(1,1)+p_dft(2,i)*d_matrix(1,2)
        p_dir(2,i) = p_dft(2,i)*d_matrix(2,2)
        p_dir(3,i) = p_dft(3,i)*d_matrix(3,3)
        if (p_dir(1,i) < 0.0) p_dir(1,i)=p_dir(1,i)+1.0d0
        if (p_dir(2,i) < 0.0) p_dir(2,i)=p_dir(2,i)+1.0d0
        if (p_dir(3,i) < 0.0) p_dir(3,i)=p_dir(3,i)+1.0d0
    end do
!    write(*,'(3f15.10)') p_dir

    ! calculate how long aimd_form needs to be:
    l=(16+2)*483
    allocate(aimd_form(3,l))

    aimd_form = 9999

    do i = 1, 483
        j = 1 + 18*(i-1)
        k = j+16
        l = k
        aimd_form(:,j:k)=au_pos1(:,1:16)
        aimd_form(:,l)=p_dir(:,i)

    end do
!    write(*,'(3f15.10)') aimd_form

    call open_for_write(30,'DFT_aimd_form2x2x4.dat')
!    write(30,'(3f15.10)') aimd_form
    close(28)
    close(29)
    close(30)

end subroutine make_dft_formate


end module emt_init_data

