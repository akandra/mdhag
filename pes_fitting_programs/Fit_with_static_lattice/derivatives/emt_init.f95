module emt_init_data
    !use EMTparms_class
    implicit none
    save

    integer :: time
    integer :: n_l, n_p
    real(8), dimension(3,6) :: celli
    logical  :: just_l
    logical :: one_p
    integer :: rep
    real(8) :: E_pseudo
    real(8), dimension(:,:), allocatable :: x_ref

contains

! COMMENT:  To improve the performance, it is better to redefine arrays with coordinates
!           in the way (coordinate, point) to escape the non-contiguous array problem
!           by passing a deferred-array to a subroutine


subroutine l_p_position(b, a_lat, rep, cell_in, control,e_aimd_max, just_l, one_p, time,  &
                                                l_aimd, n_l,n_p,celli, x_all, E_all)
!
! Purpose:
!           Reads in the gold and hydrogen positions from AIMD
!
! Date          Author          History of Revision
! ====          ======          ===================
! 01.07.2013    S. M. Janke     original
! 05.08.2013    Sascha&Svenja   fix dft read in
! 01.10.2013    Sascha&Svenja   modified celli
!
    implicit none

! Declare in and output
    real(8), intent(in)                 :: a_lat    ! lattice constant
    integer, intent(in)                 :: rep      ! Repetitions of lattice (1,2,3)
    integer, dimension(3), intent(in)   :: cell_in  ! contains cell geometry (2x2x4)
    integer, intent(out)                  :: time        ! different configuration and energies
    real(8),allocatable,dimension(:), intent(out)                 :: E_all    ! DFT energy
    real(8),dimension(3,6), intent(out)               :: celli
    integer, intent(out)                :: n_l, n_p, l_aimd  ! number of l atoms,number of p atoms, number of aimd contributions
    real(8), dimension(:,:,:), allocatable, intent(out) :: x_all
    integer , intent(inout)                        :: control
    real(8), intent(in) :: e_aimd_max
    logical :: just_l, one_p
    integer :: b


! other variables
    character(len=40)   :: position_of_l_and_p, fix_position
    character(len=40)    :: energy_l_and_p, fix_energy
    integer             :: i, j,k,u,l, m,n,q, ios, start, ende, ende2,o,npts,s,r
    real(8)             :: c11,c12,c22, c33
    real(8)             :: temp, e_max
    real(8)             :: isqrt_2
    real(8)   :: empty, empty3(3), c_matrix(3,3), d_matrix(3,3)
    character(len=35):: emptys
    ! Arrays. Structure: timestep, x,y,z of atoms
    real(8), dimension(:,:,:), allocatable :: aimd_l,  prae_read_l, read_l ! array of au after read in
    real(8), dimension(:,:),allocatable    :: aimd_p,  prae_d_p, read_p     ! array of h
    real(8), dimension(:,:,:), allocatable :: d_l, d_p     ! array after multiplying lattice image
    real(8), dimension(:,:), allocatable :: fix_l, dfix, d_ref, x_ref1 ! array of au after read in
    real(8), dimension(:,:),allocatable    :: fix_p     ! array of h
    real(8), dimension(:,:,:), allocatable:: r_l, r_p      ! Position of lattice atoms, position of particle

    real(8), dimension(:), allocatable   :: E_fix    ! energy of fixed lattice geometires
    real(8),allocatable,dimension(:)                 :: E_dft1, prae_E_dft    ! read-in-dft-energy

print *, b
    if (b == 1) then
    position_of_l_and_p = 'data/traj005/XDATCAR_005.dat'
    energy_l_and_p =      'data/traj005/analyse_005.out'
    end if


    if (b == 2) then
    position_of_l_and_p = 'data/traj010/XDATCAR_010.dat'
    energy_l_and_p =      'data/traj010/analyse_010.out'
    end if

    if (b == 3) then
    position_of_l_and_p = 'data/traj801/XDATCAR_801.dat'
    energy_l_and_p =      'data/traj801/analyse_801.out'
    end if

    if (b == 4) then
    position_of_l_and_p = 'data/traj814/XDATCAR_814.dat'
    energy_l_and_p =      'data/traj814/analyse_814.out'
    end if

    if (b == 5) then
    position_of_l_and_p = 'data/traj817/XDATCAR_817.dat'
    energy_l_and_p =      'data/traj817/analyse_817.out'
    end if

    if (b == 6) then
    position_of_l_and_p = 'data/traj818/XDATCAR_818.dat'
    energy_l_and_p =      'data/traj818/analyse_818.out'
    end if

    if (b == 7) then
    position_of_l_and_p = 'data/traj820/XDATCAR_820.dat'
    energy_l_and_p =      'data/traj820/analyse_820.out'
    end if

    if (b == 8) then
    position_of_l_and_p = 'data/traj821/XDATCAR_821.dat'
    energy_l_and_p =      'data/traj821/analyse_821.out'
    end if

    if (b == 9) then
    position_of_l_and_p = 'data/traj825/XDATCAR_825.dat'
    energy_l_and_p =      'data/traj825/analyse_825.out'
    end if

    if (b == 10) then
    position_of_l_and_p = 'data/traj831/XDATCAR_831.dat'
    energy_l_and_p =      'data/traj831/analyse_831.out'
    end if

    if (b == 11) then
    position_of_l_and_p = 'data/traj832/XDATCAR_832.dat'
    energy_l_and_p =      'data/traj832/analyse_832.out'
    end if

    if (b == 12) then
    position_of_l_and_p = 'data/traj833/XDATCAR_833.dat'
    energy_l_and_p =      'data/traj833/analyse_833.out'
    end if

    if (b == 13) then
    position_of_l_and_p = 'data/traj858/XDATCAR_858.dat'
    energy_l_and_p =      'data/traj858/analyse_858.out'
    end if
    ! 187 config, some from aimd, some from dragging Au out
    !position_of_l_and_p = 'data/trajjustau/XDATCAR_move_Au.out'
    !energy_l_and_p =      'data/trajjustau/enery_move_Au.dat'

    ! 447 config, some from aimd, some from dragging Au out
    !position_of_l_and_p = 'data/trajjustau/XDATCAR_more_drAu.out'
    !energy_l_and_p =      'data/trajjustau/energy_more_drAu.dat'



    ! 388 config, all from aimd
    !position_of_l_and_p = 'data/trajjustau/XDATCAR_more_Au.out'
    !energy_l_and_p =      'data/trajjustau/energy_more_Au.dat'


    ! 137 config, all from aimd
    !position_of_l_and_p = 'data/trajjustau/XDATCAR.out'
    !energy_l_and_p =      'data/trajjustau/energy.dat'

    fix_position = 'data/au111_2x2x4.POSCAR'
    fix_energy = 'data/energy.dat'!'data/hau111_plot.E.dat'

!------------------------------------------------------------------------------
!                       READ IN GEOMETRIES AND ENERGIES
!                       ===============================
!------------------------------------------------------------------------------

! -------------------------READ IN FIXED LATTICE-------------------------------
! read in energies
    e_max=20.0

    call open_for_read(39,fix_energy)
    i=1
    do
        read(39,*,iostat=ios) empty, empty, empty, empty, empty
        if(ios <0) exit
        if (Abs(empty+ 24.995689)<=e_max ) i=i+1
    end do
    npts = i-1

! Read in Equilibrium position geometries
    if (control == 300) then
        call open_for_read(69,'data/points_dft.dat')
        npts=1500
        allocate(fix_p(npts,3), E_fix(npts))
        E_fix = 0.0d0
        do j = 1,npts
            read(69,*) fix_p(j,:)
        end do

        control = 200
        close(69)
    else

        allocate(E_fix(npts),fix_p(npts,3))

        rewind(39)

        j=1
        do
            read(39,*,iostat=ios) empty, empty3(1), empty3(2), empty3(3), empty
            if(ios <0) exit
            !if (abs(empty)<=e_max) then
            !if (empty<=e_max .and. empty > e_max*2) then
            if (Abs(empty+ 24.995689)<=e_max ) then
                E_fix(j) = empty
                fix_p(j,:) = empty3
                j=j+1
            end if
        end do
        close(39)
    end if
! Read in the geometry of the fixed lattice
    call open_for_read(38,fix_position)
    read(38,*) emptys
    read(38,*) empty
    read(38,*) c_matrix

    d_matrix = 0.0d0
    d_matrix(1,1) = 1.0d0/c_matrix(1,1)
    d_matrix(2,2) = 1.0d0/c_matrix(2,2)
    d_matrix(3,3) = 1.0d0/c_matrix(3,3)
    d_matrix(1,2) = -d_matrix(2,2)*c_matrix(1,2)*d_matrix(1,1)

    read(38,*) k, emptys

    allocate(fix_l(3,k))
    read(38,*) fix_l

    close(38)

! Celli contains both c_ and d_matrix for pbc-procedure
    celli=0.0d0
    celli(1:2,1:2)=c_matrix(1:2,1:2)*(2*rep+1)
    celli(3,3)=c_matrix(3,3)

    celli(1:2,4:5)=d_matrix(1:2,1:2)/(2*rep+1)
    celli(3,6)=d_matrix(3,3)

    fix_p=matmul(fix_p,transpose(d_matrix))     ! Turn H-positions into direct coordinates

    do i=1,npts
        if (fix_p(i,1) < -0.001d0) fix_p(i,1)=fix_p(i,1)+1.0d0
        if (fix_p(i,2) < -0.001d0) fix_p(i,2)=fix_p(i,2)+1.0d0
        if (fix_p(i,3) < -0.001d0) fix_p(i,3)=fix_p(i,3)+1.0d0
    end do

    fix_l=matmul(d_matrix,fix_l)                ! Turn Au-positions into direct coordinates
    fix_l=Nint(fix_l*10000)/10000.0d0

    do i=1,k
!        if (fix_l(1,i) < 0.0d0) fix_l(1,i)=fix_l(1,i)+1.0d0
!        if (fix_l(2,i) < 0.0d0) fix_l(2,i)=fix_l(2,i)+1.0d0
        if (fix_l(3,i) < -0.0001d0) fix_l(3,i)=fix_l(3,i)+1.0d0
    end do

    ! Define geometry for reference energy
    allocate(d_ref(3,k+1))
    d_ref(:,1) = (/0.0d0,0.0d0,0.29618952180d0/)
    d_ref(:,2:k+1) = fix_l




!---------------------------READ IN AIMD GEOMETRIES----------------------------

! read in energies and timesteps of the AIMD trajectory
! The geometries and energies are in different files (which is why we are
! opening two here.)
    call open_for_read(17,energy_l_and_p)
    call open_for_read(18, position_of_l_and_p)
    ! The coefficients of the transformation matrix.
    read(18,'(A,/)') emptys
    read(18,*) d_matrix
    read(18,'(//)')

! According to the energy file, we define our time steps
    time = 0
    do
        read(17,*,iostat=ios)
        if(ios <0) exit
        time=time+1
    end do
    rewind(17)
    ende=time-2
!    time=602
!    ende=time-2

    allocate(E_dft1(time))
    allocate(aimd_l(time,3,k))
    allocate(dfix(3,k))
    allocate(aimd_p(time,3))
    aimd_p=0.0
    aimd_l=0.0
    E_dft1=0.0

    read(17,'(A)') empty
    read(17,*) temp, temp, E_dft1(1)
    read(18,*) aimd_l(1,:,:)
    if (just_l .eqv. .false.) read(18,*) aimd_p(1,:)
    if (just_l .eqv. .true.)  aimd_p(1,:) = (/0.0d0,0.0d0,6.0d0/)
    j=2

    do i=2,ende
        read(17,*) temp, temp, E_dft1(j)
        read(18,*) aimd_l(j,:,:)
        if (just_l .eqv. .false.) read(18,*) aimd_p(j,:)
        if (just_l .eqv. .true.)  aimd_p(1,:) = (/0.0d0,0.0d0,6.0d0/)
        temp = E_dft1(j-1)-E_dft1(j)
        if (abs(temp) >= e_aimd_max) j = j+1
    end do
    close(17)
    close(18)

    ende2 = j-1
    ! Place AIMD-atoms close to corresponding equilibrium positions
     do i=1,ende2
        dfix = aimd_l(i,:,:) - fix_l
        do q=1,k
               aimd_l(i,1,q)=aimd_l(i,1,q) - ANINT(dfix(1,q))
               aimd_l(i,2,q)=aimd_l(i,2,q) - ANINT(dfix(2,q))
               aimd_l(i,3,q)=aimd_l(i,3,q) - ANINT(dfix(3,q))
        end do
    end do
   !E_dft1=E_dft1+25.019988 ! reference energy

!------------------------------------------------------------------------------
!                   HOW MUCH AIMD CONTRIBUTION DO WE WANT?
!                   ======================================
!------------------------------------------------------------------------------

    if (control < 100) then
        o=1 + ende2/control
        allocate(prae_read_l(o,3,k))
        allocate(prae_d_p(o,3))
        allocate(prae_E_dft(o))

        q = 1
        do u=1,ende2,control
            prae_read_l(q,:,:)=aimd_l(u,:,:)
            prae_d_p(q,:) = aimd_p(u,:)
            prae_E_dft(q)=E_dft1(u)
            q=q+1
        end do
       l_aimd=o


    ! combine the DFT and the AIMD-array
        time=o+npts
        allocate(read_l(time,3,16))
        allocate(read_p(time,3))
        allocate(E_all(time))
        read_l=0
        read_p=0
        E_all=0

        do i=1,npts
            read_l(i,:,:) = fix_l
        end do
        read_l(npts+1:time,:,:) = prae_read_l
        read_p(1:npts,:) = fix_p
        read_p(npts+1:time,:) = prae_d_p
        E_all(1:npts) = E_fix
        E_all(npts+1:time)=prae_E_dft

    ! only original DFT points for fit
    else if (control == 200) then
        time = npts
        allocate(read_l(time,3,k))
        allocate(read_p(time,3))
        allocate(E_all(time))
        read_l=0
        read_p=0
        E_all=0
        do i=1,npts
            read_l(i,:,:) = fix_l
        end do
        read_p = fix_p
        l_aimd=0
        E_all=E_fix

    ! only AIMD points for fit
    ! The first geometry here will be later used as a reference energy, but
    ! to make things easier, it is calculated in this array
    else if (control == 201) then
        time=ende2
        allocate(read_l(time,3,k))
        allocate(read_p(time,3))
        allocate(E_all(time))
        read_l=0.0
        read_p=0
        E_all=0
        read_l = aimd_l
        read_p = aimd_p(1:time,:)
        l_aimd = time
        E_all=E_dft1(1:time)
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
    n_l=temp*(2*rep+1)**2
    n_p = (2*rep+1)**2
    if (one_p .eqv. .true.) n_p = 1
    !E_all = E_all * (2*rep+1)**2 ! energy needs to be increased to account for all images. From old times, when still with VASP reference energy
    E_all = E_all + 24.995689 !(E_all + 24.995689)* (2*rep+1)**2


! allocate arrays
    allocate(d_l(ende,3,n_l))
    allocate(d_p(ende,3,n_p))
    allocate(x_ref1(3,n_p+n_l))
    d_l=0.d0
    d_p = 0.0d0
    x_ref1=0.0d0


! Translation of array under consideration.
    do q=1, ende
        i = 1
        j= 1
        do r =-rep, rep
            do s=-rep, rep
                ! Set the running parameters
                d_l(q,1,i:i+temp-1) = read_l(q,1,:)+r
                d_l(q,2,i:i+temp-1) = read_l(q,2,:)+s
                d_l(q,3,i:i+temp-1) = read_l(q,3,:)
                i = i+temp

                if (one_p .eqv. .false.) then
                    d_p(q,1,j) = read_p(q,1)+r
                    d_p(q,2,j) = read_p(q,2)+s
                    d_p(q,3,j) = read_p(q,3)
                    j = j + 1
                end if

           end do
        end do
        if (one_p .eqv. .true.) d_p(q,:,1) = read_p(q,:)
    end do



        ! Translation of array for reference energy
    i=1
    j=1
    do r =-rep, rep
        do s=-rep, rep
            ! Set the running parameters
            x_ref1(1,n_p+i:n_p+n_l) = d_ref(1,2:k+1)+r
            x_ref1(2,n_p+i:n_p+n_l) = d_ref(2,2:k+1)+s
            x_ref1(3,n_p+i:n_p+n_l) = d_ref(3,2:k+1)
            i = i+temp

            if (one_p .eqv. .false.) then
                x_ref1(1,j) = d_ref(1,1)+r
                x_ref1(2,j) = d_ref(2,1)+s
                x_ref1(3,j) = d_ref(3,1)
                j=j+1
            end if

        end do
    end do
    if (one_p .eqv. .true.) x_ref1(:,1) = d_ref(:,1)


 !allocate important arrays
    allocate(r_l(ende,3,n_l))
    allocate(r_p(ende,3,n_p))

     do q=1,ende
            r_l(q,:,:) = matmul(c_matrix,d_l(q,:,:))
            r_p(q,:,:) = matmul(c_matrix,d_p(q,:,:)) ! is this right?
     end do

    ! transform reference coordinates into cartesian:
    x_ref1=matmul(c_matrix,x_ref1)


! Sascha, you have commented that you understand what we did here and that it is right.
! Dont ask questions, just ACCEPT it.
    ! Write the overall array that contains both H and Au positions


    if (just_l .eqv. .false.) then
        k=n_l+n_p
        allocate(x_all(time,3,k))
        x_all(:,:,1:n_p)=r_p
        x_all(:,:,n_p+1:k)=r_l
        x_ref=x_ref1
    else if (just_l .eqv. .true.) then
        allocate(x_all(ende,3,n_l))
        x_all=r_l
        x_ref=x_ref1(:,n_p+1:n_p+n_l)
        n_p = 0
    else
        print*,"There's something wrong with your allocation of x_all in emt_init"
    end if



    ! DON'T FORGET TO DEALLOCATE EVERYTHING!!!!!
    deallocate(d_l, r_l, r_p, E_dft1)
    deallocate(aimd_l, read_l, d_p )
    deallocate(aimd_p)
    deallocate(E_fix,fix_l,fix_p )

end subroutine l_p_position


end module emt_init_data

