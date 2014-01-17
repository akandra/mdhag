module md_init
!
! Purpose:
!    prepare the system for the md simulations
!    It should be able to do the following things:
!      1. Read in initial configuration
!      2. Construct simulations cell from the read-in configuration
!      3. Assign initial thermal velocities to atoms
!
! Date          Author          History of Revision
! ====          ======          ===================
!01.10.2013     Sascha&Svenja   Implementation of 1.
!09.10.2013     Sascha&Svenja   Implementation of 2.
!

    use atom_class
    use open_file
    use useful_things
    implicit none
    save

    integer :: start_tr     = 1     ! a trajectory to start with
    integer :: ntrajs       = 10    ! number of trajectories
    real(8) :: einc         = 5     ! incidence energy (eV)
    real(8) :: inclination  = 0     ! incidence polar angle (degree)
    real(8) :: azimuth      = 0     ! incidence azimuthal angle (degree)
    real(8) :: Tsurf        = 300   ! surface temperature (Kelvin)
    real(8) :: step         = 0.1   ! time step in fs
    integer :: nsteps       = 100   ! number of steps
    integer :: wstep        = 1     ! interval to save data
    integer :: md_algo_l    = 1     ! md propagation algorithm: 1 - verlet
                                    !                           2 - beeman
                                    !                           3 - langevin
                                    !                           4 - langevin (series)
    integer :: md_algo_p    = 0     ! 0 means no projectile
    character(len=80) :: name_p = 'Elerium'
    character(len=80) :: pot_p = 'emt'
    character(len=80) :: key_p = 'empty'
    real(8) :: mass_p = 1.0d0
    integer :: npars_p = 0
    character(len=80) :: name_l, pot_l, key_l
    integer :: npars_l
    real(8) :: mass_l

    real(8),dimension(3,3) :: cell_mat, cell_imat ! simulation cell matrix and its inverse
    real(8), dimension(:), allocatable :: pars_l, pars_p ! potential parameters


contains

subroutine simbox_init(slab, teil)
!
! Purpose:
!           Initialise the entire system:
!               1. Geometry
!               2. Interaction Potential
!               3. Velocities
!
    implicit none

    type(atoms), intent(out) :: slab, teil   ! hold r, v and f for atoms in the box

    character(len=80) :: pos_init_file, confname_file
    character(len=80) :: buffer, label
    character(len= 7) :: confname
    character(len= 1) :: coord_sys
    character(len=80) :: mdpa_name

    integer :: pos1, ios = 0, line = 0
    integer :: rep = 2  ! number of repetition layers arounf an original cell
    integer :: n_l0, n_l, n_p, n_p0=0, itemp
    integer :: i, j, k, l, s, r
    integer :: randk = 13

    integer, dimension(3) :: celldim=(/2,2,4/)  ! input cell structure

    real(8) :: cellscale, v_pdof, vinc

    real(8), dimension(3,3) :: c_matrix, d_matrix

    real(8), dimension(:,:), allocatable :: start_l, start_p, d_l, d_p
    real(8), dimension(:,:), allocatable :: pos_l, vel_l, pos_p
!    real(8), dimension(:,:), allocatable :: d_ref

    logical :: exists

    if (iargc() == 0) stop " I need an input file"
    call getarg(1, pos_init_file)

    ! seed random number generator
    call random_seed(size=randk)
    call random_seed(put=randseed)


!------------------------------------------------------------------------------
!                       READ IN INPUT FILE
!                       ==================
!------------------------------------------------------------------------------

    call open_for_read(38, pos_init_file)

    ! ios < 0: end of record condition encountered or endfile condition detected
    ! ios > 0: an error is detected
    ! ios = 0  otherwise

    do while (ios == 0)
        read(38, '(A)', iostat=ios) buffer
        if (ios == 0) then
            line = line + 1

        ! Find the first instance of whitespace.  Split label and data.
            pos1 = scan(buffer, ' ')
            label = buffer(1:pos1)
            buffer = buffer(pos1+1:)

            select case (label)
            case('start')
                read(buffer,*,iostat=ios) start_tr
            case('ntrajs')
                read(buffer,*,iostat=ios) ntrajs
            case('Einc')
                read(buffer,*,iostat=ios) einc
            case('inclination')
                read(buffer,*,iostat=ios) inclination
                inclination = inclination*deg2rad
            case('azimuth')
                read(buffer,*,iostat=ios) azimuth
                azimuth = azimuth*deg2rad
            case('Tsurf')
                read(buffer,*,iostat=ios) Tsurf
            case('step')
                read(buffer,*,iostat=ios) step
            case('nsteps')
                read(buffer,*,iostat=ios) nsteps
            case('wstep')
                read(buffer,*,iostat=ios) wstep
            case ('projectile')
                read(buffer, *, iostat=ios) name_p, mass_p, pot_p, npars_p, &
                                            key_p, mdpa_name, n_p0
                mass_p=mass_p*amu2mass
                call lower_case(mdpa_name)
                select case (mdpa_name(1:3))
                    case ('ver')
                        md_algo_p = 1
                    case ('bee')
                        md_algo_p = 2
                    case ('lan')
                        md_algo_p = 3
                    case ('sla')
                        md_algo_p = 4
                   case default
                        print *, 'algorithm ', trim(mdpa_name), ' unknown'
                        stop
                end select

            case ('lattice')
                read(buffer, *, iostat=ios) name_l, mass_l, pot_l, npars_l, &
                                            key_l, mdpa_name
                mass_l=mass_l*amu2mass
                call lower_case(mdpa_name)
                select case (mdpa_name(1:3))
                    case ('ver')
                        md_algo_l = 1
                    case ('bee')
                        md_algo_l = 2
                    case ('lan')
                        md_algo_l = 3
                     case ('sla')
                        md_algo_l = 4
                   case default
                        print *, 'algorithm ', trim(mdpa_name), ' unknown'
                        stop
                end select
            case ('celldim')
                read(buffer, *, iostat=ios) celldim
            case ('rep')
                read(buffer, *, iostat=ios) rep
            case ('conf')
                read(buffer, *, iostat=ios) confname, confname_file
!
           case default
                print *, 'Skipping invalid label at line', line, label
            end select
        end if
    end do ! ios
    close(38)

!------------------------------------------------------------------------------
!                       READ IN CONFIGURATION
!                       =====================
!------------------------------------------------------------------------------

    call open_for_read(38, confname_file)

    if (confname == 'POSCAR') then

        read(38,*) buffer
        read(38,*) cellscale
        read(38,*) c_matrix
!        read(38,*) n_l0, n_p
        read(38,'(A)') buffer
        read(buffer, *, iostat=ios) n_l0, n_p
        read(38,*) coord_sys


        ! Construct simulation cell matrix and its inverse
        d_matrix = 0.0d0
        d_matrix(1,1) = 1.0d0/c_matrix(1,1)
        d_matrix(2,2) = 1.0d0/c_matrix(2,2)
        d_matrix(3,3) = 1.0d0/c_matrix(3,3)
        d_matrix(1,2) = -d_matrix(2,2)*c_matrix(1,2)*d_matrix(1,1)

        cell_mat = 0.0d0
        cell_mat(1:2,1:2) = c_matrix(1:2,1:2)*(2*rep + 1)
        cell_mat(  3,  3) = c_matrix(3,3)

        cell_imat = 0.0d0
        cell_imat(1:2,1:2) = d_matrix(1:2,1:2)/(2*rep + 1)
        cell_imat(  3,  3) = d_matrix(3,3)

        ! Read in coordinates
        allocate(start_l(3,n_l0))
        read(38,*) start_l
        if (n_p > 0) then
            allocate(start_p(3,n_p))
            read(38,*) start_p
        end if

        ! Transform the read in coordinates into direct if they are cartesians:
        if (coord_sys == 'C' .or. coord_sys == 'c') then
            start_l=matmul(d_matrix,start_l)
            if (n_p > 0) start_p=matmul(d_matrix,start_p)
        end if

    !------------------------------------------------------------------------------
    !                    Replicate input cell
    !                    ====================
    !------------------------------------------------------------------------------
    ! 1. Allocate d_l-array:
    ! The array size is determined by the amount of repetitions of the cell image:
    ! the number of lattice atoms in the cell has to be multiplied by the number of
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

        allocate(d_l(3,n_l))
        d_l = 0.0d0

        ! Replication
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


        allocate(pos_l(3,n_l), vel_l(3,n_l))
        pos_l = matmul(c_matrix,d_l)

        ! Sample velocities of lattice atoms from thermal distribution
        ! assuming the minimum energy configuration
        v_pdof = sqrt(2.0d0*kB*Tsurf/mass_l)
        n_l0 = n_l/celldim(3)*(celldim(3)-1)    ! Exclude fixed atoms

        vel_l = 0.0d0
        do i=1,n_l0
            vel_l(1,i) = normal(0.0d0,v_pdof)
            vel_l(2,i) = normal(0.0d0,v_pdof)
            vel_l(3,i) = normal(0.0d0,v_pdof)
        enddo
        ! Set c.-of-m. velocity to zero
        vel_l(1,1:n_l0) = vel_l(1,1:n_l0) - sum(vel_l(1,1:n_l0))/n_l0
        vel_l(2,1:n_l0) = vel_l(2,1:n_l0) - sum(vel_l(2,1:n_l0))/n_l0
        vel_l(3,1:n_l0) = vel_l(3,1:n_l0) - sum(vel_l(3,1:n_l0))/n_l0

        if (md_algo_p > 0 .and. n_p == 0 .and. n_p0 == 0) then
            md_algo_p = 0
            print *, "Warning: Number of projectiles both in POSCAR and input file is 0."
            print *, "         Calculations will continue without projectile."
        end if

        itemp = n_p*(2*rep+1)**2

        ! Projectile initialization
        if (md_algo_p > 0) then    ! projectile existence justified

            if (n_p0 == 0) n_p0 = itemp

            allocate(d_p(3,n_p0))
            d_p = 0.0d0
            if (itemp < n_p0 ) then
                print *, "Warning: Number of projectiles larger than can be produced"
                print *, "         from repititions of POSCAR file."
                print *, "         All projectile positions set to zero."
            else
                i=1
                j=1
                l=1

                do r =-rep, rep
                do s =-rep, rep
                do l =   1, n_p

                    if (j > n_p0) exit
                    d_p(1,j) = start_p(1,l)+r
                    d_p(2,j) = start_p(2,l)+s
                    d_p(3,j) = start_p(3,l)

                    j=j+1

                end do
                end do
                end do
            end if

            allocate(pos_p(3,n_p0))
            pos_p = matmul(c_matrix,d_p)

        end if

        n_p = n_p0

    else

        read(38,*) cell_mat
        read(38,*) cell_imat
        read(38,*) n_l, n_p

        allocate(pos_l(3,n_l),vel_l(3,n_l),start_p(3,n_p))

        read(38,*) start_p
        read(38,*) pos_l
        read(38,*) vel_l

    endif

    close(38)

    !   Read in potential parameters

    allocate(pars_l(npars_l))

    call open_for_read(23,key_l)
    read(23,'(/)')
    do i = 1, npars_l
        read(23,*) buffer, pars_l(i)
    end do
    close(23)

    if (n_p > 0) then

        allocate(pars_p(npars_p))

        call open_for_read(23,key_p)
        read(23,'(/)')
        do i = 1, npars_p
            read(23,*) buffer, pars_p(i)
        end do
        close(23)
    end if

    ! Create slab objects
    slab = atoms(n_l)
    ! Assign slab positions and velocities
    slab%r = pos_l
    slab%v = vel_l

    ! Create projectile objects
    if (n_p > 0) then
        teil = atoms(n_p)
    else
        teil = atoms(1)
        teil%n_atoms = 0
    end if

    ! Assign projectile positions and velocities
    if (n_p > 0) then

        teil%r = pos_p(:,1:n_p)
        if (confname .ne. 'POSCAR') then
            vinc = sqrt(2.0d0*einc/mass_p)
            teil%v(1,:) =  vinc*sin(inclination)*cos(azimuth)
            teil%v(2,:) =  vinc*sin(inclination)*sin(azimuth)
            teil%v(3,:) = -vinc*cos(inclination)
        end if
        deallocate(pos_p, d_p, start_p)
    end if

    ! Create a directory for trajectory data
    inquire(file='trajs',exist=exists)
    if (.not. exists) then
        call system('mkdir trajs')
    end if

    deallocate(vel_l, pos_l, d_l, start_l)

end subroutine simbox_init

end module md_init

