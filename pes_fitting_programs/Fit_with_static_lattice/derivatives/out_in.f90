module input
    use open_file

    ! Purpose:
    !           Read in input from md_tian.inp file
    !


    implicit none
    save

    real(8), parameter          :: pi       = 3.14159265359d0
    real(8), parameter          :: kB       = 8.61733238496d-5
    real(8), parameter          :: amu2mass = 103.6382d0
    real(8), parameter          :: deg2rad  = pi/180.0d0
    real(8), parameter          :: bohr2ang = 0.529177211d0

    real(8) :: a_lat, d
    real(8) :: einc, inclination, azimuth, Tsurf, step
    real(8) :: mass_p, mass_l
    integer :: n_p
    integer :: ntraj        ! Don't confuse with ntrajs.
    character(len=80) :: directory
    real(8) :: DEbin

    contains

subroutine readin(ntrajs)

    character(len=80) :: filename, empty,directory
    character(len=80) :: file_read_in
    character(len=80) :: buffer, label
    character(len=80) :: name_p, name_l
    character(len=3)  :: pot_p, pot_l

    integer :: i,j = 0, q, k, l, b
    integer :: ios=0, itemp, line=0, pos1
    integer :: start_tr, ntrajs
    integer, dimension(2) :: wstep
    integer :: nsteps

    real(8) :: rtemp



    ! Prompt the user to enter necessary information
    write(*,*) 'Please enter the lattice constant.'
    !read(*,*) a_lat
    a_lat=4.201
    write(*,*) 'Please enter the directory name.'
    file_read_in='ad_p274_ver_6x6x4_T300_2.8_60_60'
    !read(*,*) file_read_in
    write(*,*) 'Please enter the number of finished trajectories'
    !read(*,*) ntraj
    write(*,*) 'Please enter binning interval for eloss'
    !read(*,*) DEbin
    DEbin = 0.01
    ntraj = 1000000

     ! Calculate lattice properties. Surf and a_lat need to be updated as lattice changes
    d = a_lat/sqrt(3.0d0)
    filename = trim(file_read_in)//'/md_tian.inp.1'
    directory = trim(file_read_in)//'/results/'
    call open_for_read(38,filename)

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
                ntrajs = ntrajs*400
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
                if (wstep(1)==-1) wstep(2) = nsteps + 1
            case ('projectile')
                read(buffer, *, iostat=ios) name_p, mass_p, pot_p, itemp, &
                                            empty, empty, n_p
                mass_p=mass_p*amu2mass
            case ('lattice')
                read(buffer, *, iostat=ios) name_l, mass_l, pot_l, itemp, &
                                            empty, empty, rtemp
           case default
                print *, 'Skipping unrelated information at line', line, label
            end select
        end if
    end do ! ios
    close(38)


end subroutine readin

end module input
