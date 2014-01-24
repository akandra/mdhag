module output
use open_file
use atom_class
use md_init

implicit none
save
integer :: save_counter = 1
logical :: overwrite = .true.

contains

subroutine full_conf(slab, teil, itraj)
    !
    ! Purpose:
    !           Prints out information necessary to continue simulation

    type(atoms) :: slab, teil
    integer :: ios, itraj
    character(len=8) str
    character(len=80) filename

    write(str,'(I8.8)') save_counter

    filename = 'conf/mxt_conf'//str//'.bin'

    open (753,file=filename,form="unformatted", status='replace', &
                    action='write', iostat=ios)

    ! Here comes the output

    write(753) itraj        ! Number of trajectory
    write(753) step         ! time step
    write(753) Epot
    write(753) Tsurf        ! Surface temperature
    ! number of species
    if (teil%n_atoms .ne. 0) then
        write(753) 2
    else
        write(753) 1
    end if
    ! name, number of atoms, no of fixed atoms
    ! masses, potential, number of parameters, file name parameters, paramter values
    ! propagator
    write(753) name_l, slab%n_atoms, slab%nofix,&
               mass_l, pot_l, npars_l,key_l
    write(753 )pars_l, md_algo_l

    write(753) a_lat        ! lattice constant
    write(753) cell_mat     ! Cell matrix
    write(753) cell_imat    ! inverse cell matrix
    write(753) slab%r, slab%v, slab%a, slab%dens
    if (teil%n_atoms .ne. 0) then
        write(753) name_p, teil%n_atoms, teil%nofix, &
                   mass_p, pot_p, npars_p, key_p
        write(753) pars_p, md_algo_p
        write(753) teil%r, teil%v, teil%a, teil%dens
    end if


    close(753)
    save_counter = save_counter+1

end subroutine full_conf

subroutine out_short(slab, teil,Epot, itraj, q, rmin_p, col_int)
    !
    ! Purpose :
    !           Prints out final information at end of trajectory
    !

    type(atoms) :: slab, teil
    real(8) :: Epot, Ekin_l, Ekin_p
    integer :: itraj, q
    real(8), dimension(:,:), allocatable :: rmin_p
    integer, dimension(:), allocatable   :: col_int
    character(len=8) str
    character(len=80) filename

    write(str,'(I8.8)') itraj
    filename = 'traj/mxt_fin'//str//'.dat'
    if (overwrite) then
        call open_for_write(753,filename)
    else
        call open_for_append(753,filename)
        write(753,*)
        write(753,'(A14,f15.5)') 'Time    (fs) = ', q*step
    end if

    Ekin_l = E_kin(slab,mass_l)
    Ekin_p = E_kin(teil,mass_p)
    write(753,'(A14,f15.5)') 'E_kin_p (eV) = ', Ekin_p
    write(753,'(A14,f15.5)') 'E_kin_l (eV) = ', Ekin_l
    write(753,'(A14,f15.5)') 'E_pot   (eV) = ', Epot
    write(753,'(A14,f15.5)') 'E_total (eV) = ', Epot + Ekin_l + Ekin_p
    write(753,'(A8)')'r_p (A):'
    write(753,'(3f15.5)') teil%r
    write(753,'(A11)')'v_p (A/fs):'
    write(753,'(3e15.5)') teil%v

    if (.not. overwrite) then
        write(753,'(A12)')'r_min_p (A):'
        write(753,'(3f15.5)') rmin_p
        write(753,'(A21)')'Time at surface (fs):'
        write(753,'(1000f15.5)') col_int*step
   end if

    close(753)
    overwrite = .false.

end subroutine out_short

subroutine out_detail(output_info, n, itraj)
    !
    ! Purpose :
    !           Prints out a lot of trajectory information
    !

    real(8), dimension(:,:), allocatable :: output_info
    integer :: n, i
    integer :: itraj
    character(len=8) str
    character(len=80) filename

    write(str,'(I8.8)') itraj
    filename = 'traj/mxt_trj'//str//'.dat'

    call open_for_write(753,filename)
    write(753,'(1000A15)') 'time (fs)', 'E_pot (eV)', 'E_kin_l (eV)','E_kin_p (eV)',&
                           'E_total (eV)', 'dens (A^-3)', 'r_p (A)', 'v_p (A/fs)'

    do i = 1, n
        write(753,'(10000e15.5)') i*wstep(2)*step, output_info(:,i)
    end do

    close(753)


end subroutine out_detail

end module output
