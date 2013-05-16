program testAtomClass

    use atom_class
    implicit none



    type (atom)  ,target            :: particle
    type (atoms) ,target            :: lat1
    type (atoms) ,target            :: lat2

    integer, parameter              :: n_atoms=4
    integer                         :: i,j

    real(8), dimension(3,n_atoms)   :: test1
    real(8), dimension(3,n_atoms)   :: test2
    real(8), dimension(n_atoms)     :: test3
    real(8), dimension(n_atoms)     :: test4

    real(8), pointer, dimension(:)     :: patom
    real(8), pointer, dimension(:,:)   :: patoms

    namelist / particle_list /  particle
    !namelist / lattice_list1 /  lat1   not allowed becuase has allocatable
    !namelist / lattice_list2 /  lat2   not allowed becuase has allocatable


    print *, 'particle r       v       f'
    print '(6x,I5,3F8.2,(6X,3f8.2))', particle

    print *
    print *,'output using namelist'
    write(*,nml=particle_list)

    !--------------------------------------------------------!
    !   allocate allocatable elements of atoms manually      !
    !--------------------------------------------------------!
    allocate (lat1%r(3,n_atoms))
    allocate (lat1%v(3,n_atoms))
    allocate (lat1%f(3,n_atoms))

    !--------------------------------------------------------!
    !   initialize r,v,f components of atoms structure       !
    !--------------------------------------------------------!
    do j=1,n_atoms
        do i=1,3
            lat1%r(i,j) = i     + 10*j
            lat1%v(i,j) = i + 3 + 10*j
            lat1%f(i,j) = i + 6 + 10*j
        end do
    enddo

    !-------------------------------------------------------------------------!
    !   instantiate a atoms object using defined constructor                  !
    !   allocate allocatable elements of atoms2 with user defined constructor !
    !-------------------------------------------------------------------------!
    lat2 = atoms(n_atoms)


    !-------------------------------------------------------------------------!
    !   print the results                                                     !
    !-------------------------------------------------------------------------!
    print *, '--------------------------------------------------------------------------'
    print *, '                              lattice 1                                   '
    print *, '--------------------------------------------------------------------------'
    print *, ''
    do j=1,n_atoms
        print 1000, j, "r=", (lat1%r(i,j),i=1,3), "v=", (lat1%v(i,j),i=1,3), "f=", (lat1%v(i,j),i=1,3)
    enddo

    print *, ''
    print *, '--------------------------------------------------------------------------'
    print *, '                              lattice 2                                   '
    print *, '--------------------------------------------------------------------------'
    do j=1,n_atoms
        print 1000, j, "r=", (lat2%r(i,j),i=1,3), "v=", (lat2%v(i,j),i=1,3), "f=", (lat2%f(i,j),i=1,3)
    enddo

    !-------------------------------------------------------------------------!
    !                                                                         !
    !   print locations of variable in particle, lat1, and lat2               !
    !                                                                         !
    !-------------------------------------------------------------------------!
    print *
    print *, '--------------------------------------------------------------------------'
    print *, '                        location of variables in particle                 '
    print *, '--------------------------------------------------------------------------'
    !rint 1100, 'n_atoms : ', loc(particle%n_atoms)
    print 1100, 'position: ', ((loc(particle%r(i))),i=1,3)
    print 1100, 'velocity: ', ((loc(particle%v(i))),i=1,3)
    print 1100, 'force   : ', ((loc(particle%f(i))),i=1,3)
    print *

    print *
    print *, '--------------------------------------------------------------------------'
    print *, '                        location of variables in lattice 1                '
    print *, '--------------------------------------------------------------------------'
    print 1100, 'n_atoms : ', loc(lat1%n_atoms)
    print 1100, 'position: ', (((loc(lat1%r(i,j))),i=1,3), j=1,n_atoms)
    print 1100, 'velocity: ', (((loc(lat1%v(i,j))),i=1,3), j=1,n_atoms)
    print 1100, 'force   : ', (((loc(lat1%f(i,j))),i=1,3), j=1,n_atoms)
    print *

    print *
    print *, '--------------------------------------------------------------------------'
    print *, '                        location of variables in lattice 2                '
    print *, '--------------------------------------------------------------------------'
    print 1100, 'n_atoms : ', loc(lat2%n_atoms)
    print 1100, 'position: ', (((loc(lat2%r(i,j))),i=1,3), j=1,n_atoms)
    print 1100, 'velocity: ', (((loc(lat2%v(i,j))),i=1,3), j=1,n_atoms)
    print 1100, 'force   : ', (((loc(lat2%f(i,j))),i=1,3), j=1,n_atoms)
    print *

    print *, ''
    print *, '--------------------------------------------------------------------------'
    print *, '             location of variables in declare arrays test1, test2         '
    print *, '--------------------------------------------------------------------------'
    print 1100, 'test1   : ', (((loc(test1(i,j))),i=1,3), j=1,n_atoms)
    print 1100, 'test2:    ', (((loc(test2(i,j))),i=1,3), j=1,n_atoms)
    print *

    print *, ''
    print *, '--------------------------------------------------------------------------'
    print *, '             location of variables in declare arrays test3, test4         '
    print *, '--------------------------------------------------------------------------'
    print 1100, 'test3   : ', ((loc(test3(i))),i=1,3)
    print 1100, 'test4:    ', ((loc(test4(i))),i=1,3)
    print *
    print *, 'size of test1=',size(test1)


    write(*,'(//,"use pointers to access positions")')
    patom  => particle%r
    write(*,2000) (patom(i),i=1,9)

1000 format (i5, a5, 3f5.1, 2(5x, a5, 3f5.1) )
1100 format ((a)/(3i9,2x))
2000 format (9f5.0)
end program
