program testAtomClass

    use atom_class
    implicit none

    type (atom)  ,target            :: particle     ! used to demonstrate atom object
    type (atoms) ,target            :: lat1         ! used to demonstrate manual allocation of atoms object components
    type (atoms) ,target            :: lat2         ! used to demonstrate use of construct for
    ! alloction of atoms object components
    integer, parameter              :: n_atoms=4    ! number of atoms for the atoms objects
    integer                         :: i,j

    namelist / particle_list /  particle
    !namelist / lattice_list1 /  lat1   not allowed becuase has allocatable elements.

    print '((a),9f5.1)', 'particle=', particle

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
    print *, 'atom      position            velocity            force'
    print *, '--------------------------------------------------------------------------'
    do j=1,n_atoms
        print 1000, j,((lat1%r(i,j)),i=1,3), ((lat1%v(i,j)),i=1,3), ((lat1%f(i,j)),i=1,3)
    end do
    print *, ''
    print *, '--------------------------------------------------------------------------'
    print *, '                              lattice 2                                   '
    print *, '--------------------------------------------------------------------------'
    print *, 'atom      position            velocity            force'
    print *, '--------------------------------------------------------------------------'
    do j=1,n_atoms
        print 1000, j, (lat2%r(i,j),i=1,3), (lat2%v(i,j),i=1,3), (lat2%f(i,j),i=1,3)
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


    1000 format (i5, 5x, 3f5.1, 2(5x, 3f5.1) )
    1100 format ((a)/(3i9,2x))
    2000 format (9f5.0)
end program
