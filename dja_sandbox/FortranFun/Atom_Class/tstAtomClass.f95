program testAtomClass

    use atom_class
    implicit none

    type (atom)             :: particle
    type (atoms)            :: lat1
    type (atoms2)           :: lat2
    integer                 :: n_atoms=2
    integer                 :: i,j


    print *, "particle=", particle

    !*********************************************************
    !   direct way to allocate allocatable elements of atoms *
    !*********************************************************
    allocate (lat1%r(3,n_atoms))
    allocate (lat1%v(3,n_atoms))
    allocate (lat1%f(3,n_atoms))
    !*********************************************************
    !   initialize r,v,f compents of atoms                   *
    !*********************************************************
    do j=1,n_atoms
        do i=1,3
            lat1%r(i,j) = i     + 10*j
            lat1%v(i,j) = i + 3 + 10*j
            lat1%f(i,j) = i + 6 + 10*j
        end do
    enddo

    !**************************************************************************
    !   allocate allocatable elements of atoms2 with user defined constructor *
    !**************************************************************************
    lat2 = atoms2(n_atoms)

    !*********************************************************
    !   print the results                                    *
    !*********************************************************

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

1000 format (i5, a5, 3f5.1, 2(5x, a5, 3f5.1) )
end program
