program testAtomClass

    use atom_class
    implicit none

    type (atom)              :: H, D

    print *, "H positions with using defaults ", H%r
    H = atom(0)
    print *, "H positions with using atom(0)  ", H%r
    D = H
    print *, "D positions with using D = H    ", H%r


end program
