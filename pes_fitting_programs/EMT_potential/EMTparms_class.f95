module EMTparms_class
    implicit none

    type EMTparms
        character(2)::  name    = 'Whatever'

        real(8)       ::  eta2    = 1
        real(8)       ::  kappa   = 2
        real(8)       ::  lambda  = 3

        real(8)       ::  E0      = 4
        real(8)       ::  n0      = 5
        real(8)       ::  s0      = 6
        real(8)       ::  V0      = 7

    end type EMTparms


end module
