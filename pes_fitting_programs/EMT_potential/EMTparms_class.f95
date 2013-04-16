module EMTparms_class
    implicit none

    type EMTparms
        character(2)::  name

        real        ::  eta2
        real        ::  kappa
        real        ::  lambda

        real        ::  E0
        real        ::  n0
        real        ::  s0
        real        ::  V0

    end type EMTparms


end module
