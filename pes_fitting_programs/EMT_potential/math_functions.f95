module math_functions
    !
    ! Purpose:
    !           Sports all kinds of useful functions.
    !
    implicit none

    contains

    subroutine fexp (r, a, b, c, f, df)
        ! Subroutine value of function f = a exp( b (c - r)) and its differential
        real(8), intent(in) :: r, a, b, c       ! independent variables
        real(8), intent(out) :: f               ! function value
        real(8), dimension(4), intent(out) :: df    ! array of partial derivatives

        f = a * exp( b * (c - r))

        df(1) = -b
        df(2) = 1/a
        df(3) = c-r
        df(4) = b
        df = f * df


    end subroutine fexp


end module math_functions
