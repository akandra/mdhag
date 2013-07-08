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

   subroutine fexplus (r, a, b, c, f, df)
        ! Subroutine value of function f = a exp( b (c - r)) and its differential
        ! Compared to fexp, this subroutine adds the calculated f and df to the
        ! input values
        real(8), intent(in) :: r, a, b, c       ! independent variables
        real(8), intent(inout) :: f               ! function value
        real(8), dimension(4), intent(inout) :: df    ! array of partial derivatives

        real(8), dimension(4) :: dfd            ! dummy variable
        real(8) :: fdummy

        fdummy = a * exp( b * (c - r))


        dfd(1) = -b
        dfd(2) = 1/a
        dfd(3) = c-r
        dfd(4) = b

        f = f + fdummy
        df = df + fdummy * dfd



    end subroutine fexplus


end module math_functions
