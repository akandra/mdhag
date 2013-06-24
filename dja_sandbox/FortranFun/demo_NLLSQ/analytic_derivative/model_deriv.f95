subroutine model( F, Y, X, RRR, I, JP )
    implicit none

    real(8)     Y(N), X(N), RRR(N)
    integer     I, JP
    real(8)     B, P, RES
    integer     N, M, K, IB, IP
    real(8)     A, E, E0, F, W, Y0
    real(8)     cerf, derf

    logical     first, debug
    logical     vary(20)
    integer     ii

    common/BLK1/B(20),P(20),RES,N,M,K
    COMMON/BLK5/IB(20),IP
    common/DJA1/ debug(5),first

    E = X(I)
    A = B(1)
    E0= B(2)
    W = B(3)
    Y0= B(4)
    vary = .true.

    if (IP>0) then
        do ii=1,IP
            vary(IB(ii))=.false.
        end do
    endif

    select case(jp)
        case(1)
            F = A * .5 * (erf( (E-E0)/W ) + 1.0) + Y0
            RES = Y(I) - F

        case(2)
            cerf = .5 * (erf( (E-E0)/W )+ 1.0)
            derf = 0.56419 * A * exp( -(E-E0)**2 / W**2 )

            F = A * cerf + Y0
            RES = Y(I) - F

            ! if any parameters are held constant, don't compute the partial derivatives
            ! to save computational time
            ! this is not important in this example but could be if evaluation of derivatives
            ! takes a lot of computation.
            if (vary(1)) P(1) = cerf
            if (vary(2)) P(2) = -derf / W
            if (vary(3)) P(3) = -derf * (E-E0)/W**2
            if (vary(4)) P(4) = 1

        case(4)
            F = A * .5 * (erf((E-E0)/W)+1.0) + Y0
            RRR(I)=F

        case default
            print *,"error: model called with jp=",JP
            stop 300
    end select

    RETURN

end subroutine MODEL
