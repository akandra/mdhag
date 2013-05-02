subroutine model( F, Y, X, RRR, I, JP )
    implicit none

    real(8)     Y(N), X(N,M), RRR(N)
    integer     I, JP
    real(8)     B, P, RES
    integer     N, M
    real(8)     A, E, E0, F, W, Y0

    logical     first, debug

    common/BLK1/B(20),P(20),RES,N,M
    common/DJA1/ debug(5),first

    E = X(I,1)
    A = B(1)
    E0= B(2)
    W = B(3)
    Y0= B(4)

      F = A * .5 * (erf((E-E0)/W)+1.0) + Y0

      RES = Y(I) - F

      IF (JP.EQ.2) JP = 3
      IF (JP.EQ.4) RRR(I) = F

      RETURN

end subroutine MODEL
