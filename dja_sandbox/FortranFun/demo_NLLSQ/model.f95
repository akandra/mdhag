subroutine model( F, YDAT, XDAT, RRR, I, JP )
      DIMENSION YDAT(1), XDAT(1), RRR(1)
      COMMON /BLK1/   B(20),P(20),RE,N,M,KK

      E = xdat(i)
      A = B(1)
      E0= B(2)
      W = B(3)
      Y0= B(4)

      F = A * .5 * (erf((E-E0)/W)+1.0) + Y0

      RE= YDAT(I) - F

      IF (JP.EQ.2) JP = 3
      IF (JP.EQ.4) RRR(I) = F

      RETURN

end subroutine MODEL
