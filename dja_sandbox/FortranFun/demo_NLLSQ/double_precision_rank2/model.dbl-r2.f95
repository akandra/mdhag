
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

    F = A * .5 * (erf((E-E0)/W)+1.0) + Y0 + X(I,2)
    RES= Y(I) - F
    IF (JP.EQ.4) RRR(I) = F

    if(debug(1) .and. first) then
        print '(//(a))', '-----------------MODEL FIRST CALL---------------'
        print '((a),4f10.5)', '  x1,y1=   ', X(1,1), X(1,2), Y(1)
        print '((a),4f10.5)', '  x2,y2=   ', X(2,1), X(2,2), Y(2)
        print '((a),3I10)',   '  loc(x1))=', loc(x(1,1)), loc(x(1,2))
        print '((a),3I10)',   '  loc(x2))=', loc(x(2,1)), loc(x(2,2))
        print '((a),7f7.3/4x,7f7.3/)','  B=',B(1:14)
        print '((a))', '--------------END MODEL FIRST CALL----------------'
    end if

    first=.false.

    if (JP.EQ.2) JP = 3
    if (JP.EQ.4) RRR(I)=F

    RETURN

end subroutine MODEL
