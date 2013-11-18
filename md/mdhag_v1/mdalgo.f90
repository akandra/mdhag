module mdalgo
!
!
! Purpose :
!           Contains all the algorithms needed for the md procedure
!
! Date          Author              History of Revision
! ====          ======              ===================
! 23.10.2013    Sascha & Svenja     Original

    use atom_class
    use open_file
    use force
    implicit none

contains

subroutine norm_dist(vec1,vec2,norm)

    real(8),allocatable, dimension(:,:) :: vec1, vec2
    real(8) :: norm, n1 = 0.0d0, n2 = 0.0d0
    integer :: i
    norm = 0.0d0

    do i=1,3
        norm = norm + dot_product(vec1(i,:),vec2(i,:))
        n1   = n1   + dot_product(vec1(i,:),vec1(i,:))
        n2   = n2   + dot_product(vec2(i,:),vec2(i,:))
    enddo

    n1 = max(n1,n2)
    norm=Sqrt(norm/n1)

end subroutine norm_dist


subroutine beeman_1(r,v,a,aold)
    !
    ! Purpose:
    !           1st step of Refson-Beeman algorithm,
    !           K. Refson, Physica 131B, (1985), 256.
    !           Moldy User's Manual.
    !

!    type(atom), dimension(:), allocatable :: slab, teilchen

    real(8), dimension(:,:), allocatable :: r,v,a,aold
    real(8) :: step_sq

    step_sq = step * step / 6.0d0

    r = r + step*v + step_sq*(4.0*a - aold)


end subroutine beeman_1

subroutine predict(v,vp,a,aold)

    real(8), dimension(:,:), allocatable :: r,v,vp,a,aold

    vp = v + 0.5d0*step * (3.0d0 * a - aold)



end subroutine predict

subroutine beeman_2(v,vc,a,aold,avold)

    real(8), dimension(:,:), allocatable :: v,vc,a,aold,avold

    vc = v + step/6.0d0 * (2.0d0 * a + 5.0d0 * aold - avold)



end subroutine beeman_2


subroutine newton(f, a, minv)

    real(8), dimension(:,:), allocatable :: f, a
    real(8), dimension(:), allocatable :: minv

    a(1,:) = f(1,:) * minv
    a(2,:) = f(2,:) * minv
    a(3,:) = f(3,:) * minv


end subroutine newton

subroutine fric(a, v, zetadmass)

    real(8), dimension(:,:), allocatable :: a, v
    real(8), intent(in) :: zetadmass ! zetadmass is the friction coefficient devided by the mass

    a = a - zetadmass * v


end subroutine fric


end module mdalgo
