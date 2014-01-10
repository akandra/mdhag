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

subroutine verlet_1(s)
    !
    ! Purpose:
    !           1st and 2nd steps of velocity Verlet algorithm,
    !           Allen & Tildesley, Computer Simulation of Liquids (1987).
    !

    type(atoms) :: s

    ! new positions
    s%r = s%r + step*s%v + 0.5d0*step*step*s%a
    ! half-step velocities
    s%v = s%v + 0.5d0*step*s%a

end subroutine verlet_1

subroutine verlet_2(s)
    !
    ! Purpose:
    !           3rd step of velocity Verlet algorithm,
    !           Allen & Tildesley, Computer Simulation of Liquids (1987).
    !

    type(atoms) :: s

    ! new velocities
    s%v = s%v + 0.5d0*step*s%a

end subroutine verlet_2

subroutine beeman_1(s)
    !
    ! Purpose:
    !           1st and 2nd steps of Refson-Beeman algorithm,
    !           K. Refson, Physica 131B, (1985), 256.
    !           Moldy User's Manual.
    !
    type(atoms) :: s
    real(8) :: step_sq

    step_sq = step * step / 6.0d0

    ! new positions
    s%r = s%r + step*s%v + step_sq*(4.0*s%ao - s%au)
    ! predicted velocities
    s%vp = s%v + 0.5d0*step*(3.0d0*s%ao - s%au)

end subroutine beeman_1

subroutine beeman_2(s)
    !
    ! Purpose:
    !           4th step of Refson-Beeman algorithm,
    !           K. Refson, Physica 131B, (1985), 256.
    !           Moldy User's Manual.
    !
    type(atoms) :: s

    s%vc = s%v + step/6.0d0*(2.0d0*s%a + 5.0d0*s%ao - s%au)

end subroutine beeman_2

subroutine langevin_1(s)
    !
    ! Purpose:
    !           1st and 2nd steps of Langevin Dynamics algorithm,
    !           Allen & Tildesley, Computer Simulation of Liquids (1987).
    !

    type(atoms) :: s

    ! new positions
    s%r = s%r + c1*s%v + c2*s%a + rrandom
    ! half-step velocities
    s%v = s%v + 0.5d0*step*s%a

end subroutine langevin_1

subroutine newton(s, minv)
    !
    ! Purpose:
    !           Newton equation
    !
    type(atoms) :: s
    real(8) :: minv

    s%a = s%f*minv

end subroutine newton

subroutine norm_dist(vec1, vec2, length, norm)
    !
    ! Purpose: normalised distance between 2 vectors
    !

    integer :: length
    real(8), dimension(length) :: vec1, vec2
    real(8) :: norm, n1, n2

    norm = dot_product(vec1 - vec2, vec1 - vec2)
    n1   = dot_product(vec1,vec1)
    n2   = dot_product(vec2,vec2)

    n1 = max(n1,n2)
    if (n1 .eq. 0.0d0) then
        norm = 0.0d0
    else
        norm=Sqrt(norm/n1)
    end if

end subroutine norm_dist

end module mdalgo
