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
    use useful_things
    use open_file
    use force
    implicit none

    real(8):: xi=1.0d0

contains

subroutine propagator_1(s, md_algo, imass)
    !
    ! Purpose:
    !           select algorithm and calculate 1st and 2nd step
    !

    type(atoms) :: s
    integer :: md_algo
    real(8) :: imass



    select case (md_algo)

        case (1) ! velocity Verlet Algorithm
            call verlet_1(s)

        case (2) ! Refson-Beeman Algorithm
            call beeman_1(s)

        case (3) ! Langevin
            call langevin_1(s,imass)

        case (4) ! Langevin up to 2nd order
            call langevins_1(s,imass)

   end select

end subroutine propagator_1

subroutine propagator_2(s, md_algo, imass)
    !
    ! Purpose:
    !           select algorithm and calculate 3rd and 4th step
    !

    type(atoms) :: s
    integer :: md_algo
    real(8) :: imass, norm

        select case (md_algo)

            case (1) ! velocity Verlet Algorithm
                call newton(s, imass)
                call verlet_2(s)

            case (2) ! Refson-Beeman Algorithm
                do
                    call newton(s, imass)
                    call beeman_2(s)
                    call norm_dist(s%vp, s%vc, 3*s%n_atoms, norm)
                    if (norm < 1.0e-007) then
                        s%v  = s%vc
                        s%au = s%ao
                        s%ao = s%a
                        exit
                    else
                        s%vp = s%vc
                    end if
                end do

            case (3) ! Langevin
                call ldfa(s)
                call newton(s, imass)
                call langevin_2(s,imass)

            case (4) ! Langevin up to 2nd order
                call ldfa(s)
                call newton(s, imass)
                call langevins_2(s,imass)

       end select



end subroutine propagator_2

subroutine verlet_1(s)
    !
    ! Purpose:
    !           1st and 2nd steps of velocity Verlet algorithm,
    !           Allen & Tildesley, Computer Simulation of Liquids (1987).
    !

    type(atoms) :: s

    ! half-step velocities
    s%v = s%v + 0.5d0*step*s%a
    ! new positions
    s%r = s%r + step*s%v

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

subroutine langevin_1(s, imass)
    !
    ! Purpose:
    !           1st step of Langevin Dynamics algorithm,
    !           Dellago et al. JCP 108 (1998) 1964
    !

    type(atoms) :: s
    real(8) :: imass, temp
    real(8), dimension(  s%n_atoms) :: c0, c1, c2, xidt, sigma_r, sigma_v, c_rv
    real(8), dimension(3,s%n_atoms) :: randy, cofm
    integer :: i

    xidt = (s%dens*step)
    c0 = exp(-xidt)

    c1 = (1.0d0 - c0)/s%dens
    c2 = (step - c1)/xidt

    temp = kB*Tsurf*imass
    sigma_r = sqrt(temp*(2.d0*xidt - 3.d0 + 4.d0*c0 - c0*c0))/s%dens
    sigma_v = sqrt(temp*(1.d0 - c0*c0))
    c_rv = temp*(1.d0 - c0)**2/(sigma_r*sigma_v*s%dens)

    do i =1, s%n_atoms
        randy(1,i) = normal(0.0d0,1.0d0)
        randy(2,i) = normal(0.0d0,1.0d0)
        randy(3,i) = normal(0.0d0,1.0d0)
    end do

    cofm(1,:) = sigma_r*randy(1,:)
    cofm(2,:) = sigma_r*randy(2,:)
    cofm(3,:) = sigma_r*randy(3,:)

    cofm(1,:) = cofm(1,:) - sum(cofm(1,:))/s%n_atoms
    cofm(2,:) = cofm(2,:) - sum(cofm(2,:))/s%n_atoms
    cofm(3,:) = cofm(3,:) - sum(cofm(3,:))/s%n_atoms

    ! propagate positions
    s%r(1,:) = s%r(1,:) + c1*s%v(1,:) + c2*step*s%a(1,:) + cofm(1,:)
    s%r(2,:) = s%r(2,:) + c1*s%v(2,:) + c2*step*s%a(2,:) + cofm(2,:)
    s%r(3,:) = s%r(3,:) + c1*s%v(3,:) + c2*step*s%a(3,:) + cofm(3,:)

    cofm(1,:) = sigma_v*c_rv*randy(1,:)
    cofm(2,:) = sigma_v*c_rv*randy(2,:)
    cofm(3,:) = sigma_v*c_rv*randy(3,:)

    cofm(1,:) = cofm(1,:) - sum(cofm(1,:))/s%n_atoms
    cofm(2,:) = cofm(2,:) - sum(cofm(2,:))/s%n_atoms
    cofm(3,:) = cofm(3,:) - sum(cofm(3,:))/s%n_atoms

    ! partially propagate velocities
    s%v(1,:) = c0*s%v(1,:) + (c1 - c2)*s%a(1,:) + cofm(1,:)
    s%v(2,:) = c0*s%v(2,:) + (c1 - c2)*s%a(2,:) + cofm(2,:)
    s%v(3,:) = c0*s%v(3,:) + (c1 - c2)*s%a(3,:) + cofm(3,:)

end subroutine langevin_1

subroutine langevin_2(s, imass)
    !
    ! Purpose:
    !           1st step of Langevin Dynamics algorithm,
    !           Dellago et al. JCP 108 (1998) 1964
    !

    type(atoms) :: s
    real(8) :: imass, temp
    real(8), dimension(  s%n_atoms) :: c0, c1, c2, xidt, sigma_r, sigma_v, c_rv
    real(8), dimension(3,s%n_atoms) :: randy, cofm
    integer :: i

    xidt = (s%dens*step)
    c0 = exp(-xidt)
    c1 = (1.0d0 - c0)/s%dens
    c2 = (step - c1)/xidt

    temp = kB*Tsurf*imass
    sigma_r = sqrt(temp*(2.d0*xidt - 3.d0 + 4.d0*c0 - c0*c0))/s%dens
    sigma_v = sqrt(temp*(1.d0 - c0*c0))
    c_rv = temp*(1.d0 - c0)**2/(sigma_r*sigma_v*s%dens)

    do i =1, s%n_atoms
        randy(1,i) = normal(0.0d0,1.0d0)
        randy(2,i) = normal(0.0d0,1.0d0)
        randy(3,i) = normal(0.0d0,1.0d0)
    end do

    cofm(1,:) = sigma_v*sqrt(1 - c_rv*c_rv)*randy(1,:)
    cofm(2,:) = sigma_v*sqrt(1 - c_rv*c_rv)*randy(2,:)
    cofm(3,:) = sigma_v*sqrt(1 - c_rv*c_rv)*randy(3,:)

    cofm(1,:) = cofm(1,:) - sum(cofm(1,:))/s%n_atoms
    cofm(2,:) = cofm(2,:) - sum(cofm(2,:))/s%n_atoms
    cofm(3,:) = cofm(3,:) - sum(cofm(3,:))/s%n_atoms

    ! partially propagate velocities
    s%v(1,:) = s%v(1,:) + c2*s%a(1,:) + cofm(1,:)
    s%v(2,:) = s%v(2,:) + c2*s%a(2,:) + cofm(2,:)
    s%v(3,:) = s%v(3,:) + c2*s%a(3,:) + cofm(3,:)

end subroutine langevin_2


subroutine langevins_1(s, imass)
    !
    ! Purpose:
    !           1st step of Langevin Dynamics algorithm,
    !           Allen & Tildesley, Computer Simulation of Liquids (1987).
    !           Li & Wahnström, Phys. Rev. B (1992).
    !

    type(atoms) :: s
    real(8) :: c0, c1, c2
    real(8) :: xidt, imass, sigma_r, sigma_v, c_rv,temp
    real(8), dimension(3,s%n_atoms) :: randy, cofm
    integer :: i

    xidt = (xi*step)
    c0 = 1.0d0 - xidt + 0.50d0*xidt*xidt
    c1 = (1.0d0 - 0.50d0*xidt + 2.0d0*twelfth*xidt*xidt)*step
    c2 = (0.5d0 - 2.0d0*twelfth*xidt + 0.5d0*twelfth*xidt*xidt)*step

    temp = kB*Tsurf*imass
    sigma_r = step*sqrt(temp*(8.0d0*twelfth - 0.5d0*xidt)*xidt)
    sigma_v = sqrt(temp*2.0d0*(1.0d0 - xidt)*xidt)
    c_rv = 0.5d0*sqrt3*(1.0d0 - 0.125d0*xidt)

    do i =1, s%n_atoms
        randy(1,i) = normal(0.0d0,1.0d0)
        randy(2,i) = normal(0.0d0,1.0d0)
        randy(3,i) = normal(0.0d0,1.0d0)
    end do

    cofm = sigma_r*randy
    cofm(1,:) = cofm(1,:) - sum(cofm(1,:))/s%n_atoms
    cofm(2,:) = cofm(2,:) - sum(cofm(2,:))/s%n_atoms
    cofm(3,:) = cofm(3,:) - sum(cofm(3,:))/s%n_atoms
    ! propagate positions
    s%r = s%r + c1*s%v + c2*step*s%a + cofm

    cofm = sigma_v*c_rv*randy
    cofm(1,:) = cofm(1,:) - sum(cofm(1,:))/s%n_atoms
    cofm(2,:) = cofm(2,:) - sum(cofm(2,:))/s%n_atoms
    cofm(3,:) = cofm(3,:) - sum(cofm(3,:))/s%n_atoms
    ! partially propagate velocities
    s%v = c0*s%v + (c1 - c2)*s%a + cofm
end subroutine langevins_1

subroutine langevins_2(s, imass)
    !
    ! Purpose:
    !           1st step of Langevin Dynamics algorithm,
    !           Allen & Tildesley, Computer Simulation of Liquids (1987).
    !           Li & Wahnström, Phys. Rev. B (1992).
    !Ŕ

    type(atoms) :: s
    real(8) :: c0, c1, c2
    real(8) :: xidt, ixidt, imass, sigma_r, sigma_v, c_rv, temp
    real(8), dimension(3,s%n_atoms) :: randy, cofm
    integer :: i

    xidt = (xi*step)
    ixidt = 1.0d0/xidt
    c0 = 1.0d0 - xidt + 0.50d0*xidt*xidt
    c1 = (1.0d0 - 0.50d0*xidt + 2.0d0*twelfth*xidt*xidt)*step
    c2 = (0.5d0 - 2.0d0*twelfth*xidt + 0.5d0*twelfth*xidt*xidt)*step

    temp = kB*Tsurf*imass
    sigma_r = step*sqrt(temp*(8.0d0*twelfth - 0.5d0*xidt)*xidt)
    sigma_v = sqrt(temp*2.0d0*(1.0d0 - xidt)*xidt)
    c_rv = 0.5d0*sqrt3*(1.0d0 - 0.125d0*xidt)

    do i =1, s%n_atoms
        randy(1,i) = normal(0.0d0,1.0d0)
        randy(2,i) = normal(0.0d0,1.0d0)
        randy(3,i) = normal(0.0d0,1.0d0)
    end do

    cofm = sigma_v*sqrt(1 - c_rv*c_rv)*randy

    cofm(1,:) = cofm(1,:) - sum(cofm(1,:))/s%n_atoms
    cofm(2,:) = cofm(2,:) - sum(cofm(2,:))/s%n_atoms
    cofm(3,:) = cofm(3,:) - sum(cofm(3,:))/s%n_atoms

    s%v = s%v + c2*s%a + cofm

end subroutine langevins_2

subroutine newton(s, minv)
    !
    ! Purpose:
    !           Newton equation
    !
    type(atoms) :: s
    real(8) :: minv

    s%a = s%f*minv

end subroutine newton

end module mdalgo
