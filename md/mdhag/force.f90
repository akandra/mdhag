module force
    !
    ! Purpose:
    !           Gather all potentials to calculate the forces for MD
    !           MAY THE FORCE BE WITH YOU!
    ! Definition of Zero of Energy:
    !                               The particle is at infinite distance to the
    !                               lattice.
    !
    ! Date          Author          History of Revison
    ! ====          ======          ==================
    ! 09.10.2013    Sascha&Svenja   Original
    !
    use atom_class
    use md_init

    implicit none
    save

    real(8), private, parameter                 :: beta    = 1.8093997906d0
    real(8), private, parameter                 :: twelfth= 0.0833333333333333d0
    integer, private, dimension(3), parameter   :: b       = (/12, 6, 24/)

    contains

subroutine pes (teilchen,slab,Epot)
!
! Purpose:
!       Calculates the energy according to the potential chosen
!
implicit none

!------------------------------------------------------------------------------
!                                   PREAMBLE
!                                   ========
!------------------------------------------------------------------------------

! declare variables and parameters that are passed down from the program

    type(atom), dimension(:), allocatable, intent(inout)    :: teilchen, slab
    real(8), intent(out)                                    :: Epot

! declare the variables that appear in the subroutine

    integer :: i,j, k                   ! running parameter
    real(8) :: r                        ! distance
    real(8) :: energy                   ! Batterie

! emt variables
    real(8) :: rcut, rr, acut           ! values to calculate cut-off
    real(8) :: igamma1p, igamma2p       ! inverse gamma for particle
    real(8) :: igamma1l, igamma2l      ! inverse gamma for lattice atoms
    real(8) :: voldegamma2l, vopdegamma2p
    real(8) :: theta                    ! variable for cut-off calculation
    real(8) :: chilp, chipl             ! mixing between lattice (l) and particle (p)
    real(8) :: sigma_pl                 ! mixed contribution to neutral sphere,
    real(8) :: s_p                      ! neutral sphere radius particle
    real(8), dimension(spec_l%n) :: sigma_ll ! contribution to ns, l only
    real(8), dimension(spec_l%n) :: s_l_exp  ! contribution to ns, l only
    real(8), dimension(spec_l%n) :: sigma_lp ! mixed contribution to ns
    real(8), dimension(spec_l%n) :: s_l      ! neutral sphere radius lattice atoms
    real(8) :: V_pl, V_lp, V_ll         ! Pair potential contributions
    real(8) :: vref_l, vref_p           ! reference pair pot. contrib.
    real(8) :: Ecoh                     ! cohesive energy of part & lattice atoms
    real(8), dimension(3) :: xl, xp     ! for cal. cut-off
    real(8), dimension(3) :: rnnl, rnnp ! next neighbour distance for cut-off
    real(8) :: betas0_l, betas0_p       ! beta * s0= for l and p
    real(8) :: kappadbeta_l, kappadbeta_p ! beta * kappa for l and p
    real(8), dimension(3) :: r3temp,n3temp ! temporary array variable
    real(8) :: rtemp                    ! temporary variable

    real(8) :: E_ref                    ! reference energy
    real(8) :: Ecoh_ref                  !cohesive reference energy

! For the reference energy
    real(8)                 :: rn_ltemp(spec_l%n), r3temp1(3), rtemp1
    real(8), dimension(spec_l%n) :: f_lx, f_l , s_l_ref, p_l
    real(8)                 :: vref_l_ref

! Declaration for Morse Potential
    real(8) :: expar, expar2

! Derivatives EMT
    real(8) :: gij, gip, gpi, gtemp, Qij, f_p
    real(8), dimension(3,spec_l%n*(spec_l%n-1)/2) :: nij, gij1
    real(8), dimension(3,spec_l%n) :: nip, npi, force_l, dV_lp, dV_pl
    real(8), dimension(3,spec_l%n) :: dvref_l, dvref_p
    real(8), dimension(3) ::force_p

!----------------------VALUES OF FREQUENT USE ---------------------------------
! definition of a few values that appear frequently in calculation
    energy=0.0d0
select case (spec_l%pot)

    case ('emt')

        betas0_l = beta * pars_l(7)
        kappadbeta_l = pars_l(6) / beta

    !------------------------------------------------------------------------------
    !                                  CUT-OFF
    !                                  =======
    !------------------------------------------------------------------------------
    ! We use the distance to the next-next nearest neighbours as cut-off.
    ! FOR FUTURE REVISION:
    !            cut-off should be defined via lattice constant _AND_ changeable.
    !    rcut = a_lat * sqrt3 * isqrt2
    !    rr = 4 * rcut / (sqrt_3 + 2)
    !    acut = 9.21024/(rr -rcut) ! ln(10000)

        rcut = betas0_l * sqrt3
        rr = 4.0d0 * rcut / (sqrt3 + 2.0d0)
        acut = 9.21024d0/(rr -rcut) ! ln(10000)

   ! Distances to the considered neighbours
        rnnl(1) = betas0_l
        rnnl(2) = rnnl(1) * sqrt2
        rnnl(3) = rnnl(1) * sqrt3
        xl = b * twelfth / (1.0d0 + exp(acut*(rnnl-rcut)))

    !-----------------------------------GAMMA--------------------------------------
    ! Gamma enforces the cut-off together with theta (see below)
    ! Gamma is defined as inverse.
        r3temp = rnnl-betas0_l
        igamma1l = 1.0d0 / sum(xl*exp(-pars_l(1) * r3temp))
        igamma2l = 1.0d0 /sum(xl*exp(-kappadbeta_l * r3temp))
        voldegamma2l = pars_l(5)*igamma2l
    ! in case of alloy

        if (spec_p%pot == 'emt') then

            betas0_p = beta * pars_p(7)
            kappadbeta_p = pars_p(6) / beta
            ! 'coupling' parameters between p and l
            chilp = pars_p(2) / pars_l(2)
            chipl = 1.0d0 / chilp

            rnnp(1) = betas0_p
            rnnp(2) = rnnp(1) * sqrt2
            rnnp(3) = rnnp(1) * sqrt3
            xp = b * twelfth/ (1.0d0 + exp(acut*(rnnp-rcut)))

            r3temp = rnnp-betas0_p
            igamma1p = 1.0d0 / sum(xp*exp(-pars_p(1) * r3temp))
            igamma2p = 1.0d0 / sum(xp*exp(-kappadbeta_p * r3temp))
            vopdegamma2p = pars_p(5)*igamma2p

        end if

   ! Initializing some accumulators

        sigma_ll = 0.0d0
        sigma_pl = 0.0d0
        V_ll = 0.0d0
        V_lp = 0.0d0
        V_pl = 0.0d0
        vref_p = 0.0d0
        Ecoh= 0.0d0
        force_l = 0.0d0
        force_p = 0.0d0
        dV_lp=0.0d0
        dV_pl=0.0d0
        dvref_l=0.0d0
        dvref_p=0.0d0


    case('morse')
        ! Hier gibt es nichts zu sehen. Bitte weitergehen.

    case default

        print *, 'pes subroutine: unknown lattice pes'
        stop

end select
    k = 0
    do i = 1, spec_l%n
        do j = i+1, spec_l%n

        !-----------------PERIODIC BOUNDARY CONDITIONS LATTICE-----------------
        ! Because we want them.

            r3temp(1) = slab(i)%r(1)-slab(j)%r(1)
            r3temp(2) = slab(i)%r(2)-slab(j)%r(2)
            r3temp(3) = slab(i)%r(3)-slab(j)%r(3)

            ! transform distances into direct coordinates
            r3temp= matmul(celli(1:3,4:6),r3temp)

            r3temp(1)=r3temp(1)-Anint(r3temp(1))
            r3temp(2)=r3temp(2)-Anint(r3temp(2))
            r3temp(3)=r3temp(3)-Anint(r3temp(3))

            r3temp=matmul(celli(1:3,1:3),r3temp)

            ! Length of the vector rjk, e.i.: distance between atom j and k
            r =  sqrt(sum(r3temp**2))
            ! drjk/dri; unit vector that points into the direkt of the vector
            ! between j and k.
            r3temp = r3temp / r
            nij(:, k+j-i) = r3temp


            select case (spec_l%pot)
            case('emt')
            !---------------------------THETA LATTICE------------------------------
            ! Theta enforces the cut-off together with gamma (see above). This
            ! function enacts cutoff by reducing contributions of atoms outside the
            ! cut-off to zero.

                rtemp = exp( acut * (r - rcut) )
                theta = 1.0d0 / (1.0d0 + rtemp )
                gij =   theta * acut * rtemp
                Qij =   gij + kappadbeta_l
                gij =   gij + pars_l(1)

            !----------------------------SIGMA LATTICE-----------------------------
            ! Sigma is a contribution to the neutral sphere radius.
            ! It is a list in which for each lattice atom, the contributions of the
            ! others are summed up. To enforce the cut-off, it will be later
            ! corrected by gamma.
            ! sigma_pp does not exist because there is only a single particle.

                rtemp = theta*exp(-pars_l(1) * (r - betas0_l) )

                sigma_ll(i) = sigma_ll(i) + rtemp
                sigma_ll(j) = sigma_ll(j) + rtemp

                gij = -gij * rtemp

            !-----------------------PAIR POTENTIAL LATTICE-------------------------
            ! For the lattice only.
            ! Will later be subjected to gamma to complete the cut-off.
            ! The particle has no pair potential contribution since there is only
            ! one and thus does not have a partner to interact with.

                rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
                V_ll = V_ll + rtemp

            ! Since the derivative of V_ll is one of the additive parts to the
            ! force, this derivative is written directly into the force
            ! Multiplication with vol and igamma2l outside of loop
                Qij = Qij * rtemp
                force_l(:,i) = force_l(:,i) + Qij * nij(:,k+j-i)
                force_l(:,j) = force_l(:,j) - Qij * nij(:,k+j-i)



            ! We need nij for the calculation of the reference energy.
                gij1(:,k+j-i)=nij(:,k+j-i)*gij

            case('morse')

                expar = exp(-pars_l(2)*(r-pars_l(3)))
                expar2 = expar**2
                energy = energy + pars_l(1) * (expar2 - 2*expar)

            case default

                print *, 'pes subroutine, lattice: unknown lattice pes'
                stop

            end select

        end do
        ! Here, we perform operations that we need to either go back into the
        ! previous loop (like advancing variable k)
        k=k+spec_l%n-i


    !-----------------PERIODIC BOUNDERY CONDITIONS PARTICLE--------------------

        r3temp(1) = slab(i)%r(1)-teilchen(1)%r(1)
        r3temp(2) = slab(i)%r(2)-teilchen(1)%r(2)
        r3temp(3) = slab(i)%r(3)-teilchen(1)%r(3)

    ! transform distances into direct coordinates
        r3temp= matmul(celli(1:3,4:6),r3temp)


        r3temp(1)=r3temp(1)-Anint(r3temp(1))
        r3temp(2)=r3temp(2)-Anint(r3temp(2))
        r3temp(3)=r3temp(3)-Anint(r3temp(3))

        r3temp=matmul(celli(1:3,1:3),r3temp)
        r =  sqrt(sum(r3temp**2))
        r3temp = r3temp / r
        nip(:,i) = r3temp


        select case (spec_p%pot)

        case('emt')

        !----------------------------THETA PARTICLE--------------------------------

            rtemp = exp( acut * (r - rcut) )
            theta = 1.0d0 / (1.0d0 + rtemp )
            gtemp =   theta * acut * rtemp


        !-------------------------------MIXED SIGMA--------------------------------
        ! Contributions of both particle and lattice to neutral sphere radius
        ! To fully include the cut-off, we correct them later by gamma.


            rtemp = theta*exp(-pars_p(1) * (r - betas0_p) )
            sigma_lp(i) = rtemp
            gip = -(gtemp+pars_p(1)) * rtemp

            rtemp = theta*exp(-pars_l(1) * (r - betas0_l) )
            sigma_pl = sigma_pl + rtemp
            gpi = (gtemp+pars_l(1)) * rtemp



        !--------------------MIXED PAIR POTENTIAL CONTRIUBUTION--------------------

            rtemp = theta*exp(-kappadbeta_p * (r - betas0_p))
            V_lp= V_lp + rtemp

            ! I'm not quite sure in this part. It might be a good idea to check the
            ! derivatives
            dV_lp(:,i) = (gtemp + kappadbeta_p)*rtemp*nip(:,i)

            rtemp = theta*exp(-kappadbeta_l * (r - betas0_l))
            V_pl = V_pl + rtemp
            ! I think the sum has to be added with - since the indices in nip have
            ! to be
            dV_pl(:,i) = (gtemp+kappadbeta_l)*rtemp*nip(:,i)

 !           write(*,'(3e20.10)') theta

            npi(:,i)= - nip(:,i)*gpi
            nip(:,i)=   nip(:,i)*gip

        case('morse')

            expar = exp(-pars_p(2)*(r-pars_p(3)))
            expar2 = expar**2
            energy = energy + pars_p(1) * (expar2 - 2*expar)

        case default

            print *, 'pes subroutine, particle: unknown particle pes'
            stop

        end select

    end do


select case(spec_l%pot)

    case ('emt')
    !-------------------------------CUT-OFF ENACTION-------------------------------
    ! Don't forget the gamma!

        sigma_ll = sigma_ll * igamma1l
        V_ll = V_ll * voldegamma2l
        s_l_exp = sigma_ll
        dV_lp= - (dV_lp*voldegamma2l*chilp+dV_pl*vopdegamma2p*chipl)*0.50d0
        force_l = force_l*voldegamma2l + dV_lp
        force_p(1) = sum(dV_lp(1,:))
        force_p(2) = sum(dV_lp(2,:))
        force_p(3) = sum(dV_lp(3,:))



    if (spec_p%pot == 'emt') then
        sigma_lp = sigma_lp * igamma1l
        sigma_pl = sigma_pl * igamma1p


        rtemp =chilp * pars_l(5) * igamma2l
        V_lp = V_lp * rtemp

        rtemp = pars_p(5) * igamma2p * chipl
        V_pl = V_pl * rtemp

        s_l_exp = sigma_ll+ chilp * sigma_lp

        ! Neutral Sphere Radius (s.b.)
        s_p  = -log( sigma_pl * chipl * twelfth) / ( beta * pars_p(1))


        ! Mixed Reference Pair Potential Contributions (s.b.)
        vref_p = 12.0d0 * pars_p(5) * exp( -pars_p(6) * s_p)

        ! Cohesive energy contribution due to the particle
        rtemp = exp(-pars_p(4) * s_p)* pars_p(3)
        Ecoh = ((1.0d0 + pars_p(4)*s_p) * rtemp - pars_p(3))

        f_p = (s_p * rtemp * pars_p(4)**2 + 0.5d0 * vref_p * pars_p(6)) &
              / (sigma_pl*beta*pars_p(1))

    end if


    !-----------------------------NEUTRAL SPHERE RADIUS----------------------------
    ! The neutral sphere radius is the radius in which the entire density of the
    ! atom is included
        s_l = -log( s_l_exp * twelfth ) / ( beta * pars_l(1))
        p_l=1.0d0 / (s_l_exp* beta * pars_l(1))

        f_l = s_l*p_l



   !----------------MIXED REFERENCE PAIR POTENTIAL CONTRIBUTIONS------------------
    ! These contributions have to be substracted to account for the contributions
    ! that were included twice.

    s_l_exp = exp( -pars_l(6) * s_l)
    rtemp = 12.0d0 * pars_l(5)
    vref_l =  rtemp * sum( s_l_exp )
    rtemp = -rtemp * pars_l(6)
    p_l= p_l * s_l_exp * rtemp

    !------------------------------------------------------------------------------
    !                           CALCULATING THE ENERGY
    !                           ======================
    !------------------------------------------------------------------------------


    !---------------------------COHESIVE ENERGY FUNCTION---------------------------
    ! Calculates and sums the contributions to the cohesive energy for both lattice
    ! and particle.

    s_l_exp = exp(-pars_l(4)*s_l)
    Ecoh = Ecoh &
            + sum( (1.0d0 + pars_l(4)*s_l)*s_l_exp-1.d0 )*pars_l(3)


    !-------------------------------CALCULATE Fi---------------------------------

    ! fi/f_l is the contribution to Ecoh_l that depends only on one index
    ! pj is the contribution to the reference energy that depends on one index
    f_l = pars_l(3) * pars_l(4)**2 * s_l_exp * f_l  - 0.50d0*p_l

    gij1 = gij1 * igamma1l
    nip = nip * igamma1l * chilp
    npi = npi * igamma1p

    k = 0

    do i= 1, spec_l%n
        do j = i+1, spec_l%n
            r3temp = gij1(:,k+j-i)  * (f_l(i) + f_l(j))
            force_l(:,i) = force_l(:,i) + r3temp
            force_l(:,j) = force_l(:,j) - r3temp

        end do
        k=k+spec_l%n-i

        r3temp = f_l(i)*nip(:,i) - f_p*npi(:,i)
        slab(i)%f = force_l(:,i)  + r3temp
        force_p   = force_p       - r3temp



   end do
    !-------------------------------OVERALL ENERGY---------------------------------
    ! Summation over all contributions.

    energy = energy+Ecoh - V_ll - 0.50d0 * ( V_lp + V_pl - vref_l - vref_p)
print *, energy
stop
    case('morse')
    ! EMPTY

    case default

        print *, 'pes subroutine: unknown pes'
        stop

end select

Epot = energy
teilchen(1)%f = force_p

end subroutine pes

subroutine fric_coef(r,fric_key,zeta)
    !
    ! Purpose:
    !       Here, the friction coefficient is calculated
    !       DON'T FORGET TO DEVIDE THE FRICTION THROUGH THE MASS!!!!
    !

    real(8), dimension(:,:), allocatable, intent(in) :: r
    integer, intent(in) :: fric_key
    real(8),intent(out) :: zeta ! Zeta is the friction coefficient.
    real(8) :: mass

    mass = spec_l%mass

    select case(fric_key)
    case(1)
        zeta = 1.0
    case default
        print *, 'This is not the friction you are looking for'
    end select

end subroutine fric_coef



end module force

