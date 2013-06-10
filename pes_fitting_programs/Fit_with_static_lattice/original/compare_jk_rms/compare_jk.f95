program compare_jk

    use EMTparms_class
    use open_file

    implicit none

    real(8)     :: jk_EMT_energy(1000)          ! jk EMT energy from hau11_plot.E.dat)
    real(8)     :: DFT(1000)                    ! DFT energy from hEMTfortran.dat based on f119 parameter set
    real(8)     :: r10(3,1000), r11(3)          ! arrays for H coordinate data from these two files
    integer     :: site10(1000), site11

    real(8)     :: sumsq, res                   ! for computation of rms error
    integer     :: ios10, ios11, i, j, npts

    call open_for_read(10,"hau111_plot.E.dat")  ! file with site,x,y,x,DFT energy
    call open_for_read(11,"hEMTfortran.dat")    ! file with jk EMT values

    ! *************************************************************************
    ! ***       GET THE NUMBER OF DATA POINTS BY READING TO END OF FILE     ***
    ! *************************************************************************
    do i=1,1000
        read(10,*,iostat=ios10)
        if ((ios10/=0) ) exit   !only test file 10 because it is one line shorter
    end do
    npts = i-1


    ! *************************************************************************
    ! ***       REWIND FILE AND READ THE DATA FROM BOTH FILES               ***
    ! ***       DISCARD POINTS WITH ABSOLUTE VALUE OF DFT ENERGY > 10       ***
    ! *************************************************************************
    rewind (10)
    j=1
    do i=1,npts
        read(10,*,iostat=ios10) site10(j), r10(:,j), DFT(j)
        read(11,*,iostat=ios11) site11   , r11(:)  , jk_EMT_energy(j)
        if (r11(1)/= r10(1,j) .or. r11(2)/= r10(2,j) .or. r11(3)/= r10(3,j)) then
            print *,'miscompare in coordinates j=', j
            print *,'   r10 =', r10(:,j)
            print *,'   r11 =', r11
            stop 201
        end if
        if( (abs(DFT(j))<=10) .and. (i/=npts) ) j=j+1  ! keep point only if abs(DFT)<10
    end do

    print *,'total number of points   =',npts
    print *,'points with abs(DFT)<=10 =', j
    npts=j

    ! *************************************************************************
    ! ***       CALULATE RMS ERROR FOR JK DFT AND FIT 119                   ***
    ! *************************************************************************
    sumsq=0

    call open_for_overwrite(8, "jk_compare.out")

    write(8,'((a))') '   i site, x     y     z      jk_EMT      DFT         Diff'
    do i=1,npts
        res=DFT(i)-jk_EMT_energy(i)
        write(8,'(i4,i3,3f6.2, 3f12.7)') i, site10(i), r10(:,i), DFT(i), jk_EMT_energy(i), res
        sumsq=sumsq + res**2
    end do

    write(*,*)'rms error=', sqrt(sumsq/npts)
    write(8,*)'rms error=', sqrt(sumsq/npts)

end program compare_jk
