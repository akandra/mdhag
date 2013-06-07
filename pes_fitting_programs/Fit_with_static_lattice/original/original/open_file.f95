module open_file

    interface
        subroutine open_for_read(lun, file_name)
            integer, intent(in)     :: lun
            character, intent(in)   :: file_name
        end subroutine

        subroutine open_for_write(lun, file_name)
            integer, intent(in)     :: lun
            character, intent(in)   :: file_name
        end subroutine

        subroutine open_for_overwrite(lun, file_name)
            integer, intent(in)     :: lun
            character, intent(in)   :: file_name
        end subroutine
    endinterface

contains

end module

subroutine open_for_read(lun,file_name)
    implicit none
    integer, intent(in)           :: lun
    character(len=*), intent(in)  :: file_name

    integer                       :: ios
    character(len=120)            :: error_message

    open(unit=lun, file=file_name, status='old', action='read', iostat=ios, iomsg=error_message)
    if (ios==0) return
    print '(/ "Error on open file ", (a), " for read i/o status=", i4 )', TRIM(file_name), ios
    print '( "error message=", (a) )', error_message

    STOP 101

end subroutine open_for_read

subroutine open_for_write(lun,file_name)
    implicit none
    integer, intent(in)           :: lun
    character(len=*), intent(in)  :: file_name

    integer                       :: ios
    character(120)                :: error_message
    character                     :: answer

    open(unit=lun, file=file_name, status='new', action='write', iostat=ios, iomsg=error_message)
    if (ios==0) return

    print '( /"Error on open file ", (a), " for write i/o status=", i4 )', TRIM(file_name), ios
    print '( "error message: ", (a) )', error_message
    write (*, '( "overwrite existing file (y/n)? ")',advance='no')
    read(*,*) answer
    if (answer /='y' .and. answer/='Y') STOP 102
    print '((a)/)', 'OVERWRITING EXISTING FILES'
    open(unit=lun, file=file_name, status='replace', action='write', iostat=ios, iomsg=error_message)

    if (ios==0) return
    print *, 'failed to open file ', file_name, ' for write with status=replace. i/o status =',ios
    print *, 'error message: ', error_message
    STOP 103

end subroutine open_for_write

subroutine open_for_overwrite(lun,file_name)
    implicit none
    integer, intent(in)           :: lun
    character(len=*), intent(in)  :: file_name

    integer                       :: ios
    character(120)                :: error_message
    character                     :: answer

    open(unit=lun, file=file_name, status='unknown', action='write', iostat=ios, iomsg=error_message)
    if (ios==0) return

    print '( /"Error on open file ", (a), " for write i/o status=", i4 )', TRIM(file_name), ios
    print '( "error message: ", (a) )', error_message
    write (*, '( "overwrite existing file (y/n)? ")',advance='no')
    read(*,*) answer
    if (answer /='y' .and. answer/='Y') STOP 102
    print '((a)/)', 'OVERWRITING EXISTING FILES'
    open(unit=lun, file=file_name, status='replace', action='write', iostat=ios, iomsg=error_message)

    if (ios==0) return
    print *, 'failed to open file ', file_name, ' for write with status=replace. i/o status =',ios
    print *, 'error message: ', error_message
    STOP 103

end subroutine open_for_overwrite
