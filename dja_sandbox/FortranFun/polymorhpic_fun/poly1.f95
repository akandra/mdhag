module poly1
   implicit none

   !interface
   !   subroutine fun1ir (i,r)
   !      implicit none
   !      integer :: i
   !      real    :: r
   !   end subroutine
   !
   !   subroutine fun1ri (r,i)
   !      implicit none
   !      integer :: i
   !      real    :: r
   !   end subroutine

   !end interface

   interface fun1
      module procedure fun1i
      module procedure fun1r
      module procedure fun1ir
      module procedure fun1ri
   end interface

contains
   subroutine fun1i (i)
      implicit none
      integer  :: i
      print '((a),i5)',  'i :', i
   end subroutine

   subroutine fun1r (i)
      implicit none
      real     :: i
      print '((a),f5.1)','r :', i
   end subroutine

   subroutine fun1ir (i,r)
      implicit none
      integer :: i
      real    :: r
      print '((a),i5,f5.1)','ir:', i, r
   end subroutine

   subroutine fun1ri (r,ii)
      implicit none
      integer :: ii
      real    :: r
      print '((a),f5.1,i5)','ri:', r, ii
   end subroutine



end module


