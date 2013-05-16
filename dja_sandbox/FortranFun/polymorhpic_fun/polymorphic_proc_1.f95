program main

   use poly1
   implicit none

   integer  :: i=7
   real     :: x=42
   real     :: y=41

   call fun1 (i)
   call fun1 (x)
   call fun1 (i,x)
   call fun1 (y,i)


end program main
