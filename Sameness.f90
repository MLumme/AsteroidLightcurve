!Test if two size 3 arrays are composed of the same vertices
module Sameness

contains
    function Arreq(a,b)
    implicit none
    integer,dimension(3),intent(in) :: a,b
    logical :: Arreq
    integer :: i,j,vcount
    
    vcount=0
    
    do i=1,3
       do j=1,3
          if(a(i).EQ.b(j))then
             vcount=vcount+1
          end if
       end do
    end do
    
    if(vcount.EQ.3)then
       Arreq=.TRUE.
    else
       Arreq=.FALSE.
    end if
    
  end function Arreq
end module Sameness
