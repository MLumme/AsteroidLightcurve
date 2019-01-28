module Ellipsoid_Conversion

contains
  
  !converts spherical coordinates to cartesian quaternions
  function Q_coord(sphericals,ell_param,n)
    use Quaternion_Math
    implicit none
    double precision,intent(in),dimension(3) :: ell_param
    double precision,intent(in),dimension(2,n) :: sphericals
    type(quaternion),dimension(n) :: Q_coord
    integer,intent(in) :: n
    double precision :: theta,phi,r_trig,a,b,c
    integer :: i
    
    a=ell_param(1)
    b=ell_param(2)
    c=ell_param(3)
    
    do i=1,n
       theta=sphericals(1,i)
       phi=sphericals(2,i)
       
       r_trig=(b*c)**2*(sin(theta)*cos(phi))**2+ &
            (a*c)**2*(sin(theta)*sin(phi))**2+ &  
            (a*b)**2*cos(theta)**2 
       
       Q_coord(i)=dble(0.0)
       Q_coord(i)%q(2)=a*b*c/sqrt(r_trig)*sin(theta)*cos(phi)
       Q_coord(i)%q(3)=a*b*c/sqrt(r_trig)*sin(theta)*sin(phi)
       Q_coord(i)%q(4)=a*b*c/sqrt(r_trig)*cos(theta)
       
    end do
    
  end function Q_coord
  
end module Ellipsoid_Conversion
