module Vertex_Generator

contains

  function VertGen(n,ell_param)
    integer,intent(in) :: n
    double precision,intent(in),dimension(3) :: ell_param
    double precision,dimension(2,n) :: VertGen
    integer :: i,time(8),seed(12)
    double precision :: C_max,a,b,c,pi,phi,theta,y
    double precision :: polar(2,3),C_poles(3)
    pi=4.0*atan(1.0)

    a=ell_param(1)
    b=ell_param(2)
    c=ell_param(3)
    
    !coordinates for maximum C-values at poles
    polar(:,1)=(/(0.5*pi),(0.5*pi)/)
    polar(:,2)=(/dble(0.0),dble(0.0)/)
    polar(:,3)=(/(0.5*pi),dble(0.0)/)
    
    
    !find maximum C       
    do i=1,3
       phi=polar(1,i)
       theta=polar(2,i)
       C_poles(i)=4*pi*p_func(a,b,c,theta,phi)
    end do

    C_max=maxval(abs(C_poles))

    !random number generation seed from computer clock
    call date_and_time(values=time)
    seed(1)=time(4)*(360000*time(5)+6000*time(6)+100*time(7)+time(8))
    call random_seed(put=seed)

    !keep producing andom number of points until vanted number of vertexes is 
    !reached
    i=1
    do while(i.LE.n)
       call random_number(theta)
       call random_number(phi)
       call random_number(y)

       theta=pi*theta
       phi=2*pi*phi

       !test if cordinates fulfill propabilty distribution requirements
       if(4*pi*p_func(a,b,c,theta,phi) .GE. C_max*y)then
          VertGen(1,i)=theta
          VertGen(2,i)=phi
          i=i+1
       end if
    end do
 
  end function VertGen

  !splits probabilty distribution in more easy to manage functions
  function p_func(a,b,c,theta,phi)
    double precision,intent(in) :: a,b,c,theta,phi
    double precision :: r_trig,p_func
    
    r_trig=(b*c*sin(theta)*cos(phi))**2+ &
         (a*c*sin(theta)*sin(phi))**2+ &  
         (a*b*cos(theta))**2 
    
    p_func=((a*b*c)*(r_trig)**(-0.5))**2* &
         (1+(-sin(2*theta)*((-a**2*b**2)+ &
         (a**2*c**2*cos(phi)**2)+ &
         (b**2*c**2*cos(phi)**2))/(2*r_trig))**2+ &
         (-sin(2*phi)*sin(theta)*(a**2*c**2-b**2*c**2)/ &
         (2*r_trig))**2)**0.5
    
  end function p_func
  
end module Vertex_Generator
