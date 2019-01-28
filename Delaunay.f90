program Delaunay_Triangulation
  use Vertex_Generator
  use Ellipsoid_Conversion
  use Tetrahedration
  use Quaternion_Math
  use Integrator
  implicit none

  character(len=70) :: filename
  integer :: n,i,j,rotax,ios
  type(trig_type),allocatable :: tri(:)
  double precision,allocatable :: sphericals(:,:),brightness(:,:)
  double precision,dimension(3) :: ell_param,rot
  double precision :: maxell,minell,pi,alpha,t1,t2,t3,t4,t5,t6,scale,len
  type(quaternion),allocatable :: q_c(:)
  type(quaternion) :: q_r,q_rotaxis
  
  pi=4.0*atan(1.0)
  

  !=================================================================
  !parameter input and choice of rotational axis
  write(*,*) "Give ellipsoid parameters (a,b,c)."
  do
     read(*,*) ell_param

     rotax=minloc(ell_param,1) !find shortest axis for rotation
     
     j=0                       !check if minimum axis is not unique
     do i=1,3
        if(ell_param(rotax).EQ.ell_param(i))then
           j=j+1
        end if
     end do
        
     if(any(ell_param.LE.0))then !fail if flat or negative ellipse
        write(*,*) "Ellipsoid parameter cannot be less than or equal to zero, re-enter parameters."

     else if(j.EQ.2)then !request manual choice of rotation axis if two minimum axii of same length
        write(*,*) "Two or more ellipsoid parameters are the same," 
        write(*,*) "please insert which ellipsoid axis wil be rotational axis, 1 for A, etc."
        read(*,*) rotax !override position from minloc with manual input
        exit
     else
        exit
     end if
  end do
  
  q_rotaxis=dble(0.0)
  do i=1,3
     if(i.EQ.rotax)then
        q_rotaxis%q(i+1)=dble(1.0)
     end if
  end do
  
  write(*,*) "Enter the number of vertexes."
  do  
     read(*,*) n
     if(n.LE.4)then !ensure that there is atleast one tetrahedra
        write(*,*) "Number of vertexes cant be less or equal to four, please retry."
     else
        exit
     end if
  end do
  
  write(*,*) "Give spherical coordinates of ellipsoid pole in degrees, and rotation round polar axis (lat,long,phi)."
  read(*,*) rot
  
  write(*,*) "Give phase angle in degrees."
  read(*,*) alpha
  
  write(*,*) "Give filename for result array." !get filename, attempt to open
  do
     read(*,*) filename
     open(unit=1,file=filename,status='new',iostat=ios)
     if(ios.GT.0)then  !failstate, return to beginning of do loop
        write(*,*) "Error opening file or exists allready, retry."
     else 
        exit
     end if
  end do

  
  !===========================================================
  !shift asteroid rotation and phase angle to radians
  rot(1)=rot(1)/180.*pi
  rot(2)=rot(2)/180.*pi
  rot(3)=rot(3)/180.*pi
  alpha=alpha/180*pi

  !===========================================================
  !generation of vertices to triangulate
  
  scale=maxval(abs(ell_param)) !find maximum ellipsoid radius
  ell_param(:)=ell_param(:)/scale !normalize to largest axis
  maxell=maxval(abs(ell_param))   
  minell=minval(abs(ell_param))

  allocate(sphericals(2,n),q_c(n)) !allocate room for all points
  
  sphericals=VertGen(n,ell_param)      !call vertex spherical coordinate 
                                       !generator
  q_c=q_coord(sphericals,ell_param,n) !call spherical to quaternion
                                        !converter

  !==========================================================
  !Triangulate
  call DT(n,q_c,tri,maxell,minell)              !call surface trangulation 

  !==========================================================
  !rotating surface vertices and rotational axis to inserted orientation
  t1=cos(rot(2)*0.5)
  t2=sin(rot(2)*0.5)
  t3=cos(rot(3)*0.5)
  t4=sin(rot(3)*0.5)
  t5=cos(rot(1)*0.5)
  t6=sin(rot(1)*0.5)

  q_r%q(1)=t1*t3*t5+t2*t4*t6
  q_r%q(2)=t1*t4*t5-t2*t3*t6
  q_r%q(3)=t1*t3*t6+t2*t4*t5
  q_r%q(4)=t2*t3*t5+t1*t4*t6

  len=q_len(q_r)
  q_r=q_r/len

  do i=1,n
     q_c(i)=q_prod(q_prod(q_r,q_c(i)),conjugate(q_r))
  end do

  q_rotaxis=q_prod(q_prod(q_r,q_rotaxis),conjugate(q_r))
     
  !=========================================================
  !lightcurve integration

  call Integ(q_c,tri,alpha,brightness,q_rotaxis) !call brightness integrator  
  
  n=size(brightness)/2            !write fraction of full rotation and
  do i=1,n                        !corresponding brightness on file
     write(1,*) brightness(:,i)
  end do

  !Code to print files for ocatves trimesh for triangulation performace check
  !n=size(tri)
  !do i=1,n        
  !   write(1,*) (tri(i)%i(j),j=1,3)
  !end do
  !
  !open(unit=2,file='vpoints',iostat=ios)  
  !
  !do i=1,n
  !   write(2,*) (q_c(i)%q(j),j=2,4)
  !end do

  
  
end program Delaunay_Triangulation
