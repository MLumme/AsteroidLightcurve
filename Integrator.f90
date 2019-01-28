module Integrator

contains

  subroutine Integ(q_c,tri,alpha,brightness,q_rotaxis)
    use Quaternion_Math
    implicit none
    type(quaternion),intent(in) :: q_c(:)
    type(trig_type),intent(in) :: tri(:)
    double precision :: lightdir(3,3),u(3),v(3),w(3),Dir(3),Tempcoord(3)
    double precision :: proj(3),templight(3,3),x,y,z,pi,r1,r2,alpha,theta
    double precision :: len,len2,len3,a,b,c,brightsum,ny,ny0
    double precision,allocatable :: face(:,:),brightness(:,:)
    type(quaternion) :: q_rotaxis,q_temp1,q_temp2,q_rot
    integer :: i,j,k,nn,tstep,p1,p2,p3
    pi=4.0*atan(1.0)
    
    nn=size(tri)
    tstep=200
    allocate(face(4,nn),brightness(2,tstep+1))
    face(:,:)=0
    
    lightdir(:,1)=(/dble(1),dble(0),dble(0)/)        !incoming light
    lightdir(:,2)=(/cos(alpha),sin(alpha),dble(0)/)  !outgoing light
    lightdir(:,3)=(/dble(0),dble(0),dble(1)/)        !ecliptic normal
    
    do i=1,nn !normals for all faces
       p1=tri(i)%i(1)
       p2=tri(i)%i(2)
       p3=tri(i)%i(3)
       u=q_q2v(q_c(p2)-q_c(p1))
       v=q_q2v(q_c(p3)-q_c(p1))
       w=q_q2v(q_c(p2)-q_c(p3))
       Tempcoord=q_q2v(q_c(p3))

       !calculate face normal
       Dir=(/u(2)*v(3)-u(3)*v(2),u(3)*v(1)-u(1)*v(3),u(1)*v(2)-u(2)*v(1)/)
       len=sqrt(dot_product(Dir,Dir))
       len2=q_len(q_c(p1))

       r1=acos(dot_product(Dir,Tempcoord)/(len*len2)) !angle between normal
                              !and direction from origin
       if(r1.GE.pi/2)then     !flip normal if face normal points inwards
          Dir=(-1)*Dir
       end if

       face(1,i)=Dir(1)/len   !face normal unit vector
       face(2,i)=Dir(2)/len
       face(3,i)=Dir(3)/len
       
       a=sqrt(dot_product(u,u))
       b=sqrt(dot_product(v,v))
       c=sqrt(dot_product(w,w))

       face(4,i)=sqrt((a+b+c)/2*((a+b+c)/2-a)*& !area of face
            ((a+b+c)/2-b)*((a+b+c)/2-c))
       
    end do

    do i=1,tstep+1
       
       brightsum=0     
       
       theta=2*pi/(tstep)*(i-1.0)

       q_rot%q(1)=cos(theta/2)
       q_rot%q(2)=q_rotaxis%q(2)*sin(theta/2)
       q_rot%q(3)=q_rotaxis%q(3)*sin(theta/2)
       q_rot%q(4)=q_rotaxis%q(4)*sin(theta/2)

       len=q_len(q_rot)
       q_rot=q_rot/len

       do j=1,3          
   
          !rotates vectors pointing to source and observer and ecliptic normal
          q_temp1=dble(0.0)
          q_temp1%q(2)=lightdir(1,j)
          q_temp1%q(3)=lightdir(2,j)
          q_temp1%q(4)=lightdir(3,j)
          q_temp2=q_prod(q_prod(q_rot,q_temp1),conjugate(q_rot))
          
          templight(1,j)=q_temp2%q(2)
          templight(2,j)=q_temp2%q(3)
          templight(3,j)=q_temp2%q(4)

       end do

       u=templight(:,1) !assign rotated directions and normal matrix to vectors
       v=templight(:,2)
       w=templight(:,3)
       
       do j=1,nn  !repeat to all faces
          
          Dir=face(1:3,j)                        !extract face normal
          
          proj(:)=Dir-dot_product(Dir(:),w(:))/& !project surface normal on 
               dot_product(w(:),w(:))*w(:)       !ecliptic plane
          
          len=sqrt(dot_product(proj(:),proj(:))) !projection and source/observer
          len2=sqrt(dot_product(u(:),u(:)))      !vector lengths
          len3=sqrt(dot_product(v(:),v(:)))      
          
          r1=dot_product(proj(:),u(:))/(len*len2) !angles between face normal
          r2=dot_product(proj(:),v(:))/(len*len3) !projection&observer&source
          
          if((r1.GT.0).AND.(r2.GT.0))then  !select only faces visible to both
                                           !lightsource and observer
             ny=r2
             ny0=r1
             
             brightsum=brightsum+face(4,j)*ny*ny0/(ny+ny0)  !sum of brighness 
                                                            !from all faces
          end if
          
       end do

       brightness(1,i)=1.0/tstep*(i-1) !record fraction of full circle 
       brightness(2,i)=brightsum       !and corresponding brightness
    end do
    
    x=sum(brightness(2,:))/(tstep+1) !average brighness
    brightness(2,:)=brightness(2,:)/x !normalization to average brightness
                                      !return to main      
  end subroutine Integ
  
end module Integrator
