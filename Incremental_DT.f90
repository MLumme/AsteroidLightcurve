module Tetrahedration
  
  type tetra_type
     integer :: i(4)
     double precision :: cs(4)
     logical :: l
  end type tetra_type
  
contains
  
  subroutine DT(n,q_c,tri,maxell,minell)
    use Determinant
    use Quaternion_math
    use Sameness
    implicit none
    type(tetra_type),allocatable :: tetra(:),dyntetra(:)
    type(trig_type),allocatable :: tri(:),tri_temp(:),cavity_wall(:)
    logical,allocatable :: tmask(:),dynmask(:),cavity_mask(:),tri_mask(:)
    type(quaternion),dimension(n) :: q_c
    double precision :: Dxx,Dyy,Dzz,Daa,Dcc,dist,R,maxell,minell
    double precision,dimension(4,4) :: Dx,Dy,Dz,Da,Dc
    double precision,dimension(4,n+4) :: Coord_Table
    integer,allocatable :: split(:)
    integer :: i,j,k,l,n,nn,ind,asize
    character(len=35) :: status="000.0% of points inserted in mesh."
    
    nn=1

    !Construct assist tetrahedron for boyer-watson algorithm    
    R=5*maxell    
    
    Coord_Table(1:4,1)=(/3*R**2,R,R,R/)
    Coord_Table(1:4,2)=(/3*R**2,-R,-R,R/)
    Coord_Table(1:4,3)=(/3*R**2,-R,R,-R/)
    Coord_Table(1:4,4)=(/3*R**2,R,-R,-R/)
    
    !Insert all points into coordinate table, and calculate 
    !dot product of coordinate
    do i=1,n
       Coord_Table(2:4,i+4)=q_c(i)%q(2:4)
       Coord_Table(1,i+4)=q_c(i)%q(2)**2+q_c(i)%q(3)**2+q_c(i)%q(4)**2
    end do
    
    !allocate arrays for tetrahedrons and mask for removal of split tetras
    allocate(tetra(1),tmask(1))
    tmask(:)=.TRUE.
    
    !insert assist tetras vertices into tetra-array 
    do i=1,4
       tetra%i(i)=i
    end do
    
    !mark all tetras as existing
    tetra%l=.TRUE.

    !repeat loop until all points inserted
    do i=1,n
       !repeat for all currently existing tetras
       do j=1,nn
          !calculate circumradius and -center for tetras using determinat
          if(tetra(j)%l)then
             
             Da=1.0
             Dx=1.0
             Dy=1.0
             Dz=1.0
             
             do k=1,4
                
                Da(1:3,k)=Coord_Table(2:4,tetra(j)%i(k))
                
                Dc(1:4,k)=Coord_Table(:,tetra(j)%i(k))
                
                Dx(1,k)=Coord_Table(1,tetra(j)%i(k))
                Dx(2:3,k)=Coord_Table(3:4,tetra(j)%i(k))
                
                Dy(1:2,k)=Coord_Table(1:2,tetra(j)%i(k))
                Dy(3,k)=Coord_Table(4,tetra(j)%i(k))
                
                Dz(1:3,k)=Coord_Table(1:3,tetra(j)%i(k))
                
             end do
             
             Dxx=DetGen(Dx)
             Dyy=(-1)*DetGen(Dy)
             Dzz=DetGen(Dz)
             Daa=DetGen(Da)
             Dcc=Detgen(Dc)

             tetra(j)%cs(1)=Dxx/(2*Daa)
             tetra(j)%cs(2)=Dyy/(2*Daa)
             tetra(j)%cs(3)=Dzz/(2*Daa)
             tetra(j)%cs(4)=sqrt(Dxx**2+Dyy**2+Dzz**2-4.0*Daa*Dcc)/ &
                  (2*abs(Daa))
             tetra(j)%l=.FALSE.
          end if
       end do
       
       ind=0
       allocate(split(nn))
       
       !check in which tetras the point is
       do j=1,nn
          
          dist=((tetra(j)%cs(1)-Coord_Table(2,i+4))**2+ &
               (tetra(j)%cs(2)-Coord_Table(3,i+4))**2+ &
               (tetra(j)%cs(3)-Coord_Table(4,i+4))**2)**0.5

          !mark tera/s for removal if point vithin
          if((dist-tetra(j)%cs(4)).LE.0)then
             ind=ind+1
             split(ind)=j
             tmask(j)=.FALSE.
          end if
          
       end do

       !increase the size of tetra-array
       allocate(dyntetra(nn+4*ind))
       dyntetra(1:nn)=tetra
       deallocate(tetra)
       call move_alloc(dyntetra,tetra)
       
       !mark new teras as without circumradii, add placeholder values
       do j=1,4
          tetra(nn+1:)%i(j)=0
          tetra(nn+1:)%cs(j)=0.0
       end do
       tetra(nn+1:)%l=.TRUE.
       
       !increase the size of tetramask-array
       allocate(dynmask(nn+4*ind))
       dynmask(1:nn)=tmask
       deallocate(tmask)
       call move_alloc(dynmask,tmask)
       
       tmask(nn+1:)=.FALSE.
       
       !if point in only on tetra split in 4 tetras 
       if(ind.EQ.1)then
          tetra(nn+1)%i(1)=tetra(split(1))%i(1)
          tetra(nn+1)%i(2)=tetra(split(1))%i(2)
          tetra(nn+1)%i(3)=tetra(split(1))%i(3)
          tetra(nn+1)%i(4)=i+4
          
          tetra(nn+2)%i(1)=tetra(split(1))%i(1)
          tetra(nn+2)%i(2)=tetra(split(1))%i(2)
          tetra(nn+2)%i(3)=tetra(split(1))%i(4)
          tetra(nn+2)%i(4)=i+4
          
          tetra(nn+3)%i(1)=tetra(split(1))%i(1)
          tetra(nn+3)%i(2)=tetra(split(1))%i(3)
          tetra(nn+3)%i(3)=tetra(split(1))%i(4)
          tetra(nn+3)%i(4)=i+4
          
          tetra(nn+4)%i(1)=tetra(split(1))%i(2)
          tetra(nn+4)%i(2)=tetra(split(1))%i(3)
          tetra(nn+4)%i(3)=tetra(split(1))%i(4)
          tetra(nn+4)%i(4)=i+4
          
          tmask(nn+1:nn+4)=.TRUE.
          
       !else construct cavity from group of tetras outer walls, 
       !connect point to outer walls forming new tetras
       else if(ind.GT.1)then
          allocate(cavity_wall(4*ind),cavity_mask(4*ind)) 
          
          do j=1,3
             cavity_wall(:)%i(j)=0
          end do
          
          cavity_mask=.TRUE.
          
          do j=1,ind
             k=(4*(j-1))
             l=split(j)
             
             cavity_wall(k+1)%i(:)=(/tetra(l)%i(1),tetra(l)%i(2),tetra(l)%i(3)/)
             cavity_wall(k+2)%i(:)=(/tetra(l)%i(1),tetra(l)%i(2),tetra(l)%i(4)/)
             cavity_wall(k+3)%i(:)=(/tetra(l)%i(1),tetra(l)%i(3),tetra(l)%i(4)/)
             cavity_wall(k+4)%i(:)=(/tetra(l)%i(2),tetra(l)%i(3),tetra(l)%i(4)/)
          end do

          do j=1,4*ind-1
             do k=j+1,4*ind
                if(Arreq(cavity_wall(j)%i(:),cavity_wall(k)%i(:)))then
                   cavity_mask(j)=.FALSE.
                   cavity_mask(k)=.FALSE.
                   exit
                end if
             end do
          end do

          do j=1,4*ind
             if(cavity_mask(j))then
                nn = nn + 1
                tetra(nn)%i(1) = cavity_wall(j)%i(1)
                tetra(nn)%i(2) = cavity_wall(j)%i(2)
                tetra(nn)%i(3) = cavity_wall(j)%i(3)
                tetra(nn)%i(4)=i+4
                tmask(nn)=.TRUE.
             end if
          end do

          deallocate(cavity_wall,cavity_mask)
          
       end if

       !remove old split tetras, unnecesary places
       tetra=pack(tetra,tmask)
       
       !nn=current number of tetras
       nn=size(tetra)
       
       !resize tetramask for true number of tetras
       deallocate(tmask)
       allocate(tmask(nn))
       tmask=.TRUE.
       
       !information for user
       write(unit=status(1:5),fmt="(f5.1)") float(i)/float(n)*100
       write(unit=6,fmt="(a1,a)",advance="no") char(13),adjustl(trim(status))
       
       deallocate(split)
       
    end do
    
    !mask for removal of interior triangles, triangles connected to assist verts
    tmask=.FALSE.
    
    !mark all tetras connected to assist tetras for saving
    do i=1,nn
       if(ANY(tetra(i)%i(:).EQ.1) &
            .OR.ANY(tetra(i)%i(:).EQ.2) &
            .OR.ANY(tetra(i)%i(:).EQ.3) &
            .OR.ANY(tetra(i)%i(:).EQ.4))then
          tmask(i)=.TRUE.
       end if
    end do
    
    !repack removing interior tetras
    tetra=pack(tetra,tmask)
    asize=size(tetra)
    k=4*asize 

    deallocate(tmask) 

    !allocate arrays for triangles
    allocate(tri_temp(4*asize),tri_mask(4*asize))
    tri_mask(:)=.TRUE.

    !split tetras into triangles
    do i=1,asize
       j=4*(i-1)
       tri_temp(1+j)%i(:)=(/tetra(i)%i(1),tetra(i)%i(2),tetra(i)%i(3)/)
       tri_temp(2+j)%i(:)=(/tetra(i)%i(1),tetra(i)%i(2),tetra(i)%i(4)/)
       tri_temp(3+j)%i(:)=(/tetra(i)%i(1),tetra(i)%i(3),tetra(i)%i(4)/)
       tri_temp(4+j)%i(:)=(/tetra(i)%i(2),tetra(i)%i(3),tetra(i)%i(4)/)
    end do
    
    !check if triangle connected to assist vertices, mark for removal
    do i=1,k
       if(ANY(tri_temp(i)%i(:).EQ.1) &
            .OR.ANY(tri_temp(i)%i(:).EQ.2) &
            .OR.ANY(tri_temp(i)%i(:).EQ.3) &
            .OR.ANY(tri_temp(i)%i(:).EQ.4))then
          tri_mask(i)=.FALSE.
       end if
    end do

    !return surface triangles to main, first adjust vertex numbering for lack
    !of assist tetras
    tri_temp=pack(tri_temp,tri_mask)
    
    k=size(tri_temp)

    do i=1,k
       do j=1,3
          tri_temp(i)%i(j)=tri_temp(i)%i(j)-4
       end do
    end do
    
    allocate(tri(k))
    
    tri=tri_temp
    
    write(unit=6,fmt="(A1)")
    write(6,fmt="(A23)")"Triangulation complete."
    
  end subroutine DT
  
end module Tetrahedration
