module Determinant
contains

  recursive function DetGen(M) result(D)
    implicit none
    double precision,dimension(:,:) :: M
    double precision,allocatable :: N(:,:)
    double precision :: D
    integer :: s,x
    logical,allocatable :: mask(:,:)

    s=size(M,1)
    allocate(mask(s,s),N(s-1,s-1))

    !Determinant module for square matrix of any size
    if(s==2)then	
       D=M(1,1)*M(2,2)-M(1,2)*M(2,1)
    else
       D=0
       do x=1,s
          mask=.true.
          mask(1,:)=.false.
          mask(:,x)=.false.
          N=reshape(pack(M,mask),[s-1,s-1])
          D=D+(-1.0)**(1.0+x)*M(1,x)*DetGen(N)
       end do
    end if
  end function DetGen

end module Determinant
