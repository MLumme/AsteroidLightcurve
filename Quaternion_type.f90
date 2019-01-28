!Quaternion type, required mathematical operations and conversions fom quaternions to vectors
!and vice versa. 
module Quaternion_Math

  implicit none

  type trig_type
     integer :: i(3)
  end type trig_type

  type quaternion
     double precision :: q(4)
  end type quaternion

  interface operator(+)
     module procedure q_add
  end interface operator(+)

  interface operator(-)
     module procedure q_sub
  end interface operator(-)

  interface operator(*)
     module procedure q_multi
  end interface operator(*)

  interface operator (/)
     module procedure q_div
  end interface operator (/)

  interface assignment(=)
     module procedure q_eq
  end interface assignment(=)

  interface operator (.len.)
     module procedure q_len
  end interface operator (.len.)

  interface operator(.vtoq.)
     module procedure q_v2q
  end interface operator(.vtoq.)

  interface operator(.qtov.)
     module procedure q_q2v
  end interface operator(.qtov.)

contains
  
  function conjugate(q1) result(q_conj)
    implicit none
    type(quaternion),intent(in) :: q1
    type(quaternion) :: q_conj
    integer :: i

    q_conj = q1
    do i=2,4
       q_conj%q(i) = -1*q1%q(i)
    end do
    
  end function conjugate

  function q_add(q1,q2) result(q_out)
    implicit none
    type(quaternion), intent(in) :: q1,q2
    type(quaternion) :: q_out
    
    q_out%q(1) = q1%q(1) + q2%q(1)
    q_out%q(2) = q1%q(2) + q2%q(2)
    q_out%q(3) = q1%q(3) + q2%q(3)
    q_out%q(4) = q1%q(4) + q2%q(4)
    
  end function q_add

  function q_sub(q1,q2) result(q_out)
    implicit none
    type(quaternion), intent(in) :: q1,q2
    type(quaternion) :: q_out
    
    q_out%q(1) = q1%q(1) - q2%q(1)
    q_out%q(2) = q1%q(2) - q2%q(2)
    q_out%q(3) = q1%q(3) - q2%q(3)
    q_out%q(4) = q1%q(4) - q2%q(4)

  end function q_sub
  
  function q_prod(q1,q2) result(q_out)
    implicit none
    type(quaternion), intent(in) :: q1,q2
    type(quaternion) :: q_out
    
    q_out%q(1) = q1%q(1)*q2%q(1) - q1%q(2)*q2%q(2) &
         - q1%q(3)*q2%q(3) - q1%q(4)*q2%q(4)
    q_out%q(2) = q1%q(1)*q2%q(2) + q1%q(2)*q2%q(1) &
         + q1%q(3)*q2%q(4) - q1%q(4)*q2%q(3)
    q_out%q(3) = q1%q(1)*q2%q(3) - q1%q(2)*q2%q(4) &
         + q1%q(3)*q2%q(1) + q1%q(4)*q2%q(2)
    q_out%q(4) = q1%q(1)*q2%q(4) + q1%q(2)*q2%q(3) &
         - q1%q(3)*q2%q(2) + q1%q(4)*q2%q(1)
    
  end function q_prod

  type(quaternion) function q_multi(q1,x)
    implicit none
    type(quaternion),intent(in) :: q1
    double precision,intent(in) :: x
    integer :: i

    do i=1,4
       q_multi%q(i) = q1%q(i)*x
    end do

  end function q_multi

  type(quaternion) function q_div(q1,x)
    implicit none
    type(quaternion),intent(in) :: q1
    double precision,intent(in) :: x
    integer :: i
    
    do i=1,4
       q_div%q(i) = q1%q(i)/x
    end do
    
  end function q_div

  double precision function q_len(q1)
    implicit none
    type(quaternion),intent(in) :: q1
    
    q_len = (q1%q(1)**2 + q1%q(2)**2 + q1%q(3)**2 + q1%q(4)**2)**0.5

  end function q_len

  subroutine q_eq(q1,x)
    implicit none
    type(quaternion),intent(out) :: q1
    double precision,intent(in) :: x
    integer :: i

    do i=1,4
       q1%q(i) = x
    end do
    
  end subroutine q_eq

  function q_v2q(v1) result(q1)
    implicit none
    type(quaternion) :: q1
    double precision,dimension(3),intent(in) :: v1
    
    q1%q(1) = 0
    q1%q(2) = v1(1)
    q1%q(3) = v1(2)
    q1%q(4) = v1(3)

  end function q_v2q

  function q_q2v(q1) result(v1)
    implicit none
    double precision,dimension(3) :: v1
    type(quaternion),intent(in) :: q1

    v1(1) = q1%q(2)
    v1(2) = q1%q(3)
    v1(3) = q1%q(4)
    
  end function q_q2v

end module Quaternion_Math
